#!/usr/bin/env python3
"""Prove `probe_consensus.py` recovers answers that were planted on purpose.

The probe will be run once, on a database this repository cannot see, and its
output believed. So every claim it makes is checked here first against a
synthetic database whose answers are known: a fixed number of A1 picks, a fixed
number of A2, a fixed calibration offset, a modification the calculated column
does not carry, and two load-response slopes. If the probe cannot recover those,
nothing it says about real data is worth reading.

The case that matters most is `deamidation_is_not_an_isotope`. Deamidation is
+0.9840 Da and an A1 isotope is +1.0034: 0.0193 Da apart. A tolerance wide
enough to be safe about calibration is wide enough to count every deamidated
peptide as an isotope pick, and the PAX8 run alone had 20 deamidations.

    python3 tools/test_probe_consensus.py
    python3 tools/test_probe_consensus.py -v

No pytest, no dependencies. Exit status is 1 if any case fails.
"""

import argparse
import json
import math
import os
import random
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROBE = ROOT / "tools" / "probe_consensus.py"
sys.path.insert(0, str(ROOT / "tools"))
import probe_consensus as pc  # noqa: E402

RESIDUE = pc.__dict__.get("RESIDUE") or {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259,
    "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111,
    "Y": 163.06333, "W": 186.07931,
}
WATER = 18.010565

# Planted truth. Every assertion below is stated against these numbers.
CALIB_PPM = 3.2
N_IONS, N_A1, N_A2, N_OX, N_DEAM = 500, 33, 5, 40, 20
LOADS = (25.0, 50.0, 120.0, 240.0)
SLOW, FAST = 0.25, 1.00
OX, DEAM = 15.994915, 0.984016


def plant(path, calib_ppm=CALIB_PPM):
    """A consensus-shaped database with known answers."""
    rng = random.Random(20261016)
    ions, seen = [], set()
    while len(ions) < N_IONS:
        body = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(rng.randint(6, 19)))
        seq, z = body + rng.choice("KR"), rng.choice([2, 2, 3, 3, 4])
        if (seq, z) in seen:
            continue
        seen.add((seq, z))
        ions.append((seq, z))

    a1 = set(ions[:N_A1])
    a2 = set(ions[N_A1:N_A1 + N_A2])
    ox = set(ions[N_A1 + N_A2:N_A1 + N_A2 + N_OX])
    deam = set(ions[N_A1 + N_A2 + N_OX:N_A1 + N_A2 + N_OX + N_DEAM])
    slow = set(ions[::2])

    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL, CalcMass REAL,
        PrecursorIntensity REAL)""")
    # Same shape, far fewer rows: the probe must not choose this one.
    con.execute("CREATE TABLE scratch (PeptideSequence TEXT, Charge INTEGER, "
                "PrecMz REAL, CalcMz REAL)")
    con.execute("INSERT INTO scratch VALUES ('PEPTIDEK', 2, 500.0, 500.0)")

    rid = 0
    for seq, z in ions:
        mass = WATER + sum(RESIDUE[c] for c in seq)
        calc_mz = mass / z + pc.PROTON
        amp = 10 ** rng.uniform(3.0, 6.0)
        beta = SLOW if (seq, z) in slow else FAST
        shift = (pc.C13 if (seq, z) in a1 else 2 * pc.C13 if (seq, z) in a2
                 else OX if (seq, z) in ox else DEAM if (seq, z) in deam else 0.0)
        for ng in LOADS:
            if rng.random() < 0.06:
                continue                      # ragged coverage, as real data is
            reported = (calc_mz + shift / z) * (1 + (calib_ppm + rng.gauss(0, 0.3)) / 1e6)
            rid += 1
            con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?)",
                        (rid, seq, z, "HeLa_%dng" % int(ng), reported, calc_mz, mass,
                         amp * (ng ** beta) * math.exp(rng.gauss(0, 0.05))))
    con.commit()
    con.close()


def run(*args, **kw):
    env = dict(os.environ)
    if "home" in kw:
        env["HOME"] = str(kw["home"])
    proc = subprocess.run([sys.executable, str(PROBE)] + list(args),
                          capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout + proc.stderr


def bucket(out, label):
    """(rows, ions) from the probe's offset table, or None if absent."""
    for line in out.splitlines():
        parts = line.split()
        if parts[:1] == [label] or (label.count(" ") and line.strip().startswith(label)):
            rest = line.strip()[len(label):].split()
            try:
                return int(rest[0]), int(rest[1])
            except (IndexError, ValueError):
                return None
    return None


# --- cases ------------------------------------------------------------------
# Each returns a list of failure strings; empty means the case passed.

def classifier_boundaries(_db, _out):
    bad = []
    for delta, want in ((0.0, "monoisotopic"),
                        (pc.TOL * 0.9, "monoisotopic"),
                        (pc.TOL * 2, "unexplained"),
                        (pc.C13, "A1"), (2 * pc.C13, "A2"), (3 * pc.C13, "A3"),
                        (DEAM, "deamidation"), (OX, "oxidation"),
                        (OX + pc.C13, "oxidation + A1"), (0.5, "unexplained"),
                        (1.97, "unexplained"), (-pc.C13, "below mono (-A1)")):
        got = pc.classify(delta)[0]
        if got != want:
            bad.append("classify(%+.4f) = %r, expected %r" % (delta, got, want))
    return bad


def deamidation_is_not_an_isotope(_db, _out):
    gap = pc.C13 - DEAM
    bad = []
    if pc.classify(DEAM)[0] != "deamidation":
        bad.append("a bare deamidation delta was not called deamidation")
    if pc.TOL >= gap / 2:
        bad.append("TOL %.4f is at least half the %.4f Da deamidation/A1 gap, so the "
                   "two windows overlap" % (pc.TOL, gap))
    if pc.BOUNDARY_MARGIN >= gap:
        bad.append("BOUNDARY_MARGIN %.4f is wider than the %.4f Da gap, so every A1 row "
                   "would be flagged ambiguous" % (pc.BOUNDARY_MARGIN, gap))
    return bad


def recovers_planted_isotope_picks(_db, out):
    bad = []
    for label, want_ions in (("A1", N_A1), ("A2", N_A2),
                             ("oxidation", N_OX), ("deamidation", N_DEAM),
                             ("monoisotopic", N_IONS - N_A1 - N_A2 - N_OX - N_DEAM)):
        got = bucket(out, label)
        if got is None:
            bad.append("no %r row in the offset table" % label)
        elif got[1] != want_ions:
            bad.append("%s: %d ions, planted %d" % (label, got[1], want_ions))
    if "A1 only:" not in out:
        bad.append("no A1-only summary line")
    else:
        line = [l for l in out.splitlines() if "A1 only:" in l][0]
        if "%d ions" % N_A1 not in line:
            bad.append("A1-only summary disagrees with the planted %d: %r" % (N_A1, line.strip()))
    return bad


def recovers_planted_calibration(_db, out):
    bad = []
    if "Pooled median" not in out:
        return ["no pooled calibration line"]
    line = [l for l in out.splitlines() if "Pooled median" in l][0]
    got = float(line.split("Pooled median")[1].split("ppm")[0])
    if abs(got - CALIB_PPM) > 0.2:
        bad.append("calibration %.2f ppm, planted %.1f" % (got, CALIB_PPM))
    runs = [l for l in out.splitlines() if l.strip().startswith("HeLa_")]
    if len(runs) != len(LOADS):
        bad.append("%d per-run rows, expected %d" % (len(runs), len(LOADS)))
    return bad


def calibration_resolves_the_boundary(_db, out):
    """With the offset removed, nothing should sit in the coin-flip zone."""
    if "CAUTION:" in out:
        line = [l for l in out.splitlines() if "CAUTION:" in l][0]
        return ["rows still ambiguous after correction: %r" % line.strip()]
    if "cleanly separated" not in out:
        return ["no statement either way about the deamidation/A1 boundary"]
    return []


def recovers_planted_load_slopes(_db, out):
    bad = []
    if "quartiles" not in out:
        return ["no slope quartiles line"]
    line = [l for l in out.splitlines() if "quartiles" in l][0]
    lo, hi = (float(t) for t in line.split("quartiles")[1].split("to"))
    if abs(lo - SLOW) > 0.05:
        bad.append("lower quartile %+.3f, planted %.2f" % (lo, SLOW))
    if abs(hi - FAST) > 0.05:
        bad.append("upper quartile %+.3f, planted %.2f" % (hi, FAST))
    return bad


def picks_the_richer_table(_db, out):
    if "table    observations" not in out:
        return ["chose the wrong table: %r" % out.splitlines()[2:3]]
    return []


def refuses_a_published_shaped_table(db, tmp):
    """One already-corrected m/z and nothing to compare it with: must refuse."""
    db = tmp / "published.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE psms (PeptideSequence TEXT, Charge INT, mz REAL, inten REAL)")
    con.executemany("INSERT INTO psms VALUES (?,?,?,?)",
                    [("PEPTIDEK", 2, 500.2, 1e5)] * 50)
    con.commit()
    con.close()
    rc, out = run("--db", str(db))
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1" % rc)
    if "A1" in out:
        bad.append("claimed something about isotopes with no calculated column")
    return bad


def refuses_a_file_that_is_not_a_database(db, tmp):
    bogus = tmp / "bogus.db"
    bogus.write_text("not a database")
    rc, out = run("--db", str(bogus))
    if rc != 2:
        return ["exit %d, expected 2" % rc]
    if "Could not open" not in out:
        return ["no reason given: %r" % out.strip()]
    return []


def load_test_skips_without_a_map(db, tmp):
    """Test 3 must say why it did not run, not quietly omit itself."""
    rc, out = run("--db", str(db))
    bad = []
    if "TEST 3" not in out:
        bad.append("test 3 vanished instead of reporting a skip")
    elif "SKIPPED" not in out.split("TEST 3")[1]:
        bad.append("test 3 ran without a load map")
    if rc != 0:
        bad.append("exit %d, expected 0 when no load map was asked for" % rc)
    return bad


def find_locates_the_database(db, tmp):
    """--find must pick the real one and ignore the places it should not look."""
    home = tmp / "fakehome"
    for sub in ("Documents/proteomics", "data", "Library/Caches/junk", "node_modules/pkg"):
        (home / sub).mkdir(parents=True, exist_ok=True)
    real = home / "Documents" / "proteomics" / "ARiA_Consensus.db"
    real.write_bytes(Path(db).read_bytes())
    # A real SQLite file, but tiny: must be listed and must not win.
    small = sqlite3.connect(home / "data" / "aria_old.sqlite")
    small.execute("CREATE TABLE t (a INT)")
    small.commit()
    small.close()
    # Must be ignored: wrong header, and two skipped directories.
    (home / "data" / "aria_notes.db").write_text("not a database")
    (home / "Library" / "Caches" / "junk" / "ARiA_Consensus.db").write_bytes(
        Path(db).read_bytes())
    (home / "node_modules" / "pkg" / "aria_consensus.db").write_bytes(
        Path(db).read_bytes())

    rc, out = run("--find", "--schema-only", home=home)
    bad = []
    if rc != 0:
        bad.append("exit %d, expected 0" % rc)
    if "<- using this" not in out:
        bad.append("did not say which database it chose")
    chosen = [l for l in out.splitlines() if "<- using this" in l]
    if chosen and str(real) not in chosen[0]:
        bad.append("chose %r, expected the one under Documents" % chosen[0].strip())
    if "aria_notes.db" in out:
        bad.append("offered a file that is not a SQLite database")
    if "Library" in out:
        bad.append("walked into Library, which is on the skip list")
    if "node_modules" in out:
        bad.append("walked into node_modules, which is on the skip list")
    if "aria_old.sqlite" not in out:
        bad.append("did not list the smaller candidate, so a wrong pick would be invisible")
    return bad


def find_says_so_when_nothing_matches(db, tmp):
    empty = tmp / "emptyhome"
    empty.mkdir(exist_ok=True)
    rc, out = run("--find", home=empty)
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1" % rc)
    if "found none" not in out:
        bad.append("did not say it found nothing: %r" % out.strip()[:80])
    if "--db" not in out:
        bad.append("did not point at --db as the way out")
    return bad


def neither_db_nor_find_is_an_error(db, tmp):
    rc, out = run("--schema-only")
    if rc == 0:
        return ["ran with no database at all"]
    if "--find" not in out:
        return ["the error did not mention --find: %r" % out.strip()[:80]]
    return []


def _mini_db(path, shift_da, rng_seed=7, calib_ppm=2.0):
    """A small observations table where every row carries the same mass shift."""
    rng = random.Random(rng_seed)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL)""")
    for i in range(400):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "K"
        z = rng.choice([2, 3])
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        reported = (calc_mz + shift_da / z) * (1 + (calib_ppm + rng.gauss(0, 0.3)) / 1e6)
        con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)",
                    (i + 1, seq, z, "run_A", reported, calc_mz))
    con.commit()
    con.close()


def a_mass_gated_table_is_not_evidence(db, tmp):
    """Zero picks from a table that cannot hold one must be refused as a result.

    This is the finding that came back from the real database: every row had
    passed an ingest tolerance, so an isotope pick could never have been stored,
    and the probe's zero measured the gate rather than the pipeline.
    """
    gated = tmp / "gated.db"
    _mini_db(gated, shift_da=0.0)
    rc, out = run("--db", str(gated))
    bad = []
    if rc != 0:
        bad.append("exit %d, expected 0" % rc)
    if "Any isotope:   0 rows" not in out:
        bad.append("expected zero isotope picks in a table with none")
    if "THIS RESULT IS NOT EVIDENCE" not in out:
        bad.append("read its own zero as a finding on a table that cannot hold a pick")
    for want in ("upstream", "quarantined"):
        if want not in out:
            bad.append("did not say where to ask instead (%r missing)" % want)
    return bad


def an_unfiltered_zero_carries_weight(db, tmp):
    """Zero picks where large errors DO occur is a real negative, not a vacuous one."""
    wide = tmp / "wide.db"
    _mini_db(wide, shift_da=15.994915)   # oxidation: hundreds of ppm out, not an isotope
    rc, out = run("--db", str(wide))
    bad = []
    if "THIS RESULT IS NOT EVIDENCE" in out:
        bad.append("dismissed a zero that was not excluded by construction")
    if "carries weight" not in out:
        bad.append("did not credit a zero measured on a table that allows large errors")
    if "Rows beyond" not in out:
        bad.append("did not report how many rows sit beyond the gate width")
    return bad


def limit_reads_a_spread_not_the_oldest(db, tmp):
    """--limit must not silently mean 'the oldest rows'."""
    rc, out = run("--db", str(db), "--limit", "200")
    bad = []
    if rc != 0:
        bad.append("exit %d, expected 0" % rc)
    if "spread evenly across the table" not in out:
        bad.append("did not say which rows a limited read covered")
    rc2, out2 = run("--db", str(db), "--limit", "200", "--order", "oldest")
    if "oldest" not in out2:
        bad.append("--order oldest did not say so")
    return bad


def schema_only_names_every_table(db, tmp):
    rc, out = run("--db", str(db), "--schema-only")
    bad = []
    if rc != 0:
        bad.append("exit %d, expected 0" % rc)
    for want in ("observations", "scratch", "raw_mz=PrecMz", "calc_mz=CalcMz"):
        if want not in out:
            bad.append("--schema-only did not mention %r" % want)
    return bad


NEEDS_RUN = (classifier_boundaries, deamidation_is_not_an_isotope,
             recovers_planted_isotope_picks, recovers_planted_calibration,
             calibration_resolves_the_boundary, recovers_planted_load_slopes,
             picks_the_richer_table)
STANDALONE = (refuses_a_published_shaped_table, refuses_a_file_that_is_not_a_database,
              load_test_skips_without_a_map, schema_only_names_every_table,
              find_locates_the_database, find_says_so_when_nothing_matches,
              neither_db_nor_find_is_an_error, a_mass_gated_table_is_not_evidence,
              an_unfiltered_zero_carries_weight, limit_reads_a_spread_not_the_oldest)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--verbose", action="store_true", help="print the probe's output")
    args = ap.parse_args(argv)

    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        db = tmp / "planted.db"
        plant(db)
        loads = tmp / "loads.json"
        loads.write_text(json.dumps({"%dng" % int(ng): ng for ng in LOADS}))
        rc, out = run("--db", str(db), "--load-map", str(loads))
        if rc != 0:
            print("FAIL the probe exited %d on the planted database" % rc)
            print(out)
            return 1
        if args.verbose:
            print(out)

        for case in NEEDS_RUN:
            bad = case(db, out)
            failed += bool(bad)
            print("%s %s" % ("ok  " if not bad else "FAIL", case.__name__))
            for line in bad:
                print("     %s" % line)
        for case in STANDALONE:
            bad = case(db, tmp)
            failed += bool(bad)
            print("%s %s" % ("ok  " if not bad else "FAIL", case.__name__))
            for line in bad:
                print("     %s" % line)

    total = len(NEEDS_RUN) + len(STANDALONE)
    print("\n%d/%d passed" % (total - failed, total))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
