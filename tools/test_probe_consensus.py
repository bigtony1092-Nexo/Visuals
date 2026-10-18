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
import re
import sqlite3
import statistics
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


def section(out, name):
    """Just one TEST block. Splitting on the header alone runs past it into the
    next test, which made an ar1 case fail on the word SKIPPED from TEST 5."""
    if name not in out:
        return ""
    rest = out.split(name, 1)[1]
    nxt = re.search(r"\nTEST \d", rest)
    return rest[:nxt.start()] if nxt else rest


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
    # the line now carries an IQR after the quartiles, so take the two numbers
    # that belong to the quartile pair rather than splitting the whole tail
    nums = re.findall(r"[-+]\d+\.\d+", line.split("quartiles")[1])
    if len(nums) < 2:
        return ["could not read the quartile pair from %r" % line.strip()]
    lo, hi = float(nums[0]), float(nums[1])
    if abs(lo - SLOW) > 0.05:
        bad.append("lower quartile %+.3f, planted %.2f" % (lo, SLOW))
    if abs(hi - FAST) > 0.05:
        bad.append("upper quartile %+.3f, planted %.2f" % (hi, FAST))
    return bad




def plant_volumes(path, volumes=(1.0, 2.0, 5.0), beta=1.0, break_last=None,
                  loads=None):
    """Same vial at several injection volumes, with a known response."""
    rng = random.Random(7)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL, CalcMass REAL,
        PrecursorIntensity REAL)""")
    rid, seen = 0, set()
    for _ in range(120):
        while True:
            seq = "".join(rng.choice("GASPVTLINDQEMHFYW")
                          for _ in range(rng.randint(7, 15))) + "K"
            z = rng.choice([2, 3])
            if (seq, z) not in seen:
                seen.add((seq, z))
                break
        mass = WATER + sum(RESIDUE[c] for c in seq)
        calc = mass / z + pc.PROTON
        amp = 10 ** rng.uniform(3.5, 6.0)
        if loads is not None:
            # a plain ng-named cohort, for the "one load is no axis" case
            for ng in loads:
                rid += 1
                con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?)",
                            (rid, seq, z, "HeLa_%dng" % int(ng), calc, calc, mass,
                             amp * ng * math.exp(rng.gauss(0, 0.04))))
            continue
        for v in volumes:
            rid += 1
            inten = amp * (v ** beta) * math.exp(rng.gauss(0, 0.04))
            if break_last and v == volumes[-1]:
                inten /= break_last
            con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?)",
                        (rid, seq, z, "HeLa_120ng_%guL_r1" % v, calc, calc, mass, inten))
    con.commit()
    con.close()


def an_estimated_axis_has_no_absolute_scale(db, tmp):
    """The load factor is built from the same ratios, so the unit is arbitrary.

    The planted betas are 0.25 and 1.00. On an estimated axis both come out
    multiplied by one unknown constant, so their RATIO is recoverable and neither
    value is. Reading 1.0 as "linear" here is the error the text must prevent.
    """
    rc, out = run("--db", str(db))
    sec = section(out, "TEST 3")
    bad = []
    if "ESTIMATED from the ions" not in sec:
        bad.append("did not say the axis was estimated rather than declared")
    if "NO ABSOLUTE SCALE" not in sec:
        bad.append("did not warn that the slopes have no absolute scale")
    if "NOT `linear`" not in sec:
        bad.append("did not rule out reading a slope of 1.0 as linearity")
    if "absorbed by the axis" not in sec:
        bad.append("did not say a compression shared by every ion is invisible")
    line = [l for l in sec.splitlines() if "quartiles" in l]
    if not line:
        return bad + ["no quartile line"]
    nums = re.findall(r"[-+]\d+\.\d+", line[0].split("quartiles")[1])
    lo, hi = float(nums[0]), float(nums[1])
    # the ratio survives the missing scale; the absolute values do not
    if abs(hi / lo - FAST / SLOW) > 0.4:
        bad.append("quartile ratio %.2f, planted %.2f -- the one thing an estimated"
                   " axis preserves is not preserved" % (hi / lo, FAST / SLOW))
    if abs(lo - SLOW) < 0.05 and abs(hi - FAST) < 0.05:
        bad.append("the estimated axis returned the true betas, which it cannot do"
                   " -- the fixture or the estimator is not what it claims")
    return bad


def the_normalization_cut_is_relative_not_absolute(db, tmp):
    """An absolute cut on an axis with no unit lands in an arbitrary place."""
    rc, out = run("--db", str(db))
    sec = section(out, "TEST 3")
    bad = []
    if "% of the median slope" not in sec:
        bad.append("the unfit cut is not expressed against the median slope")
    if "relative on purpose" not in sec:
        bad.append("did not say why the cut is relative")
    line = [l for l in sec.splitlines() if "must stay OUT" in l or "respond below" in l]
    if not line:
        return bad + ["no normalization-set line"]
    m = re.search(r"\(([-+]\d+\.\d+)\)", line[0])
    if not m:
        bad.append("the cut's own value is not printed: %r" % line[0].strip())
    elif abs(float(m.group(1)) - 0.85) < 0.01:
        bad.append("the cut landed at 0.85, the absolute value, not a share of the"
                   " median")
    return bad








def a_run_with_no_volume_label_says_the_axis_is_unanchored(db, tmp):
    rc, out = run("--db", str(db))
    sec = section(out, "TEST 3")
    bad = []
    if "no run name carries an injection volume" not in sec:
        bad.append("ran a volume check on run names that carry no volume")
    if "anchored outside the intensities" not in sec:
        bad.append("did not say the axis is left unanchored without it")
    return bad


def a_declared_load_map_is_still_not_circular(db, tmp):
    """Real ng values must keep working, and must say they are not circular."""
    import json as _json
    lm = tmp / "loads.json"
    lm.write_text(_json.dumps({"HeLa_%dng" % int(n): n for n in LOADS}))
    rc, out = run("--db", str(db), "--load-map", str(lm))
    sec = section(out, "TEST 3")
    bad = []
    if "DECLARED ng" not in sec:
        bad.append("did not say the axis came from the load map")
    if "not circular" not in sec:
        bad.append("did not distinguish a declared axis from an estimated one")
    if "NO ABSOLUTE SCALE" in sec:
        bad.append("warned about an arbitrary scale on a declared ng axis")
    line = [l for l in sec.splitlines() if "quartiles" in l]
    if line:
        nums = re.findall(r"[-+]\d+\.\d+", line[0].split("quartiles")[1])
        lo, hi = float(nums[0]), float(nums[1])
        if abs(lo - SLOW) > 0.05 or abs(hi - FAST) > 0.05:
            bad.append("declared axis gave %.3f/%.3f, planted %.2f/%.2f"
                       % (lo, hi, SLOW, FAST))
    return bad


def the_median_per_ion_run_is_used_not_the_max(db, tmp):
    """Max-of-replicates inflated the published load ratio to 5.89x from 5.37x."""
    import inspect
    src = inspect.getsource(pc.test_load)
    bad = []
    if "statistics.median" not in src:
        bad.append("test_load does not take a median per ion and run")
    if re.search(r"best\[key\]\.get\(.*\) < inten", src):
        bad.append("still keeping the brightest replicate per load")
    return bad




def plant_shape(path, cap=None, clamp=None, subtract=None, beta=1.0,
                loads=(25.0, 50.0, 120.0, 240.0), n=300):
    """A pure power law, optionally bent by a ceiling, a noise floor or a baseline.

    cap      soft ceiling: headroom runs out, the high end flattens
    clamp    noise floor: below it the measured value sticks, the low end flattens
    subtract a baseline removed from every intensity, which costs the low end more
    """
    rng = random.Random(11)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL, CalcMass REAL,
        PrecursorIntensity REAL)""")
    rid, seen = 0, set()
    for _ in range(n):
        while True:
            seq = "".join(rng.choice("GASPVTLINDQEMHFYW")
                          for _ in range(rng.randint(7, 15))) + "K"
            z = rng.choice([2, 3])
            if (seq, z) not in seen:
                seen.add((seq, z))
                break
        mass = WATER + sum(RESIDUE[c] for c in seq)
        calc = mass / z + pc.PROTON
        amp = 10 ** rng.uniform(3.0, 6.0)
        for ng in loads:
            v = amp * (ng ** beta) * math.exp(rng.gauss(0, 0.05))
            if cap:
                v = v * cap / (cap + v) * 2
            if clamp:
                v = max(v, clamp)
            if subtract:
                v = max(v - subtract, subtract * 0.01)
            rid += 1
            con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?)",
                        (rid, seq, z, "HeLa_%dng" % int(ng), calc, calc, mass, v))
    con.commit()
    con.close()


def _four_loads(tmp):
    import json as _json
    lm = tmp / "loads4.json"
    lm.write_text(_json.dumps({"HeLa_%dng" % int(n): n for n in (25, 50, 120, 240)}))
    return lm


def the_printed_cut_is_the_printed_median(db, tmp):
    """The shadowed-variable bug: the cut was 0.85 x the BRIGHTEST quartile.

    On real data where bright ions saturate, Q4's median is low, so the cut
    dropped and the very ions it exists to exclude got through.
    """
    rc, out = run("--db", str(db))
    sec = section(out, "TEST 3")
    bad = []
    mline = [l for l in sec.splitlines() if "median slope" in l and "quartiles" in l]
    cline = [l for l in sec.splitlines() if "respond below" in l]
    if not mline or not cline:
        return ["could not find the median and cut lines"]
    med = float(re.search(r"median slope ([-+]\d+\.\d+)", mline[0]).group(1))
    cut = float(re.search(r"\(([-+]\d+\.\d+)\)", cline[0]).group(1))
    if abs(cut - 0.85 * med) > 0.002:
        bad.append("cut %+.3f is not 0.85 x the printed median %+.3f (= %+.3f)"
                   % (cut, med, 0.85 * med))
    # the brightest quartile's median must NOT be what the cut was taken from
    q4 = [l for l in sec.splitlines() if "Q4" in l and "median slope" in l]
    if q4:
        q4med = float(re.search(r"median slope ([-+]\d+\.\d+)", q4[0]).group(1))
        if abs(cut - 0.85 * q4med) < 0.002 and abs(q4med - med) > 0.01:
            bad.append("cut equals 0.85 x the Q4 median, so `med` is still shadowed")
    return bad


def the_intensity_quartiles_carry_no_verdict(db, tmp):
    """Ranking ions by intensity ranks them partly by response."""
    rc, out = run("--db", str(db))
    sec = section(out, "TEST 3")
    bad = []
    if "DESCRIPTIVE ONLY" not in sec:
        bad.append("the quartile table is not marked descriptive")
    if "NO VERDICT IS READ FROM THAT TABLE" not in sec:
        bad.append("no statement that the table carries no verdict")
    if "122" not in sec or "104" not in sec:
        bad.append("does not quote how badly the ranking is confounded")
    for banned in ("That is saturation", "OPPOSITE of saturation",
                   "read it as a floor"):
        if banned in sec:
            bad.append("a verdict is still drawn from the quartile table: %r" % banned)
    return bad


def curvature_needs_an_axis_with_a_known_shape(db, tmp):
    """An estimated axis is warped, and the warp bends every ion the same way."""
    rc, out = run("--db", str(db))
    sec = section(out, "TEST 3")
    bad = []
    if "NOT ASSESSED" not in sec:
        bad.append("ran a curvature test against an unanchored estimated axis")
    if "no certified shape" not in sec:
        bad.append("did not say the estimated axis has no certified shape")
    if "0.000/0.181/0.391/0.577" not in sec:
        bad.append("does not quote the measured departure from proportionality")
    if "arbitrary" not in sec or "proportionality" not in sec:
        bad.append("does not separate an arbitrary scale, which is harmless, from a"
                   " departure from proportionality, which is not")
    for banned in ("A DOWNWARD BEND", "AN UPWARD BEND"):
        if banned in sec:
            bad.append("a curvature verdict was printed on an unanchored axis: %r"
                       % banned)
    return bad


def a_planted_ceiling_reads_as_a_downward_bend(db, tmp):
    p = tmp / "ceiling.db"
    plant_shape(str(p), cap=5e6)
    rc, out = run("--db", str(p), "--load-map", str(_four_loads(tmp)))
    sec = section(out, "TEST 3")
    bad = []
    if "A DOWNWARD BEND" not in sec:
        bad.append("a planted ceiling did not read as a downward bend")
    if "AN UPWARD BEND" in sec:
        bad.append("a planted ceiling also read as an upward bend")
    if "baseline subtracted" not in sec:
        bad.append("did not name the second cause of a downward bend")
    if "before calling it saturation" not in sec:
        bad.append("called it saturation without naming the ambiguity")
    return bad


def a_planted_noise_floor_reads_as_an_upward_bend(db, tmp):
    p = tmp / "clamp.db"
    plant_shape(str(p), clamp=3e6)
    rc, out = run("--db", str(p), "--load-map", str(_four_loads(tmp)))
    sec = section(out, "TEST 3")
    bad = []
    if "AN UPWARD BEND" not in sec:
        bad.append("a planted noise floor did not read as an upward bend")
    if "A DOWNWARD BEND" in sec:
        bad.append("a planted noise floor also read as a downward bend")
    if "opposite shape from saturation" not in sec:
        bad.append("did not distinguish it from saturation")
    return bad


def a_plain_power_law_shows_no_curvature(db, tmp):
    """The verdicts must fire only on their own fixture."""
    p = tmp / "neither.db"
    plant_shape(str(p))
    rc, out = run("--db", str(p), "--load-map", str(_four_loads(tmp)))
    sec = section(out, "TEST 3")
    bad = []
    for banned in ("A DOWNWARD BEND", "AN UPWARD BEND"):
        if banned in sec:
            bad.append("%r fired on a pure power law with nothing planted" % banned)
    if "Neither direction dominates" not in sec:
        bad.append("did not report the absence of curvature")
    if "still invisible" not in sec:
        bad.append("did not restate that a shared compression is invisible")
    return bad


def a_baseline_subtraction_is_not_called_saturation(db, tmp):
    """It bends the same way as a ceiling, so the verdict must not name a cause."""
    p = tmp / "subtract.db"
    plant_shape(str(p), subtract=2e5)
    rc, out = run("--db", str(p), "--load-map", str(_four_loads(tmp)))
    sec = section(out, "TEST 3")
    bad = []
    if "SATURATION:" in sec:
        bad.append("named saturation as the cause of a bend a baseline also makes")
    return bad


def curvature_needs_four_loads(db, tmp):
    p = tmp / "three.db"
    plant_shape(str(p), loads=(25.0, 50.0, 120.0))
    import json as _json
    lm = tmp / "loads3.json"
    lm.write_text(_json.dumps({"HeLa_%dng" % int(n): n for n in (25, 50, 120)}))
    rc, out = run("--db", str(p), "--load-map", str(lm))
    sec = section(out, "TEST 3")
    bad = []
    if "needs 4 distinct loads" not in sec:
        bad.append("did not say curvature needs four loads per ion")
    if "Not assessed" not in sec:
        bad.append("did not mark it unassessed")
    for banned in ("A DOWNWARD BEND", "AN UPWARD BEND"):
        if banned in sec:
            bad.append("%r fired on three loads, which cannot be split in half"
                       % banned)
    return bad


def a_half_microlitre_is_not_five(db, tmp):
    """PAX8 writes 0.5 uL as `0-5uL`; matching only the decimal read it as 5."""
    bad = []
    for text, want_v, want_label in (
            ("100ng_60_0-5uL", 0.5, "100ng_60"),
            ("100ng_60_0.5uL", 0.5, "100ng_60"),
            ("500ng_60_1uL", 1.0, "500ng_60"),
            ("500ng_60_2uL", 2.0, "500ng_60"),
            ("500ng_60_5uL", 5.0, "500ng_60"),
            ("100ng_120_1uL", 1.0, "100ng_120"),
            ("HeLa_60_2uL_20260923_1", 2.0, "HeLa_60_20260923_1")):
        v, label = pc.microlitres(text)
        if v != want_v:
            bad.append("%s read as %r uL, expected %g" % (text, v, want_v))
        if label != want_label:
            bad.append("%s grouped under %r, expected %r" % (text, label, want_label))
    if pc.microlitres("no_volume_here")[0] is not None:
        bad.append("found a volume in a name that has none")
    # and the 0.5 and 5 uL runs of one series must not share a group
    a = pc.microlitres("100ng_60_0-5uL")
    b = pc.microlitres("100ng_60_5uL")
    if a[1] == b[1] and a[0] == b[0]:
        bad.append("0-5uL and 5uL are indistinguishable")
    return bad




def the_trimmed_mean_beats_the_median_on_the_known_truth(db, tmp):
    """Tony's question, answered against plant()'s known loads.

    The load factor is a trimmed mean of the central half rather than a median.
    Scored on departure from proportionality (the warp that bends every ion) and
    on how far the axis moves under contamination.
    """
    import sqlite3 as _sq
    con = _sq.connect(str(db))
    rows = {}
    for seq, z, run, i in con.execute(
            "SELECT PeptideSequence,Charge,RunName,PrecursorIntensity FROM observations"):
        rows.setdefault((seq, z), {}).setdefault(run, []).append(i)
    con.close()
    per = {k: {r: statistics.median(v) for r, v in d.items()} for k, d in rows.items()}
    runs = sorted({r for d in per.values() for r in d},
                  key=lambda r: float(r.split("_")[1].rstrip("ng")))
    truth = {r: math.log10(float(r.split("_")[1].rstrip("ng"))) for r in runs}
    anch = {k: d for k, d in per.items() if len(d) >= 3}
    ref = {k: statistics.median(d.values()) for k, d in anch.items()}

    def warp(fn):
        fac = {r: fn([math.log10(d[r] / ref[k]) for k, d in anch.items()
                      if r in d and ref[k] > 0]) for r in runs}
        xs = [truth[r] - truth[runs[0]] for r in runs]
        ys = [fac[r] - fac[runs[0]] for r in runs]
        scale = sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)
        return max(abs(y - scale * x) for x, y in zip(xs, ys))

    w_med = warp(statistics.median)
    w_trim = warp(pc.trimmed_mean)
    bad = []
    if w_trim >= w_med:
        bad.append("the trimmed mean warps %.4f against the median's %.4f, so the"
                   " estimator change is not justified" % (w_trim, w_med))
    if w_trim > 0.015:
        bad.append("the trimmed-mean axis still departs %.4f from proportionality;"
                   " the quoted figures are stale" % w_trim)
    # robustness: contaminate 5% of ions and see how far the axis moves
    import random as _rnd
    rng = _rnd.Random(3)
    keys = list(per)
    rng.shuffle(keys)
    dirty = {k: dict(d) for k, d in per.items()}
    for k in keys[:int(0.05 * len(keys))]:
        for r in dirty[k]:
            dirty[k][r] *= 10 ** rng.choice([-2.5, 2.5])
    anch2 = {k: d for k, d in dirty.items() if len(d) >= 3}
    ref2 = {k: statistics.median(d.values()) for k, d in anch2.items()}

    def axis(src, refs, fn):
        return {r: fn([math.log10(d[r] / refs[k]) for k, d in src.items()
                       if r in d and refs[k] > 0]) for r in runs}

    for name, fn in (("trimmed mean", pc.trimmed_mean),
                     ("plain mean", lambda xs: sum(xs) / len(xs))):
        a = axis(anch, ref, fn)
        b = axis(anch2, ref2, fn)
        shift = max(abs(b[r] - a[r]) for r in runs)
        if name == "trimmed mean" and shift > 0.005:
            bad.append("the trimmed mean moved %.4f under 5%% contamination, which is"
                       " not robust" % shift)
        if name == "plain mean" and shift < 0.01:
            bad.append("the plain mean moved only %.4f under contamination, so this"
                       " fixture does not test robustness" % shift)
    return bad


def a_median_axis_would_have_false_positived(db, tmp):
    """Why the estimator change matters downstream, not just on paper."""
    import sqlite3 as _sq
    con = _sq.connect(str(db))
    rows = {}
    for seq, z, run, i in con.execute(
            "SELECT PeptideSequence,Charge,RunName,PrecursorIntensity FROM observations"):
        rows.setdefault((seq, z), {}).setdefault(run, []).append(i)
    con.close()
    per = {k: {r: statistics.median(v) for r, v in d.items()} for k, d in rows.items()}
    runs = sorted({r for d in per.values() for r in d})
    anch = {k: d for k, d in per.items() if len(d) >= 3}
    ref = {k: statistics.median(d.values()) for k, d in anch.items()}

    def bend_share(fn):
        fac = {r: fn([math.log10(d[r] / ref[k]) for k, d in anch.items()
                      if r in d and ref[k] > 0]) for r in runs}
        bends = []
        for k, d in per.items():
            pts = sorted((fac[r], math.log10(d[r])) for r in d if r in fac)
            if len({f for f, _ in pts}) < 4:
                continue
            h = (len(pts) + 1) // 2
            lo = pc.slope([f for f, _ in pts[:h]], [y for _, y in pts[:h]])
            hi = pc.slope([f for f, _ in pts[-h:]], [y for _, y in pts[-h:]])
            w = pc.slope([f for f, _ in pts], [y for _, y in pts])
            if lo is None or hi is None or not w:
                continue
            bends.append((hi - lo) / abs(w))
        down = sum(1 for b in bends if b <= -0.25)
        up = sum(1 for b in bends if b >= 0.25)
        fired = down >= 0.25 * len(bends) and down > 2 * up
        return fired, statistics.median(bends)

    med_fired, med_bend = bend_share(statistics.median)
    trim_fired, trim_bend = bend_share(pc.trimmed_mean)
    bad = []
    if not med_fired:
        bad.append("a median axis no longer false-positives on plant(), so the"
                   " reason given for the change is stale")
    if trim_fired:
        bad.append("the trimmed-mean axis still reports a bend on a fixture with"
                   " nothing planted (median bend %+.3f)" % trim_bend)
    if abs(trim_bend) >= abs(med_bend):
        bad.append("trimmed-mean bend %+.3f is no smaller than the median's %+.3f"
                   % (trim_bend, med_bend))
    return bad




def plant_pax8_shape(path):
    """The real Colorado PAX8 cohort shape, from a read-only check of the 24 .raw
    folders: one sample, 22 runs at 500 ng and 2 at 100 ng.

    The 1, 2 and 5 uL injections came from vials at 0.5, 0.25 and 0.1 ug/uL, so
    all three put 500 ng on column -- the volumes were traded against
    concentration to hold load CONSTANT. Two names also contradict their own
    records: a run named 1uL recorded 0.2 uL (100 ng), and one named
    100ng_0-5uL recorded 2.0 uL of 0.25 ug/uL (500 ng).
    """
    runs = ([("JP210_500ng_60_1uL_%d" % i, 500.0) for i in range(1, 9)]
            + [("JP210_500ng_60_2uL_%d" % i, 500.0) for i in range(1, 8)]
            + [("JP210_500ng_60_5uL_%d" % i, 500.0) for i in range(1, 8)]
            + [("JP210_100ng_60_1uL_4", 100.0), ("JP210_100ng_120_1uL_4", 100.0)])
    rng = random.Random(5)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL, CalcMass REAL,
        PrecursorIntensity REAL)""")
    rid, seen = 0, set()
    for _ in range(200):
        while True:
            seq = "".join(rng.choice("GASPVTLINDQEMHFYW")
                          for _ in range(rng.randint(7, 15))) + "K"
            z = rng.choice([2, 3])
            if (seq, z) not in seen:
                seen.add((seq, z))
                break
        mass = WATER + sum(RESIDUE[c] for c in seq)
        calc = mass / z + pc.PROTON
        amp = 10 ** rng.uniform(3.5, 6.0)
        for run, ng in runs:
            rid += 1
            con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?)",
                        (rid, seq, z, run, calc, calc, mass,
                         amp * ng * math.exp(rng.gauss(0, 0.08))))
    con.commit()
    con.close()
    return runs


def volume_is_never_converted_to_load(db, tmp):
    """Load is volume x concentration, and PAX8 traded one against the other.

    1 uL of 0.5 ug/uL, 2 uL of 0.25 and 5 uL of 0.1 are all 500 ng. Read as a
    1:2:5 ladder they span 0.699 log10 units of load; the real span is 0.000.
    """
    p = tmp / "vol.db"
    plant_volumes(str(p))
    rc, out = run("--db", str(p))
    sec = section(out, "TEST 3")
    bad = []
    if "REPORTED, not" not in sec:
        bad.append("did not say the volumes are reported rather than converted")
    if "volume x concentration" not in sec:
        bad.append("did not state that load is volume times concentration")
    if "all 500 ng on column" not in sec:
        bad.append("did not give the PAX8 counter-example")
    if "No verdict either way" not in sec:
        bad.append("still draws a verdict from the volumes")
    for banned in ("agrees to", "load axis is reading real injection amount",
                   "Same-vial is NOT confirmed"):
        if banned in sec:
            bad.append("the retired proportionality verdict is still printed: %r"
                       % banned)
    # and it must not anchor the curvature test
    if "NOT ASSESSED" not in sec:
        bad.append("a volume spread anchored the curvature test, which it cannot do")
    return bad


def volumes_moving_with_intensity_are_flagged_not_credited(db, tmp):
    """On concentration-compensated runs, a 1:2:5 intensity ladder is the red flag."""
    p = tmp / "vol2.db"
    plant_volumes(str(p))
    rc, out = run("--db", str(p))
    sec = section(out, "TEST 3")
    bad = []
    if "move roughly WITH volume" not in sec:
        bad.append("did not notice the factors tracking volume")
    if "red flag" not in sec:
        bad.append("did not call a volume-tracking ladder a red flag")
    if "confirmation of the axis" in sec and "not a confirmation" not in sec:
        bad.append("read a volume-tracking ladder as confirmation")
    return bad


def a_constant_load_across_volumes_is_recognised(db, tmp):
    """The PAX8 case: factors barely move across a 5x volume range."""
    p = tmp / "vol3.db"
    # same vial concentration-compensated: every volume puts the same load on
    plant_volumes(str(p), beta=0.0)
    rc, out = run("--db", str(p))
    sec = section(out, "TEST 3")
    bad = []
    if "constant load across volumes looks like" not in sec:
        bad.append("did not recognise constant load across a volume range")
    if "NOT a failure of the axis" not in sec:
        bad.append("read constant load across volumes as an axis failure")
    return bad


def the_real_pax8_axis_is_two_points_and_says_so(db, tmp):
    p = tmp / "pax8.db"
    runs = plant_pax8_shape(str(p))
    import json as _json
    lm = tmp / "pax8loads.json"
    lm.write_text(_json.dumps({r: ng for r, ng in runs}))
    rc, out = run("--db", str(p), "--load-map", str(lm))
    sec = section(out, "TEST 3")
    bad = []
    if "Only TWO distinct loads" not in sec:
        bad.append("did not say the axis has only two points")
    if "unchecked" not in sec:
        bad.append("did not say a two-point slope is unchecked")
    if "no residual" not in sec or "no curvature" not in sec:
        bad.append("did not say what two points cannot support")
    if "concentration-compensated" not in sec:
        bad.append("did not explain why the cohort is this shape")
    if "seen at 3+ loads" in sec:
        bad.append("still claims 3+ loads on a two-load axis")
    # the planted response is exactly linear, so a declared axis must recover 1.0
    line = [l for l in sec.splitlines() if "median slope" in l and "quartiles" in l]
    if line:
        m = float(re.search(r"median slope ([-+]\d+\.\d+)", line[0]).group(1))
        if abs(m - 1.0) > 0.05:
            bad.append("median slope %+.3f on a planted linear response at two"
                       " declared loads" % m)
    # curvature needs four loads, so it must not run here
    for banned in ("A DOWNWARD BEND", "AN UPWARD BEND"):
        if banned in sec:
            bad.append("%r fired on a two-point axis" % banned)
    return bad


def one_load_is_no_axis(db, tmp):
    p = tmp / "oneload.db"
    plant_volumes(str(p), loads=(500.0,))
    import json as _json
    lm = tmp / "oneload.json"
    lm.write_text(_json.dumps({"HeLa_500ng": 500}))
    rc, out = run("--db", str(p), "--load-map", str(lm))
    sec = section(out, "TEST 3")
    bad = []
    if "SKIPPED" not in sec:
        bad.append("did not skip when every run sits at one load")
    if "no load axis to regress against" not in sec:
        bad.append("did not say why one load is not an axis")
    return bad




PAX8_LOADS = ROOT / "examples" / "pax8-colorado-loads.json"


def plant_colorado(path, vol_effect=0.0, grad_scale=None):
    """The real 24-run Colorado cohort: names and loads from the instrument record.

    vol_effect adds a genuine trap/dilution trend per decade of injection volume,
    at CONSTANT load, which is the one thing this design can measure well.

    grad_scale multiplies every intensity at one gradient length by a constant.
    A constant factor shifts log intensity by a constant, so it cannot change
    that gradient's own load slope and must not change the other gradient's at
    all. Anything that moves under it has pooled the two.
    """
    import json as _json
    lm = {k: v for k, v in _json.loads(PAX8_LOADS.read_text()).items()
          if not k.startswith("_")}
    rng = random.Random(9)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL, CalcMass REAL,
        PrecursorIntensity REAL)""")
    rid, seen = 0, set()
    for _ in range(150):
        while True:
            seq = "".join(rng.choice("GASPVTLINDQEMHFYW")
                          for _ in range(rng.randint(7, 15))) + "K"
            z = rng.choice([2, 3])
            if (seq, z) not in seen:
                seen.add((seq, z))
                break
        mass = WATER + sum(RESIDUE[c] for c in seq)
        calc = mass / z + pc.PROTON
        amp = 10 ** rng.uniform(3.5, 6.0)
        for run, ng in lm.items():
            v, _ = pc.microlitres(run)
            trend = 1.0 + vol_effect * math.log10(v / 1.0)
            gs = 1.0
            if grad_scale:
                gs = grad_scale.get(pc.gradient_of(run), 1.0)
            rid += 1
            con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?)",
                        (rid, seq, z, run, calc, calc, mass,
                         amp * ng * trend * gs * math.exp(rng.gauss(0, 0.07))))
    con.commit()
    con.close()
    return lm


def gradients_are_never_pooled(db, tmp):
    """60 min and 120 min are different LC methods; pooling mixes method with load."""
    p = tmp / "colorado.db"
    plant_colorado(str(p))
    rc, out = run("--db", str(p), "--load-map", str(PAX8_LOADS))
    sec = section(out, "TEST 3")
    bad = []
    if "different LC" not in sec:
        bad.append("did not say the gradients are different methods")
    if "60 min (20 runs)" not in sec or "120 min (4 runs)" not in sec:
        bad.append("did not break the cohort down by gradient: %r"
                   % [l for l in sec.splitlines() if "gradient length" in l])
    if "none is combined" not in sec:
        bad.append("did not state that the gradients are not combined")
    if "mix method with load" not in sec:
        bad.append("did not say what pooling would cost")
    return bad


def rescaling_one_gradient_leaves_the_other_alone(db, tmp):
    """The behaviour behind "none is combined".

    gradients_are_never_pooled above asserts the output SAYS the gradients are
    not combined. It never checked that they are not, which is why it passed
    while the slope loop, the quartiles and the normalization set all ran over
    every run at once.

    Multiplying one gradient's intensities by a constant cannot change that
    gradient's own load slope -- a constant factor is a constant shift in log
    intensity -- and cannot touch the other gradient at all. So every number
    reported for 60 min must be identical with and without the 120 min runs
    rescaled. Pooled, they move.
    """
    a, b = tmp / "gsplit_a.db", tmp / "gsplit_b.db"
    plant_colorado(str(a))
    plant_colorado(str(b), grad_scale={120: 0.2})
    sa = section(run("--db", str(a), "--load-map", str(PAX8_LOADS))[1], "TEST 3")
    sb = section(run("--db", str(b), "--load-map", str(PAX8_LOADS))[1], "TEST 3")

    def per_gradient_medians(sec):
        """median slope lines, keyed by the gradient they belong to."""
        out, cur = {}, None
        for line in sec.splitlines():
            m = re.search(r"(\d+) min\b", line)
            if m and ("min:" in line or "GRADIENT" in line.upper()):
                cur = m.group(1)
            m2 = re.search(r"median slope ([-+]\d+\.\d+)", line)
            if m2 and "Q" not in line.split("median slope")[0][-6:]:
                out.setdefault(cur, []).append(float(m2.group(1)))
        return out

    ma, mb = per_gradient_medians(sa), per_gradient_medians(sb)
    bad = []
    if None in ma or len(ma) < 2:
        bad.append("Test 3 reports no per-gradient median slope, so the two "
                   "gradients are analysed together: %r" % sorted(ma))
        return bad
    for g in sorted(set(ma) & set(mb)):
        if g == "120":
            continue
        if ma[g] != mb[g]:
            bad.append("the %s min median slope moved from %s to %s when only "
                       "the 120 min runs were rescaled, so the gradients are "
                       "pooled" % (g, ma[g], mb[g]))
    return bad


def test_three_reports_itself_as_run(db, tmp):
    """The verdict, not just the printout.

    load_analysis ended on an expression and returned None, so the whole of
    Test 3 printed normally and the summary then listed `load response` under
    NOT RUN. Every assertion in this file reads the printed text, so 62 of them
    passed over it. The constant-load check and the volume check were dropped
    on the single-gradient path as well.
    """
    p = tmp / "ran.db"
    plant_colorado(str(p))
    rc, out = run("--db", str(p), "--load-map", str(PAX8_LOADS))
    bad = []
    tail = out.split("NOT RUN:", 1)
    if len(tail) > 1 and "load response" in tail[1].split("\n")[0]:
        bad.append("Test 3 printed its whole analysis and was still reported as"
                   " NOT RUN: %r" % tail[1].split("\n")[0].strip())
    if "CONSTANT-LOAD CHECK" not in out:
        bad.append("the constant-load check was dropped")
    return bad


def a_load_level_on_one_run_is_named_unreplicated(db, tmp):
    """Each Colorado gradient holds exactly one 100 ng run."""
    p = tmp / "colorado2.db"
    plant_colorado(str(p))
    rc, out = run("--db", str(p), "--load-map", str(PAX8_LOADS))
    sec = section(out, "TEST 3")
    bad = []
    if "SINGLE run" not in sec:
        bad.append("did not flag the load level resting on one run")
    if "not measurable" not in sec:
        bad.append("did not say the scatter cannot be measured")
    # both gradients must be flagged, not just one
    if sec.count("SINGLE run") < 2:
        bad.append("only %d gradient flagged; both have a single 100 ng run"
                   % sec.count("SINGLE run"))
    return bad


def a_disputed_load_is_declared_not_trusted(db, tmp):
    """JP213_100ng_60_0-5uL: 500 ng by the record, 125 on the other reading."""
    p = tmp / "colorado3.db"
    plant_colorado(str(p))
    rc, out = run("--db", str(p), "--load-map", str(PAX8_LOADS))
    sec = section(out, "TEST 3")
    bad = []
    if "DISPUTED LOAD" not in sec:
        bad.append("the disputed run was used without saying so")
    if "500" not in sec or "125" not in sec:
        bad.append("did not give both readings")
    if "inherits the" not in sec:
        bad.append("did not say what else the dispute affects")
    if "IN this cohort" not in sec:
        bad.append("did not say whether the disputed run is present")
    return bad


def the_constant_load_series_is_its_own_test(db, tmp):
    """At a fixed 500 ng, intensity should not track injection volume."""
    flat = tmp / "flatvol.db"
    plant_colorado(str(flat), vol_effect=0.0)
    rc, out = run("--db", str(flat), "--load-map", str(PAX8_LOADS))
    sec = section(out, "TEST 3")
    bad = []
    if "CONSTANT-LOAD CHECK" not in sec:
        bad.append("no constant-load check on a cohort built around one")
    if "should not depend on" not in sec:
        bad.append("did not state the expectation")
    # a real trend must be found and must NOT be called load
    trend = tmp / "trendvol.db"
    plant_colorado(str(trend), vol_effect=0.60)
    rc2, out2 = run("--db", str(trend), "--load-map", str(PAX8_LOADS))
    sec2 = section(out2, "TEST 3")
    if "across volumes at ONE load" not in sec2:
        bad.append("a planted volume trend at constant load was not reported")
    if "not load, because load did not change" not in sec2:
        bad.append("did not rule out reading the volume trend as load")
    if "bounds how much of any load slope" not in sec2:
        bad.append("did not say what the trend is useful for")
    return bad


def the_real_names_do_carry_volumes(db, tmp):
    """The old message said they carried none, which was false for all 24."""
    p = tmp / "colorado4.db"
    plant_colorado(str(p))
    rc, out = run("--db", str(p))          # no load map: the volume path
    sec = section(out, "TEST 3")
    bad = []
    if "no run name carries an injection volume" in sec:
        bad.append("claimed the real PAX8 names carry no volume; all 24 do")
    if "carry an injection volume, but no two share a label" not in sec:
        bad.append("did not state the real reason there is nothing to compare")
    if "vial number and replicate index stay in the label" not in sec:
        bad.append("did not say why the grouping fails")
    for n in (24,):
        if str(n) not in sec:
            bad.append("did not count how many names carry a volume")
    return bad


def the_shipped_load_map_matches_the_record(db, tmp):
    """The map is data; its shape is worth asserting so a bad edit is caught."""
    import json as _json
    doc = _json.loads(PAX8_LOADS.read_text())
    loads = {k: v for k, v in doc.items() if not k.startswith("_")}
    bad = []
    if len(loads) != 24:
        bad.append("%d runs in the map, the cohort has 24" % len(loads))
    counts = {}
    for v in loads.values():
        counts[v] = counts.get(v, 0) + 1
    if counts != {500.0: 22, 100.0: 2}:
        bad.append("load distribution is %r, the record says 22 at 500 and 2 at 100"
                   % counts)
    if "_disputed" not in doc:
        bad.append("the disputed run is not flagged in the map")
    else:
        d = doc["_disputed"]
        if len(d) != 1:
            bad.append("%d disputed entries, expected 1" % len(d))
        for run, info in d.items():
            if loads.get(run) != info.get("recorded"):
                bad.append("%s is %r in the map but the dispute records %r"
                           % (run, loads.get(run), info.get("recorded")))
            if info.get("alternative") == 100.0:
                bad.append("the alternative reading is 100 ng, which the record"
                           " rules out on both readings")
    # every key must be a real-looking run stem, and every value a positive load
    for k, v in loads.items():
        if not isinstance(v, (int, float)) or v <= 0:
            bad.append("%s has load %r" % (k, v))
        if pc.gradient_of(k) not in (60, 120):
            bad.append("%s has no readable gradient" % k)
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
    if rc not in (0, 1):
        bad.append("exit %d, expected 0 or 1" % rc)
    if "Any isotope:   0 rows" not in out:
        bad.append("expected zero isotope picks in a table with none")
    if "THIS RESULT IS NOT EVIDENCE" not in out:
        bad.append("read its own zero as a finding on a table that cannot hold a pick")
    for want in ("upstream", "discarded"):
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
    if "measured rows beyond" not in out.lower():
        bad.append("did not report how many rows sit beyond the gate width")
    return bad


def a_computed_precursor_is_not_a_measurement(db, tmp):
    """prec_mz == calc_mz means the column was computed, so no mass test applies.

    Reported from the real vendor stores: 905 of 993 Thermo samples carry a
    computed precursor m/z via the MSFragger route. A probe that read those as
    "all monoisotopic, no picks" would be reporting the arithmetic back to itself.
    """
    computed = tmp / "computed.db"
    con = sqlite3.connect(computed)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL)""")
    rng = random.Random(11)
    for i in range(300):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "R"
        z = rng.choice([2, 3])
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)",
                    (i + 1, seq, z, "run_A", calc_mz, calc_mz))   # identical, by construction
    con.commit()
    con.close()
    rc, out = run("--db", str(computed))
    bad = []
    if "computed from the sequence, not measured" not in out:
        bad.append("read a computed precursor column as if it were measured")
    if "100%" not in out:
        bad.append("did not state what share of rows were computed")
    if "Nothing measured remains" not in out:
        bad.append("did not say the source is untestable when nothing measured is left")
    # A bucket row is the label followed by two counts; prose is not.
    if any(re.match(r"monoisotopic\s+\d+\s+\d+", l.strip()) for l in out.splitlines()):
        bad.append("still classified computed rows into the monoisotopic bucket")
    return bad


def a_measured_precursor_is_not_called_computed(db, tmp):
    """The planted database has real scatter; the computed warning must stay quiet."""
    rc, out = run("--db", str(db))
    if "computed from the sequence, not measured" in out:
        return ["called a measured column computed"]
    return []


def a_measured_subset_survives_a_mostly_computed_table(db, tmp):
    """The real database is 97% computed BY ROW but holds 48 measured samples.

    Pooling the two buries the testable subset: a computed row has delta exactly
    zero, so it lands in the monoisotopic bucket, holds the widest-error figure
    down, and turns a usable result into "not evidence". The measured rows must
    be reported on their own.
    """
    mixed = tmp / "mixed.db"
    con = sqlite3.connect(mixed)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL)""")
    rng = random.Random(23)
    rid = 0
    # 3000 computed rows in one run, 120 measured rows in another, 20 of which
    # sit on A1. By row the table is 96% computed; the answer is in the other 4%.
    for i in range(3000):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "R"
        z = 2
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        rid += 1
        con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)",
                    (rid, seq, z, "computed_run", calc_mz, calc_mz))
    for i in range(120):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(13)) + "K"
        z = 2
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        shift = pc.C13 if i < 20 else 0.0
        rid += 1
        con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)",
                    (rid, seq, z, "measured_run",
                     (calc_mz + shift / z) * (1 + 2.0 / 1e6), calc_mz))
    con.commit()
    con.close()
    rc, out = run("--db", str(mixed))
    bad = []
    if "1 measured throughout, 1 computed throughout, 0 mixed" not in out:
        bad.append("did not classify the runs as one measured and one computed")
    if "the 120 measured rows" not in out:
        bad.append("did not restrict the result to the measured rows")
    a1 = [l for l in out.splitlines() if l.strip().startswith("A1 only:")]
    if not a1 or "20 rows" not in a1[0]:
        bad.append("lost the 20 planted A1 rows: %r" % (a1[0].strip() if a1 else "no line"))
    if "THIS RESULT IS NOT EVIDENCE" in out:
        bad.append("dismissed a result it actually found")
    mono = [l for l in out.splitlines() if l.strip().startswith("monoisotopic")]
    if mono and "3100" in mono[0]:
        bad.append("counted the computed rows as monoisotopic")
    return bad


def _ppm_db(path, sigma, hard_cut, n=4000, seed=5):
    """Measured rows whose mass error is gaussian, optionally cut at hard_cut ppm."""
    rng = random.Random(seed)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL)""")
    rid = 0
    while rid < n:
        ppm = rng.gauss(0, sigma)
        if hard_cut is not None and abs(ppm) > hard_cut:
            continue                      # rejected before storage, as a gate would
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "K"
        z = 2
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        rid += 1
        con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)",
                    (rid, seq, z, "run_A", calc_mz * (1 + ppm / 1e6), calc_mz))
    con.commit()
    con.close()


def a_cut_distribution_is_called_cut(db, tmp):
    """Errors rejected past 10 ppm leave density at the edge; say so."""
    cutdb = tmp / "cut.db"
    _ppm_db(cutdb, sigma=8.0, hard_cut=10.0)
    rc, out = run("--db", str(cutdb))
    bad = []
    if "was CUT at about" not in out:
        bad.append("did not notice a distribution that stops dead at its cutoff")
    if "In the last 5%" not in out:
        bad.append("did not report the density at the edge")
    return bad


def a_tapered_distribution_is_not_called_cut(db, tmp):
    """Unfiltered scatter thins out on its own; do not call that a cutoff."""
    taper = tmp / "taper.db"
    _ppm_db(taper, sigma=3.0, hard_cut=None)
    rc, out = run("--db", str(taper))
    bad = []
    if "was CUT at about" in out:
        bad.append("called ordinary scatter a cutoff")
    if "looks like it ended rather than being cut" not in out:
        bad.append("did not say the distribution tapered")
    return bad


def a_taper_does_not_rescue_a_zero(db, tmp):
    """The two diagnostics must not read as contradicting each other.

    The real run came back tapered — 0.02% of rows in the top 5% of the range —
    which rules out a cutoff at 28 ppm. It does not establish that a row 358 ppm
    out would have been stored, because nothing in the sample goes near that. The
    zero stays uninformative, and what is missing is a positive control.
    """
    taper = tmp / "taper2.db"
    _ppm_db(taper, sigma=3.0, hard_cut=None, n=4000, seed=31)
    rc, out = run("--db", str(taper))
    bad = []
    if "was CUT at about" in out:
        bad.append("called ordinary scatter a cutoff")
    if "it looks like it ended rather than being cut" not in out:
        bad.append("did not report the taper")
    if "rules out a cutoff AT" not in out:
        bad.append("let the taper be read as covering cutoffs higher up")
    if "THIS RESULT IS NOT EVIDENCE" not in out:
        bad.append("credited a zero on a sample holding no error near a pick's size")
    if "positive control" not in out:
        bad.append("did not name the control that would make the zero mean something")
    return bad


def a_small_clean_sample_is_not_called_cut(db, tmp):
    """The defect in a max/p99 ratio: it flags small clean samples. This must not."""
    small = tmp / "small.db"
    _ppm_db(small, sigma=3.0, hard_cut=None, n=300, seed=9)
    rc, out = run("--db", str(small))
    if "was CUT at about" in out:
        return ["called a 300-row clean sample cut"]
    return []


def _ar1_db(path, normalized, seed=41, n_prec=200, frags=8):
    """Precursors with fragments, either absolute or normalized within spectrum."""
    rng = random.Random(seed)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, peptide_sequence TEXT, charge_state INTEGER,
        sample_name TEXT, precursor_mz REAL, calc_mz REAL,
        precursor_intensity REAL, intensity REAL, frag_type TEXT)""")
    rid = 0
    for p_i in range(n_prec):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "K"
        z = 2
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        prec = 10 ** rng.uniform(4, 7)
        # Absolute fragment intensities: each a fixed share of the precursor, so
        # ar1 is meaningful and ar1_sum is well below 1.
        raw = [prec * rng.uniform(0.01, 0.12) for _ in range(frags)]
        vals = [v / sum(raw) for v in raw] if normalized else raw
        for j, v in enumerate(vals):
            rid += 1
            con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?,?)",
                        (rid, seq, z, "run_A", calc_mz * (1 + 2.0 / 1e6), calc_mz,
                         prec, v, ("y%d" % (j + 1)) if j % 2 else ("b%d" % (j + 1))))
    con.commit()
    con.close()


def ar1_is_computed_from_absolute_intensities(db, tmp):
    """Each fragment is 1-12% of its precursor, so ar1 must land in that band."""
    adb = tmp / "ar1.db"
    _ar1_db(adb, normalized=False)
    rc, out = run("--db", str(adb))
    bad = []
    if "TEST 4" not in out:
        return ["test 4 did not appear"]
    t4 = section(out, "TEST 4")
    if "SKIPPED" in t4:
        bad.append("skipped with both intensity columns present: %r" % t4.strip()[:90])
    if "NORMALIZED" in t4:
        bad.append("called absolute intensities normalized")
    med = re.search(r"ar1 per row\s+median ([\d.]+)", t4)
    if not med:
        bad.append("no median ar1 line")
    elif not 0.01 <= float(med.group(1)) <= 0.12:
        bad.append("median ar1 %s outside the planted 0.01-0.12 band" % med.group(1))
    if "By fragment type" not in t4:
        bad.append("did not break ar1 down by fragment type")
    if "Loading invariance is NOT tested" not in t4:
        bad.append("did not say loading invariance was untested")
    return bad


def ar1_refuses_normalized_fragments(db, tmp):
    """Fragments stored as shares of their own spectrum carry no precursor scale.

    This is the defect the published PAX8 fragment tables have: every spectrum
    sums to 1.0000. Dividing a share by the precursor intensity returns the share
    scaled, not a product-to-precursor ratio.
    """
    ndb = tmp / "ar1norm.db"
    _ar1_db(ndb, normalized=True)
    rc, out = run("--db", str(ndb))
    bad = []
    t4 = section(out, "TEST 4")
    if "THE PRODUCT INTENSITIES ARE NORMALIZED" not in t4:
        bad.append("computed ar1 from fragments that sum to 1 per spectrum")
    if "100%" not in t4:
        bad.append("did not state what share of precursors summed to 1")
    if re.search(r"ar1 per row\s+median", t4):
        bad.append("reported an ar1 distribution anyway")
    return bad


def ar1_skips_without_a_product_column(db, tmp):
    """The planted database has a precursor intensity only; say so, do not guess."""
    rc, out = run("--db", str(db))
    t4 = section(out, "TEST 4")
    bad = []
    if "SKIPPED" not in t4:
        bad.append("ran ar1 with no product intensity column")
    if "product intensity" not in t4:
        bad.append("did not name the missing column")
    return bad


def _rejects_db(path, per_row, with_picks=True):
    """A store with a quarantine table, holding either rejected rows or tallies."""
    rng = random.Random(77)
    con = sqlite3.connect(path)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, peptide_sequence TEXT, charge_state INTEGER,
        sample_name TEXT, precursor_mz REAL, calc_mz REAL)""")
    if per_row:
        con.execute("""CREATE TABLE ingest_row_rejects (
            id INTEGER PRIMARY KEY, peptide_sequence TEXT, charge_state INTEGER,
            sample_name TEXT, precursor_mz REAL, calc_mz REAL, reason TEXT)""")
    else:
        con.execute("""CREATE TABLE ingest_row_rejects (
            sample_name TEXT, reason TEXT, n_rows INTEGER)""")
        con.execute("INSERT INTO ingest_row_rejects VALUES ('run_A','precursor_ppm',412)")
        con.execute("INSERT INTO ingest_row_rejects VALUES ('run_A','product_ppm',98)")

    rid = 0
    for i in range(200):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "K"
        z = 2
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        rid += 1
        # Accepted rows: inside a 10 ppm gate, so test 1 can find nothing.
        con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)",
                    (rid, seq, z, "run_A", calc_mz * (1 + rng.gauss(0, 3) / 1e6), calc_mz))
        if per_row:
            # Rejected rows: half of them picks, the rest ordinary large errors.
            shift = pc.C13 if (with_picks and i % 2 == 0) else 0.0
            ppm = 0.0 if shift else rng.uniform(40, 300) * rng.choice([1, -1])
            con.execute("INSERT INTO ingest_row_rejects VALUES (?,?,?,?,?,?,?)",
                        (rid, seq, z, "run_A",
                         (calc_mz + shift / z) * (1 + ppm / 1e6), calc_mz, "precursor_ppm"))
    con.commit()
    con.close()


def rejected_picks_are_the_positive_control(db, tmp):
    """Picks in the quarantine prove the gate is what removed them."""
    rdb = tmp / "rejects.db"
    _rejects_db(rdb, per_row=True)
    rc, out = run("--db", str(rdb))
    bad = []
    if "TEST 5" not in out:
        return ["test 5 did not appear"]
    t5 = section(out, "TEST 5")
    if "ingest_row_rejects" not in t5:
        bad.append("did not find the quarantine table")
    if "ARE ISOTOPE PICKS" not in t5:
        bad.append("did not report picks among the rejected rows")
    if "positive control" not in t5:
        bad.append("did not say what the rejected picks establish")
    if "measures the gate and not the pipeline" not in t5:
        bad.append("did not connect the rejects back to test 1's zero")
    return bad


def a_tally_only_quarantine_is_not_a_control(db, tmp):
    """Counts by reason record that rows were dropped, not what they were."""
    tdb = tmp / "tally.db"
    _rejects_db(tdb, per_row=False)
    rc, out = run("--db", str(tdb))
    t5 = section(out, "TEST 5")
    bad = []
    if "TALLIES ONLY" not in t5:
        bad.append("treated a counts-by-reason table as classifiable")
    if "ARE ISOTOPE PICKS" in t5:
        bad.append("claimed picks from a table holding no rows")
    if "keep the rejected row" not in t5:
        bad.append("did not say what the table would need to hold")
    return bad


def an_empty_quarantine_is_a_real_negative(db, tmp):
    """Rows rejected but none of them picks: that IS informative, unlike test 1's zero."""
    edb = tmp / "nopicks.db"
    _rejects_db(edb, per_row=True, with_picks=False)
    rc, out = run("--db", str(edb))
    t5 = section(out, "TEST 5")
    bad = []
    if "ARE ISOTOPE PICKS" in t5:
        bad.append("found picks where none were planted")
    if "real negative" not in t5:
        bad.append("did not credit a negative drawn from classifiable rejected rows")
    return bad


def no_quarantine_means_no_control(db, tmp):
    """The planted database has no rejects table; say the control is missing."""
    rc, out = run("--db", str(db))
    t5 = section(out, "TEST 5")
    bad = []
    if "SKIPPED" not in t5:
        bad.append("claimed to read a quarantine that does not exist")
    if "no positive control" not in t5:
        bad.append("did not say that leaves no positive control")
    return bad


def normalized_fragments_are_caught_ACROSS_RUNS(db, tmp):
    """The guard must key on a spectrum, not an ion.

    This is the case that slipped past: fragments normalized per spectrum, read
    across three runs, sum to 3.0 per (sequence, charge) rather than 1.0, so a
    guard keyed on the ion never fires. The original test used one run and
    passed while the guard was broken for every realistic input.
    """
    multi = tmp / "ar1norm_multirun.db"
    rng = random.Random(77)
    con = sqlite3.connect(multi)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, peptide_sequence TEXT, charge_state INTEGER,
        sample_name TEXT, precursor_mz REAL, calc_mz REAL,
        precursor_intensity REAL, intensity REAL, frag_type TEXT)""")
    rid = 0
    for p_i in range(80):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "K"
        z = 2
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        for sample in ("run_A", "run_B", "run_C"):
            prec = 10 ** rng.uniform(4, 7)
            raw = [rng.uniform(0.01, 0.12) for _ in range(8)]
            tot = sum(raw)
            for j, v in enumerate(raw):
                rid += 1
                con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?,?,?,?)",
                            (rid, seq, z, sample, calc_mz * (1 + 2.0 / 1e6), calc_mz,
                             prec, v / tot, "y%d" % (j + 1)))
    con.commit()
    con.close()
    rc, out = run("--db", str(multi))
    t4 = section(out, "TEST 4")
    bad = []
    if "THE PRODUCT INTENSITIES ARE NORMALIZED" not in t4:
        bad.append("missed per-spectrum normalization spanning 3 runs")
    if re.search(r"ar1 per row\s+median", t4):
        bad.append("reported an ar1 distribution from normalized fragments")
    if "spectra (run, sequence, charge)" not in t4:
        bad.append("did not say it keyed on the spectrum")
    return bad


def spread_sampling_does_not_duplicate_rows(db, tmp):
    """Each block must be bounded above, or slices overrun into each other."""
    rc, out = run("--db", str(db), "--limit", "5000", "--order", "spread")
    bad = []
    m = re.search(r"(\d+) rows carried both a reported and a calculated", out)
    if not m:
        return ["no row count reported"]
    read = int(m.group(1))
    total = 1893          # the planted table's size
    if read > total:
        bad.append("read %d rows from a %d-row table: slices overlapped" % (read, total))
    a1 = [l for l in out.splitlines() if l.strip().startswith("A1 only:")]
    if a1 and "33 ions" not in a1[0]:
        bad.append("ion count distorted by duplication: %r" % a1[0].strip())
    return bad


def a_cut_far_beyond_the_histogram_is_still_called_cut(db, tmp):
    """Picks sit at ~1000 ppm, past the 400 ppm histogram limit.

    Taking 95% of an uncapped maximum against capped bins puts the edge window
    past every bin, finds zero rows, and calls a hard-cut distribution "ended" —
    in exactly the case the tool was built for.
    """
    wide = tmp / "wide_ppm.db"
    rng = random.Random(4)
    con = sqlite3.connect(wide)
    con.execute("""CREATE TABLE observations (
        id INTEGER PRIMARY KEY, PeptideSequence TEXT, Charge INTEGER,
        RunName TEXT, PrecMz REAL, CalcMz REAL)""")
    for i in range(600):
        seq = "".join(rng.choice("GASPVTLINDQEMHFYW") for _ in range(12)) + "K"
        z = 2
        calc_mz = (WATER + sum(RESIDUE[c] for c in seq)) / z + pc.PROTON
        # A sixth of the rows sit at a pick's distance, the rest near zero:
        # appreciable density right at the edge, which is the definition of cut.
        shift = pc.C13 if i % 6 == 0 else 0.0
        con.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)",
                    (i + 1, seq, z, "run_A",
                     (calc_mz + shift / z) * (1 + rng.gauss(0, 1) / 1e6), calc_mz))
    con.commit()
    con.close()
    rc, out = run("--db", str(wide))
    t1 = section(out, "TEST 1")
    bad = []
    if "looks like it ended rather than being cut" in t1:
        bad.append("called a distribution with 17% of rows at the maximum 'ended'")
    if "+" not in t1.split("Mass error")[1][:160]:
        bad.append("did not mark the capped percentiles as lower bounds")
    if "histogram limit" not in t1:
        bad.append("did not say the percentiles hit the histogram limit")
    return bad


def a_skipped_test_is_not_a_pass(db, tmp):
    """Exit status must reflect every test, not just the load response."""
    rc, out = run("--db", str(db))
    bad = []
    if rc != 1:
        bad.append("exit %d with ar1 and rejects both skipped" % rc)
    if "NOT RUN:" not in out:
        bad.append("did not name the tests that did not run")
    if "skipped test is not a passed one" not in out:
        bad.append("did not say why the exit status is 1")
    return bad


def limit_reads_a_spread_not_the_oldest(db, tmp):
    """--limit must not silently mean 'the oldest rows'."""
    rc, out = run("--db", str(db), "--limit", "200")
    bad = []
    if rc not in (0, 1):
        bad.append("exit %d, expected 0 or 1" % rc)
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
              an_estimated_axis_has_no_absolute_scale,
              the_printed_cut_is_the_printed_median,
              the_intensity_quartiles_carry_no_verdict,
              curvature_needs_an_axis_with_a_known_shape,
              a_planted_ceiling_reads_as_a_downward_bend,
              a_planted_noise_floor_reads_as_an_upward_bend,
              a_plain_power_law_shows_no_curvature,
              a_baseline_subtraction_is_not_called_saturation,
              curvature_needs_four_loads,
              a_half_microlitre_is_not_five,
              the_trimmed_mean_beats_the_median_on_the_known_truth,
              a_median_axis_would_have_false_positived,
              the_normalization_cut_is_relative_not_absolute,
              volume_is_never_converted_to_load,
              volumes_moving_with_intensity_are_flagged_not_credited,
              a_constant_load_across_volumes_is_recognised,
              the_real_pax8_axis_is_two_points_and_says_so,
              one_load_is_no_axis,
              gradients_are_never_pooled,
         rescaling_one_gradient_leaves_the_other_alone,
         test_three_reports_itself_as_run,
              a_load_level_on_one_run_is_named_unreplicated,
              a_disputed_load_is_declared_not_trusted,
              the_constant_load_series_is_its_own_test,
              the_real_names_do_carry_volumes,
              the_shipped_load_map_matches_the_record,
              a_run_with_no_volume_label_says_the_axis_is_unanchored,
              a_declared_load_map_is_still_not_circular,
              the_median_per_ion_run_is_used_not_the_max,
              schema_only_names_every_table,
              find_locates_the_database, find_says_so_when_nothing_matches,
              neither_db_nor_find_is_an_error, a_mass_gated_table_is_not_evidence,
              an_unfiltered_zero_carries_weight, limit_reads_a_spread_not_the_oldest,
              a_computed_precursor_is_not_a_measurement,
              a_measured_precursor_is_not_called_computed,
              a_measured_subset_survives_a_mostly_computed_table,
              a_cut_distribution_is_called_cut,
              a_tapered_distribution_is_not_called_cut,
              a_taper_does_not_rescue_a_zero,
              a_small_clean_sample_is_not_called_cut,
              ar1_is_computed_from_absolute_intensities,
              ar1_refuses_normalized_fragments,
              ar1_skips_without_a_product_column,
              rejected_picks_are_the_positive_control,
              a_tally_only_quarantine_is_not_a_control,
              an_empty_quarantine_is_a_real_negative,
              no_quarantine_means_no_control,
              normalized_fragments_are_caught_ACROSS_RUNS,
              spread_sampling_does_not_duplicate_rows,
              a_cut_far_beyond_the_histogram_is_still_called_cut,
              a_skipped_test_is_not_a_pass)


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
        # 1 is correct here: the planted table has no product-intensity column
        # and no rejects table, so ar1 and rejects skip, and a skipped test is
        # not a passed one. Only a crash or a usage error (2) is a setup failure.
        if rc not in (0, 1) or "TEST 1" not in out:
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
