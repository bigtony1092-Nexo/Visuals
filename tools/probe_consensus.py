#!/usr/bin/env python3
"""Ask ARiA_Consensus.db the three questions the published data cannot answer.

The PAX8 page carries one already-corrected m/z per PSM. Testing whether the
pipeline ever names an isotope peak instead of the monoisotopic one needs the
precursor m/z *as the instrument reported it*, before correction, and that
survives only in the consensus database.

    python3 tools/probe_consensus.py --find --schema-only
    python3 tools/probe_consensus.py --db /path/to/ARiA_Consensus.db
    python3 tools/probe_consensus.py --db ... --load-map loads.json

--find searches the usual places rather than being told the path. It opens
nothing but the 16-byte SQLite header, reports every candidate it saw, and uses
the largest.

Column names are discovered, not assumed: any table holding a peptide sequence,
a charge, a reported m/z and a calculated m/z or mass is a candidate, and the
one with the most rows wins. --table/--raw-mz/--calc-mz/--seq/--charge override
that when the guess is wrong, and --schema-only prints what it saw so you can
pass the right ones.

Three tests. Each is skipped with a stated reason rather than faked if the
columns it needs are absent:

  1. ISOTOPE PICK   (reported - calculated) x z, clustered. A cluster at
     +1.0034 is the pipeline putting a name on the A1 peak. The published ion
     keys yielded 1 such pick; the study reported 33. Run second, on deltas
     corrected by test 2, because an uncorrected offset drags deamidation
     (+0.9840) up toward A1 (+1.0034) and the two are 0.0193 Da apart.
     It also refuses two vacuous cases: a reported m/z identical to the
     calculated one (computed, not measured), and a sample whose widest error
     is far short of what a pick would be (filtered before storage).
  2. CALIBRATION    the ppm offset of the monoisotopic cluster, per run. This
     is the correction applied before publication, measured rather than
     assumed. Reported first because test 1 depends on it.
  3. LOAD RESPONSE  log-log slope of precursor intensity against ng loaded, per
     ion. Should sit near 1. Recomputed from the published data it came out at
     +0.239, which turned out to be the load axis rather than the detector.
     Needs --load-map.
  4. AR1            product intensity over precursor intensity, on the same row.
     Refuses first if the fragments were normalized within their own spectrum,
     which destroys the scale ar1 depends on.

Exit status is 0 if every requested test ran, 1 if any was skipped for want of
a column, 2 if the database could not be opened.
"""

import argparse
import json
import math
import os
import re
import sqlite3
import statistics
import sys
from collections import Counter, defaultdict

PROTON = 1.007276
C13 = 1.003355          # monoisotopic spacing of the isotope envelope
TOL = 0.008             # Da, on the charge-multiplied delta

# Deamidation is +0.984016 Da and an A1 isotope is +1.003355: 0.019339 Da apart.
# A tolerance of 0.02 counts every deamidated peptide as an isotope pick, and the
# PAX8 run alone had 20 deamidations. TOL must therefore stay below half that gap
# (0.00967) or the two windows overlap and a delta sitting between them is
# assigned by table order rather than by distance.
#
# 0.008 Da is roomy for mass accuracy on a corrected m/z: the calibrated
# residuals here scatter by about 0.3 ppm, some 0.0003 Da. It is tight for an
# UNcorrected one, where a 4 ppm offset on a 2000 Da peptide is already 0.008 —
# which is why the calibration is measured first. Past that point rows fall into
# "unexplained" rather than being mislabelled, so the failure is loud.
#
# A delta is tested against known modifications before isotopes: a calculated
# m/z that omits a peptide's modification produces an offset that says nothing
# about isotope selection.
MODS = (
    ("deamidation", 0.984016),
    ("oxidation", 15.994915),
    ("carbamidomethyl", 57.021464),
    ("acetyl", 42.010565),
    ("phospho", 79.966331),
    ("ubiquitin remnant (GG)", 114.042927),
    ("pyroglutamate", -17.026549),
    ("water loss", -18.010565),
    ("Met loss", -131.040485),
)
DEAMIDATION = 0.984016
# An A1 call is a coin flip only near the midpoint between deamidation and A1.
# A band wider than their 0.0193 Da separation would flag every A1 row, which
# says nothing; this one flags the rows that genuinely could go either way.
BOUNDARY = (DEAMIDATION + C13) / 2
BOUNDARY_MARGIN = 0.004

# An isotope pick is off by 1.003/z Da, which is hundreds of ppm: 1000 ppm at
# m/z 500 and z=2, 250 ppm at m/z 1000 and z=4. So a table whose widest observed
# error is a few tens of ppm cannot contain one, and reporting zero picks from
# such a table says nothing about the pipeline. It says the rows were filtered
# on mass before they were stored. GATE_PPM is the width below which the probe
# refuses to read its own zero as a result.
GATE_PPM = 50.0

# Some pipelines write a COMPUTED precursor m/z into the reported column — the
# MSFragger route does this for most Thermo samples in the vendor stores, where
# prec_mz equals calc_mz exactly. Such a row cannot show any mass error, so a
# mass test over them is not a lenient test, it is no test. Real measured values
# essentially never land exactly on the calculated one.
IDENTICAL_TOL = 1e-9

# A fixed width threshold cannot tell a filtered distribution from a clean one.
# The giveaway for a cutoff is density right at the edge: a distribution that was
# cut still has appreciable weight in its last few ppm and then nothing, while one
# that ended naturally has thinned to a handful of rows by the time it reaches its
# own maximum. A ratio of maximum to 99th percentile will not do, because for
# ordinary scatter that ratio depends on how many rows were read — about 1.4 at a
# few hundred rows and 1.9 at a hundred thousand — so it would call a small clean
# sample cut. Density at the edge does not have that defect.
EDGE_FRACTION = 0.05    # the top 5% of the observed range is "the edge"
EDGE_SHARE = 0.005      # density there above this share of all rows means cut
EDGE_MIN = 5            # and at least this many rows, so small samples stay quiet
PPM_BIN = 0.1          # histogram resolution, ppm
PPM_BINS = 4000        # up to 400 ppm, then an overflow bin

# Candidate column names, normalized (lowercased, non-alphanumerics dropped) and
# in priority order: the first one present wins.
HINTS = {
    "product_intensity": ("intensity", "productintensity", "fragmentintensity",
                          "prodintensity", "fragintensity"),
    "frag_type": ("fragtype", "fragmenttype", "iontype", "fragion", "annotation"),
    "raw_mz": ("precmz", "precursormz", "reportedmz", "observedmz", "obsmz",
               "rawmz", "mzreported", "mzobserved", "measuredmz", "mzprec"),
    "calc_mz": ("calcmz", "calculatedmz", "theoreticalmz", "theormz",
                "expectedmz", "mzcalc", "mztheoretical"),
    "calc_mass": ("calcmass", "calculatedmass", "theoreticalmass", "theormass",
                  "monoisotopicmass", "neutralmass", "masscalc", "peptidemass"),
    "seq": ("peptidesequence", "sequence", "peptide", "strippedsequence",
            "modifiedsequence", "pepseq"),
    "charge": ("charge", "precursorcharge", "chargestate", "z"),
    "run": ("run", "runname", "file", "filename", "rawfile", "sourcefile",
            "sample", "samplename", "runid"),
    "intensity": ("precursorintensity", "matchedintensity", "apexintensity",
                  "integratedarea", "intensity", "area", "sumintensity"),
}


# Where a consensus database plausibly sits, and what it is plausibly called.
# --find walks these so the person running the probe does not have to know the
# path. Directories that are large and never hold one are skipped outright.
FIND_ROOTS = ("~", "~/Documents", "~/Desktop", "~/Downloads", "~/data", "~/ARiA",
              "/data", "/opt", "/srv", "/mnt", "/media", "/Volumes")
FIND_NAMES = ("consensus", "aria")
FIND_SKIP = {".git", "node_modules", "__pycache__", "Library", "AppData",
             ".cache", ".Trash", "venv", ".venv", "site-packages"}
FIND_MAX_DEPTH = 5


def looks_like_sqlite(path):
    """Cheap check: SQLite files start with a fixed 16-byte header."""
    try:
        with open(path, "rb") as fh:
            return fh.read(16) == b"SQLite format 3\x00"
    except OSError:
        return False


def find_databases():
    """Candidate consensus databases, largest first.

    Matching is on the filename, so a database named something else entirely is
    missed; that is what --db is for. Nothing is opened except to read the
    16-byte header, and nothing outside the listed roots is touched.
    """
    found, seen_dirs = {}, set()
    for root in FIND_ROOTS:
        base = os.path.expanduser(root)
        if not os.path.isdir(base):
            continue
        base_depth = base.rstrip(os.sep).count(os.sep)
        for dirpath, dirnames, filenames in os.walk(base, followlinks=False):
            real = os.path.realpath(dirpath)
            if real in seen_dirs:
                dirnames[:] = []
                continue
            seen_dirs.add(real)
            if dirpath.count(os.sep) - base_depth >= FIND_MAX_DEPTH:
                dirnames[:] = []
            dirnames[:] = [d for d in dirnames if d not in FIND_SKIP
                           and not d.startswith(".")]
            for name in filenames:
                low = name.lower()
                if not low.endswith((".db", ".sqlite", ".sqlite3")):
                    continue
                if not any(h in low for h in FIND_NAMES):
                    continue
                full = os.path.join(dirpath, name)
                if full in found or not looks_like_sqlite(full):
                    continue
                try:
                    found[full] = os.path.getsize(full)
                except OSError:
                    pass
    return sorted(found.items(), key=lambda kv: kv[1], reverse=True)


def norm(name):
    return re.sub(r"[^a-z0-9]", "", name.lower())


def pick(cols, role):
    by = {norm(c): c for c in cols}
    for hint in HINTS[role]:
        if hint in by:
            return by[hint]
    return None


def tables(con):
    rows = con.execute(
        "SELECT name FROM sqlite_master WHERE type IN ('table','view') "
        "AND name NOT LIKE 'sqlite_%' ORDER BY name")
    return [r[0] for r in rows]


def table_columns(con, table):
    return [r[1] for r in con.execute('PRAGMA table_info("%s")' % table)]


def row_count(con, table):
    try:
        return con.execute('SELECT COUNT(*) FROM "%s"' % table).fetchone()[0]
    except sqlite3.Error:
        return 0


def survey(con):
    """What every table offers, so a wrong guess is visible rather than silent."""
    out = []
    for t in tables(con):
        cols = table_columns(con, t)
        roles = {role: pick(cols, role) for role in HINTS}
        out.append({"table": t, "rows": row_count(con, t), "columns": cols, "roles": roles})
    return out


def score(entry):
    """Rank a table by how much of the test it can support."""
    r = entry["roles"]
    if not (r["seq"] and r["charge"] and r["raw_mz"]):
        return (0, 0)
    have_truth = bool(r["calc_mz"] or r["calc_mass"])
    return (2 if have_truth else 1, entry["rows"])


def classify(delta_neutral):
    """Name the offset between the reported precursor and the calculated one.

    Returns (label, isotope_steps). isotope_steps is the number of C13 spacings
    to subtract before the residual can be read as a mass calibration, or None
    when the offset is not understood well enough to calibrate from.

    Modifications are tested before isotopes, because a calculated m/z that does
    not carry a peptide's modification produces a delta that has nothing to do
    with isotope selection, and because deamidation sits only 0.0193 Da from A1.
    """
    if abs(delta_neutral) <= TOL:
        return "monoisotopic", 0
    for name, mass in MODS:
        for k in (0, 1, 2):
            if abs(delta_neutral - mass - k * C13) <= TOL:
                if k:
                    return "%s + A%d" % (name, k), None
                return name, None
    for k in (1, 2, 3, 4):
        if abs(delta_neutral - k * C13) <= TOL:
            return "A%d" % k, k
        if abs(delta_neutral + k * C13) <= TOL:
            return "below mono (-A%d)" % k, -k
    return "unexplained", None


def rowid_span(con, table):
    """(min, max) rowid, or None for a WITHOUT ROWID table."""
    try:
        lo, hi = con.execute('SELECT MIN(rowid), MAX(rowid) FROM "%s"' % table).fetchone()
    except sqlite3.Error:
        return None
    if lo is None or hi is None or hi <= lo:
        return None
    return lo, hi


def fetch(con, cfg, limit, order="oldest", blocks=64):
    r = cfg["roles"]
    want = [r["seq"], r["charge"], r["raw_mz"]]
    truth = r["calc_mz"] or r["calc_mass"]
    truth_is_mass = not r["calc_mz"]
    want.append(truth)
    extra = [c for c in (r["run"], r["intensity"]) if c]
    cols = ", ".join('"%s"' % c for c in want + extra)
    where = 'WHERE "%s" IS NOT NULL AND "%s" IS NOT NULL' % (r["raw_mz"], truth)
    base = 'SELECT %s FROM "%s" ' % (cols, cfg["table"])

    queries = []
    span = rowid_span(con, cfg["table"]) if (limit and order == "spread") else None
    if span and limit >= blocks:
        # Read equal slices from evenly spaced points across the rowid range, so a
        # limited read is not a read of the oldest rows only.
        lo, hi = span
        per = max(1, limit // blocks)
        step = (hi - lo) / blocks
        for b in range(blocks):
            start = int(lo + b * step)
            queries.append((base + where + ' AND rowid >= %d LIMIT %d' % (start, per), ()))
    else:
        sql = base + where
        if limit:
            if order == "newest":
                sql += " ORDER BY rowid DESC"
            sql += " LIMIT %d" % limit
        queries.append((sql, ()))

    for sql, params in queries:
        for row in con.execute(sql, params):
            seq, z, raw, calc = row[0], row[1], row[2], row[3]
            rest = dict(zip(extra, row[4:]))
            try:
                z = int(z)
                raw = float(raw)
                calc = float(calc)
            except (TypeError, ValueError):
                continue
            if z < 1 or raw <= 0 or calc <= 0:
                continue
            calc_mz = calc / z + PROTON if truth_is_mass else calc
            yield seq, z, raw, calc_mz, rest


def measure_calibration(records, cfg, out):
    """Pass 1: the ppm offset of the monoisotopic cluster, per run.

    Only rows already within TOL of the calculated m/z are used, so this needs no
    prior knowledge of which rows are isotope picks. The offset matters twice: it
    is the correction that was applied before publication, and leaving it in place
    pushes every delta upward, which drives deamidation (+0.9840) toward A1
    (+1.0034) and makes the two indistinguishable.
    """
    run_col = cfg["roles"]["run"]
    ppm = defaultdict(list)
    for seq, z, raw, calc_mz, rest in records:
        if abs(raw - calc_mz) < IDENTICAL_TOL:
            continue        # computed, not measured: it would fix the offset at zero
        if abs((raw - calc_mz) * z) <= TOL:
            run = str(rest.get(run_col)) if run_col else "(all runs)"
            ppm[run].append((raw - calc_mz) / calc_mz * 1e6)

    out.append("TEST 2 - MASS CALIBRATION (measured first; test 1 uses it)")
    if not ppm:
        out.append("  No measured row sat on its calculated m/z, so there is no cluster to")
        out.append("  calibrate from. Rows whose reported m/z equals the calculated one are")
        out.append("  excluded here: they are computed, and would peg the offset at zero.")
        out.append("  Test 1 runs on uncorrected deltas below.")
        return {}
    if not run_col:
        out.append("  No run column found, so this is pooled across every run in the table.")
    out.append("  %-34s %7s %10s %10s" % ("run", "n", "median ppm", "spread"))
    offsets = {}
    for run in sorted(ppm):
        vals = ppm[run]
        offsets[run] = statistics.median(vals)
        spread = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        out.append("  %-34s %7d %10.2f %10.2f" % (run[-34:], len(vals), offsets[run], spread))
    allv = [v for vals in ppm.values() for v in vals]
    out.append("  Pooled median %.2f ppm over %d measurements." % (
        statistics.median(allv), len(allv)))
    out.append("  This is the correction that was applied before publication.")
    out.append("  Rows whose reported m/z equals the calculated one are excluded; they")
    out.append("  would peg every offset at zero.")
    return offsets


def test_isotope(records, cfg, offsets, out):
    """Pass 2: name every delta, on calibration-corrected m/z."""
    buckets = Counter()
    distinct = defaultdict(set)
    examples = {}
    unexplained = Counter()
    unexplained_eg = {}
    run_col = cfg["roles"]["run"]
    total = 0
    measured = 0
    ambiguous = 0
    widest_ppm = 0.0
    beyond_gate = 0
    identical = 0
    by_run = defaultdict(lambda: [0, 0])   # run -> [measured, computed]
    hist = [0] * (PPM_BINS + 1)            # |ppm| distribution of measured rows

    for seq, z, raw, calc_mz, rest in records:
        total += 1
        run = str(rest.get(run_col)) if run_col else "(all runs)"
        ppm = offsets.get(run, 0.0)
        corrected = raw / (1 + ppm / 1e6)
        delta = (corrected - calc_mz) * z
        # How far from its calculated m/z this row sits, in ppm. An isotope pick
        # is hundreds of ppm out, so the widest value seen bounds what this table
        # is capable of showing.
        if abs(raw - calc_mz) < IDENTICAL_TOL:
            # Computed from the sequence, not measured. It can show no error of any
            # size, so counting it as "monoisotopic" would be counting arithmetic.
            identical += 1
            by_run[run][1] += 1
            continue
        measured += 1
        by_run[run][0] += 1
        err_ppm = abs(corrected - calc_mz) / calc_mz * 1e6
        widest_ppm = max(widest_ppm, err_ppm)
        hist[min(int(err_ppm / PPM_BIN), PPM_BINS)] += 1
        if err_ppm > GATE_PPM:
            beyond_gate += 1
        label, k = classify(delta)
        buckets[label] += 1
        if label in ("A1", "deamidation") and abs(delta - BOUNDARY) <= BOUNDARY_MARGIN:
            ambiguous += 1
        key = (seq, z)
        distinct[label].add(key)
        examples.setdefault(label, (seq, z, corrected, calc_mz, delta))
        if label == "unexplained":
            step = round(delta, 2)
            unexplained[step] += 1
            unexplained_eg.setdefault(step, (seq, z, corrected, calc_mz))

    out.append("")
    out.append("TEST 1 - ISOTOPE PICK")
    out.append("  %d rows carried both a reported and a calculated precursor." % total)
    if identical:
        meas_runs = sum(1 for m, c in by_run.values() if m and not c)
        comp_runs = sum(1 for m, c in by_run.values() if c and not m)
        mixed_runs = sum(1 for m, c in by_run.values() if m and c)
        out.append("  %d of them (%.0f%%) have a reported m/z identical to the calculated"
                   % (identical, 100.0 * identical / total if total else 0))
        out.append("  one — computed from the sequence, not measured. A computed row can show")
        out.append("  no error of any size, so it is EXCLUDED below rather than being counted")
        out.append("  as a monoisotopic hit. What follows is the %d measured rows." % measured)
        out.append("  By run: %d measured throughout, %d computed throughout, %d mixed."
                   % (meas_runs, comp_runs, mixed_runs))
        if not measured:
            out.append("")
            out.append("  Nothing measured remains. This source cannot be tested for mass error.")
            return {"a1_ions": 0, "any_isotope_ions": 0, "rows": total, "measured": 0}
    if offsets:
        out.append("  Deltas below are after removing each run's calibration offset.")
    else:
        out.append("  No calibration was available, so deltas below are uncorrected.")
    if not total:
        out.append("  Nothing to cluster. The two m/z columns never co-occur on a row.")
        return {"a1_ions": 0, "any_isotope_ions": 0, "rows": 0}
    out.append("")
    out.append("  %-24s %8s %8s   %s" % ("offset", "rows", "ions", "example"))
    order = ["monoisotopic", "A1", "A2", "A3", "A4"]
    order += sorted(b for b in buckets if b not in order and b != "unexplained")
    order += ["unexplained"]
    for label in order:
        if label not in buckets:
            continue
        seq, z, raw, calc_mz, delta = examples[label]
        out.append("  %-24s %8d %8d   %s %d+ %.4f vs %.4f (%+.4f)" % (
            label, buckets[label], len(distinct[label]),
            seq[:24] + ("..." if len(seq) > 24 else ""), z, raw, calc_mz, delta))

    def tally(pred):
        rows = sum(n for b, n in buckets.items() if pred(b))
        ions = set()
        for label, keys in distinct.items():
            if pred(label):
                ions |= keys
        return rows, len(ions)

    a1_rows, a1_ions = tally(lambda b: b == "A1")
    any_rows, any_ions = tally(lambda b: b.startswith("A"))
    out.append("")
    out.append("  A1 only:       %d rows, %d ions." % (a1_rows, a1_ions))
    out.append("  Any isotope:   %d rows, %d ions." % (any_rows, any_ions))
    out.append("  The study's 33 is one of these two; which depends on whether it counted")
    out.append("  A1 alone. Compare ions, not rows: one ion appears once per run.")
    def pct(q):
        """The |ppm| value at quantile q, from the histogram."""
        want = q * measured
        seen = 0
        for i, n in enumerate(hist):
            seen += n
            if seen >= want:
                return i * PPM_BIN
        return widest_ppm

    out.append("")
    out.append("  Mass error of the measured rows, |ppm|:")
    out.append("    median %.1f   p95 %.1f   p99 %.1f   max %.1f"
               % (pct(0.5), pct(0.95), pct(0.99), widest_ppm))
    out.append("  Measured rows beyond %.0f ppm: %d of %d."
               % (GATE_PPM, beyond_gate, measured))
    edge_from = widest_ppm * (1 - EDGE_FRACTION)
    edge_rows = sum(n for i, n in enumerate(hist) if i * PPM_BIN >= edge_from)
    cut = (measured and edge_rows >= EDGE_MIN
           and edge_rows / measured >= EDGE_SHARE)
    if measured:
        out.append("  In the last 5%% of that range (%.1f to %.1f ppm): %d rows, %.2f%% of all."
                   % (edge_from, widest_ppm, edge_rows, 100.0 * edge_rows / measured))
    if cut:
        out.append("  That is appreciable density right at the edge, then nothing. The")
        out.append("  distribution was CUT at about %.0f ppm rather than ending there, so"
                   % widest_ppm)
        out.append("  what lies past the cut was never stored and cannot be counted here.")
    elif measured:
        out.append("  The distribution has thinned to almost nothing by its own maximum, so")
        out.append("  it looks like it ended rather than being cut.")
        out.append("  Read that narrowly: it rules out a cutoff AT %.0f ppm. It says nothing"
                   % widest_ppm)
        out.append("  about a cutoff higher up, because no row in this sample goes near one.")
    if any_ions == 0 and beyond_gate == 0:
        out.append("")
        out.append("  THIS RESULT IS NOT EVIDENCE. An isotope pick is 1.003/z Da from the")
        out.append("  calculated m/z, which is hundreds of ppm: about 1000 ppm at m/z 500")
        out.append("  and z=2. No row in this sample is even %.0f ppm out, so a pick could" % GATE_PPM)
        out.append("  not have been stored here and could not have been found. Zero is a")
        out.append("  property of this table, not a measurement of the pipeline.")
        out.append("")
        out.append("  Something bounded these errors before they were written — an ingest")
        out.append("  tolerance, or a search engine that only emits matches inside its own")
        out.append("  tolerance. Ask upstream of it: the search output (DIA-NN report,")
        out.append("  FragPipe psm.tsv, vendor spectrum files) or whatever it discarded.")
        out.append("")
        out.append("")
        out.append("  WHAT WOULD MAKE A ZERO HERE MEAN SOMETHING: a positive control. One row,")
        out.append("  anywhere in this source, carrying a mass error of a pick's magnitude —")
        out.append("  a few hundred ppm — that the pipeline nevertheless stored. An")
        out.append("  entrapment sequence, a deliberate misassignment, a quarantine record.")
        out.append("  Without one, 'no picks' and 'picks cannot get in' produce identical")
        out.append("  output, and a smooth tail at %.0f ppm does not separate them."
                   % widest_ppm)
        out.append("")
        out.append("  AND CHECK ONE THING THERE FIRST. This test assumes the calculated m/z")
        out.append("  was derived from the sequence alone. If it was instead chosen to")
        out.append("  explain the observation — the modification or isotope that leaves the")
        out.append("  smallest residual — then the difference between reported and")
        out.append("  calculated is small by construction, for a pick as much as for a")
        out.append("  monoisotope, and no mass test over those two columns can ever find")
        out.append("  anything. The pick is recorded in the assignment, not the residual.")
    elif any_ions == 0:
        out.append("")
        out.append("  Zero picks, and %d rows do sit beyond %.0f ppm, so a pick was not"
                   % (beyond_gate, GATE_PPM))
        out.append("  excluded by construction. This zero carries weight.")
    if ambiguous:
        out.append("")
        out.append("  CAUTION: %d rows fell within %.3f Da of %.4f, the midpoint between"
                   % (ambiguous, BOUNDARY_MARGIN, BOUNDARY))
        out.append("  deamidation (%+.4f) and A1 (%+.4f), which are only %.4f Da apart."
                   % (DEAMIDATION, C13, C13 - DEAMIDATION))
        out.append("  Those calls went to whichever was nearer and should not be trusted")
        out.append("  individually. Subtract them from both counts to get a safe range.")
    else:
        out.append("")
        out.append("  No row landed near the deamidation/A1 boundary (%.4f +- %.3f), so the"
                   % (BOUNDARY, BOUNDARY_MARGIN))
        out.append("  two are cleanly separated in this data.")

    if unexplained:
        out.append("")
        out.append("  Deltas matching no isotope and no modification in the table above.")
        out.append("  Each recurring value is probably one modification worth naming:")
        for step, n in unexplained.most_common(10):
            seq, z, raw, calc_mz = unexplained_eg[step]
            out.append("    %+9.2f Da  x%-6d  %s %d+" % (
                step, n, seq[:30] + ("..." if len(seq) > 30 else ""), z))
    return {"a1_ions": a1_ions, "any_isotope_ions": any_ions, "rows": total,
            "measured": measured, "computed": identical}


def slope(xs, ys):
    """Least-squares slope, or None if x never varies."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx <= 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def test_ar1(con, cfg, limit, out):
    """Test 4: ar1, the product intensity over the precursor intensity on the same row.

    The quantity asked for is per-row and needs no fitting. What it needs is for
    the product intensity to be ABSOLUTE. A pipeline that normalizes each
    spectrum — dividing every fragment by the spectrum's total, or by its base
    peak — destroys the scale that makes ar1 mean anything, and the published
    PAX8 fragment tables do exactly that, summing to 1.0000 per precursor. So the
    first thing this does is check for it, by summing ar1 across each precursor:
    sums piling up at 1.0 mean the fragments are shares of their own spectrum,
    and the ratio carries no information about the precursor at all.
    """
    out.append("")
    out.append("TEST 4 - AR1, PRODUCT OVER PRECURSOR")
    r = cfg["roles"]
    prod = r.get("product_intensity")
    prec = r.get("intensity")
    if prod and prec and prod == prec:
        prod = None
    if not prod or not prec:
        missing = "a product intensity" if not prod else "a precursor intensity"
        out.append("  SKIPPED: no column for %s in %s. ar1 needs both on the same row;"
                   % (missing, cfg["table"]))
        out.append("  pass --product-intensity / --intensity if they are named differently.")
        return False

    cols = [r["seq"], r["charge"], prod, prec]
    if r.get("frag_type"):
        cols.append(r["frag_type"])
    sql = 'SELECT %s FROM "%s" WHERE "%s" > 0 AND "%s" > 0' % (
        ", ".join('"%s"' % c for c in cols), cfg["table"], prod, prec)
    if limit:
        sql += " LIMIT %d" % limit

    ratios = []
    per_precursor = defaultdict(float)     # ar1 summed, for reporting
    prod_sum = defaultdict(float)          # raw product intensity summed, for the guard
    prod_max = defaultdict(float)
    by_type = defaultdict(list)
    for row in con.execute(sql):
        seq, z, ip, iprec = row[0], row[1], row[2], row[3]
        ftype = row[4] if len(row) > 4 else None
        try:
            ar1 = float(ip) / float(iprec)
        except (TypeError, ValueError, ZeroDivisionError):
            continue
        if ar1 <= 0:
            continue
        ratios.append(ar1)
        key = (seq, z)
        per_precursor[key] += ar1
        prod_sum[key] += float(ip)
        prod_max[key] = max(prod_max[key], float(ip))
        if ftype:
            by_type[str(ftype)].append(ar1)

    if len(ratios) < 20:
        out.append("  SKIPPED: only %d rows carried both intensities." % len(ratios))
        return False

    sums = sorted(per_precursor.values())
    n_prec = len(prod_sum)
    out.append("  %d rows, %d precursors." % (len(ratios), n_prec))

    # The guard is on the RAW product intensities, not on ar1. A spectrum
    # normalized to its total has sum(product) == 1 per precursor; one normalized
    # to its base peak has max(product) == 1. Summing ar1 instead would test
    # 1/precursor, which is not a constant and catches neither.
    by_total = sum(1 for v in prod_sum.values() if abs(v - 1.0) < 0.01)
    by_base = sum(1 for v in prod_max.values() if abs(v - 1.0) < 0.01)
    for count, how, what in ((by_total, "its own total", "sum"),
                             (by_base, "its base peak", "largest fragment")):
        if n_prec and count / n_prec >= 0.5:
            out.append("")
            out.append("  THE PRODUCT INTENSITIES ARE NORMALIZED. %d of %d precursors (%.0f%%)"
                       % (count, n_prec, 100.0 * count / n_prec))
            out.append("  have a %s of exactly 1.000, so each spectrum was divided by"
                       % what)
            out.append("  %s before storage. That removes the scale ar1 depends on:" % how)
            out.append("  dividing a share by a precursor intensity returns the share,")
            out.append("  scaled, and says nothing about the product-to-precursor ratio.")
            out.append("  ar1 needs absolute fragment intensities. This is the defect the")
            out.append("  published PAX8 fragment tables carry — they sum to 1.0000.")
            return False

    ratios.sort()
    def q(f):
        return ratios[min(int(f * len(ratios)), len(ratios) - 1)]
    out.append("  ar1 per row     median %.4f   p05 %.4f   p95 %.4f   max %.3f"
               % (q(0.5), q(0.05), q(0.95), ratios[-1]))
    out.append("  ar1 > 1 (product above precursor): %d rows, %.1f%%"
               % (sum(1 for v in ratios if v > 1), 100.0 * sum(1 for v in ratios if v > 1) / len(ratios)))
    out.append("  ar1_sum per precursor   median %.4f   p95 %.4f"
               % (sums[len(sums)//2], sums[min(int(0.95*len(sums)), len(sums)-1)]))
    if by_type:
        out.append("")
        out.append("  By fragment type (median ar1, rows):")
        for t in sorted(by_type, key=lambda t: -len(by_type[t]))[:8]:
            v = sorted(by_type[t])
            out.append("    %-12s %.4f   %d" % (t[:12], v[len(v)//2], len(v)))
    out.append("")
    out.append("  Loading invariance is NOT tested here: it needs the same ion at two known")
    out.append("  loads, which needs a load map and a run column. ar1 is only expected to")
    out.append("  hold steady inside the linear range, so a stated range must come first.")
    return True


def test_load(con, cfg, load_map, limit, out):
    """Test 3: does precursor intensity track how much was loaded?"""
    out.append("")
    out.append("TEST 3 - LOAD RESPONSE")
    r = cfg["roles"]
    if not r["intensity"]:
        out.append("  SKIPPED: no intensity column in %s. Pass --table/--intensity if one exists"
                   " elsewhere." % cfg["table"])
        return False
    if not r["run"]:
        out.append("  SKIPPED: no run column, so intensities cannot be tied to a load.")
        return False
    if not load_map:
        out.append("  SKIPPED: needs --load-map, a JSON object mapping run name (or any")
        out.append("  substring of it) to nanograms loaded, e.g. {\"HeLa_120\": 120}.")
        return False

    def ng_for(run):
        run = str(run)
        if run in load_map:
            return load_map[run]
        for frag, ng in load_map.items():
            if frag in run:
                return ng
        return None

    sql = 'SELECT "%s", "%s", "%s", "%s" FROM "%s" WHERE "%s" > 0' % (
        r["seq"], r["charge"], r["run"], r["intensity"], cfg["table"], r["intensity"])
    if limit:
        sql += " LIMIT %d" % limit
    best = defaultdict(dict)   # (seq, z) -> ng -> max intensity
    unmatched = Counter()
    for seq, z, run, inten in con.execute(sql):
        ng = ng_for(run)
        if ng is None:
            unmatched[str(run)] += 1
            continue
        try:
            z, inten, ng = int(z), float(inten), float(ng)
        except (TypeError, ValueError):
            continue
        if inten <= 0 or ng <= 0:
            continue
        key = (seq, z)
        if best[key].get(ng, 0) < inten:
            best[key][ng] = inten

    slopes = []
    for key, per_load in best.items():
        if len(per_load) < 3:
            continue
        loads = sorted(per_load)
        s = slope([math.log10(l) for l in loads],
                  [math.log10(per_load[l]) for l in loads])
        if s is not None:
            slopes.append(s)

    if unmatched:
        out.append("  %d runs matched no load-map entry (%d rows dropped): %s" % (
            len(unmatched), sum(unmatched.values()),
            ", ".join(list(unmatched)[:4])))
    if len(slopes) < 5:
        out.append("  SKIPPED: only %d ions appear at 3 or more distinct loads. A slope needs"
                   " at least three points." % len(slopes))
        return False
    slopes.sort()
    q = statistics.quantiles(slopes, n=4) if len(slopes) >= 4 else [None] * 3
    out.append("  %d ions seen at 3+ loads." % len(slopes))
    out.append("  median slope %+.3f   quartiles %+.3f to %+.3f" % (
        statistics.median(slopes), q[0], q[2]))
    out.append("  %d of %d above 0.5 (%.0f%%)" % (
        sum(1 for s in slopes if s > 0.5), len(slopes),
        100.0 * sum(1 for s in slopes if s > 0.5) / len(slopes)))
    out.append("  Proportional response is 1.0. The published data gave +0.239; if this")
    out.append("  column reproduces that, the compression is real and belongs in the")
    out.append("  contract. If it lands near 1.0, the published intensity was the problem.")
    return True


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", help="path to ARiA_Consensus.db")
    p.add_argument("--find", action="store_true",
                   help="search the usual places for a consensus database instead of"
                        " being told the path; the largest match wins")
    p.add_argument("--schema-only", action="store_true",
                   help="print every table, its row count and which role each column filled")
    p.add_argument("--table", help="override the table the tests read")
    p.add_argument("--raw-mz", help="column holding the precursor m/z as reported")
    p.add_argument("--calc-mz", help="column holding the calculated precursor m/z")
    p.add_argument("--calc-mass", help="column holding the calculated neutral mass")
    p.add_argument("--seq", help="peptide sequence column")
    p.add_argument("--charge", help="precursor charge column")
    p.add_argument("--run", help="run or file column")
    p.add_argument("--intensity", help="precursor intensity column")
    p.add_argument("--product-intensity", dest="product_intensity",
                   help="product (fragment) intensity column, for ar1")
    p.add_argument("--frag-type", dest="frag_type", help="fragment type column, for ar1")
    p.add_argument("--load-map", help="JSON file mapping run name or substring to ng loaded")
    p.add_argument("--limit", type=int, default=0, help="read at most N rows per test")
    p.add_argument("--order", choices=("spread", "oldest", "newest"), default="spread",
                   help="which rows a --limit reads: spread evenly across the table"
                        " (default), the first rows SQLite reaches, or the last")
    args = p.parse_args(argv)

    if not args.db:
        if not args.find:
            p.error("pass --db PATH, or --find to search for it")
        hits = find_databases()
        if not hits:
            print("Searched for a consensus database and found none.")
            print("Looked under: %s" % ", ".join(FIND_ROOTS))
            print("for a .db, .sqlite or .sqlite3 file whose name contains"
                  " 'consensus' or 'aria'.")
            print("If it is somewhere else, or called something else, pass --db PATH.")
            return 1
        args.db = hits[0][0]
        print("Found %d candidate database%s; using the largest:" % (
            len(hits), "" if len(hits) == 1 else "s"))
        for path, size in hits[:5]:
            print("  %6.1f MB  %s%s" % (size / 1e6, path,
                                        "   <- using this" if path == args.db else ""))
        print()

    try:
        con = sqlite3.connect("file:%s?mode=ro" % args.db, uri=True)
        con.execute("SELECT 1 FROM sqlite_master LIMIT 1")
    except sqlite3.Error as exc:
        print("Could not open %s: %s" % (args.db, exc), file=sys.stderr)
        return 2

    entries = survey(con)
    if args.schema_only:
        for e in entries:
            print("%s  (%d rows)" % (e["table"], e["rows"]))
            filled = {r: c for r, c in e["roles"].items() if c}
            print("  columns: %s" % ", ".join(e["columns"]))
            print("  recognized: %s" % (", ".join("%s=%s" % (r, c) for r, c in filled.items())
                                        or "nothing"))
            print()
        return 0

    chosen = None
    if args.table:
        for e in entries:
            if e["table"] == args.table:
                chosen = e
        if chosen is None:
            print("No table named %s. Run --schema-only to see what is there." % args.table,
                  file=sys.stderr)
            return 2
    else:
        ranked = sorted(entries, key=score, reverse=True)
        if ranked and score(ranked[0])[0] > 0:
            chosen = ranked[0]

    if chosen is None:
        print("No table carries a sequence, a charge and a reported m/z together.")
        print("Run --schema-only and pass the columns explicitly.")
        return 1

    roles = dict(chosen["roles"])
    for role, flag in (("raw_mz", args.raw_mz), ("calc_mz", args.calc_mz),
                       ("calc_mass", args.calc_mass), ("seq", args.seq),
                       ("charge", args.charge), ("run", args.run),
                       ("intensity", args.intensity),
                       ("product_intensity", args.product_intensity),
                       ("frag_type", args.frag_type)):
        if flag:
            roles[role] = flag
    cfg = {"table": chosen["table"], "roles": roles}

    out = []
    out.append("ARiA_Consensus.db probe")
    out.append("  database %s" % args.db)
    out.append("  table    %s (%d rows)" % (chosen["table"], chosen["rows"]))
    out.append("  columns  " + ", ".join("%s=%s" % (r, c) for r, c in roles.items() if c))
    missing = [r for r in ("seq", "charge", "raw_mz") if not roles.get(r)]
    if missing:
        out.append("")
        out.append("STOPPED: no column found for %s." % ", ".join(missing))
        print("\n".join(out))
        return 1
    if not (roles.get("calc_mz") or roles.get("calc_mass")):
        out.append("")
        out.append("STOPPED: found a reported m/z but nothing calculated to compare it with.")
        out.append("Pass --calc-mz or --calc-mass. Without one, an isotope pick and a")
        out.append("modification are indistinguishable and the test would be a guess.")
        print("\n".join(out))
        return 1
    if roles.get("calc_mz") and roles.get("calc_mass"):
        out.append("  (two truth columns present; using %s)" % roles["calc_mz"])
    out.append("")

    if args.limit:
        out.append("  sample   %d rows, %s" % (args.limit, {
            "spread": "spread evenly across the table",
            "oldest": "the first rows reached, which on an append-only table are the oldest",
            "newest": "the last rows by rowid",
        }[args.order]))
        out.append("")
    offsets = measure_calibration(fetch(con, cfg, args.limit, args.order), cfg, out)
    test_isotope(fetch(con, cfg, args.limit, args.order), cfg, offsets, out)

    load_map = None
    if args.load_map:
        with open(args.load_map) as fh:
            load_map = json.load(fh)
    ran_load = test_load(con, cfg, load_map, args.limit, out)
    test_ar1(con, cfg, args.limit, out)

    print("\n".join(out))
    return 0 if (ran_load or not args.load_map) else 1


if __name__ == "__main__":
    sys.exit(main())
