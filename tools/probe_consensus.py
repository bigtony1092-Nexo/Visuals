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
  3. LOAD RESPONSE  log-log slope of precursor intensity against the run's load,
     per ion. With --load-map the axis is declared nanograms and 1.0 means a
     proportional response; recomputed from the published data it came out at
     +0.239, which turned out to be the load axis rather than the detector.
     WITHOUT a load map the axis is estimated from the ions themselves and has
     no absolute unit: 1.0 then means "responds like the typical ion", ratios
     between ions survive, and linearity cannot be read at all. Saturation is
     measured inside each ion, because ranking ions by intensity ranks them
     partly by response.
  4. AR1            product intensity over precursor intensity, on the same row.
     Refuses first if the fragments were normalized within their own spectrum,
     which destroys the scale ar1 depends on.
  5. REJECTS        the same classifier run over whatever the ingest gate threw
     away. This is the positive control the other tests cannot supply: rejected
     rows at +1.003/z Da show that picks occur and the gate removes them. A
     rejects table holding only tallies by reason is reported as unusable, not
     read as a negative.

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
REJECT_NAME = re.compile(r"reject|quarantin|discard|excluded|failed", re.I)
REJECT_COUNT_HINTS = ("count", "n", "nrows", "rows", "num", "total", "reason", "rule")

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
            # Bound each slice above as well as below. Without the upper bound a
            # block whose range holds fewer than `per` rows runs on into the next
            # block's rows, so a small table returns the same row many times and
            # every count downstream is inflated.
            stop = int(lo + (b + 1) * step) if b + 1 < blocks else hi + 1
            if stop <= start:
                continue
            queries.append((base + where
                            + ' AND rowid >= %d AND rowid < %d LIMIT %d'
                            % (start, stop, per), ()))
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
    HIST_MAX = PPM_BINS * PPM_BIN      # everything past this lands in one bin

    def pct(q):
        """The |ppm| value at quantile q, from the histogram.

        Returns (value, capped). The last bin is an overflow holding everything
        beyond HIST_MAX, so a quantile landing in it is only a lower bound and
        must say so rather than reporting the cap as a measurement.
        """
        want = q * measured
        seen = 0
        for i, n in enumerate(hist):
            seen += n
            if seen >= want:
                return (i * PPM_BIN, i >= PPM_BINS)
        return (widest_ppm, False)

    qs = [pct(q) for q in (0.5, 0.95, 0.99)]
    out.append("")
    out.append("  Mass error of the measured rows, |ppm|:")
    out.append("    median %s   p95 %s   p99 %s   max %.1f"
               % tuple(["%.1f%s" % (v, "+" if cap else "") for v, cap in qs]
                       + [widest_ppm]))
    if any(cap for _, cap in qs):
        out.append("    (+ means at or past the %.0f ppm histogram limit, so a lower bound)"
                   % HIST_MAX)
    out.append("  Measured rows beyond %.0f ppm: %d of %d."
               % (GATE_PPM, beyond_gate, measured))
    # The edge test must run on the same scale the histogram can express. Taking
    # 95% of an uncapped maximum against capped bins puts the window past every
    # bin, finds zero rows, and calls a cut distribution "ended".
    top = min(widest_ppm, HIST_MAX)
    edge_from = top * (1 - EDGE_FRACTION)
    edge_rows = sum(n for i, n in enumerate(hist) if i * PPM_BIN >= edge_from)
    cut = (measured and edge_rows >= EDGE_MIN
           and edge_rows / measured >= EDGE_SHARE)
    if measured:
        out.append("  In the last 5%% of that range (%.1f to %.1f ppm): %d rows, %.2f%% of all."
                   % (edge_from, top, edge_rows, 100.0 * edge_rows / measured))
    if cut:
        out.append("  That is appreciable density right at the edge, then nothing. The")
        out.append("  distribution was CUT at about %.0f ppm rather than ending there, so"
                   % top)
        out.append("  what lies past the cut was never stored and cannot be counted here.")
    elif measured:
        out.append("  The distribution has thinned to almost nothing by its own maximum, so")
        out.append("  it looks like it ended rather than being cut.")
        out.append("  Read that narrowly: it rules out a cutoff AT %.0f ppm. It says nothing"
                   % top)
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
    # The normalization guard has to key on a SPECTRUM, not on an ion. Keyed on
    # (sequence, charge) alone, a source normalized per spectrum and read across
    # three runs sums to 3.0 rather than 1.0 and the refusal never fires — which
    # is the defect this guard exists to catch, in the guard itself.
    run_col = r.get("run")
    if run_col:
        cols.append(run_col)
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
        rest = list(row[4:])
        run = rest.pop(0) if run_col else None
        ftype = rest.pop(0) if (r.get("frag_type") and rest) else None
        try:
            ar1 = float(ip) / float(iprec)
        except (TypeError, ValueError, ZeroDivisionError):
            continue
        if ar1 <= 0:
            continue
        ratios.append(ar1)
        key = (run, seq, z) if run_col else (seq, z)
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
    out.append("  %d rows, %d %s." % (len(ratios), n_prec,
               "spectra (run, sequence, charge)" if run_col
               else "precursors (NO run column: a per-spectrum normalization "
                    "spanning runs cannot be detected)"))

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


def test_rejects(con, entries, limit, out):
    """Test 5: what the ingest gate threw away.

    This is the positive control the other tests cannot supply. Test 1 reports
    zero isotope picks from tables where a pick could not have been stored, and
    no tail shape distinguishes "none occurred" from "none could get in". The
    quarantine does: if rows rejected on mass carry their reported and calculated
    m/z, the same classifier can be run over them, and a cluster at +1.003/z Da
    is both the control AND, plausibly, where the picks went.

    A rejects table that stores only tallies by reason cannot serve. It records
    that something was dropped, not what it was, so it is reported as unusable
    rather than treated as a negative.
    """
    out.append("")
    out.append("TEST 5 - WHAT THE GATE REJECTED")
    found = [e for e in entries if REJECT_NAME.search(e["table"])]
    if not found:
        out.append("  SKIPPED: no table whose name suggests rejected or quarantined rows.")
        out.append("  Without one there is no positive control: nothing in this database")
        out.append("  shows that a pick-sized error would have been visible had it occurred.")
        return False

    usable = False
    for e in found:
        r = e["roles"]
        have = all(r.get(k) for k in ("seq", "charge", "raw_mz")) and (
            r.get("calc_mz") or r.get("calc_mass"))
        out.append("")
        out.append("  %s (%d rows)" % (e["table"], e["rows"]))
        if not have:
            counts = [c for c in e["columns"]
                      if norm(c) in REJECT_COUNT_HINTS or "count" in norm(c)]
            out.append("    columns: %s" % ", ".join(e["columns"]))
            if counts:
                out.append("    TALLIES ONLY. This records that rows were dropped and why, not")
                out.append("    what they were. A count cannot be classified, so it cannot show")
                out.append("    that a pick-sized error would have been caught. To serve as a")
                out.append("    control it would have to keep the rejected row's reported and")
                out.append("    calculated m/z, sequence and charge.")
            else:
                out.append("    No reported and calculated m/z together, so the rejected rows")
                out.append("    cannot be classified. Not usable as a control.")
            continue

        usable = True
        cfg = {"table": e["table"], "roles": r}
        buckets = Counter()
        picks = 0
        total = 0
        widest = 0.0
        example = None
        for seq, z, raw, calc_mz, _rest in fetch(con, cfg, limit, "spread"):
            total += 1
            delta = (raw - calc_mz) * z
            label, _k = classify(delta)
            buckets[label] += 1
            widest = max(widest, abs(raw - calc_mz) / calc_mz * 1e6)
            if label.startswith("A"):
                picks += 1
                if example is None:
                    example = (seq, z, raw, calc_mz, delta)
        if not total:
            out.append("    No row carried both a reported and a calculated m/z.")
            continue
        out.append("    %d rejected rows classified; widest error %.0f ppm." % (total, widest))
        for label, n in buckets.most_common(8):
            out.append("      %-24s %8d" % (label, n))
        out.append("")
        if picks:
            seq, z, raw, calc_mz, delta = example
            out.append("    %d REJECTED ROWS ARE ISOTOPE PICKS. This is the positive control:"
                       % picks)
            out.append("    the pipeline does produce them, and the gate is what removes them,")
            out.append("    so a zero in test 1 measures the gate and not the pipeline.")
            out.append("    These rows are also the place to count picks from.")
            out.append("    e.g. %s %d+ %.4f vs %.4f (%+.4f)"
                       % (seq[:30], z, raw, calc_mz, delta))
        else:
            out.append("    No rejected row is a pick either. Combined with test 1 that is a")
            out.append("    real negative for this source — errors of a pick's magnitude were")
            out.append("    neither stored nor rejected, so none arrived. Check the gate")
            out.append("    actually ran on these files before leaning on it.")
    return usable


# Test 3's load axis, when no ng values exist. A run's load is estimated from the
# ions themselves: per run, the median over ions of log10(this run's intensity /
# that ion's reference intensity), where the reference is the ion's median across
# the runs it appears in.
#
# The axis so built has NO ABSOLUTE SCALE, and that is the single most important
# thing about the result. The factor is the median of the very ratios the slopes
# are computed from, so every slope comes out divided by the cohort's typical
# response. On this file's own planted data -- half the ions at a true 0.25, half
# at a true 1.00 -- the estimated axis reports 0.51 and 2.06. Both are wrong by
# the same ~2.06x, and their ratio, 4.04 against a true 4.00, is right.
#
# So on an estimated axis:
#   - a slope of 1.0 does NOT mean a linear response. It means this ion responds
#     like the typical ion in the cohort, whatever that response is.
#   - RATIOS between ions are preserved, so ions can be ranked and compared, and
#     a set that tracks load can be told from one that does not.
#   - any compression shared by every ion is invisible. The factor absorbs it
#     exactly, and no amount of data here will show it.
#   - only the injection-volume check can put a scale back on the axis, because
#     it is the one input that does not come from the intensities.
#
# Every threshold below is therefore stated RELATIVE to the cohort's median
# slope, never as an absolute number: an absolute cut on an arbitrary scale is a
# cut at an arbitrary place.
MIN_RUNS_FOR_REFERENCE = 3     # an ion needs this many runs to anchor a reference
MIN_RUN_SHARE = 0.6            # and must appear in this share of runs to vote on a load
# The load factor is a trimmed mean of the central half, not a median. Measured
# against the known truth on this file's own two-mode fixture (half the ions at a
# true 0.25, half at 1.00), worst departure from proportionality in log10 units:
#
#   median              0.0284      plain mean          0.0061
#   trimmed 25%         0.0075
#
# and under a 5% contamination of 300x wild ions, the worst run's factor moves:
#
#   median              0.0009      plain mean          0.0445
#   trimmed 25%         0.0007
#
# So the trimmed mean is four times less warped than the median and sixty times
# less outlier-sensitive than the mean. It matters downstream, not just on paper:
# on the clean fixture a median axis made the curvature test report SATURATION
# where nothing was planted, and a trimmed-mean axis does not.
TRIM_FRACTION = 0.25
# A quartile gap is read as a fraction of the cohort's median slope, because the
# estimated axis fixes no unit for an absolute gap to be measured in.
NORMALIZATION_FRACTION = 0.85  # an ion must reach this share of the median to normalize with
# Curvature is measured inside one ion, as a share of that ion's own slope, so
# the axis's arbitrary unit cancels. Splitting an ion's loads in half needs four.
MIN_LOADS_FOR_BEND = 4
BEND_FRACTION = 0.25
# A run name's gradient length, e.g. the 60 in PAX8_JP210_500ng_60_1uL_..._1.
# Taken as the only bare integer token in a plausible gradient range, which skips
# the acquisition date, the replicate index and the ng and uL tokens (those carry
# units). Reported, never guessed at: if the rule finds nothing the runs are not
# split and the output says so.
GRADIENT_MIN, GRADIENT_MAX = 10, 600
TOKEN = re.compile(r"[^0-9A-Za-z]+")


def gradient_of(run):
    """The gradient length in a run name, or None."""
    hits = []
    for tok in TOKEN.split(str(run)):
        if tok.isdigit() and GRADIENT_MIN <= int(tok) <= GRADIENT_MAX:
            hits.append(int(tok))
    return hits[0] if len(hits) == 1 else None


# A volume may be written 0.5uL or 0-5uL; the second form is what the PAX8 run
# names use. Matching only the plain decimal read "0-5uL" as FIVE microlitres and
# filed it under a label ending "_0", which is a 10x error on the one axis that
# is supposed to anchor the scale.
UL = re.compile(r"(\d+(?:[.-]\d+)?)\s*(?:u|\u00b5)l", re.I)


def microlitres(text):
    """The volume in a run name, or None. Accepts 0.5uL and 0-5uL alike."""
    m = UL.search(str(text))
    if not m:
        return None, None
    try:
        v = float(m.group(1).replace("-", "."))
    except ValueError:
        return None, None
    if v <= 0:
        return None, None
    label = UL.sub("_", str(text))
    label = re.sub(r"[_\-. ]{2,}", "_", label).strip("_- .")
    return v, label


def trimmed_mean(values, frac=TRIM_FRACTION):
    """Mean of the central share. Median-robust, mean-linear; see TRIM_FRACTION.

    int(n * frac) trims nothing at n = 3 and one from each end at n = 5, which
    keeps 100% and 60% where the name says 50%. The share actually kept is
    returned alongside so a caller can print what happened rather than the
    nominal fraction.
    """
    xs = sorted(values)
    k = int(len(xs) * frac)
    core = xs[k:len(xs) - k] or xs
    return sum(core) / len(core)


def trimmed_kept(n, frac=TRIM_FRACTION):
    """The share of n values trimmed_mean actually averages."""
    k = int(n * frac)
    kept = n - 2 * k
    return (kept if kept > 0 else n) / float(n)


def estimate_load_factors(per_ion_run, out):
    """Each run's load in log10 units, estimated from the ions it shares.

    Returns (factors, anchors) where factors maps run -> log10 load relative to
    the cohort, and anchors is the number of ions that voted.
    """
    runs = sorted({run for per_run in per_ion_run.values() for run in per_run})
    need = max(MIN_RUNS_FOR_REFERENCE, int(math.ceil(MIN_RUN_SHARE * len(runs))))
    anchors = {key: per_run for key, per_run in per_ion_run.items()
               if len(per_run) >= need}
    if not anchors:
        out.append("  no ion appears in %d of the %d runs, so no run's load can be"
                   " estimated from shared ions." % (need, len(runs)))
        return None, 0
    ref = {key: statistics.median(per_run.values()) for key, per_run in anchors.items()}
    factors = {}
    for run in runs:
        ratios = [math.log10(per_run[run] / ref[key])
                  for key, per_run in anchors.items()
                  if run in per_run and ref[key] > 0]
        if len(ratios) >= MIN_RUNS_FOR_REFERENCE:
            factors[run] = trimmed_mean(ratios)
    out.append("  load estimated from %d ion(s) present in >=%d of %d runs,"
               " trimmed mean of the central %d%%."
               % (len(anchors), need, len(runs), 100 - 2 * int(100 * TRIM_FRACTION)))
    if len(factors) < len(runs):
        out.append("  %d run(s) had too few shared ions to place on the load axis."
                   % (len(runs) - len(factors)))
    return factors, len(anchors)


def volume_check(factors, out):
    """What the injection volumes in run names do and do not say about load.

    They do not say load. Load is volume x concentration, and on the Colorado
    PAX8 runs the two were traded off against each other deliberately: the 1, 2
    and 5 uL injections came from vials at 0.5, 0.25 and 0.1 ug/uL, so all three
    put 500 ng on column. Read as a 1:2:5 ladder they span 0.699 log10 units of
    load; the real span is 0.000.

    Two of those names are wrong in the other direction as well: a run named
    "1uL" recorded 0.2 uL of 0.5 ug/uL (100 ng), and one named "100ng_0-5uL"
    recorded 2.0 uL of 0.25 ug/uL (500 ng).

    So this reports the volumes and refuses to convert them. It returns None
    always: a volume spread can never anchor the load axis. Only a declared load
    per run can, which is what --load-map is for.
    """
    groups = defaultdict(list)
    for run, f in factors.items():
        v, label = microlitres(run)
        if v is not None:
            groups[label].append((v, run, f))
    usable = {k: sorted(v) for k, v in groups.items() if len(v) >= 2}
    found = sum(len(v) for v in groups.values())
    if not found:
        out.append("  no run name carries an injection volume. Nothing here is"
                   " anchored outside the intensities; only a declared load per"
                   " run (--load-map) can do that.")
        return None
    if not usable:
        # The real PAX8 names hit this: all 24 carry a volume, and no two share a
        # label because the vial number and replicate index stay in it. Printing
        # "no run name carries a volume" here was simply false.
        out.append("  %d run name(s) carry an injection volume, but no two share a"
                   " label once the volume is removed, so there is nothing to" % found)
        out.append("  compare. The vial number and replicate index stay in the label:")
        for label in sorted(groups)[:4]:
            out.append("    %s  <- %s" % (label, groups[label][0][1]))
        out.append("  Volume could not tell you load here even if they did group:"
                   " load is")
        out.append("  volume x concentration, and this cohort varied both. Use"
                   " --load-map.")
        return None
    out.append("  injection volumes appear in the run names. They are REPORTED, not")
    out.append("  converted to load: load is volume x concentration, and on the")
    out.append("  Colorado PAX8 runs the two were traded off deliberately -- 1, 2 and")
    out.append("  5 uL drawn from 0.5, 0.25 and 0.1 ug/uL vials, all 500 ng on column.")
    out.append("  Reading them as a 1:2:5 ladder would invent 0.699 log10 units of load")
    out.append("  where the real span is 0.000. Two of those names are wrong the other")
    out.append("  way too: one \"1uL\" run recorded 0.2 uL (100 ng), and one named")
    out.append("  \"100ng_0-5uL\" recorded 2.0 uL of 0.25 ug/uL (500 ng).")
    for label, members in sorted(usable.items()):
        out.append("  %s:" % (label or "(unlabelled)"))
        base_v, base_run, base_f = members[0]
        for v, run, f in members:
            out.append("    %-28s %4g uL   estimated factor %+.3f (vs %+.3f at %g uL)"
                       % (str(run)[:28], v, f, base_f, base_v))
        spread = max(f for _, _, f in members) - min(f for _, _, f in members)
        ratio = math.log10(max(v for v, _, _ in members)
                           / min(v for v, _, _ in members))
        out.append("    estimated spread %.3f log10 units over a %.3f volume ratio."
                   % (spread, ratio))
        if spread < 0.05 and ratio > 0.1:
            out.append("    The factors barely move across a %.1fx volume range, which is"
                       % (10 ** ratio))
            out.append("    what a constant load across volumes looks like. Consistent with")
            out.append("    concentration compensating volume. NOT a failure of the axis.")
        elif spread > 0.5 * ratio:
            out.append("    The factors move roughly WITH volume. If these runs were")
            out.append("    concentration-compensated to one load, that is a red flag in")
            out.append("    the intensities, not a confirmation of the axis.")
    out.append("  No verdict either way: volume cannot anchor load without the")
    out.append("  concentration, and the concentration is not in these names. Supply")
    out.append("  --load-map with the recorded ng per run.")
    return None


def constant_load_check(per_ion_run, factors, out):
    """At a FIXED load, does intensity depend on injection volume?

    The Colorado design has replicates in exactly one place: 500 ng delivered at
    1, 2 and 5 uL from vials at 0.5, 0.25 and 0.10 ug/uL, three runs each. The
    load is constant, so the expectation is no change. A trend is trap loading or
    a dilution effect, and it is measurable here precisely because the load is
    not varying.
    """
    by_load_vol = defaultdict(lambda: defaultdict(list))
    for run, f in factors.items():
        v, _ = microlitres(run)
        g = gradient_of(run)
        if v is not None:
            by_load_vol[(round(f, 9), g)][v].append(run)
    series = {k: v for k, v in by_load_vol.items() if len(v) >= 2}
    if not series:
        return
    out.append("  CONSTANT-LOAD CHECK: at one load, intensity should not depend on")
    out.append("  injection volume. A trend is trap loading or dilution, not load.")
    for (load, g), vols in sorted(series.items()):
        out.append("    load factor %+.3f%s, %d volume(s):"
                   % (load, ", %s min" % g if g else "", len(vols)))
        meds, reps = [], []
        for v in sorted(vols):
            runs_at = vols[v]
            # median over ions of the median intensity at this volume, in log10
            per_ion = []
            for key, per_run in per_ion_run.items():
                here = [per_run[r] for r in runs_at if r in per_run]
                if here:
                    per_ion.append(math.log10(statistics.median(here)))
            if not per_ion:
                continue
            meds.append((v, statistics.median(per_ion), len(runs_at)))
            reps.append(len(runs_at))
        if len(meds) < 2:
            continue
        base = meds[0][1]
        for v, m, n in meds:
            out.append("      %4g uL  n=%d  median log10 intensity %+.3f  (%+.3f)%s"
                       % (v, (n if m else 0), m, m - base,
                          "   <- single run, not judgeable" if n < 2 else ""))
        # Judge the trend on the REPLICATED volumes only. A singleton column has
        # no scatter of its own, and letting it suppress the verdict loses the
        # finding -- the Colorado series has three replicated volumes and one
        # singleton, and the singleton is the disputed run.
        rep = [(v, m, n) for v, m, n in meds if n >= 2]
        singles = [(v, m) for v, m, n in meds if n < 2]
        if len(rep) < 2:
            spread = max(m for _, m, _ in meds) - min(m for _, m, _ in meds)
            out.append("      fewer than two replicated volumes, so the %.3f log10 spread"
                       % spread)
            out.append("      has no scatter to be judged against.")
            continue
        spread = max(m for _, m, _ in rep) - min(m for _, m, _ in rep)
        out.append("      over the %d replicated volume(s): %.3f log10 units."
                   % (len(rep), spread))
        for v, m in singles:
            near = min(rep, key=lambda t: abs(math.log10(t[0] / v)))
            want = math.log10(v / near[0])
            out.append("      the %g uL run is unreplicated. It sits %+.3f from the %g uL"
                       % (v, m - near[1], near[0]))
            out.append("      median; if its load matched, it would sit near 0, and if it"
                       " carried")
            out.append("      a different load that offset is evidence about which. Where"
                       " a run's")
            out.append("      load is disputed, this is the measurement that speaks to it.")
        if spread < 0.05:
            out.append("      %.3f log10 units across the volumes: flat, as a constant"
                       % spread)
            out.append("      load should be. The dilution and the trap are not showing.")
        else:
            out.append("      %.3f log10 units (%.2fx) across volumes at ONE load. That is"
                       % (spread, 10 ** spread))
            out.append("      not load, because load did not change. Trap loading or a")
            out.append("      dilution effect, and it bounds how much of any load slope")
            out.append("      elsewhere is really the injection.")


def report_normalization_set(slopes, vals, cohort_med, declared, out):
    """Which ions are unfit to normalize with.

    Readable on an estimated axis as well as a declared one, because it is a
    comparison BETWEEN ions and ratios are what that axis preserves. The cut is a
    share of the cohort median, never an absolute slope: on an axis with no unit
    an absolute cut lands in an arbitrary place.
    """
    cut = NORMALIZATION_FRACTION * cohort_med
    unfit = [key for s, _, key in slopes if s < cut]
    out.append("  %d of %d ion(s) respond below %.0f%% of the median slope (%+.3f) and"
               % (len(unfit), len(vals), 100 * NORMALIZATION_FRACTION, cut))
    out.append("  must stay OUT of any within-run normalization set: an ion that does not")
    out.append("  track load cannot be used to cancel it. See docs/top3-estimator.md.")
    if not declared:
        out.append("  That cut is relative on purpose. An absolute one (`slope below 0.85`)")
        out.append("  would be a cut at an arbitrary point on an axis with no unit. It is")
        out.append("  readable here because it compares ions with each other, and ratios")
        out.append("  are what an estimated axis keeps.")


def load_analysis(per_ion_run, factors, declared, out, label=""):
    """The slope, quartile, curvature and normalization analysis for ONE set
    of runs.

    It is a function because it must be applied per gradient. 60 min and 120
    min are different LC methods, and a pooled slope mixes method with load --
    which the output above says plainly. For a while it said that while the
    slope loop, the quartiles and the normalization set all ran over every run
    at once, so the sentence was true of the printed breakdown and false of
    every number under it. Scaling one gradient's intensities by 0.2, which
    cannot change either gradient's own slope, moved the reported median from
    +1.004 to +1.368.
    """
    pre = ("  [%s] " % label) if label else "  "
    n_loads = len({round(f, 9) for f in factors.values()})
    if n_loads < 2:
        out.append("  SKIPPED: every run sits at the same load (%d distinct value),"
                   " so there is no load axis to regress against." % n_loads)
        return False
    if n_loads == 2:
        out.append("  Only TWO distinct loads. Two points determine a slope exactly, so")
        out.append("  every ion's slope below is unchecked: there is no third point to")
        out.append("  disagree with it, no residual, and no curvature. The Colorado PAX8")
        out.append("  cohort is this shape -- 100 ng on two runs and 500 ng on twenty-two,")
        out.append("  because the 1, 2 and 5 uL injections were concentration-compensated")
        out.append("  to one load. Treat each slope as a single measurement, not a fit.")

    slopes = []
    for key, per_run in per_ion_run.items():
        pts = [(factors[run], per_run[run]) for run in per_run if run in factors]
        if len({f for f, _ in pts}) < min(3, n_loads):
            continue
        s = slope([f for f, _ in pts], [math.log10(i) for _, i in pts])
        if s is not None:
            slopes.append((s, statistics.median(per_run.values()), key))
    if len(slopes) < 5:
        out.append("  SKIPPED: only %d ion(s) appear at %d or more distinct loads."
                   % (len(slopes), min(3, n_loads)))
        return False

    # Only a DECLARED load per run anchors the axis. Injection volumes cannot:
    # load is volume x concentration, and the Colorado PAX8 runs traded the two
    # off against each other on purpose. See volume_check.
    anchored = declared
    volume_lines = []
    if not declared:
        volume_check(factors, volume_lines)

    vals = sorted(s for s, _, _ in slopes)
    q = statistics.quantiles(vals, n=4) if len(vals) >= 4 else [None] * 3
    out.append("  %d ions seen at %d+ loads." % (len(vals), min(3, n_loads)))
    out.append("  median slope %+.3f   quartiles %+.3f to %+.3f   IQR %.3f" % (
        statistics.median(vals), q[0], q[2], q[2] - q[0]))

    cohort_med = statistics.median(vals)
    if declared:
        out.append("  Proportional response is 1.0. The published data gave +0.239; if this")
        out.append("  column reproduces that, the compression is real and belongs in the")
        out.append("  contract. If it lands near 1.0, the published intensity was the problem.")
    else:
        out.append("  THE NUMBERS ABOVE HAVE NO ABSOLUTE SCALE. The axis is built from the")
        out.append("  same ratios, so every slope is divided by the cohort's typical")
        out.append("  response: 1.0 means `responds like the typical ion`, NOT `linear`.")
        out.append("  On this tool's planted data a true 0.25 and a true 1.00 come out as")
        out.append("  0.51 and 2.06 -- both wrong by one factor, their 4.0x ratio intact.")
        out.append("  So: ions can be ranked against each other, and a set that tracks load")
        out.append("  can be told from one that does not. Linearity cannot be read at all,")
        out.append("  and a compression shared by every ion is absorbed by the axis and")
        out.append("  invisible. Only the volume check below can put a scale back on this.")

    # Saturation: does slope depend on the ion's own intensity? The construction
    # cannot pin this, because it fixes one centre for all ions at once.
    by_inten = sorted(slopes, key=lambda t: t[1])
    n = len(by_inten)
    quartiles = [by_inten[i * n // 4:(i + 1) * n // 4] for i in range(4)]
    out.append("  slope by the ion's own intensity quartile, DESCRIPTIVE ONLY:")
    meds = []
    for i, group in enumerate(quartiles):
        if not group:
            meds.append(None)
            continue
        qmed = statistics.median([s for s, _, _ in group])
        meds.append(qmed)
        out.append("    Q%d %s  n=%-5d median slope %+.3f" % (
            i + 1, ("dimmest" if i == 0 else "brightest" if i == 3 else "       "),
            len(group), qmed))
    if meds[3] is not None and all(m is not None for m in meds[:3]) and cohort_med:
        rest = statistics.median(meds[:3])
        gap = (rest - meds[3]) / abs(cohort_med)
        out.append("    Q4 is %+.0f%% of the median slope away from the other three."
                   % (100 * gap))
    out.append("  NO VERDICT IS READ FROM THAT TABLE. Ranking ions by intensity ranks")
    out.append("  them partly by response: an ion that gains more per unit load is")
    out.append("  brighter at high load for that reason alone. On this tool's own planted")
    out.append("  data, which contains no saturation and no floor, that ranking sorts 122")
    out.append("  of 125 ions into the right response class -- so the table would `find`")
    out.append("  an effect in data built without one. Ranking by the lowest load instead")
    out.append("  still sorts 104 of 125. The quartiles stay as description.")

    # Saturation measured INSIDE each ion, where amplitude cannot confound it.
    # It still needs an ANCHORED axis. An estimated axis has no certified shape:
    # its scale is arbitrary by design, but it can also depart from
    # PROPORTIONALITY, and that is what bends every ion's curve the same way. On
    # this tool's two-mode fixture the trimmed-mean axis runs
    # 0.000/0.181/0.391/0.577 where the truth is 0.000/0.301/0.681/0.982: a
    # constant rescale would be harmless, and the departure from one is not.
    # Curvature read against such an axis reads the axis rather than the ions --
    # the same compression-shared-by-every-ion that is invisible everywhere else
    # in this test, reappearing as a false verdict. The median estimator made
    # exactly that call on a fixture with nothing planted.
    if not anchored:
        out.append("  Within-ion curvature NOT ASSESSED: it needs an axis whose SHAPE is")
        out.append("  known, and an estimated one has no certified shape. An arbitrary")
        out.append("  scale would be harmless here; a departure from proportionality is")
        out.append("  not, because it bends every ion the same way. On this tool's two-mode")
        out.append("  fixture the axis runs 0.000/0.181/0.391/0.577 against a true")
        out.append("  0.000/0.301/0.681/0.982. With a median estimator that was enough to")
        out.append("  report a downward bend on data built with none. Declare ng with")
        out.append("  --load-map, or get the volume check to agree, and this runs.")
        report_normalization_set(slopes, vals, cohort_med, declared, out)
        out.extend(volume_lines)
        return True
    bends = []
    for key, per_run in per_ion_run.items():
        pts = sorted((factors[run], math.log10(per_run[run]))
                     for run in per_run if run in factors)
        pts = [(f, y) for f, y in pts]
        if len({f for f, _ in pts}) < MIN_LOADS_FOR_BEND:
            continue
        half = (len(pts) + 1) // 2
        lo = slope([f for f, _ in pts[:half]], [y for _, y in pts[:half]])
        hi = slope([f for f, _ in pts[-half:]], [y for _, y in pts[-half:]])
        whole = slope([f for f, _ in pts], [y for _, y in pts])
        if lo is None or hi is None or not whole:
            continue
        # As a share of the ion's own slope, so the arbitrary axis unit cancels.
        bends.append(((hi - lo) / abs(whole), statistics.median(per_run.values()), key))
    if not bends:
        out.append("  Within-ion curvature needs %d distinct loads per ion and no ion has"
                   % MIN_LOADS_FOR_BEND)
        out.append("  that many, so saturation cannot be separated from a spread in")
        out.append("  response at all in this cohort. Not assessed.")
    else:
        bv = sorted(b for b, _, _ in bends)
        down = sum(1 for b in bv if b <= -BEND_FRACTION)
        up = sum(1 for b in bv if b >= BEND_FRACTION)
        out.append("  within-ion curvature over %d ion(s) with %d+ loads: median %+.3f"
                   % (len(bv), MIN_LOADS_FOR_BEND, statistics.median(bv)))
        out.append("    %d (%.0f%%) bend DOWN at high load, %d (%.0f%%) bend UP"
                   % (down, 100.0 * down / len(bv), up, 100.0 * up / len(bv)))
        out.append("  This is each ion against itself, so amplitude cannot confound it.")
        if down >= 0.25 * len(bv) and down > 2 * up:
            bright = [b for b, i, _ in bends
                      if i >= statistics.median([i for _, i, _ in bends])]
            out.append("  A DOWNWARD BEND: a quarter or more of ions flatten at their own")
            out.append("  high end, and far more bend down than up. Median bend among the")
            out.append("  brighter half is %+.3f." % statistics.median(bright))
            out.append("  Two things do this and curvature cannot separate them: headroom")
            out.append("  running out (saturation), and a baseline subtracted from every")
            out.append("  intensity, which costs the low end proportionally more. Check how")
            out.append("  the intensity column was computed before calling it saturation.")
        elif up >= 0.25 * len(bv) and up > 2 * down:
            out.append("  AN UPWARD BEND: ions steepen at their high end. That is a noise")
            out.append("  floor -- below it the measured intensity sticks at the floor")
            out.append("  instead of falling, so the low end flattens.")
            out.append("  This is the opposite shape from saturation, and must not be read")
            out.append("  as a spread in response.")
        else:
            out.append("  Neither direction dominates, so no curvature shows. A compression")
            out.append("  shared by every ion at every load is still invisible: it is in the")
            out.append("  axis, not in the residual.")

    report_normalization_set(slopes, vals, cohort_med, declared, out)


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

    # load_map is None when --load-map was not given, which is the usual case now.
    given = load_map or {}
    notes = {k: v for k, v in given.items() if str(k).startswith("_")}
    runs_map = {k: v for k, v in given.items() if not str(k).startswith("_")}
    disputed = notes.get("_disputed") or {}

    def ng_for(run):
        run = str(run)
        if run in runs_map:
            return runs_map[run]
        for frag, ng in runs_map.items():
            if frag in run:
                return ng
        return None

    sql = 'SELECT "%s", "%s", "%s", "%s" FROM "%s" WHERE "%s" > 0' % (
        r["seq"], r["charge"], r["run"], r["intensity"], cfg["table"], r["intensity"])
    if limit:
        sql += " LIMIT %d" % limit
    # Median per (ion, run), not max: taking the brightest replicate inflated the
    # published load ratio to 5.89x where the median gives 5.37x (top3-estimator.md).
    rows = defaultdict(lambda: defaultdict(list))
    for seq, z, run, inten in con.execute(sql):
        try:
            z, inten = int(z), float(inten)
        except (TypeError, ValueError):
            continue
        if inten <= 0:
            continue
        rows[(seq, z)][str(run)].append(inten)
    per_ion_run = {key: {run: statistics.median(v) for run, v in per_run.items()}
                   for key, per_run in rows.items()}
    if not per_ion_run:
        out.append("  SKIPPED: no usable intensity rows.")
        return False

    declared = bool(runs_map)
    if given and not runs_map:
        out.append("  the load map holds only comment keys (leading underscore) and no"
                   " run entries, so the axis is estimated, not declared.")
    if declared:
        factors, unmatched = {}, Counter()
        for run in {run for per_run in per_ion_run.values() for run in per_run}:
            ng = ng_for(run)
            if ng is None:
                unmatched[run] += 1
            else:
                try:
                    if float(ng) > 0:
                        factors[run] = math.log10(float(ng))
                except (TypeError, ValueError):
                    unmatched[run] += 1
        if unmatched:
            out.append("  %d run(s) matched no load-map entry and are excluded: %s"
                       % (len(unmatched), ", ".join(list(unmatched)[:4])))
        out.append("  load axis is DECLARED ng from --load-map, so the regression below"
                   " is not circular.")
        anchors = None
    else:
        out.append("  no --load-map given, so the load axis is ESTIMATED from the ions.")
        factors, anchors = estimate_load_factors(per_ion_run, out)
        if not factors:
            out.append("  SKIPPED: no load axis could be built.")
            return False

    # 60 min and 120 min are different LC methods. A load slope pooled across them
    # mixes method with load, and on the Colorado cohort each method holds exactly
    # one 100 ng run, so pooling is also the only way to make the low end look
    # replicated when it is not.
    grads = {}
    for run in factors:
        grads.setdefault(gradient_of(run), []).append(run)
    if len(grads) > 1 and None not in grads:
        out.append("  %d gradient length(s) in this cohort: %s. These are different LC"
                   % (len(grads), ", ".join("%s min (%d run%s)"
                                            % (g, len(r), "" if len(r) == 1 else "s")
                                            for g, r in sorted(grads.items()))))
        out.append("  methods, so a pooled load slope would mix method with load. Each is")
        out.append("  reported separately below; none is combined.")
        for g, runs_in in sorted(grads.items()):
            by_load = {}
            for run in runs_in:
                by_load.setdefault(round(factors[run], 9), []).append(run)
            desc = ", ".join("%d run%s" % (len(v), "" if len(v) == 1 else "s")
                             for _, v in sorted(by_load.items()))
            out.append("    %s min: %d load level(s), %s"
                       % (g, len(by_load), desc))
            single = [k for k, v in by_load.items() if len(v) == 1]
            if single and len(by_load) >= 2:
                out.append("      %d level(s) rest on a SINGLE run, so the slope at this"
                           % len(single))
                out.append("      gradient is anchored by one unreplicated point. Its"
                           " scatter is")
                out.append("      not measurable and nothing checks it.")
    elif len(grads) > 1:
        out.append("  some run names carry a gradient length and some do not, so the runs")
        out.append("  were NOT split by method. If this cohort mixes gradients, the slope")
        out.append("  below mixes method with load.")

    if declared and disputed:
        for run, info in disputed.items():
            hit = [r for r in factors if str(run) in str(r) or str(r) in str(run)]
            out.append("  DISPUTED LOAD on %s: the map says %s ng, the alternative"
                       " reading is %s ng." % (str(run)[:44], info.get("recorded"),
                                               info.get("alternative")))
            out.append("    %s" % info.get("why", "no reason recorded"))
            out.append("    %s this cohort. Every number below that uses it inherits the"
                       % ("IN" if hit else "NOT in"))
            out.append("    dispute; re-run with the other value to see what moves.")

    # Per gradient, never pooled. The breakdown printed above says these are
    # different LC methods; running one analysis over both would contradict it.
    if len(grads) > 1 and None not in grads:
        any_ok = False
        for g, runs_in in sorted(grads.items()):
            sub_f = {k: v for k, v in factors.items() if k in set(runs_in)}
            sub_p = {}
            for key, per_run in per_ion_run.items():
                kept = {rn: iv for rn, iv in per_run.items() if rn in sub_f}
                if kept:
                    sub_p[key] = kept
            out.append("")
            out.append("  GRADIENT %s min - %d run(s), analysed alone"
                       % (g, len(runs_in)))
            if load_analysis(sub_p, sub_f, declared, out, label="%s min" % g):
                any_ok = True
        if declared:
            constant_load_check(per_ion_run, factors, out)
            volume_check(factors, out)
        return any_ok

    ok = load_analysis(per_ion_run, factors, declared, out)
    if not ok:
        return False
    if declared:
        constant_load_check(per_ion_run, factors, out)
        volume_check(factors, out)
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
    ran_ar1 = test_ar1(con, cfg, args.limit, out)
    ran_rej = test_rejects(con, entries, args.limit, out)

    skipped = []
    if args.load_map and not ran_load:
        skipped.append("load response")
    if not ran_ar1:
        skipped.append("ar1")
    if not ran_rej:
        skipped.append("rejects")
    if skipped:
        out.append("")
        out.append("NOT RUN: %s. Each says above why. Exit status is 1 because a"
                   % ", ".join(skipped))
        out.append("skipped test is not a passed one.")
    print("\n".join(out))
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
