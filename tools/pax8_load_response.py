#!/usr/bin/env python3
"""Load response and the constant-load series for the Colorado PAX8 cohort.

This is NOT the consensus probe. The 24 PAX8 runs are not in ARiA_Consensus.db
-- checked against all 863 distinct observations.sample_name values and the
69-row sample_registry: no PAX8 or JP21x name, no Colorado acquisition date. So
Test 3 run against that database is about other runs entirely, and this cohort
needs a tool that reads where it actually lives.

    python3 tools/pax8_load_response.py --apex E:/Colorado_ApexRT \
        --search E:/Colorado_search
    python3 tools/pax8_load_response.py --apex DIR --search DIR --schema-only

THE INTENSITY HAS TO BE JOINED, BECAUSE FRAGPIPE HAS NONE

psm.tsv and ion.tsv carry an Intensity column that is 0 in every row of all 21
searched runs: IonQuant is in the workflow but cannot quantify MGF input. So an
identification and its intensity live in different files and must be joined:

  1  psm.tsv                     peptide, charge, Calculated M/Z, Qvalue, and a
                                 Spectrum string whose NNNNN is an index
  2  idts_out/<run>_plainFP_table.csv
                                 that index -> rt_min, drift, pepmass, charge
  3  Summary/<run>_ApexRT_temp1.csv
                                 Function 1 MS1 peak groups: Mz, ScanTimeMin,
                                 Drift, Intensity, IsSaturated, PeakGroupIndex

The PSM's calculated m/z, and the rt and drift from (2), locate a peak group in
(3); its summed Intensity is the only measured precursor intensity these runs
have. If this tool is pointed at psm.tsv alone it refuses rather than reporting
a column of zeros as abundance.

UNITS, BOTH OF WHICH ARE EASY TO GET WRONG

  psm.tsv Retention is in SECONDS. plainFP_table rt_min is in MINUTES. The join
  uses (2), and cross-checks Retention/60 against it rather than trusting either.

  ApexRT Drift is in TWIMS BINS on this Waters instrument, not 1/K0 and not
  milliseconds. Every drift tolerance here is in bins and is labelled as such.

WHAT THIS COHORT CAN AND CANNOT ANSWER

It is a CONSTANT-LOAD design. 1.0 uL x 0.5 ug/uL, 2.0 x 0.25 and 5.0 x 0.10 all
put 500 ng on column, so the 1:2:5 volumes do not span load. Reading them as a
load ladder invents 0.699 log10 units of load where the real span is 0.000.

  LOAD RESPONSE is two points, 100 ng and 500 ng, with exactly ONE 100 ng run
  per gradient. A ratio, not a fit: I(500)/I(100) against an expected 5.0,
  resting on a single unreplicated low point. No curvature verdict -- two levels
  cannot support one.

  SATURATION does not have to be inferred here. ApexRT sets IsSaturated per
  detection, so the saturated share is reported per load level directly.

  THE CONSTANT-LOAD SERIES is the well-designed part: at 500 ng, three volumes,
  two vials x three replicates each, 60 min only. Intensity should not depend on
  volume. A trend is trap loading or a dilution effect, and it bounds how much of
  any load ratio is really the injection.

  60 min and 120 min are different LC methods and are never pooled.

ApexRT exits 3221225477 (0xC0000005, an access violation) on these runs. Its
outputs were checked read-only and are not truncated: for all 96 files named in
apexrt_done.json the current size equals the size recorded at completion, each
ends with a newline, and each last line has the header's column count. Row counts
are consistent within each group, and the two 100 ng runs are lower than the 500
ng runs as the load would predict.

Exit 0 when every test that had inputs ran, 1 on a finding, 2 when an input could
not be read.
"""

import argparse
import bisect
import csv
import glob
import json
import math
import os
import random
import re
import statistics
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ApexRT temp1/temp2 columns, in order. The Summary variants append the last two.
APEX_COLS = ("Function", "ScanIndex", "ScanTimeMin", "Mz", "MzFwhm", "Drift",
             "DriftFwhm", "QuadMass", "Counts", "Intensity", "IsSaturated",
             "PeakGroupIndex", "LinkCount")

# Join windows. The m/z window is FragPipe's own precursor tolerance; the rt and
# drift windows are the reader's and are stated on every run.
JOIN_PPM = 10.0
JOIN_RT_MIN = 0.2
JOIN_DRIFT_BINS = 2.0
QVALUE_CUT = 0.01

BOOTSTRAP = 2000        # resamples for the ratio's interval
EXPECTED_RATIO = 5.0           # 500 ng over 100 ng
FLAT_TOL = 0.05                # log10 units: below this a volume series is flat
PUBLISHED_RATIO = 5.37         # the figure to recheck, from docs/top3-estimator.md






def read_plainfp(run, apex_dir):
    """idx -> (rt_min, drift_bins, pepmass, charge) from the iDTs spectrum table.

    Three runs name this file <run>.plain_table.csv instead of
    <run>_plainFP_table.csv, so both spellings are looked for and the one found
    is reported.
    """
    base = os.path.join(apex_dir, run, "idts_out")
    for name in ("%s_plainFP_table.csv" % run, "%s.plain_table.csv" % run):
        path = os.path.join(base, name)
        if os.path.exists(path):
            break
    else:
        return None, None, "no plainFP table in %s" % base
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                idx = int(str(row["idx"]).strip())
                rt = float(row["rt_min"])
                drift = float(row["drift"])
            except (KeyError, TypeError, ValueError):
                continue
            z = str(row.get("charge", "")).strip().rstrip("+")
            try:
                z = int(z)
            except ValueError:
                z = None
            try:
                pepmass = float(row.get("pepmass"))
            except (TypeError, ValueError):
                pepmass = None
            out[idx] = (rt, drift, pepmass, z)
    return out, os.path.basename(path), None


SPECTRUM_IDX = re.compile(r"\.(\d+)\.(\d+)\.(\d+)$")


def read_psms(path, qcut=QVALUE_CUT):
    """Confident PSMs: (idx, peptide, charge, calc_mz, qvalue).

    Rejects the Intensity column outright. It is 0 in every row of every searched
    run because IonQuant cannot quantify MGF input, so reading it as abundance
    would report zeros as measurements.
    """
    psms, dropped, decoys, zero_intensity, total = [], 0, 0, 0, 0
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            total += 1
            # Counted over EVERY row, before any filter. The claim being checked
            # is "this column is zero throughout the file", and a count taken
            # after filtering could never equal the row count and so could never
            # support it.
            try:
                if float(row.get("Intensity", 0) or 0) == 0.0:
                    zero_intensity += 1
            except (TypeError, ValueError):
                pass
            if str(row.get("Is Decoy", "")).strip().lower() in ("true", "1", "yes"):
                decoys += 1
                continue
            try:
                q = float(row["Qvalue"])
            except (KeyError, TypeError, ValueError):
                dropped += 1
                continue
            if q > qcut:
                dropped += 1
                continue
            m = SPECTRUM_IDX.search(str(row.get("Spectrum", "")))
            if not m:
                dropped += 1
                continue
            try:
                calc_mz = float(row["Calculated M/Z"])
                z = int(row["Charge"])
            except (KeyError, TypeError, ValueError):
                dropped += 1
                continue
            psms.append((int(m.group(1)), str(row.get("Peptide", "")).strip(),
                         z, calc_mz, q))
    return psms, {"total": total, "dropped": dropped, "decoys": decoys,
                  "zero_intensity": zero_intensity}


def read_apex_groups(path):
    """Function-1 MS1 peak groups: mz, rt, drift, summed intensity, saturated.

    Streams the file and aggregates by PeakGroupIndex, because these tables run
    to 148k-179k rows per run and the whole cohort is millions.
    """
    acc = {}
    rows = 0
    feat = flagged = 0
    flagged_i, clean_i = [], []
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        first = next(reader, None)
        if first is None:
            return [], 0, "empty", {}
        # The file may or may not carry a header; the first field of a data row
        # is the function number.
        header = not str(first[0]).strip().isdigit()
        if not header:
            reader = _chain(first, reader)
        for rec in reader:
            rows += 1
            if len(rec) < 13:
                continue
            try:
                if int(float(rec[0])) != 1:
                    continue
                rt = float(rec[2])
                mz = float(rec[3])
                drift = float(rec[5])
                inten = float(rec[9])
                sat = int(float(rec[10]))
                group = int(float(rec[11]))
            except (TypeError, ValueError):
                continue
            if inten <= 0:
                continue
            feat += 1
            if sat:
                flagged += 1
                flagged_i.append(inten)
            else:
                clean_i.append(inten)
            cur = acc.get(group)
            if cur is None:
                acc[group] = [inten * mz, inten * rt, inten * drift, inten, sat, 1]
            else:
                cur[0] += inten * mz
                cur[1] += inten * rt
                cur[2] += inten * drift
                cur[3] += inten
                cur[4] = max(cur[4], sat)
                cur[5] += 1
    groups = []
    for g, (wmz, wrt, wdr, tot, sat, n) in acc.items():
        groups.append((wmz / tot, wrt / tot, wdr / tot, tot, sat, n))
    groups.sort()
    stats = {"features": feat, "flagged": flagged,
             "lowest_flagged": min(flagged_i) if flagged_i else None,
             "highest_clean": max(clean_i) if clean_i else None}
    return groups, rows, None, stats


def _chain(first, reader):
    yield first
    for rec in reader:
        yield rec


def join_run(run, apex_dir, search_dir, out, ppm=JOIN_PPM, rt_tol=JOIN_RT_MIN,
             drift_tol=JOIN_DRIFT_BINS, qcut=QVALUE_CUT):
    """One run's identified ions with a measured intensity, or a reason there are none."""
    hits = sorted(glob.glob(os.path.join(search_dir, "*%s*" % run, "**", "psm.tsv"),
                            recursive=True))
    psm_path = hits[0] if hits else None
    if not psm_path:
        return None, "no psm.tsv: this run was not searched, or the search produced none"
    inner = os.path.basename(os.path.dirname(psm_path))
    sanitised = inner != "%s_apexrt_idts" % run
    summary = os.path.join(apex_dir, run, "Summary", "%s_ApexRT_temp1.csv" % run)
    if not os.path.exists(summary):
        return None, "no Summary/%s_ApexRT_temp1.csv, so there is no measured intensity" % run

    spectra, table_name, err = read_plainfp(run, apex_dir)
    if err:
        return None, err
    psms, stats = read_psms(psm_path, qcut)
    if not psms:
        return None, ("no PSM at Qvalue <= %g (%d row(s), %d decoy(s))"
                      % (qcut, stats["total"], stats["decoys"]))
    groups, apex_rows, gerr, sat_stats = read_apex_groups(summary)
    if gerr or not groups:
        return None, "could not read peak groups from %s (%s)" % (summary, gerr)

    mzs = [g[0] for g in groups]
    ions, unmatched, saturated, rt_mismatch = {}, 0, 0, 0
    for idx, pep, z, calc_mz, q in psms:
        meta = spectra.get(idx)
        if meta is None:
            unmatched += 1
            continue
        rt, drift = meta[0], meta[1]
        tol = calc_mz * ppm / 1e6
        lo = bisect.bisect_left(mzs, calc_mz - tol)
        hi = bisect.bisect_right(mzs, calc_mz + tol)
        best = None
        for g in groups[lo:hi]:
            if abs(g[1] - rt) > rt_tol or abs(g[2] - drift) > drift_tol:
                continue
            if best is None or g[3] > best[3]:
                best = g
        if best is None:
            unmatched += 1
            continue
        if best[4]:
            saturated += 1
        key = (pep, z)
        prev = ions.get(key)
        # Several PSMs can land on one ion; keep the strongest match, and note
        # that this is a per-ion summary rather than a per-spectrum one.
        if prev is None or best[3] > prev[0]:
            ions[key] = (best[3], bool(best[4]), q)
    info = {"sat": sat_stats, "psm_path": psm_path,
            "inner": inner if sanitised else None,
            "psm_rows": stats["total"], "psm_kept": len(psms),
            "decoys": stats["decoys"], "zero_intensity": stats["zero_intensity"],
            "apex_rows": apex_rows, "groups": len(groups),
            "spectra": len(spectra), "table": table_name,
            "ions": len(ions), "unmatched": unmatched, "saturated": saturated}
    return (ions, info), None


def load_record(path):
    doc = json.load(open(path))
    runs = doc.get("runs") or []
    if not runs:
        raise ValueError("%s has no `runs` array" % path)
    return doc, {r["run"]: r for r in runs}


def match_run(name, record):
    """A feature table's run string against the recorded stems."""
    name = str(name)
    if name in record:
        return name
    base = os.path.splitext(os.path.basename(name))[0]
    if base in record:
        return base
    for stem in record:
        if stem in name or name in stem:
            return stem
    return None



def name_vs_record(record, out):
    """Every run whose file name disagrees with the recorded load.

    A cross-check only. It can never be used to infer load: load is volume times
    concentration, and this cohort varied both.
    """
    bad = [r for r in record.values() if r.get("name_disagrees")]
    out.append("")
    out.append("NAME VERSUS RECORD")
    if not bad:
        out.append("  every file name agrees with the recorded load.")
        return 0
    out.append("  %d of %d run name(s) disagree with the instrument record:"
               % (len(bad), len(record)))
    for r in sorted(bad, key=lambda r: r["run"]):
        out.append("    %s" % r["run"])
        out.append("      recorded %g uL x %g ug/uL = %g ng"
                   % (r["recorded_uL"], r["conc_ug_per_uL"], r["load_ng"]))
        if r.get("note"):
            out.append("      %s" % r["note"])
    out.append("  The names are not used for anything here. Load is volume times")
    out.append("  concentration and this cohort varied both, so a name can never")
    out.append("  establish load.")
    return len(bad)


def bootstrap_median(values, n=BOOTSTRAP, seed=20261017):
    """A percentile interval on the median, because ~30 ions is thin.

    Resamples IONS. The single unreplicated low-load run is one measurement in
    every resample, so its own error is not in this interval and nothing in the
    cohort can put it there.
    """
    if len(values) < 5:
        return None, None
    rng = random.Random(seed)
    k = len(values)
    meds = []
    for _ in range(n):
        meds.append(statistics.median([values[rng.randrange(k)] for _ in range(k)]))
    meds.sort()
    return meds[int(0.025 * n)], meds[min(int(0.975 * n), n - 1)]


def volume_confound(record, out):
    """Is the load contrast also a volume contrast, and by how much?

    On this cohort it is exactly that. The 100 ng runs are 0.2 uL of the 0.5
    ug/uL vial and the matched 500 ng runs are 1.0 uL of the SAME vial, so the
    5x in load is delivered entirely as 5x in volume. The constant-load series
    spans the same 5x in volume at a fixed 500 ng, which is what makes the two
    separable at all: whatever that series shows per decade of volume can be
    subtracted from the load ratio.

    Returns (volume_decades, matched_conc) or (None, None).
    """
    levels = defaultdict(set)
    for r in record.values():
        levels[r["load_ng"]].add((r["conc_ug_per_uL"], r["recorded_uL"]))
    if len(levels) < 2:
        return None, None
    lo, hi = min(levels), max(levels)
    shared = {c for c, _ in levels[lo]} & {c for c, _ in levels[hi]}
    if not shared:
        out.append("  the two load levels share no vial concentration, so the load")
        out.append("  contrast is not a pure volume contrast.")
        return None, None
    conc = sorted(shared)[0]
    v_lo = sorted(v for c, v in levels[lo] if c == conc)[0]
    v_hi = sorted(v for c, v in levels[hi] if c == conc)[0]
    decades = math.log10(v_hi / v_lo)
    out.append("  THE LOAD CONTRAST IS A VOLUME CONTRAST. %g ng is %g uL and %g ng is"
               % (lo, v_lo, hi))
    out.append("  %g uL, both from the same %g ug/uL vial, so the %.1fx in load is"
               % (v_hi, conc, hi / lo))
    out.append("  delivered entirely as %.1fx in volume. Load and injection volume are"
               % (v_hi / v_lo))
    out.append("  the same variable here and no ratio below separates them on its own.")
    out.append("  The constant-load series spans the same volume range at a fixed load,")
    out.append("  which is what makes the subtraction below possible.")
    return decades, conc


def load_response(ions, record, out, volume_slope=None, const_volumes=()):
    """I(500)/I(100) per ion, per gradient. A ratio, not a fit."""
    out.append("")
    out.append("LOAD RESPONSE  (two points, per gradient)")
    by_grad = defaultdict(lambda: defaultdict(list))
    for stem, r in record.items():
        by_grad[r["gradient"]][r["load_ng"]].append(stem)
    decades, conc = volume_confound(record, out)
    findings = 0
    for grad in sorted(by_grad):
        levels = by_grad[grad]
        out.append("  %s min: %s" % (grad, ", ".join(
            "%g ng x%d" % (ng, len(v)) for ng, v in sorted(levels.items()))))
        if len(levels) < 2:
            out.append("    one load level, so there is no response to measure.")
            continue
        lo_ng, hi_ng = min(levels), max(levels)
        lo_runs, hi_runs = levels[lo_ng], levels[hi_ng]
        # Restrict the high side to the runs that share the low side's vial
        # concentration. Then the only thing that differs is volume, the
        # contrast is a matched pair, and the volume subtraction below is exact
        # rather than averaged over whatever volumes happen to be in the group.
        lo_concs = {record[s]["conc_ug_per_uL"] for s in lo_runs}
        matched = [s for s in hi_runs
                   if record[s]["conc_ug_per_uL"] in lo_concs]
        if matched and len(matched) < len(hi_runs):
            out.append("    the %g ng side is restricted to the %s ug/uL vial(s), %d of"
                       % (hi_ng, "/".join("%g" % c for c in sorted(lo_concs)),
                          len(matched)))
            out.append("    %d run(s), so only the VOLUME differs across the contrast."
                       % len(hi_runs))
            out.append("    The other %d run(s) reach %g ng at other volumes, which would"
                       % (len(hi_runs) - len(matched), hi_ng))
            out.append("    mix the injection into the ratio.")
            hi_runs = matched
        elif not matched:
            out.append("    no %g ng run shares a vial concentration with the %g ng"
                       % (hi_ng, lo_ng))
            out.append("    side, so this ratio mixes load with volume irreducibly.")
        if len(lo_runs) == 1:
            out.append("    the %g ng level is a SINGLE run (%s)."
                       % (lo_ng, lo_runs[0]))
            out.append("    Every ratio below rests on that one measurement. Its scatter")
            out.append("    is not knowable and nothing checks it.")
        ratios = []
        for key, per_run in ions.items():
            lo = [per_run[r] for r in lo_runs if r in per_run]
            hi = [per_run[r] for r in hi_runs if r in per_run]
            if lo and hi:
                ratios.append(statistics.median(hi) / statistics.median(lo))
        if len(ratios) < 5:
            out.append("    only %d ion(s) appear at both levels; too few to report."
                       % len(ratios))
            findings += 1
            continue
        ratios.sort()
        expected = hi_ng / lo_ng
        med = statistics.median(ratios)
        q = statistics.quantiles(ratios, n=4) if len(ratios) >= 4 else [None] * 3
        out.append("    %d ion(s) at both levels. median I(%g)/I(%g) = %.3f"
                   % (len(ratios), hi_ng, lo_ng, med))
        out.append("    quartiles %.3f to %.3f.  expected %.3f if the response is"
                   " proportional" % (q[0], q[2], expected))
        lo_ci, hi_ci = bootstrap_median(ratios)
        if lo_ci is not None:
            out.append("    bootstrap 95%% interval on the median: %.3f to %.3f"
                       " over %d resample(s)." % (lo_ci, hi_ci, BOOTSTRAP))
            out.append("    That interval is over IONS only. It does not cover the one")
            out.append("    unreplicated %g ng run, whose own error is not in it and"
                       % lo_ng)
            out.append("    cannot be estimated from this cohort.")
            if lo_ci <= expected <= hi_ci:
                out.append("    %.3f is inside the interval, so this cohort does not"
                           " distinguish" % expected)
                out.append("    the response from proportional.")
        shortfall = med / expected
        out.append("    that is %.0f%% of proportional." % (100 * shortfall))
        used_lo = {record[s]["recorded_uL"] for s in lo_runs}
        used_hi = {record[s]["recorded_uL"] for s in hi_runs}
        if len(used_lo) == 1 and len(used_hi) == 1:
            span = math.log10(list(used_hi)[0] / list(used_lo)[0])
        else:
            span = None
        if span is not None and volume_slope is not None:
            # Subtract the volume effect measured at constant load. Both
            # contrasts span the same volume range, which is the only reason
            # this is possible.
            adj = med / (10 ** (volume_slope * span))
            out.append("    the constant-load series gives %+.3f per decade of volume,"
                       % volume_slope)
            out.append("    and this contrast carries %.3f decade(s) of volume (%g to %g"
                       % (span, list(used_lo)[0], list(used_hi)[0]))
            out.append("    uL), so %.3f of the ratio is the injection."
                       % (10 ** (volume_slope * span)))
            out.append("    LOAD ONLY: %.3f, which is %.0f%% of proportional."
                       % (adj, 100 * adj / expected))
            out.append("    The subtraction extrapolates a slope measured between %g and"
                       % min(const_volumes or [1.0]))
            out.append("    %g uL down to %g uL, where it was not measured. If the"
                       % (max(const_volumes or [5.0]), list(used_lo)[0]))
            out.append("    injection behaves differently at that volume this correction")
            out.append("    is wrong, and nothing here can tell.")
            out.append("    AND THE LOW VOLUME IS THE WEAK POINT. %g uL in a 10 uL loop"
                       % list(used_lo)[0])
            out.append("    under partial-loop injection is far below the range that loop")
            out.append("    is comfortable in. Pickup error there is probably the largest")
            out.append("    uncertainty in every number above, and nothing in this cohort")
            out.append("    measures it: the whole load contrast rests on that one volume.")
        elif span is None and volume_slope is not None:
            out.append("    the volume part was NOT subtracted: the runs on one side of")
            out.append("    this contrast sit at several volumes, so there is no single")
            out.append("    volume span to correct by.")
        elif decades:
            out.append("    the volume part of this ratio was NOT subtracted: the")
            out.append("    constant-load series produced no usable slope, so load and")
            out.append("    injection stay confounded in the number above.")
        if shortfall < 0.8:
            out.append("    The response is compressed. With two levels this cannot be")
            out.append("    separated into detector, ion path or injection -- and the")
            out.append("    constant-load series below bounds the injection part.")
            findings += 1
        out.append("    NO curvature or saturation verdict: two load levels cannot")
        out.append("    support one, whatever the spread above looks like.")
    return findings


def constant_load(ions, record, out, slope_out=None):
    """At a fixed load, does intensity track injection volume?

    slope_out, when given a list, receives the median per-ion slope so the load
    response can subtract the injection's share of its own ratio.
    """
    out.append("")
    out.append("CONSTANT-LOAD SERIES  (one load, several volumes)")
    groups = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for stem, r in record.items():
        groups[(r["gradient"], r["load_ng"])][r["recorded_uL"]][r.get("vial")].append(stem)
    findings = 0
    for (grad, ng), vols in sorted(groups.items()):
        if len(vols) < 2:
            continue
        out.append("  %s min at %g ng, %d volume(s):" % (grad, ng, len(vols)))
        per_vol = {}
        for uL in sorted(vols):
            runs_here = [s for v in vols[uL].values() for s in v]
            vals = []
            for key, per_run in ions.items():
                here = [per_run[s] for s in runs_here if s in per_run]
                if here:
                    vals.append(math.log10(statistics.median(here)))
            if vals:
                per_vol[uL] = (statistics.median(vals), len(runs_here),
                               sorted(vols[uL]))
        if len(per_vol) < 2:
            out.append("    no ion appears at two volumes.")
            continue
        base = per_vol[min(per_vol)][0]
        for uL in sorted(per_vol):
            m, n, vials = per_vol[uL]
            out.append("    %4g uL  n=%d  vials %s  median log10 I %+.3f  (%+.3f)%s"
                       % (uL, n, "/".join(str(v) for v in vials), m, m - base,
                          "   <- single run" if n < 2 else ""))
        rep = {u: t for u, t in per_vol.items() if t[1] >= 2}
        if len(rep) < 2:
            out.append("    fewer than two replicated volumes, so nothing here has a")
            out.append("    scatter to be judged against.")
            continue
        # Per-ion slope of log intensity against log volume, which is the quantity
        # asked for: zero is the expectation at constant load.
        slopes = []
        for key, per_run in ions.items():
            pts = []
            for uL, t in sorted(rep.items()):
                runs_here = [s for v in vols[uL].values() for s in v]
                here = [per_run[s] for s in runs_here if s in per_run]
                if here:
                    pts.append((math.log10(uL), math.log10(statistics.median(here))))
            if len(pts) >= 2 and len({x for x, _ in pts}) >= 2:
                mx = sum(x for x, _ in pts) / len(pts)
                my = sum(y for _, y in pts) / len(pts)
                den = sum((x - mx) ** 2 for x, _ in pts)
                if den > 0:
                    slopes.append(sum((x - mx) * (y - my) for x, y in pts) / den)
        if not slopes:
            out.append("    no ion spans two replicated volumes.")
            continue
        slopes.sort()
        med = statistics.median(slopes)
        q = statistics.quantiles(slopes, n=4) if len(slopes) >= 4 else [None] * 3
        out.append("    per-ion slope of log I against log volume, over %d ion(s):"
                   % len(slopes))
        out.append("      median %+.3f   quartiles %+.3f to %+.3f" % (med, q[0], q[2]))
        if slope_out is not None and ng == max(
                r["load_ng"] for r in record.values()):
            slope_out.append((med, sorted(rep)))
        out.append("      expected 0.000, because the load does not change across these")
        out.append("      volumes. A non-zero slope is the injection, not the load.")
        # vial-to-vial spread at one volume is the replicate noise this must clear
        spreads = []
        for uL, t in rep.items():
            if len(t[2]) >= 2:
                per_vial = []
                for vial, stems in vols[uL].items():
                    vals = [math.log10(statistics.median(
                        [per_run[s] for s in stems if s in per_run]))
                        for key, per_run in ions.items()
                        if any(s in per_run for s in stems)]
                    if vals:
                        per_vial.append(statistics.median(vals))
                if len(per_vial) >= 2:
                    spreads.append(max(per_vial) - min(per_vial))
        if spreads:
            noise = max(spreads)
            out.append("      vial-to-vial spread at a fixed volume: up to %.3f log10"
                       " units." % noise)
            swing = abs(med) * (math.log10(max(rep)) - math.log10(min(rep)))
            out.append("      the volume trend moves %.3f over the volume range."
                       % swing)
            if swing <= noise:
                out.append("      That does not clear the replicate noise, so NO verdict:")
                out.append("      a trend this size is indistinguishable from the vials.")
            else:
                out.append("      That clears the replicate noise. Trap loading or a")
                out.append("      dilution effect, and it bounds how much of the load")
                out.append("      ratio above is really the injection.")
                findings += 1
        else:
            out.append("      only one vial per volume, so there is no replicate noise")
            out.append("      to judge the trend against. No verdict.")
    return findings


def recheck_published(ions, record, out, published=PUBLISHED_RATIO):
    """Recompute the published load ratio with the recorded loads.

    The question behind it: was PAX8_JP213_100ng_60_0-5uL counted as 100 ng? Its
    record says 2.0 uL of 0.25 ug/uL, so 500 ng. If it sat in the 100 ng group a
    500 ng run was averaged into the low end and the ratio came out too small.
    """
    out.append("")
    out.append("RECHECK OF THE PUBLISHED %.2fx" % published)
    sixty = {s: r for s, r in record.items() if r["gradient"] == 60}
    suspect = [s for s, r in sixty.items()
               if r.get("name_disagrees") and "100ng" in s and r["load_ng"] != 100.0]
    if not suspect:
        out.append("  no 60 min run is named 100ng while recording another load, so the")
        out.append("  mis-grouping this checks for cannot have happened.")
        return 0
    lo = [s for s, r in sixty.items() if r["load_ng"] == min(
        r2["load_ng"] for r2 in sixty.values())]
    hi = [s for s, r in sixty.items() if r["load_ng"] == max(
        r2["load_ng"] for r2 in sixty.values())]
    out.append("  %d run(s) named 100ng at 60 min record a different load: %s"
               % (len(suspect), ", ".join(suspect)))
    out.append("  By the record the 60 min low group is %d run(s) and the high group"
               " %d." % (len(lo), len(hi)))

    def ratio(low_runs):
        vals = []
        for key, per_run in ions.items():
            a = [per_run[s] for s in low_runs if s in per_run]
            b = [per_run[s] for s in hi if s in per_run]
            if a and b:
                vals.append(statistics.median(b) / statistics.median(a))
        return (statistics.median(vals), len(vals)) if vals else (None, 0)

    correct, n1 = ratio(lo)
    wrong, n2 = ratio(sorted(set(lo) | set(suspect)))
    if correct is None:
        out.append("  no ion appears in both groups, so nothing can be recomputed.")
        return 1
    out.append("  with the recorded loads:            %.3f  (%d ions)" % (correct, n1))
    if wrong is not None and n2:
        out.append("  with the named-100ng runs grouped:  %.3f  (%d ions)" % (wrong, n2))
        out.append("  The published figure was %.2f." % published)
        near_wrong = abs(wrong - published) < abs(correct - published)
        if near_wrong:
            out.append("  The published value sits closer to the MIS-GROUPED figure, which")
            out.append("  is consistent with a 500 ng run having been counted as 100 ng and")
            out.append("  dragging the low end up, so the ratio came out too small.")
            out.append("  Consistent with, not proof of: the published number was computed")
            out.append("  elsewhere and this only shows the grouping would do it.")
        else:
            out.append("  The published value sits closer to the correctly grouped figure,")
            out.append("  so the mis-grouping does not explain it. Look elsewhere.")
    return 1


def saturation_report(flags, infos, record, out):
    """ApexRT's IsSaturated flag: per run, per gradient, and where it starts.

    THE INSTRUMENT'S FLAG, NOT GROUND TRUTH. ApexRT sets it; nothing here
    verifies it, and a detector can clip without setting it or set it on a
    detection that was fine. It is evidence, and better evidence than an
    inference from two load levels, which is all this cohort could otherwise
    offer.
    """
    out.append("")
    out.append("SATURATION  (ApexRT IsSaturated)")
    out.append("  This is the instrument's own flag. Nothing here checks it: a detector")
    out.append("  can clip without setting it, or set it on a detection that was fine.")
    out.append("  It is evidence, not ground truth -- and it is the only direct")
    out.append("  saturation evidence this cohort has, because two load levels cannot")
    out.append("  support a curvature verdict.")

    out.append("  per run, over ALL function-1 features (not just the matched ions):")
    by_grad = defaultdict(lambda: [0, 0, [], []])
    any_feat = False
    for run in sorted(infos):
        sat = infos[run].get("sat") or {}
        rec = record.get(run)
        feat, flg = sat.get("features") or 0, sat.get("flagged") or 0
        if not feat:
            continue
        any_feat = True
        out.append("    %-44s %8d feature(s), %6d flagged (%5.2f%%)"
                   % (run[:44], feat, flg, 100.0 * flg / feat))
        if rec:
            g = by_grad[rec["gradient"]]
            g[0] += feat
            g[1] += flg
            if sat.get("lowest_flagged") is not None:
                g[2].append(sat["lowest_flagged"])
            if sat.get("highest_clean") is not None:
                g[3].append(sat["highest_clean"])
    if not any_feat:
        out.append("    no run reported a feature count.")
        return 0
    out.append("  per gradient, with the intensity bracket the flag starts in:")
    for grad in sorted(by_grad):
        feat, flg, lows, highs = by_grad[grad]
        out.append("    %s min: %d of %d feature(s) flagged (%.2f%%)"
                   % (grad, flg, feat, 100.0 * flg / feat if feat else 0.0))
        if lows and highs:
            out.append("      lowest flagged intensity %s, highest unflagged %s."
                       % (_num(min(lows)), _num(max(highs))))
            out.append("      The flag therefore turns on somewhere between those two;")
            out.append("      they overlap when it is not a simple threshold.")
        elif lows:
            out.append("      lowest flagged intensity %s; nothing is unflagged."
                       % _num(min(lows)))
        else:
            out.append("      nothing flagged at this gradient.")

    by_level = defaultdict(lambda: [0, 0])
    for run, per_ion in flags.items():
        rec = record.get(run)
        if not rec:
            continue
        key = (rec["gradient"], rec["load_ng"])
        for _, sat in per_ion.items():
            by_level[key][0] += 1
            by_level[key][1] += 1 if sat else 0
    if not by_level:
        return 0
    out.append("  among the matched, identified ions:")
    findings = 0
    for (grad, ng), (n, sat) in sorted(by_level.items()):
        out.append("    %s min, %g ng: %d of %d ion(s) flagged (%.1f%%)"
                   % (grad, ng, sat, n, 100.0 * sat / n if n else 0.0))
    levels = sorted(by_level)
    if len(levels) >= 2:
        lo, hi = by_level[levels[0]], by_level[levels[-1]]
        lo_share = lo[1] / lo[0] if lo[0] else 0.0
        hi_share = hi[1] / hi[0] if hi[0] else 0.0
        if hi_share > lo_share + 0.02:
            out.append("  The higher load flags more (%.1f%% against %.1f%%)."
                       % (100 * hi_share, 100 * lo_share))
            out.append("  Any load ratio over flagged ions is a FLOOR on the true ratio,")
            out.append("  because the high end was clipped. The ratio test below splits")
            out.append("  on the flag for exactly this reason.")
            findings += 1
        else:
            out.append("  The flagged share does not rise with load (%.1f%% against"
                       " %.1f%%)," % (100 * hi_share, 100 * lo_share))
            out.append("  so clipping is not depressing the high end of the ratio.")
    out.append("  Flagged ions must stay OUT of any within-run normalization set"
               " (item 5):")
    out.append("  an ion whose high end was clipped cannot be used to cancel load.")
    return findings


def _num(v):
    if v is None:
        return "-"
    return "%.3g" % v


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--apex", help="ApexRT/iDTs root, holding one directory per run")
    p.add_argument("--search", help="FragPipe root, holding <run>_apexrt_idts per run")
    p.add_argument("--loads", default=os.path.join(ROOT, "examples",
                                                   "pax8_recorded_loads.json"),
                   help="recorded load record (default: examples/pax8_recorded_loads.json)")
    p.add_argument("--schema-only", action="store_true",
                   help="report what was found per run and stop, measuring nothing")
    p.add_argument("--ppm", type=float, default=JOIN_PPM,
                   help="precursor window for the join (default %g)" % JOIN_PPM)
    p.add_argument("--rt-tol", type=float, default=JOIN_RT_MIN,
                   help="retention window in MINUTES (default %g)" % JOIN_RT_MIN)
    p.add_argument("--drift-tol", type=float, default=JOIN_DRIFT_BINS,
                   help="drift window in TWIMS BINS, not 1/K0 (default %g)"
                        % JOIN_DRIFT_BINS)
    p.add_argument("--qvalue", type=float, default=QVALUE_CUT,
                   help="PSM Qvalue cut (default %g)" % QVALUE_CUT)
    p.add_argument("--psm", help="a single psm.tsv, which this tool will REFUSE:"
                                 " FragPipe carries no intensity for these runs")
    args = p.parse_args(argv)

    out = []
    try:
        doc, record = load_record(args.loads)
    except (OSError, ValueError, KeyError) as exc:
        print("could not read the load record %s: %s" % (args.loads, exc),
              file=sys.stderr)
        return 2
    counts = {}
    for r in record.values():
        counts[r["load_ng"]] = counts.get(r["load_ng"], 0) + 1
    out.append("load record  %s" % args.loads)
    out.append("  %d run(s), %s" % (len(record), ", ".join(
        "%g ng x%d" % (ng, n) for ng, n in sorted(counts.items()))))
    for run, info in (doc.get("_disputed") or {}).items():
        out.append("  DISPUTED: %s is %s ng by the record, %s ng on the other reading."
                   % (run, info.get("recorded"), info.get("alternative")))
        out.append("    %s" % info.get("why", ""))
        out.append("    Every number below that uses this run inherits the dispute.")

    findings = name_vs_record(record, out)

    if args.psm and not args.apex:
        out.append("")
        out.append("REFUSED: --psm alone cannot measure anything.")
        out.append("  psm.tsv and ion.tsv carry an Intensity column that is 0 in every")
        out.append("  row of all 21 searched runs: IonQuant is in the workflow but")
        out.append("  cannot quantify MGF input. Reading it would report zeros as")
        out.append("  abundance. The intensity is in ApexRT's MS1 peak groups and has")
        out.append("  to be joined: pass --apex as well.")
        print("\n".join(out))
        return 2

    if not (args.apex and args.search):
        out.append("")
        out.append("NO INPUT DIRECTORIES, so nothing was measured.")
        out.append("  --apex   ApexRT/iDTs root  (Scott-PC: E:/Colorado_ApexRT)")
        out.append("  --search FragPipe root     (Scott-PC: E:/Colorado_search)")
        out.append("  Intensity has to come from ApexRT: FragPipe's own Intensity column")
        out.append("  is zero throughout, so an identification and its abundance live in")
        out.append("  different files and are joined on m/z, retention and drift.")
        print("\n".join(out))
        return 1 if findings else 0

    out.append("")
    out.append("JOIN  psm.tsv -> plainFP_table -> ApexRT Summary temp1 (Function 1)")
    out.append("  windows: %.1f ppm, %.2f min, %.1f drift BINS (Waters TWIMS, not 1/K0)"
               % (args.ppm, args.rt_tol, args.drift_tol))
    out.append("  PSMs at Qvalue <= %g, decoys dropped." % args.qvalue)
    ions_by_run, flags, infos, skipped = {}, {}, {}, []
    for run in sorted(record):
        got, why = join_run(run, args.apex, args.search, out, args.ppm,
                            args.rt_tol, args.drift_tol, args.qvalue)
        if got is None:
            skipped.append((run, why))
            continue
        per_ion, info = got
        ions_by_run[run] = {k: v[0] for k, v in per_ion.items()}
        flags[run] = {k: v[1] for k, v in per_ion.items()}
        infos[run] = info
        out.append("    %-44s %5d ion(s) of %4d PSM(s); %d unmatched, %d saturated"
                   % (run[:44], info["ions"], info["psm_kept"], info["unmatched"],
                      info["saturated"]))
        if info["zero_intensity"] == info["psm_rows"] and info["psm_rows"]:
            out.append("      (every psm.tsv Intensity in this run is 0, as expected;"
                       " the intensity above is ApexRT's)")
        if info["table"] and not info["table"].endswith("_plainFP_table.csv"):
            out.append("      spectrum table is named %s, not <run>_plainFP_table.csv"
                       % info["table"])
        if info.get("inner"):
            out.append("      FragPipe's inner folder is %s, not the run stem: the"
                       % info["inner"])
            out.append("      search sanitised the name. Found by searching rather than")
            out.append("      by building the path, which would have called this run")
            out.append("      unsearched.")
    for run, why in skipped:
        out.append("    %-44s SKIPPED: %s" % (run[:44], why))
    if skipped:
        out.append("  %d of %d run(s) contributed nothing. Every number below is over"
                   % (len(skipped), len(record)))
        out.append("  the rest, and a load level can lose its only run this way.")
        findings += 1

    if args.schema_only:
        print("\n".join(out))
        return 0
    if not ions_by_run:
        out.append("  CANNOT RUN: no run yielded an identified ion with an intensity.")
        print("\n".join(out))
        return 2

    # pivot to ion -> run -> intensity
    ions = defaultdict(dict)
    for run, per_ion in ions_by_run.items():
        for key, inten in per_ion.items():
            ions[key][run] = inten
    out.append("  %d distinct ion(s) across %d run(s)" % (len(ions), len(ions_by_run)))

    findings += saturation_report(flags, infos, record, out)
    # The constant-load series runs first: its slope is what lets the load ratio
    # be split into load and injection, and the two contrasts are confounded
    # without it.
    slope_out = []
    const_lines = []
    findings += constant_load(ions, record, const_lines, slope_out)
    findings += load_response(ions, record, out,
                              slope_out[0][0] if slope_out else None,
                              slope_out[0][1] if slope_out else ())
    out.extend(const_lines)
    findings += recheck_published(ions, record, out)
    print("\n".join(out))
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
