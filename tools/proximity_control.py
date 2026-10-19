#!/usr/bin/env python3
"""What share of coordinates that CANNOT be right still find a peak?

`bare_proximity` is the tier an ion falls to when nothing decoy-controlled
scored it: a ppm window, a charge, and enough consecutive scans. The tier has
always carried a borrowed sentence -- that the same bare mechanism passed
reversed-tryptic decoys at 33-70% -- attributed to work elsewhere in this
project and describing the retired cube pathway, not the matcher that runs now.
So the false-positive rate of the code that actually runs has never been
measured. This measures it.

Reversed-tryptic decoys cannot do this job, which is worth stating because the
borrowed number came from them. Reversing a peptide permutes its residues and
so PRESERVES its monoisotopic mass exactly; a reversed decoy has the same m/z
as the peptide it decoys, at every charge. Queried through an m/z window it
matches 100% by construction. Whatever the 33-70% measured, it was not this
window. Entrapment -- real peptides from an organism that is not in the sample --
is the construction that gives a different mass on the same peptide lattice.

Four arms, because one number cannot hold them and each answers a different
question:

  target       real targets, observed m/z. The comparator. Not the rate
               aria_run already prints: that loop is PEAK-direction (for each
               peak, scan targets) and this is PROBE-direction (for each probe,
               find a peak). Different denominators, so the comparison is only
               meaningful computed in the same direction.
  entrapment   peptides that cannot be in the sample, their own calculated m/z,
               a real target's charge and coordinates. The primary arm.
  calc_target  the real targets with m/z RECOMPUTED from their sequences. This
               is the diagnostic for the one bias that would otherwise sink the
               whole measurement -- see below.
  swap         target A's m/z with target B's retention time and drift. The only
               arm that asks whether rt and drift contribute anything beyond
               m/z. Emitted only when one of them is actually compared.

Three ways this could have been a check that cannot fail, and what is done
about each:

CALCULATED VERSUS OBSERVED m/z. `make_targets.py` joins a search to a feature
table precisely so a target carries an OBSERVED m/z, because a search reports a
calculated one. An entrapment probe has no observation, so its m/z must be
calculated. This project has a documented, live calibration offset -- 1.570 ppm
over 24 runs, 24 of 24 below zero (docs/atlas-convergence.md). Against a 10 ppm
window that offset sits the probe window off-centre while the target window
stays centred, and the probe rate comes out low for a reason that has nothing
to do with the window being good. So: the offset is MEASURED from the targets
that carry a sequence (observed minus calculated), applied to every probe m/z,
and reported. The `calc_target` arm then bounds what is left -- if the target
rate falls when the same targets are queried on calculated m/z, the correction
is incomplete and the entrapment rate is a floor, not an estimate.

AXES THAT ARE CARRIED BUT NEVER MATCHED. `--rt-tol` has no default and the
shipped PAX8 target list carries `rt_lo` and `rt_hi`, neither of which is an
`rt`. In both cases a probe given "a real target's retention time" is given a
field nothing compares, and the measured quantity collapses to a global
m/z-crowding rate. Every arm therefore reports `axes_matched`, and no arm
claims a coordinate it did not match on.

A DENOMINATOR PADDED WITH IMPOSSIBLE PROBES. A probe whose m/z falls outside
the acquired range, or whose charge no peak carries, contributes one to n and
zero to k and flatters the rate. Those are excluded and counted, not warned
about.

And a rate near zero is not readable on its own: on a thin peak table every
arm reports zero. The Katz ratio of target to probe is what makes the number
mean anything, so it is reported for every arm and propagated as null -- never
as infinite discrimination -- when either numerator is zero.

    python3 tools/proximity_control.py peaks.jsonl targets.csv \\
        --entrapment ecoli_tryptic.txt --precursor-ppm 10 --rt-tol 0.5

Exit 2 if the control could not be run, 1 if it ran and the window did not
discriminate, 0 if it ran and discriminated.
"""

import argparse
import bisect
import json
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from aria_run import (PROTON, field, match_one, axes_participating,  # noqa: E402
                      load_targets, read_peaks, ppm)
from ladder_stats import wilson, katz_ratio  # noqa: E402
# The modification table and its tolerance live in probe_consensus, which chose
# TOL = 0.008 because 0.02 is wide enough to merge deamidation (+0.984016) with
# an A1 isotope pick (+1.003355), only 0.019339 Da apart. Reusing them rather
# than picking a tolerance here, which is how that mistake gets made twice.
from probe_consensus import MODS, TOL  # noqa: E402
# MIN_CONTROL_PROBES lives in validate_export because that is where it is
# enforced, and this module already imports that one for the residue table. A
# second copy of a shared number is how WATER came to be defined twice here.
from validate_export import peptide_mass, MIN_CONTROL_PROBES  # noqa: E402

# Beyond this, observed-minus-calculated is not a calibration offset. A real
# instrument offset is single-digit ppm; this project's documented one is 1.570
# ppm. A target 16 Da from its own sequence is carrying a modification the target
# list does not record, and averaging it into the calibration measures the
# modification. On the shipped PAX8 list 8 of 45 targets are in this state --
# seven oxidised methionines and one deamidation -- and leaving them in moved the
# measured offset from -1.401 to -0.787 ppm, a 44% error against the quantity
# being corrected for.
MAX_CALIBRATION_PPM = 50.0
# A probe whose m/z is this close to its paired target's is testing the target.
SELF_PPM_FACTOR = 1.0


def load_sequences(path):
    """Peptide sequences from a plain list, a CSV/JSONL with a peptide column,
    or a FASTA, which is digested at K and R without cutting before P."""
    text = open(path).read()
    if text.lstrip().startswith(">"):
        return digest_fasta(text)
    rows = None
    try:
        rows = read_peaks(path)
    except Exception:
        rows = None
    seqs = []
    if rows:
        for r in rows:
            s = field(r, "peptide")
            if s:
                seqs.append(str(s).strip().upper())
    if not seqs:
        for line in text.splitlines():
            s = line.strip().upper()
            # Two independent guards, which is deliberate: this one skips the
            # provenance header make_entrapment.py writes, and the mass filter
            # below independently drops anything that is not a peptide. Either
            # alone would do it; a header line reaching the probe list would
            # become a trial that can never match and pad the denominator.
            if s and all(c.isalpha() for c in s):
                seqs.append(s)
    return [s for s in dict.fromkeys(seqs) if s and peptide_mass(s) is not None]


def digest_fasta(text):
    """Fully tryptic peptides, no missed cleavages, 7-30 residues."""
    out, seq = [], []
    def flush(chunk):
        s = "".join(chunk)
        start = 0
        for i, c in enumerate(s):
            if c in "KR" and not (i + 1 < len(s) and s[i + 1] == "P"):
                pep = s[start:i + 1]
                if 7 <= len(pep) <= 30:
                    out.append(pep)
                start = i + 1
        tail = s[start:]
        if 7 <= len(tail) <= 30:
            out.append(tail)
    for line in text.splitlines():
        if line.startswith(">"):
            flush(seq)
            seq = []
        else:
            seq.append(line.strip().upper())
    flush(seq)
    return list(dict.fromkeys(out))


def mz_of(seq, z):
    m = peptide_mass(seq)
    if m is None or not z:
        return None
    return (m + z * PROTON) / z


def modification_like(delta):
    """The name of the modification this neutral-mass difference matches.

    Tested against modifications before isotopes, as probe_consensus does: a
    calculated m/z that omits a peptide's modification produces an offset that
    says nothing about isotope selection.
    """
    for name, mass in MODS:
        if abs(delta - mass) <= TOL:
            return name
    return None


def calibration_offset(targets):
    """Observed minus calculated m/z, in ppm, from targets that carry a sequence.

    This is the quantity that would otherwise bias every probe downward. It is
    measured as a median so one mis-annotated target cannot set it.

    Targets whose observed mass is further from their sequence than any
    calibration could explain are excluded first, because the median is not
    robust enough on its own: on the shipped PAX8 list 8 of 45 carry an
    unrecorded modification, and their presence moved the median by 0.6 ppm --
    44% of the value being measured. A median survives a few wild outliers in
    VALUE but still counts them in RANK, which is the part that bit here.

    Returns (offset_ppm, n_used, spread_ppm, excluded), where excluded is a list
    of (sequence, charge, neutral delta, modification name or None).
    """
    errs, excluded = [], []
    for t in targets:
        seq = t.get("peptide")
        if not seq or not t.get("charge"):
            continue
        seq = str(seq).strip().upper()
        calc = mz_of(seq, t["charge"])
        if calc is None:
            continue
        err = (t["mz"] - calc) / calc * 1e6
        if abs(err) > MAX_CALIBRATION_PPM:
            delta = (t["mz"] - calc) * t["charge"]
            excluded.append((seq, t["charge"], delta, modification_like(delta)))
            continue
        errs.append(err)
    if len(errs) < 3:
        return None, len(errs), None, excluded
    spread = st.median([abs(e - st.median(errs)) for e in errs]) * 1.4826
    return st.median(errs), len(errs), spread, excluded


def calibration_from_matches(hits):
    """Observed PEAK m/z minus sequence-calculated, from this run's own matches.

    The quantity that matters is where this run's peaks sit relative to a
    calculated mass, because that is the gap a calculated probe m/z has to
    cross. Measuring it from the target list instead gives the offset of the
    CONSENSUS, which cannot vary between runs -- on 20 real PAX8 runs it
    returned the identical -1.401 ppm twenty times, labelled as each run's own.
    The real per-run means range from -2.307 to -0.537 ppm with a between-run SD
    of 0.502, so a consensus figure mis-centres the probes by up to a ppm in the
    worst run.

    Returns (offset_ppm, n_used, spread_ppm, excluded) like its fallback.
    """
    errs, excluded = [], []
    for probe, peak in hits:
        seq = probe.get("peptide")
        if not seq or not probe.get("charge"):
            continue
        seq = str(seq).strip().upper()
        calc = mz_of(seq, probe["charge"])
        if calc is None:
            continue
        err = (peak["mz"] - calc) / calc * 1e6
        if abs(err) > MAX_CALIBRATION_PPM:
            delta = (peak["mz"] - calc) * probe["charge"]
            excluded.append((seq, probe["charge"], delta, modification_like(delta)))
            continue
        errs.append(err)
    if len(errs) < 3:
        return None, len(errs), None, excluded
    spread = st.median([abs(e - st.median(errs)) for e in errs]) * 1.4826
    return st.median(errs), len(errs), spread, excluded


def peak_rows(peaks):
    """Peaks reduced to what the window compares, with the scan count kept.

    A peak with no determined charge is dropped here for the same reason the
    real loop leaves it unassigned: without z there is no mass, so a match would
    rest on m/z coincidence alone.
    """
    out = []
    for i, row in enumerate(peaks):
        mz = field(row, "mz")
        zf = field(row, "charge")
        rt = field(row, "rt")
        if mz is None or zf is None or rt is None:
            continue
        dr = field(row, "drift")
        ns = field(row, "n_scans")
        out.append({"i": i, "mz": float(mz), "z": int(zf), "rt": float(rt),
                    "drift": None if dr is None else float(dr),
                    "n_scans": None if ns is None else int(ns)})
    return out


def index_peaks(prows):
    """Peaks bucketed by charge and sorted by m/z, for candidate lookup."""
    by_z = {}
    for pk in prows:
        by_z.setdefault(pk["z"], []).append(pk)
    out = {}
    for z, rows in by_z.items():
        rows.sort(key=lambda p: p["mz"])
        out[z] = (rows, [p["mz"] for p in rows])
    return out


def candidates(index, probe, ppm_tol):
    """Every peak whose m/z could satisfy the window, and no fewer.

    A CANDIDATE filter, not a second window: match_one still decides, including
    charge, retention time and drift. The bound is exact rather than
    approximate, which is the only reason this is safe -- match_one's m/z test
    divides by the TARGET's m/z, and in the probe direction the probe IS the
    target, so a peak matches only if

        |mz_peak - mz_probe| <= ppm_tol * mz_probe / 1e6

    and every such peak is returned. A loose filter here would silently narrow
    the window being measured, which is the one thing this control cannot do; a
    test asserts this path returns exactly what a brute-force scan returns.

    Without it the control is O(probes x peaks): 2,000 probes against 50,000
    peaks took 44 seconds per arm when measured, and a real export carries far
    more peaks than that.
    """
    z = probe.get("charge")
    if z:
        buckets = [index[z]] if z in index else []
    else:
        # match_one skips the charge test when the target carries no charge, so
        # such a probe can match a peak at any charge.
        buckets = list(index.values())
    half = ppm_tol * probe["mz"] / 1e6
    lo, hi = probe["mz"] - half, probe["mz"] + half
    out = []
    for rows, mzs in buckets:
        i = bisect.bisect_left(mzs, lo)
        while i < len(mzs) and mzs[i] <= hi:
            out.append(rows[i])
            i += 1
    return out


def probe_direction(probes, prows, ppm_tol, rt_tol, drift_tol, use_drift,
                    min_scans, index=None):
    """For each probe, the peak in its window nearest in ppm.

    Calls match_one -- the same function the annotation loop calls -- with the
    probe as a one-element target list. The window is therefore identical by
    construction rather than by inspection.

    Nearest in ppm, not first seen: that is the annotation loop's own tie-break,
    and it makes the answer independent of the order peaks arrive in. Taking the
    first made `distinct_peaks_matched` depend on input order once the candidate
    index reordered them.
    """
    if index is None:
        index = index_peaks(prows)
    k = 0
    matched_peaks = set()
    hits = []
    for p in probes:
        best = None
        for pk in candidates(index, p, ppm_tol):
            if pk["n_scans"] is not None and pk["n_scans"] < min_scans:
                continue          # would have been rejected as a real ion too
            if not match_one(pk["mz"], pk["z"], pk["rt"], pk["drift"], [p],
                             ppm_tol, rt_tol, drift_tol, use_drift):
                continue
            d = ppm(pk["mz"], p["mz"])
            if best is None or d < best[0]:
                best = (d, pk)
        if best is not None:
            k += 1
            matched_peaks.add(best[1]["i"])
            hits.append((p, best[1]))
    return k, len(matched_peaks), hits


def build_entrapment(seqs, targets, offset_ppm, mz_lo, mz_hi, ppm_tol, charges):
    """Each entrapment peptide at its own m/z, on a real target's coordinates.

    The charge is taken from the paired target rather than drawn, so the probe
    charge distribution matches the target charge distribution exactly and a
    mismatch cannot depress the rate. Probes that could not match for a reason
    other than the window are excluded and counted.
    """
    probes, excluded = [], {"out_of_range": 0, "unseen_charge": 0, "self_match": 0,
                            "degenerate_with_a_target": 0, "unmassable": 0}
    usable = [t for t in targets if t.get("charge")]
    if not usable:
        return probes, excluded
    tmz = [(t["mz"], t["charge"]) for t in usable]
    for n, seq in enumerate(seqs):
        t = usable[n % len(usable)]
        z = t["charge"]
        mz = mz_of(seq, z)
        if mz is None:
            excluded["unmassable"] += 1
            continue
        if offset_ppm is not None:
            mz *= (1.0 + offset_ppm / 1e6)
        if not (mz_lo <= mz <= mz_hi):
            excluded["out_of_range"] += 1
            continue
        if z not in charges:
            excluded["unseen_charge"] += 1
            continue
        if ppm(mz, t["mz"]) <= ppm_tol * SELF_PPM_FACTOR:
            excluded["self_match"] += 1
            continue
        if any(zz == z and ppm(mz, m) <= ppm_tol for m, zz in tmz):
            excluded["degenerate_with_a_target"] += 1
            continue
        probes.append({"mz": mz, "charge": z, "rt": t.get("rt"),
                       "drift": t.get("drift"), "peptide": seq,
                       "paired_with": t.get("peptide") or ("%.4f" % t["mz"])})
    return probes, excluded


def build_calc_target(targets, mz_lo, mz_hi, charges, skip=()):
    """The real targets, m/z recomputed from their own sequences.

    Everything else about these probes is the real target. If this arm's rate is
    materially below the observed-m/z target rate, the gap is the penalty for
    using a calculated m/z, and it bounds how much of the entrapment arm's
    shortfall is calibration rather than window.
    """
    probes = []
    excluded = {"no_sequence": 0, "out_of_range": 0, "unseen_charge": 0,
                "unrecorded_modification": 0}
    skip = set(skip)
    for t in targets:
        seq = t.get("peptide")
        if not seq or not t.get("charge"):
            excluded["no_sequence"] += 1
            continue
        seq = str(seq).strip().upper()
        if (seq, t["charge"]) in skip:
            # Its sequence-derived mass is wrong for a composition reason, so it
            # measures a missing modification, not the price of a calculated
            # coordinate. Left in, these drag this arm down and the report calls
            # the entrapment rate a floor for the wrong reason.
            excluded["unrecorded_modification"] += 1
            continue
        mz = mz_of(seq, t["charge"])
        if mz is None:
            excluded["no_sequence"] += 1
            continue
        if not (mz_lo <= mz <= mz_hi):
            excluded["out_of_range"] += 1
            continue
        if t["charge"] not in charges:
            excluded["unseen_charge"] += 1
            continue
        probes.append({"mz": mz, "charge": t["charge"], "rt": t.get("rt"),
                       "drift": t.get("drift"), "peptide": seq})
    return probes, excluded


def build_swap(targets, rt_tol, drift_tol, use_drift):
    """Target A's m/z on target B's coordinates, B far enough from A to be a
    different address. Answers whether rt and drift do any work at all."""
    probes, excluded = [], {"no_distant_partner": 0, "no_coordinate": 0}
    usable = [t for t in targets if t.get("charge")]
    for a in usable:
        partner = None
        for b in usable:
            if b is a:
                continue
            if rt_tol is not None and a.get("rt") is not None and b.get("rt") is not None:
                if abs(a["rt"] - b["rt"]) > 2 * rt_tol:
                    partner = b
                    break
            elif use_drift and a.get("drift") is not None and b.get("drift") is not None:
                if abs(a["drift"] - b["drift"]) > 2 * drift_tol:
                    partner = b
                    break
        if partner is None:
            excluded["no_distant_partner"] += 1
            continue
        probes.append({"mz": a["mz"], "charge": a["charge"], "rt": partner.get("rt"),
                       "drift": partner.get("drift"),
                       "peptide": a.get("peptide"), "paired_with": partner.get("peptide")})
    return probes, excluded


def strata(probes, hits, n_strata=4):
    """Per-m/z-quartile rates, so an ion can carry a number that is its own.

    A single cohort rate stamped onto every ion is the same constant repeated,
    which is decoration and invites exactly the misreading the validator
    forbids. m/z crowding is not uniform across the range, so the quartile an
    ion sits in is a real distinction.
    """
    if not probes:
        return []
    matched = {id(p) for p, _ in hits}
    ordered = sorted(probes, key=lambda p: p["mz"])
    out = []
    size = max(1, len(ordered) // n_strata)
    for s in range(n_strata):
        chunk = ordered[s * size:(s + 1) * size] if s < n_strata - 1 else ordered[s * size:]
        if not chunk:
            continue
        k = sum(1 for p in chunk if id(p) in matched)
        lo, hi = wilson(k, len(chunk))
        out.append({"stratum": s, "mz_lo": chunk[0]["mz"], "mz_hi": chunk[-1]["mz"],
                    "queried": len(chunk), "matched": k,
                    "rate": k / len(chunk), "ci": [lo, hi]})
    return out


def expected_collision_rate(prows, probes, ppm_tol, rt_tol, axes):
    """How often a probe should find a peak by chance alone, from geometry.

    This replaces comparing the measured rate against the 33-70% the project
    recorded elsewhere. That figure came from reversed-tryptic decoys, and
    reversal permutes a peptide's residues and so preserves its mass exactly --
    a reversed decoy matches an m/z window 100% by construction, so whatever the
    33-70% measured it was not this window. Disagreeing with it would prove
    nothing. This expectation is derived from THIS run instead.

    Each peak covers a box: 2 * ppm_tol wide in m/z, and 2 * rt_tol in retention
    time when retention time is actually compared. Summing those boxes over the
    peaks at a probe's charge, as a fraction of the acquired volume, is the
    chance a coordinate drawn at random from that volume lands on one.

    It is reported as CONTEXT, not enforced as a gate. The estimate assumes the
    peak windows are spread uniformly over the acquired volume, and real peaks
    are not: they cluster. A measured rate below this expectation is therefore
    just as likely to mean the peaks cluster away from the probes as that the
    probes are suppressed, and no fixture distinguishes the two without building
    in the clustering by hand. What it does tell a reader honestly is whether
    clearing the bar means anything at all: when chance alone is 0.001%, a
    measured 0% is unremarkable and the target arm beside it carries the weight.

    Two further caveats in opposite directions: real peptide masses cluster on
    the mass-defect lattice, so a peptide-plausible probe should land MORE often
    than a uniformly random m/z; and overlapping boxes are double-counted here,
    which inflates it.
    """
    if not prows or not probes:
        return None
    mz_lo = min(p["mz"] for p in prows)
    mz_hi = max(p["mz"] for p in prows)
    mz_span = mz_hi - mz_lo
    if mz_span <= 0:
        return None
    use_rt = "rt" in axes and rt_tol
    rt_lo = min(p["rt"] for p in prows)
    rt_hi = max(p["rt"] for p in prows)
    rt_span = rt_hi - rt_lo
    if use_rt and rt_span <= 0:
        use_rt = False
    by_z = {}
    for pk in prows:
        frac = 2.0 * ppm_tol * pk["mz"] / 1e6 / mz_span
        if use_rt:
            frac *= min(1.0, 2.0 * rt_tol / rt_span)
        by_z[pk["z"]] = by_z.get(pk["z"], 0.0) + frac
    tot = 0.0
    for pr in probes:
        tot += min(1.0, by_z.get(pr["charge"], 0.0))
    return tot / len(probes)


def control_for_ion(rep, mz):
    """The control block an ion at this m/z carries.

    This ion's own m/z stratum, not the cohort rate. The cohort rate is one
    number for the whole run; writing it onto every ion would repeat a constant
    thousands of times and invite reading a per-ion field as per-ion evidence,
    which is exactly what the tier is not. An ion whose m/z falls outside every
    probe stratum takes the nearest, and says so: its mz_lo and mz_hi will not
    contain its own m/z.
    """
    e = rep["arms"]["entrapment"]
    block = {"arm": "entrapment", "axes_matched": list(rep["axes_matched"])}
    st = e.get("strata") or []
    pick = None
    for cand in st:
        if cand["mz_lo"] <= mz <= cand["mz_hi"]:
            pick = cand
            break
    if pick is None and st:
        pick = min(st, key=lambda c: min(abs(mz - c["mz_lo"]), abs(mz - c["mz_hi"])))
    if pick is None:
        block.update({"queried": 0, "matched": 0, "stratum_false_rate": None,
                      "stratum_ci": None, "mz_lo": None, "mz_hi": None})
        return block
    block.update({"queried": pick["queried"], "matched": pick["matched"],
                  "stratum_false_rate": pick["rate"],
                  "stratum_ci": [pick["ci"][0], pick["ci"][1]],
                  "mz_lo": pick["mz_lo"], "mz_hi": pick["mz_hi"]})
    return block


def run_control(peaks, targets, seqs, ppm_tol, rt_tol, drift_tol, use_drift, min_scans):
    """Every arm, with the counts each rate was computed from."""
    prows = peak_rows(peaks)
    if not prows:
        return None, "no peak carried an m/z, a charge and a retention time, so" \
                     " nothing could be queried"
    if not targets:
        return None, "no targets, so there is no window to measure"

    mz_lo = min(p["mz"] for p in prows)
    mz_hi = max(p["mz"] for p in prows)
    charges = {p["z"] for p in prows}
    axes = axes_participating(targets, rt_tol, drift_tol, use_drift)

    rep = {"axes_matched": axes, "n_peaks_queryable": len(prows),
           "mz_range": [mz_lo, mz_hi], "charges_seen": sorted(charges),
           "min_scans": min_scans, "precursor_ppm": ppm_tol,
           "rt_tol": rt_tol, "drift_tol": drift_tol if use_drift else None,
           "arms": {}}

    index = index_peaks(prows)

    def arm(name, probes, excluded, note=None):
        if not probes:
            rep["arms"][name] = {"arm": name, "queried": 0, "matched": 0,
                                 "rate": None, "ci": None, "excluded": excluded,
                                 "axes_matched": axes, "note": note or
                                 "no probe survived construction"}
            return
        k, ndist, hits = probe_direction(probes, prows, ppm_tol, rt_tol,
                                         drift_tol, use_drift, min_scans, index)
        lo, hi = wilson(k, len(probes))
        rep["arms"][name] = {"arm": name, "queried": len(probes), "matched": k,
                             "rate": k / len(probes), "ci": [lo, hi],
                             "distinct_peaks_matched": ndist,
                             "excluded": excluded, "axes_matched": axes,
                             "strata": strata(probes, hits), "note": note,
                             # Dropped before the report is written; the
                             # calibration needs the matched pairs.
                             "_hits": hits}

    arm("target", [dict(t) for t in targets if t.get("charge")],
        {"no_charge": sum(1 for t in targets if not t.get("charge"))})

    # The calibration comes from the target arm's own matches, so it is this
    # run's offset and not the consensus list's. Falling back to the target list
    # only when nothing matched, which is a state the report names.
    # Two separate facts, and both are wanted.
    #
    # Which targets contradict their own sequence is a property of the LIST, so
    # it is audited against the list and reported whether or not those targets
    # matched anything here. Deriving it from the matches instead would hide
    # every bad target that found no peak -- which it did, until three cases
    # written for the audit failed.
    _, _, _, mod_excluded = calibration_offset(targets)
    # The offset actually applied to the probes is a property of THIS RUN, so it
    # comes from this run's own matched peaks.
    t_hits = rep["arms"]["target"].get("_hits") or []
    offset, n_off, spread, _ = calibration_from_matches(t_hits)
    cal_source = "this run's matched peaks"
    if offset is None:
        offset, n_off, spread, _ = calibration_offset(targets)
        cal_source = ("the target list's own coordinates, because fewer than 3"
                      " targets with a sequence matched a peak in this run --"
                      " a consensus figure, identical for every run")
    rep["calibration"] = {"offset_ppm": offset, "n_targets_used": n_off,
                          "spread_ppm": spread, "measured_from": cal_source,
                          "excluded_unrecorded_modification": [
                              {"peptide": seq, "charge": z, "delta_da": d,
                               "looks_like": name}
                              for seq, z, d, name in mod_excluded]}

    ep, eex = build_entrapment(seqs or [], targets, offset, mz_lo, mz_hi, ppm_tol, charges)
    arm("entrapment", ep, eex,
        None if seqs else "no entrapment sequences were supplied")
    cp, cex = build_calc_target(targets, mz_lo, mz_hi, charges,
                                skip=[(seq, z) for seq, z, _, _ in mod_excluded])
    arm("calc_target", cp, cex)
    # The same targets on their OBSERVED m/z, so the calculated-coordinate
    # penalty is a difference over one set of probes. The target arm's own rate
    # cannot serve: it is computed over every target with a charge, while this
    # arm drops those with no sequence, outside the range, or contradicting
    # their own mass -- on 20 real PAX8 runs that was 542 probes against 900,
    # and comparing the two rates compared different populations.
    if cp:
        obs = []
        for pr in cp:
            t = next((t for t in targets
                      if (t.get("peptide") or "").strip().upper() == pr["peptide"]
                      and t.get("charge") == pr["charge"]), None)
            if t:
                obs.append({"mz": t["mz"], "charge": t["charge"], "rt": t.get("rt"),
                            "drift": t.get("drift"), "peptide": pr["peptide"]})
        if obs:
            k_o, _, _ = probe_direction(obs, prows, ppm_tol, rt_tol, drift_tol,
                                        use_drift, min_scans, index)
            lo_o, hi_o = wilson(k_o, len(obs))
            rep["arms"]["calc_target"]["paired_observed"] = {
                "queried": len(obs), "matched": k_o, "rate": k_o / len(obs),
                "ci": [lo_o, hi_o]}
    if "rt" in axes or "drift" in axes:
        sp, sex = build_swap(targets, rt_tol, drift_tol, use_drift)
        arm("swap", sp, sex)
    else:
        rep["arms"]["swap"] = {"arm": "swap", "queried": 0, "matched": 0,
                               "rate": None, "ci": None, "axes_matched": axes,
                               "note": "neither rt nor drift is compared, so there is"
                                       " no coordinate to swap and no second axis"
                                       " whose contribution could be measured"}

    t = rep["arms"]["target"]
    for name, a in rep["arms"].items():
        if name == "target" or not a.get("queried"):
            a["ratio_vs_target"] = None
            continue
        a["ratio_vs_target"] = katz_ratio(t["matched"], t["queried"],
                                          a["matched"], a["queried"])
    ent_probes = ep if ep else []
    rep["expected_by_chance"] = expected_collision_rate(
        prows, ent_probes, ppm_tol, rt_tol, axes)
    rep["cohort_false_rate"] = rep["arms"]["entrapment"]["rate"]
    rep["findings"] = findings(rep)
    for a in rep["arms"].values():
        a.pop("_hits", None)
    return rep, None


def findings(rep):
    """Everything about this report that argues against trusting it."""
    out = []
    e = rep["arms"]["entrapment"]
    t = rep["arms"]["target"]
    c = rep["arms"]["calc_target"]

    if "rt" not in rep["axes_matched"] and "drift" not in rep["axes_matched"]:
        out.append("Only m/z and charge were compared, so every rate here is a GLOBAL"
                   " m/z-crowding rate. A probe given a target's retention time was"
                   " given a field nothing looked at. Pass --rt-tol, and give the"
                   " targets an `rt` column, before reading any of this as a rate"
                   " at an occupied address.")
    if rep["calibration"]["offset_ppm"] is None:
        out.append("No target carried a sequence, so the observed-minus-calculated"
                   " offset could not be measured and the entrapment m/z is"
                   " UNCORRECTED. This project has a documented 1.57 ppm offset;"
                   " against a %.0f ppm window that makes the entrapment rate a"
                   " floor, not an estimate." % rep["precursor_ppm"])
    elif abs(rep["calibration"]["offset_ppm"]) > 0.25 * rep["precursor_ppm"]:
        out.append("The observed-minus-calculated offset is %.2f ppm, which is %.0f%%"
                   " of the %.0f ppm window. It has been applied to every probe, but"
                   " an offset this large relative to the window means the correction"
                   " itself matters to the answer."
                   % (rep["calibration"]["offset_ppm"],
                      100 * abs(rep["calibration"]["offset_ppm"]) / rep["precursor_ppm"],
                      rep["precursor_ppm"]))
    mx = rep["calibration"].get("excluded_unrecorded_modification") or []
    n_seq = (rep["calibration"]["n_targets_used"] or 0) + len(mx)
    if mx and n_seq and len(mx) / n_seq > 0.10:
        out.append("%d of %d target(s) with a sequence carry a modification the target"
                   " list does not record (%s). They are excluded from the calibration"
                   " and from the calc_target arm, but a list this incomplete about its"
                   " own modifications is worth fixing at the source: every one of them"
                   " is an ion whose stated sequence does not match its stated mass."
                   % (len(mx), n_seq,
                      ", ".join(sorted({e["looks_like"] or "unexplained" for e in mx}))))
    # Compared against the SAME probes on their observed m/z, never against the
    # target arm, whose population differs.
    po = (c or {}).get("paired_observed") or {}
    if c.get("queried") and po.get("queried") and po.get("rate"):
        if c["rate"] < 0.75 * po["rate"]:
            out.append("The same %d target(s) match at %.1f%% on CALCULATED m/z against"
                       " %.1f%% on their observed m/z. That gap is the price of a"
                       " calculated coordinate, not a property of the window, and the"
                       " entrapment arm pays it too: its %s is a floor."
                       % (po["queried"], 100 * c["rate"], 100 * po["rate"],
                          "rate" if e["rate"] is None else "%.1f%%" % (100 * e["rate"])))
    if e.get("queried") and e["queried"] < MIN_CONTROL_PROBES:
        lo, hi = e["ci"]
        out.append("Only %d entrapment probes survived construction. The 95%% interval"
                   " is %.1f%% to %.1f%%, which is too wide to certify anything."
                   % (e["queried"], 100 * lo, 100 * hi))
    exp = rep.get("expected_by_chance")
    if e.get("queried") and exp is not None and exp < 0.01:
        out.append("Chance alone would have put only %.3g%% of these probes on a peak, from"
                   " this run's peak-window geometry. A low measured rate therefore says"
                   " very little: on a peak table this sparse almost any probe set clears"
                   " the bar. The number to read is the target arm beside it, not this one."
                   % (100 * exp))
    if e.get("queried") and e.get("distinct_peaks_matched") is not None \
            and e["matched"] > max(1, e["distinct_peaks_matched"]):
        out.append("%d entrapment probes matched only %d distinct peaks, so the trials"
                   " are clustered and the Wilson interval is narrower than the data"
                   " supports." % (e["matched"], e["distinct_peaks_matched"]))
    if t.get("queried") and not t.get("matched"):
        out.append("No real target matched either, so a low probe rate says the peak"
                   " table is thin, not that the window is good.")
    return out


def discriminated(rep):
    """Did the window tell targets from probes? Two-sided on purpose.

    A rule of the form "fail when the rate exceeds a cut" passes most easily
    when the control is broken, because every way of breaking it drives the rate
    toward zero. So a near-zero rate is only accepted alongside a target rate it
    is actually below, enough probes to support the bound, and no finding that
    says the number is a floor.
    """
    e = rep["arms"]["entrapment"]
    t = rep["arms"]["target"]
    if not e.get("queried") or not t.get("queried"):
        return False, "an arm had no probes"
    if e["queried"] < MIN_CONTROL_PROBES:
        return False, "%d probes is too few for the bound to mean anything" % e["queried"]
    if not t["matched"]:
        return False, "no real target matched, so there is nothing to discriminate from"
    if e["ci"][1] >= t["ci"][0]:
        return False, ("the entrapment interval (%.1f%%-%.1f%%) overlaps the target"
                       " interval (%.1f%%-%.1f%%): the window did not separate"
                       " coordinates that can be right from ones that cannot"
                       % (100 * e["ci"][0], 100 * e["ci"][1],
                          100 * t["ci"][0], 100 * t["ci"][1]))
    return True, "entrapment %.1f%% against target %.1f%%, intervals disjoint" % (
        100 * e["rate"], 100 * t["rate"])


def render(rep):
    out = ["PROXIMITY CONTROL  the false-positive rate of the bare_proximity window"]
    out.append("  axes actually compared: %s" % ", ".join(rep["axes_matched"]))
    out.append("  %d queryable peaks, m/z %.3f-%.3f, charges %s, >=%d scans"
               % (rep["n_peaks_queryable"], rep["mz_range"][0], rep["mz_range"][1],
                  "/".join(str(c) for c in rep["charges_seen"]), rep["min_scans"]))
    cal = rep["calibration"]
    mx = cal.get("excluded_unrecorded_modification") or []
    if mx:
        named = {}
        for e in mx:
            named[e["looks_like"] or "unexplained"] = named.get(
                e["looks_like"] or "unexplained", 0) + 1
        out.append("  %d target(s) excluded from the calibration: observed mass further"
                   " from the sequence" % len(mx))
        out.append("  than any calibration explains, so the target list does not record"
                   " their")
        out.append("  modification -- %s. Averaging those in measures the modification."
                   % ", ".join("%d %s" % (n, k) for k, n in sorted(named.items())))
    if cal["offset_ppm"] is None:
        out.append("  calibration: NOT MEASURED (%d target(s) carried a usable sequence)"
                   % cal["n_targets_used"])
    else:
        out.append("  calibration: observed - calculated = %+.2f ppm (n=%d, spread %.2f),"
                   " applied to every probe,"
                   % (cal["offset_ppm"], cal["n_targets_used"], cal["spread_ppm"]))
        out.append("  measured from %s" % cal.get("measured_from", "unstated"))
    out.append("")
    out.append("  %-12s %7s %7s %8s %18s %18s" %
               ("arm", "probes", "hit", "rate", "95% CI", "target:probe"))
    for name in ("target", "entrapment", "calc_target", "swap"):
        a = rep["arms"].get(name)
        if not a:
            continue
        if not a.get("queried"):
            out.append("  %-12s %7s %7s %8s   %s" % (name, "-", "-", "-",
                                                     a.get("note") or "not run"))
            continue
        r = a["ratio_vs_target"]
        ci = "%.1f%%-%.1f%%" % (100 * a["ci"][0], 100 * a["ci"][1])
        out.append("  %-12s %7d %7d %7.1f%% %18s %18s"
                   % (name, a["queried"], a["matched"], 100 * a["rate"], ci,
                      "-" if r is None else "%.1f-%.1f" % r))
    exp = rep.get("expected_by_chance")
    if exp is not None:
        out.append("")
        if exp <= 0:
            shown = "0 (no peak window overlaps any probe's charge)"
        elif exp < 0.001:
            shown = "%.4g%% -- about 1 probe in %s" % (100 * exp, "{:,}".format(int(round(1 / exp))))
        else:
            shown = "%.2f%%" % (100 * exp)
        out.append("  Chance alone, from peak-window geometry on this run: %s." % shown)
        out.append("  Context, not a gate: it assumes the peak windows are spread evenly")
        out.append("  and real peaks cluster, so a rate below it is as likely to mean the")
        out.append("  peaks sit away from the probes as that the probes are suppressed.")
        out.append("  What it does say is whether clearing the bar means anything: when")
        out.append("  chance alone is this small, the target arm beside it carries the")
        out.append("  weight, not the probe rate.")
    out.append("")
    out.append("  A ratio of `-` means one arm matched nothing, which has no upper")
    out.append("  bound and is reported as null rather than as perfect separation.")
    ok, why = discriminated(rep)
    out.append("")
    out.append("  DISCRIMINATED: %s -- %s" % ("yes" if ok else "NO", why))
    if rep["findings"]:
        out.append("")
        out.append("  What argues against trusting these numbers:")
        for f in rep["findings"]:
            out.append("    - " + f)
    return "\n".join(out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("peaks")
    p.add_argument("targets")
    p.add_argument("--entrapment", help="peptides that cannot be in this sample:"
                                        " a FASTA, a sequence list, or a table with"
                                        " a peptide column")
    p.add_argument("--precursor-ppm", type=float, default=10.0)
    p.add_argument("--rt-tol", type=float)
    p.add_argument("--drift-tol", type=float)
    p.add_argument("--min-scans", type=int, default=2)
    p.add_argument("--json", help="write the full report here")
    args = p.parse_args(argv)

    peaks = read_peaks(args.peaks)
    targets = load_targets(args.targets)
    seqs = load_sequences(args.entrapment) if args.entrapment else []
    if args.entrapment and not seqs:
        print("No usable peptide sequence in %s. The entrapment arm is the"
              " measurement; without it there is nothing to report."
              % args.entrapment, file=sys.stderr)
        return 2
    use_drift = args.drift_tol is not None

    rep, err = run_control(peaks, targets, seqs, args.precursor_ppm, args.rt_tol,
                           args.drift_tol, use_drift, args.min_scans)
    if rep is None:
        print("Control could not be run: %s" % err, file=sys.stderr)
        return 2
    if not seqs:
        print("Control could not be run: no --entrapment list, so the rate the"
              " bare_proximity tier needs was not measured. The other arms are"
              " printed because they are diagnostics, not the measurement.",
              file=sys.stderr)
        print(render(rep))
        return 2
    print(render(rep))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rep, fh, indent=2, sort_keys=True)
    ok, _ = discriminated(rep)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
