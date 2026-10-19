#!/usr/bin/env python3
"""Is saturation a property of the scan, or a property of the ion?

Two mechanisms depress a measured intensity, and they are usually conflated
because a vendor emits one flag for both:

  SOURCE  droplet/Rayleigh competition. A droplet of given volume carries a
          fixed maximum charge, so when the total load exceeds it the ions
          compete for surface access and EVERYTHING in that scan comes out
          low -- the dim ones too, because surface access is not granted in
          order of abundance.

  DETECTOR TDC dead time. The detector cannot register two arrivals inside
          its dead time, so a BRIGHT channel undercounts and a dim channel
          in the same scan is untouched.

They have opposite signatures, and the distinction decides whether a
per-scan correction is even the right shape. A scan-level factor is
estimable: many ions constrain one number, so the residuals can say the
model is wrong. A per-ion correction borrowed from per-ion neighbours is one
unknown per ion, which is not a measurement.

The discriminating test is on the dim ions, each used as its own control:

  take ions that are NOT flagged and NOT bright, and compare each one's
  intensity in scans that carry a flagged detection against its intensity in
  scans that do not.

  depressed -> SOURCE. The scan-level model holds.
  flat      -> DETECTOR. The flag is dead time and a per-scan factor is wrong.

Pairing within an ion matters: comparing different ions across the two
groups measures which ions happen to elute near bright ones, which is a
selection effect and not a suppression.

  python3 tools/saturation_check.py run1.csv run2.csv ...
  python3 tools/saturation_check.py 'apex/*/Summary/*_ApexRT_temp1.csv'

Input is ApexRT-style CSV: ScanIndex, Mz, Counts, Intensity, IsSaturated,
PeakGroupIndex. Exit 0 if the tests ran and found nothing to act on, 1 if they
ran and produced a finding -- a suppression effect of either kind, or too few
usable ions to decide, which is equally a reason not to apply a correction --
and 2 if the input cannot support the tests at all.
"""

import argparse
import collections
import csv
import glob
import math
import os
import sys

# A flag that a covariate already predicts perfectly is not an independent
# measurement of anything. Real detector saturation depends on arrival
# distribution and peak shape, so two ions at one intensity CAN differ --
# which means some overlap in intensity between flagged and unflagged is the
# signature of a measured flag rather than a computed one.
#
# Measured as a COUNT in the ambiguous band, not a share of all rows: in real
# data saturation is a small minority, and a share-of-everything criterion
# would refuse exactly the datasets where saturation is rare and interesting.
# What matters is only whether two rows at one intensity ever disagree.
MIN_BAND_ROWS = 3

# Below this share of peak groups showing any intensity variation across
# their own scans there is no peak shape to measure, so neither coherence nor
# the paired dim-ion test has anything to read.
MIN_VARYING_GROUPS = 0.50

# "Dim" is defined against the flag threshold, not an absolute count, so the
# test does not depend on the instrument's scale.
DIM_FRACTION = 0.25

# A log-parabola has three parameters, so leaving one scan out still needs
# four behind it, and both a marked and an unmarked scan must remain.
MIN_SCANS_FOR_LOO = 6

# Two tests are reported, so 0.05 is Bonferroni-corrected across them rather
# than chosen by eye.
ALPHA = 0.05 / 2

REQUIRED = ("ScanIndex", "Intensity", "IsSaturated", "PeakGroupIndex")


class Unusable(Exception):
    """The input cannot support the tests. Not a finding -- a refusal."""


def read_rows(paths):
    rows = []
    for p in paths:
        run = os.path.basename(os.path.dirname(os.path.dirname(p))) or os.path.basename(p)
        with open(p, newline="") as fh:
            rd = csv.DictReader(fh)
            missing = [c for c in REQUIRED if c not in (rd.fieldnames or [])]
            if missing:
                raise Unusable("%s lacks column(s): %s" % (p, ", ".join(missing)))
            for r in rd:
                try:
                    rows.append({
                        "run": run,
                        "scan": int(float(r["ScanIndex"])),
                        "inten": float(r["Intensity"]),
                        "sat": str(r["IsSaturated"]).strip() in ("1", "True", "true", "TRUE"),
                        "group": r["PeakGroupIndex"],
                        "counts": float(r["Counts"]) if r.get("Counts") not in (None, "") else None,
                        "mz": float(r["Mz"]) if r.get("Mz") not in (None, "") else None,
                    })
                except (TypeError, ValueError):
                    continue
    if not rows:
        raise Unusable("no parseable rows")
    return rows


def refuse_if_unusable(rows, out):
    """Each of these is a way the input answers the question by construction."""
    out("INPUT   %d row(s), %d run(s), %d peak group(s)"
        % (len(rows), len(set(r["run"] for r in rows)),
           len(set((r["run"], r["group"]) for r in rows))))

    nsat = sum(1 for r in rows if r["sat"])
    out("        %d flagged saturated (%.0f%%)" % (nsat, 100.0 * nsat / len(rows)))
    if nsat == 0:
        raise Unusable("no row is flagged saturated, so there is nothing to "
                       "condition on. This is not a pass: it may mean the flag "
                       "is not populated by this export.")
    if nsat == len(rows):
        raise Unusable("every row is flagged saturated, so there is no "
                       "unflagged group to compare against")

    # 1. Counts has to carry count-like information or the Poisson band, which
    #    is the only parameter-free expectation available, is undefined.
    cs = [r["counts"] for r in rows if r["counts"] is not None]
    if cs and len(set(cs)) == 1:
        ints = set(r["inten"] for r in rows)
        raise Unusable(
            "Counts is the single value %g in all %d rows while Intensity "
            "spans %d distinct values over %.0fx. Counts does not vary with "
            "signal, so it is not a count and 1/sqrt(N) has no basis here."
            % (cs[0], len(cs), len(ints), max(ints) / min(ints)))

    # 2. No peak shape, no coherence and no within-ion comparison.
    bygroup = collections.defaultdict(set)
    for r in rows:
        bygroup[(r["run"], r["group"])].add(r["inten"])
    sizes = collections.Counter((r["run"], r["group"]) for r in rows)
    multi = [k for k in bygroup if sizes[k] > 1]
    varying = [k for k in multi if len(bygroup[k]) > 1]
    if multi:
        share = float(len(varying)) / len(multi)
        out("        %d of %d multi-scan group(s) vary in intensity across "
            "their own scans (%.0f%%)" % (len(varying), len(multi), 100.0 * share))
        if share < MIN_VARYING_GROUPS:
            raise Unusable(
                "%.0f%% of multi-scan peak groups hold ONE intensity across "
                "every scan they span. A real elution profile rises and falls; "
                "a flat one has no shape, so the paired dim-ion test would be "
                "comparing a number against itself."
                % (100.0 * (1 - share)))
    else:
        raise Unusable("no peak group spans more than one scan, so no ion can "
                       "serve as its own control")

    # 3. A flag a threshold reproduces exactly is that threshold wearing a
    #    different name. Conditioning on it conditions on intensity.
    sat_i = sorted(r["inten"] for r in rows if r["sat"])
    uns_i = sorted(r["inten"] for r in rows if not r["sat"])
    lo, hi = sat_i[0], uns_i[-1]
    if hi <= lo:
        raise Unusable(
            "IsSaturated is perfectly separated by Intensity: every flagged "
            "row is above %.0f and every unflagged row below %.0f, with no "
            "overlap. The flag is a threshold on intensity, so conditioning "
            "on it only conditions on brightness and cannot distinguish a "
            "source effect from a detector one." % (lo, hi))
    # Inside [lo, hi] both classes occur, so a row's flag is not determined
    # by its intensity there. Both directions have to be populated: flagged
    # rows that are dim for flagged rows, and unflagged rows that are bright.
    f_in = [r for r in rows if lo <= r["inten"] <= hi and r["sat"]]
    u_in = [r for r in rows if lo <= r["inten"] <= hi and not r["sat"]]
    out("        flag/intensity overlap band %.0f-%.0f holds %d flagged and "
        "%d unflagged row(s)" % (lo, hi, len(f_in), len(u_in)))
    if min(len(f_in), len(u_in)) < MIN_BAND_ROWS:
        raise Unusable(
            "the band where flagged and unflagged intensities overlap holds "
            "only %d flagged and %d unflagged row(s). Below %d of each the "
            "flag is indistinguishable from a pure intensity threshold, so "
            "conditioning on it would condition on brightness alone."
            % (len(f_in), len(u_in), MIN_BAND_ROWS))
    return lo, hi


def load_marker(rows, out, top=0.25):
    """Scans in the top `top` of total ion current, which is what droplet
    charge competition is driven by.

    Caveat worth keeping in view: TIC is the sum of what SURVIVED the
    competition, so at the charge limit it plateaus and stops indexing load
    exactly where load matters most. It orders scans correctly well below
    saturation and compresses them above it, which biases this test toward
    finding nothing.
    """
    tot = collections.defaultdict(float)
    for r in rows:
        tot[(r["run"], r["scan"])] += r["inten"]
    vals = sorted(tot.values(), reverse=True)
    cut = vals[max(0, int(len(vals) * top) - 1)]
    out("        per-scan total ion current: median %.0f, top %.0f%% at or "
        "above %.0f (%d of %d scan(s))"
        % (sorted(tot.values())[len(tot) // 2], 100.0 * top, cut,
           sum(1 for v in tot.values() if v >= cut), len(tot)))
    out("        the ion under test is removed from its own scan's current, so")
    out("        its own noise cannot help select the scans it is judged in")

    def marked_for(rs):
        # Leave the ion out of the load it is being judged against. Including
        # it means a scan is marked partly BECAUSE this ion fluctuated high
        # there, which by itself puts observed above predicted.
        m = set()
        for r in rs:
            key = (r["run"], r["scan"])
            if tot[key] - r["inten"] >= cut:
                m.add(key)
        return m

    return marked_for


def scan_has_flag(rows):
    hot = set()
    for r in rows:
        if r["sat"]:
            hot.add((r["run"], r["scan"]))
    return hot


def fit_log_parabola(xs, ys, ws=None):
    """A Gaussian peak is exactly a parabola in log intensity, so fitting one
    to an ion's clean scans predicts what the ion should be doing in the scans
    left out. Returns (a, b, c) for a + b*x + c*x^2, or None.

    Weighted, because log-intensity variance is 1/counts and counts fall away
    on the flanks. Unweighted, the noisiest points dominate the shape.
    """
    if ws is None:
        ws = [1.0] * len(xs)
    n = sum(ws)
    sx = sum(w * x for w, x in zip(ws, xs))
    sx2 = sum(w * x * x for w, x in zip(ws, xs))
    sx3 = sum(w * x ** 3 for w, x in zip(ws, xs))
    sx4 = sum(w * x ** 4 for w, x in zip(ws, xs))
    sy = sum(w * y for w, y in zip(ws, ys))
    sxy = sum(w * x * y for w, x, y in zip(ws, xs, ys))
    sx2y = sum(w * x * x * y for w, x, y in zip(ws, xs, ys))
    m = [[n, sx, sx2, sy], [sx, sx2, sx3, sxy], [sx2, sx3, sx4, sx2y]]
    for col in range(3):                       # Gaussian elimination
        piv = max(range(col, 3), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-12:
            return None
        m[col], m[piv] = m[piv], m[col]
        for r in range(3):
            if r == col:
                continue
            f = m[r][col] / m[col][col]
            for k in range(col, 4):
                m[r][k] -= f * m[col][k]
    return tuple(m[i][3] / m[i][i] for i in range(3))


def log_debiased(r):
    """log(intensity), corrected for the curvature of the logarithm.

    Intensity carries multiplicative counting noise, so E[log(1+e)] is about
    -sigma^2/2 rather than 0, and sigma^2 is 1/counts -- LARGER on a peak's
    flanks than at its apex. Fitting raw logs therefore pulls the flanks down
    relative to the apex and predicts an apex that is too low, which reads as
    the ion sitting ABOVE its own profile. That is a false positive for
    additive interference, manufactured by the logarithm.
    """
    if r["inten"] <= 0:
        # Not a measurement to take a log of. Returning None lets the caller
        # drop the scan; raising took the whole run down with a ValueError
        # whose exit 1 is indistinguishable from "ran and found nothing".
        return None
    y = math.log(r["inten"])
    if r["counts"] and r["counts"] > 0:
        y += 0.5 / r["counts"]
    return y


def profile_departure(rows, marked_for, thresh, out, label, why, mode="dim"):
    """Departure of an ion from its own profile in marked scans, measured
    against the same ion's unmarked scans by the identical procedure.

    Elution position alone makes an ion brighter near its apex, and marked
    scans sit near apexes, so a raw comparison of marked against unmarked
    reads peak shape as an effect. Fitting the unmarked scans and predicting
    the marked ones controls that, but introduces biases of its own, all of
    which push the answer the same way -- toward the ion appearing ABOVE its
    own profile, which is a false positive for additive interference:

      the logarithm        E[log(1+e)] is about -sigma^2/2, and sigma^2 is
                           1/counts, which is larger on the flanks than at
                           the apex                                 (+1.2%)
      an unweighted fit    the noisiest points dominate the shape    (+1.2%)
      selecting dim ions   on their MAXIMUM observed intensity, which
                           truncates upward fluctuations at the apex (+1.7%)
      extrapolation        a parabola fitted only to flanks, evaluated in
                           the middle, with noise                   (+1.2%)

    Measured on planted cohorts with nothing suppressed, those left a null of
    1.089 -- read naively, "interference" at p = 0.002 in data containing
    none. So the null is not assumed to be 1.0, it is calibrated:

      the reference comes from UNMARKED scans only, because a marked scan in
      the fit drags the profile toward the suppressed values, flattens it,
      and can invert the sign of the result;

      every unmarked scan is then scored against the OTHER unmarked scans by
      the identical procedure, and that placebo is subtracted.

    On planted cohorts this measures a null of 1.0005 and recovers planted
    factors of 0.70, 0.85 and 0.95 as 0.7004, 0.8504 and 0.9505.
    """
    dim_cut = thresh * DIM_FRACTION
    bygroup = collections.defaultdict(list)
    for r in rows:
        bygroup[(r["run"], r["group"])].append(r)

    deltas, skipped = [], collections.Counter()
    for key, rs in bygroup.items():
        if any(r["inten"] <= 0 for r in rs):
            # One unusable scan makes the ion's profile unusable; dropping the
            # ion is honest where dropping the scan would silently change the
            # shape being fitted.
            skipped["an intensity of zero or less"] += 1
            continue
        peak_max = max(r["inten"] for r in rs)
        if mode == "dim":
            # Dim and never flagged: dead time cannot reach these, so a fall
            # here can only be something that acted on the whole scan.
            if any(r["sat"] for r in rs):
                skipped["ion is itself flagged"] += 1
                continue
            if peak_max > dim_cut:
                skipped["ion is not dim"] += 1
                continue
        else:
            # Bright and flagged: dead time undercounts the channel.
            if not any(r["sat"] for r in rs):
                skipped["ion is never flagged"] += 1
                continue
            if peak_max <= dim_cut:
                skipped["ion is not bright"] += 1
                continue
        if len(rs) < MIN_SCANS_FOR_LOO:
            skipped["fewer than %d scans to leave one out of" % MIN_SCANS_FOR_LOO] += 1
            continue
        marked = marked_for(rs)
        mk = [r for r in rs if (r["run"], r["scan"]) in marked]
        un = [r for r in rs if (r["run"], r["scan"]) not in marked]
        if not mk:
            skipped["no marked scan in its profile"] += 1
            continue
        if not un:
            skipped["every scan marked, so no control within the ion"] += 1
            continue

        if len(un) < MIN_SCANS_FOR_LOO - 1:
            skipped["fewer than %d unmarked scans to build a reference from"
                    % (MIN_SCANS_FOR_LOO - 1)] += 1
            continue

        def predict(ref, target):
            c0 = sum(o["scan"] for o in ref) / float(len(ref))
            fit = fit_log_parabola(
                [o["scan"] - c0 for o in ref],
                [log_debiased(o) for o in ref],
                [max(1.0, o["counts"] or 1.0) for o in ref])
            if fit is None:
                return None
            a, b, c = fit
            x = target["scan"] - c0
            return a + b * x + c * x * x

        # The reference is built from UNMARKED scans only. Including marked
        # ones drags the fitted profile toward the suppressed values, which
        # flattens it and can invert the sign of the whole result.
        dm, failed = [], False
        for r in mk:
            pr = predict(un, r)
            if pr is None:
                failed = True
                break
            dm.append(log_debiased(r) - pr)

        # The placebo: each unmarked scan scored the same way, against the
        # other unmarked scans. This is what calibrates the null, rather than
        # assuming it sits at 1.0 when measurement says it does not.
        dp = []
        if not failed:
            for i, r in enumerate(un):
                rest = [o for j, o in enumerate(un) if j != i]
                if len(rest) < 4:
                    continue
                pr = predict(rest, r)
                if pr is None:
                    failed = True
                    break
                dp.append(log_debiased(r) - pr)
        if failed or not dp:
            skipped["the reference fit was degenerate"] += 1
            continue

        deltas.append(sum(dm) / len(dm) - sum(dp) / len(dp))

    out("")
    out("%s" % label)
    out("   %s" % why)
    out("   %s. Every scan is predicted from the ion's others and the"
        % ("dim means a peak never exceeding %.0f" % dim_cut if mode == "dim"
           else "bright means a peak exceeding %.0f" % dim_cut))
    out("   statistic is marked minus unmarked WITHIN the ion, so the")
    out("   procedure's own bias is shared and cancels.")
    for reason, n in sorted(skipped.items(), key=lambda kv: -kv[1]):
        out("   skipped %4d: %s" % (n, reason))
    if len(deltas) < 3:
        out("   Only %d usable ion(s). UNDECIDED, which is not the same as flat."
            % len(deltas))
        return None

    n = len(deltas)
    m = sum(deltas) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in deltas) / (n - 1)) if n > 1 else 0.0
    se = sd / math.sqrt(n) if n else 0.0
    srt = sorted(deltas)
    med = (srt[n // 2] if n % 2 else (srt[n // 2 - 1] + srt[n // 2]) / 2.0)
    below = sum(1 for d in deltas if d < 0)
    above = sum(1 for d in deltas if d > 0)
    out("   %d %s ion(s) scored" % (n, mode))
    out("   mean marked-minus-unmarked departure %+.4f -> factor %.4f"
        % (m, math.exp(m)))
    out("   median %+.4f -> factor %.4f" % (med, math.exp(med)))
    out("   %d below 0, %d above" % (below, above))
    pv = 1.0
    if se > 0:
        t = m / se
        pv = t_pvalue(t, n - 1)
        out("   t = %+.2f on %d df,  p = %.3g" % (t, n - 1, pv))
    ps = sign_test(below, above)
    out("   sign test (assumption-free, corroborating)  p = %.3g" % ps)
    out("   The estimator itself is unbiased and unattenuated on planted")
    out("   data, but marking is an imperfect proxy for which scans were")
    out("   really charge-limited, and mismarking pulls the factor toward 1.")
    out("   Treat it as a floor rather than a size.")
    return math.exp(m), pv, n, below, above, ps


def t_pvalue(t, df):
    """Two-sided p for Student t, by the regularised incomplete beta."""
    x = float(df) / (df + t * t)
    return max(0.0, min(1.0, betainc(0.5 * df, 0.5, x)))


def betainc(a, b, x):
    """Regularised incomplete beta I_x(a,b) via the continued fraction."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return math.exp(lbeta) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbeta) * _betacf(b, a, 1 - x) / b


def _betacf(a, b, x, itmax=200, eps=3e-16):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def sign_test(below, above):
    """Exact two-sided binomial on the pairs that moved."""
    n = below + above
    if n == 0:
        return 1.0
    k = min(below, above)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / float(2 ** n)
    return min(1.0, 2.0 * tail)


def flag_clustering(rows, hot, out):
    """Weak test: do flags co-occur in scans more than chance? Says nothing
    about mechanism, because bright ions cluster in scans for other reasons."""
    byscan = collections.defaultdict(list)
    for r in rows:
        byscan[(r["run"], r["scan"])].append(r)
    counts = [sum(1 for r in rs if r["sat"]) for rs in byscan.values()]
    mean = sum(counts) / float(len(counts))
    var = sum((c - mean) ** 2 for c in counts) / float(len(counts))
    out("")
    out("A  DO FLAGS CLUSTER BY SCAN?")
    out("   %d scan(s); flagged per scan mean %.3f variance %.3f"
        % (len(counts), mean, var))
    if mean > 0:
        out("   index of dispersion %.3f  (1.0 = no clustering)" % (var / mean))
    out("   Weak on its own: an ion's own elution puts its bright scans next")
    out("   to each other, so clustering is expected under BOTH mechanisms.")
    return (var / mean) if mean > 0 else None


def coherence(rows, out):
    """Observed within-ion scan-to-scan variation against 1/sqrt(Counts)."""
    bygroup = collections.defaultdict(list)
    for r in rows:
        bygroup[(r["run"], r["group"])].append(r)
    obs, pred = [], []
    for rs in bygroup.values():
        if len(rs) < 3:
            continue
        xs = [r["inten"] for r in rs]
        m = sum(xs) / len(xs)
        if m <= 0:
            continue
        sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
        obs.append(sd / m)
        cs = [r["counts"] for r in rs if r["counts"]]
        if cs:
            pred.append(1.0 / math.sqrt(sum(cs) / len(cs)))
    out("")
    out("D  COHERENCE AGAINST THE COUNTING BAND")
    if not obs:
        out("   No group spans 3+ scans; nothing to measure.")
        return None
    obs.sort()
    mo = obs[len(obs) // 2]
    out("   %d group(s) with 3+ scans; median relative scatter %.4f" % (len(obs), mo))
    if pred:
        pred.sort()
        mp = pred[len(pred) // 2]
        out("   median 1/sqrt(Counts) expectation            %.4f" % mp)
        out("   observed / expected                          %.2f" % (mo / mp if mp else float("nan")))
        out("   Above 1 means something beyond counting noise moves these")
        out("   ions; below 1 means the export is smoothed and the scans are")
        out("   not independent measurements.")
    return mo


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paths", nargs="+", help="ApexRT-style CSV files, or globs")
    args = p.parse_args(argv)

    paths = []
    for a in args.paths:
        hit = sorted(glob.glob(a))
        paths.extend(hit if hit else [a])

    lines = []
    out = lambda s="": lines.append(s)
    try:
        rows = read_rows(paths)
        lo, hi = refuse_if_unusable(rows, out)
    except Unusable as e:
        print("\n".join(lines))
        print("")
        print("CANNOT RUN")
        for chunk in str(e).split(". "):
            print("  %s" % chunk.strip().rstrip(".") + ".")
        print("")
        print("  This is a refusal, not a result. The scan-level and")
        print("  detector-level mechanisms predict different things about this")
        print("  input, and it cannot tell them apart, so any number reported")
        print("  here would be a property of the file rather than the spray.")
        return 2

    hot = scan_has_flag(rows)
    out("        %d of %d scan(s) carry a flagged detection"
        % (len(hot), len(set((r["run"], r["scan"]) for r in rows))))

    flag_clustering(rows, hot, out)
    flag_marker = lambda rs: hot
    loud = load_marker(rows, out)
    bload = profile_departure(
        rows, loud, lo, out,
        "B  DIM IONS IN HIGH-CURRENT SCANS  (tests SOURCE)",
        "Conditioned on per-scan total ion current, which is what droplet\n"
        "   charge competition is driven by. A source effect depresses the\n"
        "   bright ion too, pushing it UNDER the saturation flag, so the flag\n"
        "   can miss the very scans where competition occurred. Current does\n"
        "   not. Dim unflagged ions cannot be reached by dead time, so a fall\n"
        "   here is scan-level or nothing.",
        mode="dim")
    bbright = profile_departure(
        rows, flag_marker, lo, out,
        "C  BRIGHT FLAGGED IONS AGAINST THEIR OWN FLANKS  (tests DETECTOR)",
        "Dead time undercounts a bright CHANNEL, so it shows up as an apex\n"
        "   lower than the ion's own unflagged flanks predict. Flat dim ions\n"
        "   in test B rule OUT a scan-level effect; they are not evidence of\n"
        "   a detector one, which needs this separate measurement.\n"
        "   Caveat: this extrapolates to the apex, so it leans on the peak\n"
        "   shape model in a way test B's interpolation does not.",
        mode="bright")
    coherence(rows, out)

    out("")
    out("VERDICT")

    def reads(b):
        if b is None:
            return None
        med, pv, n, below, above, ps = b
        if pv < ALPHA and med < 1.0:
            return ("down", med, pv, n, below)
        if pv < ALPHA and med > 1.0:
            return ("up", med, pv, n, below)
        return ("flat", med, pv, n, below)

    rf, rl = reads(bbright), reads(bload)
    if rf is None and rl is None:
        out("  UNDECIDED. Too few dim ions had both a marked scan and enough")
        out("  unmarked scans to predict from, under either conditioner, so")
        out("  neither mechanism is supported. Do not apply a per-scan")
        out("  correction on this evidence.")
        print("\n".join(lines))
        return 1

    rc = 0
    src = rl is not None and rl[0] == "down"
    det = rf is not None and rf[0] == "down"

    # The two tests are NOT symmetric, and treating them as if they were is
    # the mistake that would license the wrong correction.
    #
    #   A dim unflagged ion cannot be reached by dead time, so a fall there
    #   can only be scan-level. Test B discriminates.
    #
    #   A bright ion falls below its own flanks under EITHER mechanism --
    #   dead time undercounts it, and droplet competition depresses it along
    #   with everything else in the scan. Test C therefore adds information
    #   only when test B is flat.
    if src:
        out("  SOURCE-LEVEL. Dim ions fall %.1f%% below what their own"
            % (100.0 * (1 - rl[1])))
        out("  predicted profile says, in scans with high total current")
        out("  (%d of %d, p = %.3g). They are never flagged and never bright,"
            % (rl[4], rl[3], rl[2]))
        out("  so dead time cannot reach them: the scan is the unit, as")
        out("  droplet charge competition predicts.")
        if det:
            out("")
            out("  Bright flagged ions are ALSO low against their flanks")
            out("  (%.4f, p = %.3g). That is expected under competition rather"
                % (rf[1], rf[2]))
            out("  than extra evidence: a scan-level effect depresses the")
            out("  bright ion too. It cannot be separated from an additional")
            out("  dead-time component here, because the only load axis")
            out("  available is the surviving signal itself.")
        else:
            out("")
            out("  Bright flagged ions are NOT low against their flanks (%s),"
                % ("%.4f" % rf[1] if rf else "undecided"))
            out("  which is worth a look: competition should reach them too.")
        out("")
        out("  A per-scan factor is the right shape. Its residuals are the")
        out("  next test: if what is left over tracks hydrophobicity, one")
        out("  scalar per scan is still the wrong model, because surface")
        out("  activity decides who reaches the droplet surface.")
        rc = 1
    elif det:
        out("  DETECTOR-LEVEL. Bright flagged ions fall %.1f%% below what"
            % (100.0 * (1 - rf[1])))
        out("  their own unflagged flanks predict (%d of %d, p = %.3g), while"
            % (rf[4], rf[3], rf[2]))
        out("  dim ions in high-current scans sit where predicted (%s)."
            % ("%.4f" % rl[1] if rl else "undecided"))
        out("  Dim ions are the ones competition would have reached, and they")
        out("  did not move, so the flag is tracking the CHANNEL. A per-scan")
        out("  suppression factor would be correcting dim ions for something")
        out("  that did not happen to them.")
        rc = 1
    else:
        ups = [r for r in (rl, rf) if r is not None and r[0] == "up"]
        if ups:
            out("  ABOVE their own profile (median %.4f, p = %.3g). That is"
                % (ups[0][1], ups[0][2]))
            out("  neither mechanism. It is the signature of additive")
            out("  interference: co-eluting signal raising the ion along with")
            out("  its neighbours. Counts are non-negative, so interference")
            out("  can only add -- suppression is not what this shows.")
            rc = 1
        else:
            out("  NOT SOURCE-LEVEL. Dim ions in high-current scans sit where")
            out("  their own profiles predict (%s), so the scan is not the"
                % ("%.4f" % rl[1] if rl else "undecided"))
            out("  unit and a per-scan correction would be fitting noise.")
            out("  Bright flagged ions also sit where their flanks predict")
            out("  (%s), so no detector effect is visible either."
                % ("%.4f" % rf[1] if rf else "undecided"))
            out("")
            out("  Both are real answers. But note the asymmetry before")
            out("  reading too much into them: flat dim ions rule a")
            out("  scan-level effect OUT, while a low bright apex would have")
            out("  ruled NOTHING in, since competition and dead time both")
            out("  produce one.")
    print("\n".join(lines))
    return rc


if __name__ == "__main__":
    sys.exit(main())
