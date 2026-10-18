#!/usr/bin/env python3
"""Is scan-to-scan variation counting noise, and is any of it shared?

Counting gives the one parameter-free prediction available in a mass
spectrum: if a channel's value is a number of detected ions, its variance
equals its mean, so the Fano factor var/mean is 1 and relative scatter is
1/sqrt(N). Anything above that is variation the detector did not have to
produce; anything below means the scans are not independent measurements.

That matters because it decides the shape of a coherence test. If noise is
Poisson, a tolerance band must scale as 1/sqrt(N) and a fixed percentage is
wrong. If noise is multiplicative, relative scatter is flat with level and a
fixed percentage is right. The answer is not a matter of taste and it is not
the same at all intensities.

  A  COUNTING  Fano factor by level, from second differences across the whole
     series, against a Poisson null carrying the same trend.
  B  SHARED    with two channels where one is contained in the other (a base
     peak inside a total current), whether the part rises or falls as the
     whole moves. Rising together is common-mode delivery; falling as the
     whole rises is competition for a fixed resource.

Second differences rather than hand-picked flat stretches, because selecting
windows by how little they vary selects on the quantity being measured. A
range or variance criterion makes the Fano factor come out low by
construction. Second differences remove any local linear trend without
choosing which scans to keep, and the Poisson null absorbs what curvature
remains.

  python3 tools/counting_noise.py scans.jsonl
  python3 tools/counting_noise.py scans.jsonl --part bpi --whole tic

Input is one row per scan, JSONL or CSV, carrying a scan index and at least
one numeric channel; a run column separates acquisitions. Exit 0 if the tests
ran, 1 if they ran and found something to act on, 2 if the input cannot
support them.
"""

import argparse
import collections
import csv
import glob
import json
import math
import os
import sys

# A Poisson null should land near 1. Further than this and the estimator,
# not the data, is deciding the answer.
NULL_TOLERANCE = (0.5, 2.0)

# Below this many usable scans in a run there is nothing to bin by level.
MIN_SCANS = 60

# A channel whose neighbouring values repeat this often has been smoothed or
# forward-filled, so consecutive scans are not independent draws.
MAX_REPEAT_SHARE = 0.20

# Levels are split into this many bins, lowest to highest.
NBINS = 6

# A Fano ratio at or below this counts as "counting describes it". Three
# rather than one, because the null itself carries trend curvature and a
# factor of two either way is not evidence of a different mechanism.
COUNTING_RATIO = 3.0

# Smoothing half-width for the trend the null is generated from. Wide enough
# that the trend carries little of the noise it is meant to exclude.
TREND_HALFWIDTH = 15

NON_CHANNEL = ("scan", "run", "rt", "time", "scantimemin", "scanindex")


class Unusable(Exception):
    """The input cannot support the tests. Not a finding -- a refusal."""


def read_rows(paths):
    rows = []
    for p in paths:
        try:
            fh = open(p, newline="")
        except OSError as e:
            raise Unusable("cannot open %s: %s" % (p, e))
        with fh:
            head = fh.read(2048)
            fh.seek(0)
            if head.lstrip().startswith("{"):
                for line in fh:
                    line = line.strip()
                    if line:
                        try:
                            rows.append(json.loads(line))
                        except ValueError:
                            continue
            else:
                rows.extend(list(csv.DictReader(fh)))
    if not rows:
        raise Unusable("no parseable rows in %s" % ", ".join(paths))
    return rows


def channels_of(rows):
    names = []
    for k in rows[0]:
        if k.lower() in NON_CHANNEL:
            continue
        try:
            float(rows[0][k])
        except (TypeError, ValueError):
            continue
        names.append(k)
    if not names:
        raise Unusable("no numeric channel found; columns are %s"
                       % ", ".join(sorted(rows[0].keys())))
    return names


def scan_key(r):
    for k in ("scan", "ScanIndex", "scanindex", "ScanIndex "):
        if k in r:
            return float(r[k])
    for k in r:
        if k.lower() in ("scan", "scanindex"):
            return float(r[k])
    raise Unusable("no scan index column; one of scan / ScanIndex is required")


def series_of(rows, chan):
    by = collections.defaultdict(list)
    for r in rows:
        run = r.get("run") or r.get("Run") or "all"
        try:
            v = float(r[chan])
            s = scan_key(r)
        except (KeyError, TypeError, ValueError):
            continue
        by[run].append((s, v))
    out = {}
    for run, pts in by.items():
        pts.sort()
        out[run] = [v for _, v in pts]
    return out


def refuse_if_unusable(y, chan, run, out):
    if len(y) < MIN_SCANS:
        raise Unusable("%s/%s has %d scan(s), below the %d needed to bin by "
                       "level" % (run, chan, len(y), MIN_SCANS))
    if len(set(y)) == 1:
        raise Unusable("%s/%s holds the single value %g in all %d scans, so "
                       "there is no variation to attribute"
                       % (run, chan, y[0], len(y)))
    if min(y) < 0:
        raise Unusable("%s/%s holds negative values, which cannot be counts "
                       "and have no Poisson expectation" % (run, chan))
    rep = sum(1 for a, b in zip(y, y[1:]) if a == b)
    share = float(rep) / max(1, len(y) - 1)
    if share > MAX_REPEAT_SHARE:
        raise Unusable(
            "%s/%s repeats its previous value in %.0f%% of scans. A smoothed "
            "or forward-filled series has had the scan-to-scan variation "
            "removed, so a Fano factor computed from it measures the "
            "smoothing." % (run, chan, 100.0 * share))
    out("        %-6s %5d scans, %d distinct, %.0f%% repeat previous"
        % (chan, len(y), len(set(y)), 100.0 * share))


def second_differences(y):
    """d_i = y_i - 2 y_{i+1} + y_{i+2}. Removes any local linear trend, so no
    scan has to be selected or discarded. For independent noise of variance
    s^2 the result has variance 6 s^2."""
    return [(y[i] - 2.0 * y[i + 1] + y[i + 2],
             (y[i] + y[i + 1] + y[i + 2]) / 3.0)
            for i in range(len(y) - 2)]


def fano_by_level(y, nbins=NBINS):
    pairs = [(m, d) for d, m in second_differences(y) if m > 0]
    pairs.sort()
    out = []
    per = max(25, len(pairs) // nbins)
    for s in range(0, len(pairs), per):
        ch = pairs[s:s + per]
        if len(ch) < 25:
            continue
        mm = sum(m for m, _ in ch) / len(ch)
        ds = [d for _, d in ch]
        md = sum(ds) / len(ds)
        var = sum((d - md) ** 2 for d in ds) / (len(ds) - 1) / 6.0
        out.append((mm, var / mm if mm > 0 else float("nan"), len(ch)))
    return out


def smooth(y, k=TREND_HALFWIDTH):
    n = len(y)
    return [sum(y[max(0, i - k):min(n, i + k + 1)])
            / (min(n, i + k + 1) - max(0, i - k)) for i in range(n)]


def poisson(lam, rng):
    if lam <= 0:
        return 0.0
    if lam < 30:
        L = math.exp(-lam)
        k, p = 0, 1.0
        while True:
            p *= rng.random()
            if p <= L:
                return float(k)
            k += 1
            if k > 400:
                return float(k)
    return max(0.0, lam + math.sqrt(lam) * rng.gauss(0, 1))


def loo_residuals(y, w=3):
    """Each scan predicted from w neighbours each side by a local line. Used
    for the shared-term test, where a deviation per scan is needed rather
    than a variance per level."""
    res = {}
    for i in range(w, len(y) - w):
        xs = [j - i for j in range(i - w, i + w + 1) if j != i]
        ys = [y[j] for j in range(i - w, i + w + 1) if j != i]
        m = len(xs)
        sx = sum(xs)
        sxx = sum(x * x for x in xs)
        sy = sum(ys)
        sxy = sum(x * v for x, v in zip(xs, ys))
        den = m * sxx - sx * sx
        if den == 0:
            continue
        b = (m * sxy - sx * sy) / den
        a = (sy - b * sx) / m
        if a > 0:
            res[i] = (y[i] - a) / a
    return res


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan"), 1.0, 0, 0
    r = sxy / math.sqrt(sxx * syy)
    pos = sum(1 for a, b in zip(xs, ys) if (a - mx) * (b - my) > 0)
    neg = sum(1 for a, b in zip(xs, ys) if (a - mx) * (b - my) < 0)
    m = pos + neg
    if m == 0:
        return r, 1.0, pos, neg
    z = (min(pos, neg) - m / 2.0) / math.sqrt(m / 4.0)
    return r, math.erfc(abs(z) / math.sqrt(2)), pos, neg


def main(argv=None):
    import random
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paths", nargs="+", help="per-scan tables, JSONL or CSV")
    p.add_argument("--channel", action="append", default=[],
                   help="channel to test; repeatable. Default: every numeric one")
    p.add_argument("--part", help="channel contained in --whole, e.g. a base peak")
    p.add_argument("--whole", help="channel containing --part, e.g. a total current")
    p.add_argument("--seed", type=int, default=4)
    args = p.parse_args(argv)

    paths = []
    for a in args.paths:
        hit = sorted(glob.glob(a))
        paths.extend(hit if hit else [a])

    lines = []
    out = lambda s="": lines.append(s)
    rng = random.Random(args.seed)
    rc = 0
    try:
        rows = read_rows(paths)
        chans = args.channel or channels_of(rows)
        missing = [c for c in chans if c not in rows[0]]
        if missing:
            raise Unusable("no such channel(s): %s" % ", ".join(missing))
        if (args.part is None) != (args.whole is None):
            raise Unusable("--part and --whole go together; one alone cannot "
                           "say what contains what")

        out("INPUT   %d scan row(s) from %d file(s)" % (len(rows), len(paths)))
        series = {c: series_of(rows, c) for c in chans}
        runs = sorted(set(r for c in chans for r in series[c]))
        out("        %d run(s): %s" % (len(runs), ", ".join(runs[:4])
                                       + (" ..." if len(runs) > 4 else "")))
        for run in runs:
            for c in chans:
                if run in series[c]:
                    refuse_if_unusable(series[c][run], c, run, out)

        # --- A. counting ------------------------------------------------
        out("")
        out("A  IS IT COUNTING NOISE?")
        out("   Fano factor var/mean from second differences, against a")
        out("   Poisson null carrying the same smoothed trend. The null")
        out("   absorbs whatever trend curvature the differences leave, so")
        out("   the ratio is what matters, not the raw factor.")
        out("")
        out("   %-10s %-6s %12s %10s %9s %9s  %s"
            % ("run", "chan", "mean level", "F_obs", "F_null", "ratio", "rel SD"))
        verdicts = []
        for run in runs:
            for c in chans:
                y = series[c].get(run)
                if not y:
                    continue
                obs = fano_by_level(y)
                trend = smooth(y)
                nul = fano_by_level([poisson(t, rng) for t in trend])
                if not obs or not nul:
                    continue
                bad_null = [fn for _, fn, _ in nul
                            if not (NULL_TOLERANCE[0] <= fn <= NULL_TOLERANCE[1])]
                for (m, f, n), (_, fn, _) in zip(obs, nul):
                    ratio = f / fn if fn else float("nan")
                    out("   %-10s %-6s %12.0f %10.2f %9.2f %9.1f  %6.2f%%"
                        % (run[:10], c, m, f, fn, ratio,
                           100.0 * math.sqrt(max(0.0, f * m)) / m))
                    verdicts.append((run, c, m, ratio))
                if bad_null:
                    out("   %-10s %-6s NULL OFF 1.0 in %d of %d bin(s) "
                        "(%.2f to %.2f): read the ratios with care, the"
                        % (run[:10], c, len(bad_null), len(nul),
                           min(bad_null), max(bad_null)))
                    out("   %-10s %-6s estimator is contributing." % ("", ""))

        # where does counting stop describing it?
        out("")
        byc = collections.defaultdict(list)
        for run, c, m, ratio in verdicts:
            byc[c].append((m, ratio))
        for c, pts in sorted(byc.items()):
            pts.sort()
            near = [m for m, r in pts if r <= COUNTING_RATIO]
            far = [m for m, r in pts if r > COUNTING_RATIO]
            if near and far:
                hi_near, lo_far = max(near), min(far)
                if lo_far > hi_near:
                    out("   %s: counting describes it at and below %.0f, and"
                        % (c, hi_near))
                    out("   not at or above %.0f." % lo_far)
                else:
                    # Runs are binned by their own levels, so the two sets can
                    # interleave. Saying "up to X and not above Y" with Y < X
                    # would be nonsense; the honest statement is the overlap.
                    out("   %s: the crossing is not clean -- levels up to %.0f"
                        % (c, hi_near))
                    out("   still look like counting while levels from %.0f do"
                        % lo_far)
                    out("   not, so the two overlap between %.0f and %.0f."
                        % (lo_far, hi_near))
                out("   A tolerance band has to change form across it:")
                out("   1/sqrt(N) below, a flat percentage above.")
                rc = 1
            elif far:
                out("   %s: nowhere does counting describe it (ratio %.0f to"
                    % (c, min(r for _, r in pts)))
                out("   %.0f). Scatter is multiplicative here, so a FLAT"
                    % max(r for _, r in pts))
                out("   percentage is the right form of band and 1/sqrt(N)")
                out("   would be far too tight.")
                rc = 1
            else:
                out("   %s: counting describes it at every level measured" % c)
                out("   (ratio %.1f to %.1f). A band scaling as 1/sqrt(N) is"
                    % (min(r for _, r in pts), max(r for _, r in pts)))
                out("   the right form, and a flat percentage would be too")
                out("   loose at the top.")

        # --- B. shared term ---------------------------------------------
        if args.part and args.whole:
            out("")
            out("B  IS ANY OF IT SHARED, AND WHICH WAY?")
            pr, wr = series.get(args.part, {}), series.get(args.whole, {})
            tot_xs, tot_ys = [], []
            for run in runs:
                a, b = pr.get(run), wr.get(run)
                if not a or not b or len(a) != len(b):
                    continue
                over = sum(1 for x, w in zip(a, b) if x > w)
                if over:
                    raise Unusable(
                        "%s exceeds %s in %d scan(s), so it is not contained "
                        "in it and removing it would be meaningless"
                        % (args.part, args.whole, over))
                shares = [x / w for x, w in zip(a, b) if w > 0]
                sh = sum(shares) / len(shares)
                # The part is removed from the whole. Left in, the part is
                # correlated with itself and the sign of this test is forced.
                rest = [w - x for x, w in zip(a, b)]
                if sh > 0.5:
                    raise Unusable(
                        "%s is %.0f%% of %s on average; once it is removed "
                        "too little independent signal is left for the "
                        "remainder to say anything"
                        % (args.part, 100.0 * sh, args.whole))
                ra, rb = loo_residuals(a), loo_residuals(rest)
                common = sorted(set(ra) & set(rb))
                tot_xs += [rb[i] for i in common]
                tot_ys += [ra[i] for i in common]
                out("   %-10s %s is %.1f%% of %s; %d paired scan(s)"
                    % (run[:10], args.part, 100.0 * sh, args.whole, len(common)))
            if len(tot_xs) < 50:
                out("   Only %d paired scan(s). UNDECIDED." % len(tot_xs))
            else:
                r, pv, pos, neg = pearson(tot_xs, tot_ys)
                out("   %d scan(s) pooled;  r = %+.3f,  sign test p = %.3g"
                    % (len(tot_xs), r, pv))
                out("   (%d deviations same sign, %d opposite)" % (pos, neg))
                out("   %s removed from %s first, so the part cannot be"
                    % (args.part, args.whole))
                out("   correlated with itself.")
                out("")
                if r > 0.1 and pv < 0.01:
                    out("   COMMON-MODE. %s rises as the rest of %s rises"
                        % (args.part, args.whole))
                    out("   (r = %+.3f). Everything in a scan moves together," % r)
                    out("   which is delivery varying scan to scan, not")
                    out("   competition for a fixed resource. It also means a")
                    out("   shared per-scan factor exists and is estimable:")
                    out("   many ions constrain one number.")
                    out("   It explains %.1f%% of the variance in %s, so most"
                        % (100.0 * r * r, args.part))
                    out("   of the scatter is still per-ion.")
                    rc = 1
                elif r < -0.1 and pv < 0.01:
                    out("   COMPETITION. %s FALLS as the rest of %s rises"
                        % (args.part, args.whole))
                    out("   (r = %+.3f). That is the signature of a fixed" % r)
                    out("   resource being shared out -- a droplet carrying a")
                    out("   fixed maximum charge would do this.")
                    rc = 1
                else:
                    out("   NO SHARED TERM FOUND. r = %+.3f at p = %.3g, so"
                        % (r, pv))
                    out("   %s does not track the rest of %s either way."
                        % (args.part, args.whole))
                    out("   Scatter is per-channel here, and a per-scan")
                    out("   correction would have nothing to correct.")
            out("")
            out("   What this cannot settle: a positive r can also come from")
            out("   real signal co-varying inside the window the local fit")
            out("   treats as flat, and an acquisition that never approaches")
            out("   its charge limit will show no competition whether or not")
            out("   the mechanism is real. Absence here is not absence of the")
            out("   mechanism.")
    except Unusable as e:
        print("\n".join(lines))
        print("")
        print("CANNOT RUN")
        print("  %s." % str(e).rstrip("."))
        print("")
        print("  This is a refusal, not a result. Counting noise has one")
        print("  prediction and this input cannot be held against it, so any")
        print("  factor printed here would describe the file rather than the")
        print("  measurement.")
        return 2

    print("\n".join(lines))
    return rc


if __name__ == "__main__":
    sys.exit(main())
