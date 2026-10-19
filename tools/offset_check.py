#!/usr/bin/env python3
"""Ask whether a shared offset is hiding behind the Atlas's controls.

docs/atlas-convergence.md established the problem and prescribed three checks.
None of them was in any tool, and nothing in the chain ran them. This is them.

The problem, in one line: scatter averages away and a shared offset does not, so
every control that measures spread reports improving health as an offset
persists. On the 24 Colorado PAX8 runs the per-run means sit at -1.570 ppm with
24 of 24 below zero, and the between-run SD of 0.643 ppm -- the number the
premise rests on -- is unchanged by that offset. At 1,000 runs the coordinate
would be known to 0.02 ppm and still be 1.57 ppm wrong: a confidence interval
excluding the truth by about 77 standard errors.

    python3 tools/offset_check.py ppm.jsonl
    python3 tools/offset_check.py ppm.jsonl --anchor lockmass.jsonl

Input is JSON lines or CSV, one row per observation:

    {"run": "CH163", "uid": "LSIDSQSSSSGPR_2", "ppm": -1.84}

`ppm` is observed minus calculated, in parts per million. One row per PSM; the
tool does the run-level aggregation itself, because a t test on pooled PSMs
answers a different question from one on run means and the pooled version is
the one that looks significant for the wrong reason.

WHAT IT DOES THAT A RESTATEMENT WOULD NOT

It injects offsets of a known size and reports which controls move. A control
whose value is identical at 0 and 20 ppm of shared offset is not insensitive to
offsets, it is BLIND to them by construction, and the output says which of the
project's controls are in that class. That turns "spread cannot protect the
Atlas" from an argument into a measurement.

THE CIRCULARITY IT REFUSES

If the stored m/z was computed from the sequence rather than measured, then
observed minus calculated is zero by arithmetic and a calibration fitted to it
returns 0.00 ppm over however many rows there are. That happened in this project
over 1,994,456 of them. Rows whose ppm is zero to machine precision are counted
and excluded, and if they are the majority the tool refuses outright.

Exit 0 when no offset is detected, 1 when one is, 2 when the input cannot carry
the question.
"""

import argparse
import csv
import json
import math
import os
import statistics
import sys
from collections import defaultdict

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
from ladder_stats import t_critical  # noqa: E402

# Offsets injected to measure each control's sensitivity. A control that returns
# the same value at 0.0 and at 20.0 ppm is blind, not merely tolerant.
PROBE_OFFSETS = (0.0, 1.0, 5.0, 20.0)
IDENTICAL = 1e-12        # a ppm this small is arithmetic, not a measurement
BLIND_TOL = 1e-9         # a control that moves less than this is blind
SIGN_ALPHA = 0.01


def read_rows(path):
    text = open(path).read()
    stripped = text.lstrip()
    rows = []
    if stripped.startswith("["):
        rows = json.loads(text)
    elif stripped.startswith("{"):
        for n, line in enumerate(text.splitlines(), 1):
            line = line.strip()
            if line:
                try:
                    rows.append(json.loads(line))
                except ValueError as exc:
                    raise ValueError("%s line %d: %s" % (path, n, exc))
    else:
        rows = list(csv.DictReader(text.splitlines()))
    out, bad = [], 0
    for r in rows:
        run = r.get("run") or r.get("run_name") or r.get("sample") or r.get("sample_id")
        val = r.get("ppm", r.get("ppm_error", r.get("error_ppm")))
        try:
            val = float(val)
        except (TypeError, ValueError):
            bad += 1
            continue
        if run is None:
            bad += 1
            continue
        out.append((str(run), str(r.get("uid") or r.get("ion_uid") or ""), val))
    return out, bad


def run_means(rows):
    by_run = defaultdict(list)
    for run, _uid, ppm in rows:
        by_run[run].append(ppm)
    return {run: statistics.mean(v) for run, v in by_run.items()}, by_run


def one_sample_t(values):
    """Mean, 95% interval, t and df against a null of zero."""
    n = len(values)
    if n < 2:
        return None
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    if sd == 0:
        return {"mean": mean, "sd": 0.0, "n": n, "se": 0.0, "t": None,
                "df": n - 1, "lo": mean, "hi": mean,
                "note": "every value identical, so there is no scatter to test against"}
    se = sd / math.sqrt(n)
    tc = t_critical(n - 1)
    return {"mean": mean, "sd": sd, "n": n, "se": se, "t": mean / se,
            "df": n - 1, "lo": mean - tc * se, "hi": mean + tc * se, "note": None}


def sign_test(values):
    """How many share a sign, and the two-sided probability under a fair coin."""
    pos = sum(1 for v in values if v > 0)
    neg = sum(1 for v in values if v < 0)
    n = pos + neg
    if n == 0:
        return None
    k = max(pos, neg)
    # two-sided exact binomial tail at p = 0.5
    tail = sum(math.comb(n, i) for i in range(k, n + 1)) / (2 ** n)
    return {"pos": pos, "neg": neg, "n": n, "k": k, "p": min(1.0, 2 * tail)}


# --- the controls whose sensitivity is being measured ------------------------
# Each takes the per-run means and returns one number. The point of the table is
# which of them MOVE when a shared offset is added.

def c_between_sd(means):
    return statistics.stdev(means) if len(means) > 1 else 0.0


def c_pooled_within_sd(by_run):
    num = den = 0.0
    for v in by_run.values():
        if len(v) > 1:
            num += (len(v) - 1) * statistics.variance(v)
            den += len(v) - 1
    return math.sqrt(num / den) if den else 0.0


def c_range(means):
    return max(means) - min(means)


def c_iqr(means):
    if len(means) < 4:
        return 0.0
    q = statistics.quantiles(means, n=4)
    return q[2] - q[0]


def c_se_of_mean(means):
    return (statistics.stdev(means) / math.sqrt(len(means))) if len(means) > 1 else 0.0


def c_cv_of_spread(means):
    """Spread over |centre| -- the one spread-based control an offset DOES move,
    and it moves the wrong way: a bigger offset makes it look better."""
    m = statistics.mean(means)
    sd = statistics.stdev(means) if len(means) > 1 else 0.0
    return (sd / abs(m)) if m else float("inf")


def c_mean_of_means(means):
    return statistics.mean(means)


def c_share_one_sign(means):
    s = sign_test(means)
    return (s["k"] / s["n"]) if s and s["n"] else 0.0


CONTROLS = [
    ("between-run SD", c_between_sd, "spread"),
    ("pooled within-run SD", None, "spread"),      # needs by_run, handled apart
    ("range of run means", c_range, "spread"),
    ("IQR of run means", c_iqr, "spread"),
    ("SE of the Atlas mean", c_se_of_mean, "spread"),
    ("spread / |centre|", c_cv_of_spread, "spread"),
    ("mean of run means", c_mean_of_means, "offset"),
    ("share sharing a sign", c_share_one_sign, "offset"),
]


def sensitivity(means, by_run, out):
    """Which controls move when a known offset is added, and which do not."""
    out.append("")
    out.append("CONTROL SENSITIVITY TO A SHARED OFFSET")
    out.append("  Each control recomputed with a known offset added to every run.")
    out.append("  A control whose value is identical at 0 and at %g ppm is not"
               % PROBE_OFFSETS[-1])
    out.append("  tolerant of offsets, it is BLIND to them.")
    out.append("")
    head = "  %-22s %10s" + " %10s" * (len(PROBE_OFFSETS) - 1) + "   %s"
    out.append(head % (("control", "+0 ppm") + tuple("+%g" % o for o in PROBE_OFFSETS[1:])
                       + ("verdict",)))
    blind, seeing = [], []
    for name, fn, kind in CONTROLS:
        vals = []
        for off in PROBE_OFFSETS:
            shifted = [m + off for m in means]
            if fn is None:
                shifted_by_run = {r: [v + off for v in vs] for r, vs in by_run.items()}
                vals.append(c_pooled_within_sd(shifted_by_run))
            else:
                vals.append(fn(shifted))
        moved = max(abs(v - vals[0]) for v in vals[1:]) if len(vals) > 1 else 0.0
        if moved <= BLIND_TOL:
            verdict = "BLIND"
            blind.append(name)
        else:
            verdict = "sees it"
            seeing.append(name)
        cells = tuple("%10.4g" % v for v in vals)
        out.append("  %-22s" % name + "".join(" %10s" % c for c in cells)
                   + "   %s" % verdict)
    out.append("")
    out.append("  %d of %d control(s) are blind: %s"
               % (len(blind), len(CONTROLS), ", ".join(blind)))
    out.append("  That is not a defect in them. A spread is a difference between")
    out.append("  runs, and adding the same number to every run changes no")
    out.append("  difference. They are measuring the term that averages away.")
    if "spread / |centre|" in seeing:
        vals = [c_cv_of_spread([m + o for m in means]) for o in PROBE_OFFSETS]
        out.append("")
        out.append("  Note `spread / |centre|`: it moves, and it is the worst of the")
        out.append("  eight. Not because it falls -- it is NOT MONOTONE. The centre can")
        out.append("  cross zero, so the ratio blows up at the offset that happens to")
        out.append("  land the centre near zero and shrinks on either side of it: %s."
                   % ", ".join("%+g ppm -> %.3f" % (o, v)
                               for o, v in zip(PROBE_OFFSETS, vals)))
        out.append("  So it reports health as a U-shape in the thing it should be")
        out.append("  monotone in, and a large offset makes the cohort look MORE")
        out.append("  consistent than no offset at all. A control that is wrong in both")
        out.append("  directions is worse than one that is blind, because blindness at")
        out.append("  least does not argue.")
    return blind, seeing


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("ppm", help="per-PSM ppm errors (observed minus calculated)")
    p.add_argument("--anchor", help="an independent reference's ppm errors, if one exists")
    p.add_argument("--target", type=float, default=0.10,
                   help="the ppm the Atlas is aiming for (default 0.10)")
    args = p.parse_args(argv)

    out = []
    try:
        rows, bad = read_rows(args.ppm)
    except (OSError, ValueError) as exc:
        print("could not read %s: %s" % (args.ppm, exc), file=sys.stderr)
        return 2
    if bad:
        out.append("%d row(s) had no run or no numeric ppm and were skipped" % bad)
    if not rows:
        print("no usable rows in %s" % args.ppm, file=sys.stderr)
        return 2

    # The circularity that returned 0.00 ppm over 1,994,456 rows in this project.
    exact = [r for r in rows if abs(r[2]) < IDENTICAL]
    if exact:
        share = 100.0 * len(exact) / len(rows)
        out.append("")
        out.append("CIRCULARITY CHECK")
        out.append("  %d of %d row(s) (%.1f%%) have a ppm of zero to machine precision."
                   % (len(exact), len(rows), share))
        out.append("  That is arithmetic, not a measurement: the stored m/z was computed")
        out.append("  from the sequence it is being compared against. This project once")
        out.append("  reported a calibration of 0.00 ppm over 1,994,456 such rows.")
        if share > 50.0:
            out.append("  REFUSED: the majority of this input is its own answer. There is")
            out.append("  no mass error here to test for an offset.")
            print("\n".join(out))
            return 2
        out.append("  Excluded from everything below.")
        rows = [r for r in rows if abs(r[2]) >= IDENTICAL]

    means_by_run, by_run = run_means(rows)
    runs = sorted(means_by_run)
    means = [means_by_run[r] for r in runs]
    out.append("")
    out.append("INPUT  %d observation(s) over %d run(s)" % (len(rows), len(runs)))
    if len(runs) < 2:
        out.append("  One run, so between-run and shared cannot be separated at all:")
        out.append("  a single run's mean IS its offset and nothing distinguishes them.")
        print("\n".join(out))
        return 2

    within = c_pooled_within_sd(by_run)
    between = c_between_sd(means)
    out.append("  pooled within-run SD  %7.3f ppm" % within)
    out.append("  between-run SD        %7.3f ppm" % between)
    if within:
        out.append("  between / within      %7.3f" % (between / within))
        tot = between ** 2 + within ** 2
        out.append("  between-run share of variance %.1f%%" % (100 * between ** 2 / tot))

    # 1. Is the mean of the run means zero?
    out.append("")
    out.append("1  IS THE MEAN OF THE RUN MEANS ZERO?")
    t = one_sample_t(means)
    if t is None:
        out.append("  fewer than two runs.")
        return 2
    out.append("  mean of run means  %+.3f ppm   95%% CI %+.3f to %+.3f"
               % (t["mean"], t["lo"], t["hi"]))
    if t["t"] is None:
        out.append("  %s" % t["note"])
        offset_found = abs(t["mean"]) > 0
    else:
        out.append("  t = %+.2f on %d df,  SE %.3f" % (t["t"], t["df"], t["se"]))
        excludes = not (t["lo"] <= 0.0 <= t["hi"])
        offset_found = excludes
        if excludes:
            out.append("  The interval EXCLUDES zero. This is a shift, not scatter")
            out.append("  around a true value.")
        else:
            out.append("  The interval includes zero, so no shared offset is detected")
            out.append("  at this sample size. That is not the same as none existing:")
            out.append("  the detectable size is about %.3f ppm here." % (2 * t["se"]))

    # 2. Do the per-run means share a sign?
    out.append("")
    out.append("2  DO THE RUN MEANS SHARE A SIGN?")
    s = sign_test(means)
    if s is None or s["n"] == 0:
        out.append("  every run mean is exactly zero.")
    else:
        out.append("  %d of %d on one side (%d above, %d below),  p = %.3g"
                   % (s["k"], s["n"], s["pos"], s["neg"], s["p"]))
        if s["p"] < SIGN_ALPHA:
            out.append("  Assumption-free, and it catches exactly the failure that")
            out.append("  averaging hides. A fair coin does not do this.")
            offset_found = True
        else:
            out.append("  Consistent with scatter about zero.")

    # what more samples would and would not do
    if offset_found:
        out.append("")
        out.append("WHAT MORE SAMPLES WOULD DO")
        # A target of zero is not reachable at any n, and the old arithmetic
        # printed "0 run(s) reach the 0.00 ppm target" and then skipped the
        # explanation entirely -- on the exact command this repository's own
        # documentation gives.
        if args.target <= 0:
            need = 0
            out.append("  A target of %.2f ppm is not a target: no n makes"
                       " SE = %.3f/sqrt(n) reach zero." % (args.target, between))
            out.append("  Pass a positive --target to see how many runs it would")
            out.append("  take. The shift below does not shrink with n either way.")
        else:
            need = int(math.ceil((between / args.target) ** 2))
            out.append("  SE is %.3f/sqrt(n), so %d run(s) reach the %.2f ppm target."
                       % (between, need, args.target))
        if need > 1:
            se_at = between / math.sqrt(need)
            sigmas = abs(t["mean"]) / se_at if se_at else float("inf")
            out.append("  At %d runs the coordinate would be known to %.3f ppm and still"
                       % (need, se_at))
            out.append("  sit %+.3f ppm off -- an interval excluding the truth by about"
                       % t["mean"])
            out.append("  %.0f standard errors." % sigmas)
        out.append("  sqrt(n) shrinks the interval AROUND the shift, never toward zero.")
        out.append("  More samples make a shared offset more confidently wrong.")

    blind, seeing = sensitivity(means, by_run, out)

    # 3. Does an independent reference agree?
    out.append("")
    out.append("3  DOES AN INDEPENDENT REFERENCE AGREE?")
    if not args.anchor:
        out.append("  No --anchor given, so the SIZE of any offset above is unanchored.")
        out.append("  Checks 1 and 2 detect that an offset EXISTS. Neither can say how")
        out.append("  big it is, because both are computed from the same calculated")
        out.append("  m/z the offset is measured against: if that calculation is itself")
        out.append("  off by a constant, every number here moves with it and nothing")
        out.append("  internal notices.")
        out.append("  An anchor has to be independent of the mass computation -- a")
        out.append("  known-amount spike, a lock mass measured separately, or the same")
        out.append("  peptides on an instrument calibrated another way.")
    else:
        try:
            arows, abad = read_rows(args.anchor)
        except (OSError, ValueError) as exc:
            out.append("  could not read %s: %s" % (args.anchor, exc))
            arows = []
        if not arows:
            out.append("  the anchor file carried no usable rows.")
        else:
            ameans_by_run, aby_run = run_means(arows)
            ameans = [ameans_by_run[r] for r in sorted(ameans_by_run)]
            at = one_sample_t(ameans)
            shared = sorted(set(means_by_run) & set(ameans_by_run))
            out.append("  anchor: %d observation(s) over %d run(s), %d run(s) shared"
                       % (len(arows), len(ameans), len(shared)))
            if at:
                out.append("  anchor mean of run means %+.3f ppm  95%% CI %+.3f to %+.3f"
                           % (at["mean"], at["lo"], at["hi"]))
                gap = t["mean"] - at["mean"]
                out.append("  difference from the test set: %+.3f ppm" % gap)
                if not shared:
                    out.append("  NO RUN IS IN BOTH, so this compares two cohorts and not")
                    out.append("  one measurement against a reference. The difference")
                    out.append("  above mixes any cohort difference into the offset.")
                elif abs(gap) <= 2 * math.sqrt(t["se"] ** 2 + at["se"] ** 2):
                    out.append("  The two agree within their joint interval, so the offset")
                    out.append("  is in the mass computation or shared by both, not in one")
                    out.append("  of them. An anchor that shares the offset cannot find it.")
                else:
                    out.append("  They disagree by more than their joint interval, so the")
                    out.append("  offset is in one of the two and this pair locates it.")

    out.append("")
    out.append("VERDICT")
    if offset_found:
        out.append("  A shared offset is present. Every spread-based control above")
        out.append("  reports it as health, and %d of %d controls cannot see it at all."
                   % (len(blind), len(CONTROLS)))
        out.append("  Do not read a tightening interval on this dimension as"
                   " convergence.")
    else:
        out.append("  No shared offset detected at this sample size, by either check.")
        out.append("  The blind controls above are still blind; they simply have")
        out.append("  nothing to miss right now.")
    print("\n".join(out))
    return 1 if offset_found else 0


if __name__ == "__main__":
    sys.exit(main())
