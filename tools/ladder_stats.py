#!/usr/bin/env python3
"""Confidence intervals and t-tests for the real-vs-reversed ladder comparison.

Two jobs, kept strictly apart.

ONE: put intervals on the replication numbers, which have never carried any. A
found-again rate of 4/93 is reported as 4.30%; its Wilson 95% interval is
1.69-10.54%, and the 13.5x reduction it feeds has an interval of 5.0-36.4x. Four
events cannot support three significant figures, and a slide that shows the point
estimate alone reads as far more precise than the data is.

TWO: test a candidate variable WITHOUT letting it score anything. A new per-ion
quantity (--covariate) is run through a paired t-test between each real ion and
its own reversed counterpart, with the confidence interval on the mean difference
and an effect size. It is never thresholded, never summed into a score, and
never used to include or exclude an ion. This tool has no code path that can
score with it: see `score_with_covariate`, which raises.

    python3 tools/ladder_stats.py counts.json
    python3 tools/ladder_stats.py counts.json --covariate var1
    python3 tools/ladder_stats.py counts.json --check-symmetry

Input is one JSON object:

    {"gates": [{"name": "5+ fragments",
                "real_total": 643, "real_again": 373,
                "rev_total": 93,  "rev_again": 4}],
     "ions": [{"uid": "...", "sequence": "...", "charge": 2,
               "real_ladder": 21, "rev_ladder": 19,
               "real": {"var1": 0.42}, "rev": {"var1": 0.11}}]}

`gates` drives job one. `ions` drives the symmetry check and the t-tests; each
needs `real` and `rev` sub-objects holding the same covariate names.

Exit status is 1 if the ladders are asymmetric, because that biases every
count-gated comparison above and must be fixed before the intervals mean
anything.
"""

import argparse
import json
import math
import sys

Z = 1.959964          # two-sided 95%


def wilson(k, n, z=Z):
    """Score interval on a proportion. Correct at k=0 and k=n, unlike Wald."""
    if n <= 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    centre = p + z * z / (2 * n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    lo, hi = (centre - half) / d, (centre + half) / d
    # At k=0 and k=n the exact bounds are 0 and 1, but the arithmetic lands a
    # few ulps off and a bound printed as 0.9999999999999999 reads as an error.
    return (0.0 if lo < 1e-12 else lo), (1.0 if hi > 1 - 1e-12 else hi)


def katz_ratio(k1, n1, k2, n2, z=Z):
    """95% interval on a ratio of two proportions, on the log scale.

    Undefined when either numerator is zero, which is the honest answer: a ratio
    against zero events has no upper bound.
    """
    if k1 == 0 or k2 == 0 or n1 == 0 or n2 == 0:
        return None
    p1, p2 = k1 / n1, k2 / n2
    se = math.sqrt((1 - p1) / k1 + (1 - p2) / k2)
    centre = math.log(p1 / p2)
    return math.exp(centre - z * se), math.exp(centre + z * se)


def t_critical(df):
    """Two-sided 95% t, good to about 3 decimals over the range that matters."""
    if df <= 0:
        return float("nan")
    table = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
             7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 12: 2.179, 15: 2.131,
             20: 2.086, 25: 2.060, 30: 2.042, 40: 2.021, 60: 2.000, 120: 1.980}
    if df in table:
        return table[df]
    keys = sorted(table)
    if df > keys[-1]:
        return Z
    lo = max(k for k in keys if k < df)
    hi = min(k for k in keys if k > df)
    w = (df - lo) / (hi - lo)
    return table[lo] + w * (table[hi] - table[lo])


def paired_t(diffs):
    """Mean difference, its 95% interval, t, df and Cohen's d_z."""
    n = len(diffs)
    if n < 2:
        return None
    mean = sum(diffs) / n
    var = sum((x - mean) ** 2 for x in diffs) / (n - 1)
    sd = math.sqrt(var)
    if sd == 0:
        return {"n": n, "mean": mean, "lo": mean, "hi": mean, "t": float("inf"),
                "df": n - 1, "d": float("inf"), "sd": 0.0}
    se = sd / math.sqrt(n)
    tc = t_critical(n - 1)
    return {"n": n, "mean": mean, "lo": mean - tc * se, "hi": mean + tc * se,
            "t": mean / se, "df": n - 1, "d": mean / sd, "sd": sd}


def score_with_covariate(*_args, **_kwargs):
    """Deliberately absent.

    A covariate is tested here and scored nowhere. Keeping the refusal as code
    rather than as a comment means a later caller that tries to score with it
    fails loudly instead of quietly acquiring a new term.
    """
    raise NotImplementedError(
        "a covariate is tested here, never scored. Establish its interval and "
        "effect size first, decide separately whether it earns a place in "
        "scoring, and implement that where scoring lives.")


def report_gates(gates, out):
    out.append("INTERVALS ON THE FOUND-AGAIN RATES")
    out.append("  Wilson 95% on each rate; the ratio's interval on the log scale.")
    out.append("")
    out.append("  %-16s %-24s %-24s %s" % ("gate", "real found again",
                                           "reversed found again", "ratio (95% CI)"))
    for g in gates:
        rn, rk = g["real_total"], g["real_again"]
        dn, dk = g["rev_total"], g["rev_again"]
        rlo, rhi = wilson(rk, rn)
        dlo, dhi = wilson(dk, dn)
        ci = katz_ratio(rk, rn, dk, dn)
        ratio = (rk / rn) / (dk / dn) if dk and dn and rn else float("inf")
        cis = "%.1f-%.1f" % ci if ci else "no upper bound (zero events)"
        out.append("  %-16s %4d/%-5d %5.1f%% (%.1f-%.1f)  %4d/%-5d %5.2f%% (%.2f-%.2f)  %5.2fx (%s)"
                   % (g.get("name", "?"), rk, rn, 100 * rk / rn, 100 * rlo, 100 * rhi,
                      dk, dn, 100 * dk / dn, 100 * dlo, 100 * dhi, ratio, cis))
        if dk and dk < 10:
            out.append("      %d reversed events is too few for three significant figures."
                       % dk)
            out.append("      Quote the interval, not the point estimate.")
    out.append("")


def check_symmetry(ions, out):
    """The reversed ladder must be the same length as the real one.

    A reversed ladder that is systematically shorter loses ions at any gate
    defined as an absolute fragment count, and the loss concentrates exactly at
    the gate, so a few per cent of trimming costs far more than a few per cent of
    decoys. That makes every rate above a floor rather than an estimate.
    """
    pairs = [(i["real_ladder"], i["rev_ladder"]) for i in ions
             if isinstance(i.get("real_ladder"), int) and isinstance(i.get("rev_ladder"), int)
             and i["real_ladder"] > 0]
    out.append("LADDER SYMMETRY")
    if not pairs:
        out.append("  SKIPPED: no ion carries both real_ladder and rev_ladder counts.")
        return False
    short = sum(1 for r, v in pairs if v < r)
    longer = sum(1 for r, v in pairs if v > r)
    ratios = sorted(v / r for r, v in pairs)
    med = ratios[len(ratios) // 2]
    mean = sum(ratios) / len(ratios)
    out.append("  %d ions. reversed/real length: median %.3f, mean %.3f"
               % (len(pairs), med, mean))
    out.append("  reversed shorter in %d, longer in %d, equal in %d"
               % (short, longer, len(pairs) - short - longer))
    if abs(1 - med) < 0.005 and short + longer < 0.1 * len(pairs):
        out.append("  SYMMETRIC. Both ladders are built the same way, so a count gate")
        out.append("  treats them alike and the rates above are estimates, not floors.")
        return True
    out.append("")
    out.append("  ASYMMETRIC — the reversed ladder is %.1f%% %s than the real one."
               % (abs(1 - med) * 100, "shorter" if med < 1 else "longer"))
    out.append("  This biases every gate defined as an absolute fragment count, and the")
    out.append("  loss lands exactly at the gate: ions sitting one or two fragments above")
    out.append("  it fall below. So a few per cent of trimming costs far more than a few")
    out.append("  per cent of decoys, and every rate above is a %s rather than an estimate."
               % ("floor" if med < 1 else "ceiling"))
    out.append("")
    out.append("  THE FIX, in the ladder generator: build both ladders by the same rule.")
    out.append("  Either apply every trim to both (reversed b(L-1), exact coincidences,")
    out.append("  whatever else) or apply none to either. Partial trimming of one side")
    out.append("  only is the one option that cannot be defended, because it hands the")
    out.append("  decoy a handicap and then reports the margin as evidence.")
    return False


def report_covariate(ions, name, out):
    """Paired t-test on one candidate variable. Tested, never scored."""
    diffs, have = [], 0
    for i in ions:
        a = (i.get("real") or {}).get(name)
        b = (i.get("rev") or {}).get(name)
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            diffs.append(a - b)
            have += 1
    out.append("CANDIDATE VARIABLE %r — TESTED, NOT SCORED" % name)
    if have < 2:
        out.append("  SKIPPED: %d ions carry %r on both the real and reversed side."
                   % (have, name))
        out.append("  Each ion needs real.%s and rev.%s to be paired." % (name, name))
        return
    r = paired_t(diffs)
    out.append("  paired over %d ions (real minus its own reversed counterpart)" % r["n"])
    out.append("  mean difference  %+.4f   95%% CI %+.4f to %+.4f"
               % (r["mean"], r["lo"], r["hi"]))
    out.append("  t(%d) = %+.3f   Cohen's d_z = %+.3f   sd of differences %.4f"
               % (r["df"], r["t"], r["d"], r["sd"]))
    crosses = r["lo"] <= 0 <= r["hi"]
    out.append("")
    if crosses:
        out.append("  The interval includes zero, so this variable does not separate real")
        out.append("  from reversed in this data. It has not earned a place in scoring.")
    else:
        out.append("  The interval excludes zero, so it does separate real from reversed")
        out.append("  here. That is necessary and not sufficient: before it scores anything,")
        out.append("  establish that the separation is not an artefact of the comparison")
        out.append("  (see the four flipped badges in docs/verification-run.md), and that it")
        out.append("  adds something the existing terms do not already carry.")
    out.append("  Nothing in this tool applies a threshold to %r or sums it into a score."
               % name)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("counts", help="JSON file with `gates` and/or `ions`")
    p.add_argument("--covariate", action="append", default=[],
                   help="a per-ion variable to t-test; repeatable; never scored")
    p.add_argument("--check-symmetry", action="store_true",
                   help="only run the ladder symmetry check")
    args = p.parse_args(argv)

    with open(args.counts) as fh:
        data = json.load(fh)
    out = []
    symmetric = True
    ions = data.get("ions") or []

    if args.check_symmetry:
        symmetric = check_symmetry(ions, out)
        print("\n".join(out))
        return 0 if symmetric else 1

    if data.get("gates"):
        report_gates(data["gates"], out)
    if ions:
        symmetric = check_symmetry(ions, out)
        out.append("")
    for name in args.covariate:
        report_covariate(ions, name, out)
        out.append("")
    if not data.get("gates") and not ions:
        print("Nothing to do: the file has neither `gates` nor `ions`.", file=sys.stderr)
        return 1
    print("\n".join(out).rstrip())
    return 0 if symmetric else 1


if __name__ == "__main__":
    sys.exit(main())
