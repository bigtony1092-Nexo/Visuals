#!/usr/bin/env python3
"""Prove the intervals are right and the covariate cannot score.

The intervals are checked against values computed by hand, because a wrong
confidence interval is worse than none: it looks like rigour. The covariate path
is checked to be a dead end by design.

    python3 tools/test_ladder_stats.py
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "ladder_stats.py"
sys.path.insert(0, str(ROOT / "tools"))
import ladder_stats as ls  # noqa: E402

REPLICATION = {"gates": [
    {"name": "3+ fragments", "real_total": 2190, "real_again": 922,
     "rev_total": 1181, "rev_again": 152},
    {"name": "5+ fragments", "real_total": 643, "real_again": 373,
     "rev_total": 93, "rev_again": 4}]}


def run(path, *extra):
    p = subprocess.run([sys.executable, str(TOOL), str(path)] + list(extra),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def write(tmp, obj, name="counts.json"):
    f = tmp / name
    f.write_text(json.dumps(obj))
    return f


def wilson_matches_known_values(tmp):
    """Checked against the standard score interval, including the edges."""
    bad = []
    for k, n, want_lo, want_hi in ((4, 93, 0.0169, 0.1054),
                                   (373, 643, 0.5417, 0.6176),
                                   (93, 643, 0.1196, 0.1738)):
        lo, hi = ls.wilson(k, n)
        if abs(lo - want_lo) > 0.001 or abs(hi - want_hi) > 0.001:
            bad.append("wilson(%d,%d) = %.4f-%.4f, expected %.4f-%.4f"
                       % (k, n, lo, hi, want_lo, want_hi))
    # Wald would give a negative lower bound at k=0; Wilson must not.
    lo, hi = ls.wilson(0, 50)
    if lo != 0.0:
        bad.append("wilson(0,50) lower bound %.4f, must be 0" % lo)
    if not 0.0 < hi < 0.15:
        bad.append("wilson(0,50) upper bound %.4f is implausible" % hi)
    lo, hi = ls.wilson(50, 50)
    if hi != 1.0:
        bad.append("wilson(50,50) upper bound %.4f, must be 1" % hi)
    return bad


def ratio_interval_is_wide_at_four_events(tmp):
    """The whole point: 13.5x is not 13.5x."""
    ci = ls.katz_ratio(93, 643, 4, 373)
    bad = []
    if ci is None:
        return ["no interval returned for 93/643 against 4/373"]
    lo, hi = ci
    if not (4.5 < lo < 5.5):
        bad.append("lower bound %.2f, expected about 5.0" % lo)
    if not (33 < hi < 40):
        bad.append("upper bound %.2f, expected about 36" % hi)
    if hi / lo < 5:
        bad.append("interval spans only %.1fx; four events cannot be that precise"
                   % (hi / lo))
    # A ratio against zero events has no upper bound and must say so.
    if ls.katz_ratio(93, 643, 0, 373) is not None:
        bad.append("returned a finite interval for a ratio against zero events")
    return bad


def the_slide_numbers_get_intervals(tmp):
    rc, out = run(write(tmp, REPLICATION))
    bad = []
    if rc != 0:
        bad.append("exit %d with no ions to check" % rc)
    if "13.49x" not in out and "13.5" not in out:
        bad.append("did not report the point ratio")
    if "5.2-35.3" not in out and "5.0-36.4" not in out and "(5." not in out:
        bad.append("did not report the ratio interval: %r" % out[:200])
    if "too few for three significant figures" not in out:
        bad.append("did not warn that 4 events cannot carry the point estimate")
    return bad


def asymmetric_ladders_fail_the_run(tmp):
    """A 6.5% trim on the reversed side must fail, and name the fix."""
    ions = [{"uid": "i%d" % n, "real_ladder": 20,
             "rev_ladder": 19 if n % 2 else 18} for n in range(40)]
    rc, out = run(write(tmp, {"ions": ions}, "asym.json"))
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1 for asymmetric ladders" % rc)
    if "ASYMMETRIC" not in out:
        bad.append("did not flag the asymmetry")
    if "floor" not in out:
        bad.append("did not say the rates become floors")
    if "same rule" not in out and "both" not in out:
        bad.append("did not state the fix")
    return bad


def symmetric_ladders_pass(tmp):
    ions = [{"uid": "i%d" % n, "real_ladder": 20, "rev_ladder": 20} for n in range(40)]
    rc, out = run(write(tmp, {"ions": ions}, "sym.json"))
    bad = []
    if rc != 0:
        bad.append("exit %d, expected 0 for symmetric ladders" % rc)
    if "SYMMETRIC" not in out:
        bad.append("did not credit symmetry")
    if "estimates, not floors" not in out:
        bad.append("did not say the rates become estimates")
    return bad


def a_covariate_that_separates_is_reported_not_scored(tmp):
    ions = [{"uid": "i%d" % n, "real_ladder": 20, "rev_ladder": 20,
             "real": {"var1": 0.6 + 0.01 * n}, "rev": {"var1": 0.2 + 0.01 * n}}
            for n in range(30)]
    rc, out = run(write(tmp, {"ions": ions}, "cov.json"), "--covariate", "var1")
    bad = []
    if "TESTED, NOT SCORED" not in out:
        bad.append("did not label the covariate as tested rather than scored")
    if "interval excludes zero" not in out:
        bad.append("did not detect a clear separation")
    if "necessary and not sufficient" not in out:
        bad.append("treated separation alone as enough to justify scoring")
    if "Nothing in this tool applies a threshold" not in out:
        bad.append("did not state that it never thresholds the covariate")
    return bad


def a_covariate_that_does_not_separate_says_so(tmp):
    rng = [0.3, -0.2, 0.1, -0.1, 0.25, -0.3, 0.05, -0.05, 0.2, -0.15]
    ions = [{"uid": "i%d" % n, "real_ladder": 20, "rev_ladder": 20,
             "real": {"var1": 0.5 + d}, "rev": {"var1": 0.5}}
            for n, d in enumerate(rng)]
    rc, out = run(write(tmp, {"ions": ions}, "cov0.json"), "--covariate", "var1")
    bad = []
    if "interval includes zero" not in out:
        bad.append("claimed separation where the interval crosses zero")
    if "has not earned a place in scoring" not in out:
        bad.append("did not say it has not earned scoring")
    return bad


def scoring_with_a_covariate_raises(tmp):
    """The refusal is code, not a comment, so a later caller fails loudly."""
    try:
        ls.score_with_covariate({"var1": 1.0})
    except NotImplementedError as exc:
        if "never scored" not in str(exc):
            return ["raised, but without saying why: %s" % exc]
        return []
    return ["score_with_covariate did not raise"]


def a_missing_covariate_is_named(tmp):
    ions = [{"uid": "i1", "real_ladder": 20, "rev_ladder": 20,
             "real": {"other": 1.0}, "rev": {"other": 2.0}}]
    rc, out = run(write(tmp, {"ions": ions}, "miss.json"), "--covariate", "var1")
    bad = []
    if "SKIPPED" not in out:
        bad.append("did not skip a covariate absent from the data")
    if "real.var1" not in out:
        bad.append("did not name the fields it needs")
    return bad


CASES = [wilson_matches_known_values, ratio_interval_is_wide_at_four_events,
         the_slide_numbers_get_intervals, asymmetric_ladders_fail_the_run,
         symmetric_ladders_pass, a_covariate_that_separates_is_reported_not_scored,
         a_covariate_that_does_not_separate_says_so, scoring_with_a_covariate_raises,
         a_missing_covariate_is_named]


def main():
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for case in CASES:
            bad = case(tmp)
            failed += bool(bad)
            print("%s %s" % ("ok  " if not bad else "FAIL", case.__name__))
            for line in bad:
                print("     %s" % line)
    print("\n%d/%d passed" % (len(CASES) - failed, len(CASES)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
