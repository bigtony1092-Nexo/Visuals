#!/usr/bin/env python3
"""Prove the offset checks fire on a shift and stay quiet without one.

The claim being protected is that a control measuring spread cannot see a shared
offset. These cases plant offsets of known size and assert which controls move,
so the claim is a measurement rather than an argument -- and assert the reverse
too, because a check that fires on everything is no better than one that fires
on nothing.

    python3 tools/test_offset_check.py
"""

import json
import math
import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "offset_check.py"
sys.path.insert(0, str(ROOT / "tools"))
import offset_check as O  # noqa: E402

# the documented Colorado shape
OFFSET, BETWEEN, WITHIN, RUNS, PER_RUN = -1.570, 0.643, 2.277, 24, 53


def plant(path, offset=OFFSET, between=BETWEEN, within=WITHIN, runs=RUNS,
          per_run=PER_RUN, seed=42, computed=0):
    """A cohort with a known shared offset. `computed` rows are exact zeros."""
    rng = random.Random(seed)
    rows = []
    for i in range(runs):
        run = "CH%03d" % (163 + i)
        mean = offset + rng.gauss(0, between)
        for k in range(per_run):
            rows.append({"run": run, "uid": "ION_%03d_2" % k,
                         "ppm": mean + rng.gauss(0, within)})
    for k in range(computed):
        rows.append({"run": "CH163", "uid": "CALC_%d" % k, "ppm": 0.0})
    Path(path).write_text("".join(json.dumps(r) + "\n" for r in rows))
    return rows


def run(*args):
    p = subprocess.run([sys.executable, str(TOOL)] + [str(a) for a in args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def table(out):
    """The sensitivity table as {control: [values]}."""
    rows = {}
    started = False
    for line in out.splitlines():
        if "control" in line and "+0 ppm" in line:
            started = True
            continue
        if not started:
            continue
        if not line.strip() or line.strip().startswith(("5 of", "That is", "Note")):
            if rows:
                break
            continue
        # The verdict is "BLIND" (one token) or "sees it" (two), so splitting
        # a fixed number of fields from the right lands differently per row.
        # Cut the verdict off by name first.
        for verdict in ("BLIND", "sees it"):
            if line.rstrip().endswith(verdict):
                body = line.rstrip()[:-len(verdict)]
                break
        else:
            continue
        parts = body.rsplit(None, 4)
        if len(parts) != 5:
            continue
        name = parts[0].strip()
        try:
            rows[name] = [float(v) for v in parts[1:5]]
        except ValueError:
            continue
    return rows


def case_a_planted_offset_is_found(tmp):
    f = tmp / "offset.jsonl"
    plant(f)
    rc, out = run(f)
    bad = []
    if rc != 1:
        bad.append("exit %d on a cohort with a planted offset" % rc)
    if "interval EXCLUDES zero" not in out:
        bad.append("the t test did not exclude zero")
    line = [l for l in out.splitlines() if "mean of run means" in l and "CI" in l]
    if not line:
        return bad + ["no mean-of-means line"]
    got = float(line[0].split()[4])
    if abs(got - OFFSET) > 0.5:
        bad.append("recovered %+.3f ppm, planted %+.3f" % (got, OFFSET))
    if "24 of 24 on one side" not in out:
        bad.append("the sign test did not report 24 of 24: %r"
                   % [l.strip() for l in out.splitlines() if "one side" in l])
    if "1.19e-07" not in out:
        bad.append("sign-test p is not the documented 1.19e-07")
    return bad


def case_no_offset_is_not_reported_as_one(tmp):
    """The check must be able to stay quiet, or it says nothing when it fires."""
    f = tmp / "zero.jsonl"
    plant(f, offset=0.0, seed=7)
    rc, out = run(f)
    bad = []
    if rc != 0:
        bad.append("exit %d on a cohort centred at zero" % rc)
    if "interval EXCLUDES zero" in out:
        bad.append("found a shift in a cohort planted without one")
    if "No shared offset detected" not in out:
        bad.append("did not say that no offset was detected")
    if "still blind" not in out:
        bad.append("did not note the blind controls remain blind with nothing to miss")
    if "detectable size is about" not in out:
        bad.append("did not state what size it could have detected")
    return bad


def case_the_spread_controls_are_bit_identical(tmp):
    """The load-bearing measurement: adding a constant changes no difference."""
    f = tmp / "sens.jsonl"
    plant(f)
    rc, out = run(f)
    rows = table(out)
    bad = []
    expect_blind = ("between-run SD", "pooled within-run SD", "range of run means",
                    "IQR of run means", "SE of the Atlas mean")
    for name in expect_blind:
        if name not in rows:
            bad.append("%r missing from the sensitivity table" % name)
            continue
        vals = rows[name]
        if max(abs(v - vals[0]) for v in vals[1:]) > 1e-9:
            bad.append("%r moved under a 20 ppm offset: %r -- the claim that spread"
                       " cannot see an offset would be false" % (name, vals))
    for name in ("mean of run means", "share sharing a sign"):
        if name in rows and max(abs(v - rows[name][0]) for v in rows[name][1:]) <= 1e-9:
            bad.append("%r did NOT move under a 20 ppm offset, so it cannot be the"
                       " control that catches one" % name)
    if "5 of 8 control(s) are blind" not in out:
        bad.append("the blind count is not 5 of 8: %r"
                   % [l.strip() for l in out.splitlines() if "are blind" in l])
    if "changes no" not in out:
        bad.append("did not explain why a spread cannot see a shared offset")
    return bad


def case_the_ratio_control_is_non_monotone(tmp):
    """spread/|centre| is wrong in both directions, which is worse than blind."""
    f = tmp / "cv.jsonl"
    plant(f)
    rc, out = run(f)
    rows = table(out)
    bad = []
    vals = rows.get("spread / |centre|")
    if not vals:
        return ["spread / |centre| missing from the table"]
    # it must rise then fall: the centre crosses zero between +0 and +5
    if not (vals[1] > vals[0] and vals[3] < vals[0]):
        bad.append("values %r are monotone; the non-monotonicity claim would be"
                   " wrong" % vals)
    if "NOT MONOTONE" not in out:
        bad.append("did not name the non-monotonicity")
    if "worse than one that is blind" not in out:
        bad.append("did not rank it against the blind controls")
    return bad


def case_computed_rows_are_refused(tmp):
    """Observed minus calculated is zero when the m/z came from the sequence."""
    f = tmp / "circular.jsonl"
    # a minority: excluded with a note
    plant(f, computed=100)
    rc, out = run(f)
    bad = []
    if "CIRCULARITY CHECK" not in out:
        bad.append("exact zeros were not flagged")
    if "1,994,456" not in out:
        bad.append("did not cite the project's own instance of this")
    if "Excluded from everything below" not in out:
        bad.append("did not say the rows were excluded")
    # a majority: refused outright
    g = tmp / "mostly_circular.jsonl"
    plant(g, runs=4, per_run=5, computed=500)
    rc2, out2 = run(g)
    if rc2 != 2:
        bad.append("exit %d when most rows are their own answer" % rc2)
    if "REFUSED" not in out2:
        bad.append("did not refuse an input that is mostly arithmetic")
    if "no mass error here" not in out2:
        bad.append("did not say why a majority is fatal")
    return bad


def case_one_run_cannot_separate_offset_from_mean(tmp):
    f = tmp / "onerun.jsonl"
    plant(f, runs=1, per_run=40)
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on a single run" % rc)
    if "a single run's mean IS its offset" not in out:
        bad.append("did not say why one run cannot answer this")
    return bad


def case_more_samples_are_shown_to_make_it_worse(tmp):
    f = tmp / "conv.jsonl"
    plant(f)
    rc, out = run(f, "--target", "0.10")
    bad = []
    if "WHAT MORE SAMPLES WOULD DO" not in out:
        bad.append("did not project what sample count buys")
    if "standard errors" not in out:
        bad.append("did not express the error in standard errors at the target")
    if "never toward zero" not in out:
        bad.append("did not say sqrt(n) shrinks around the shift")
    if "more confidently wrong" not in out:
        bad.append("did not state the conclusion plainly")
    line = [l for l in out.splitlines() if "run(s) reach the" in l]
    if line:
        m = re.search(r"so (\d+) run\(s\) reach", line[0])
        n = int(m.group(1)) if m else 0
        # between-run SD near 0.87 here, so (0.87/0.1)^2 is about 75
        if not (40 <= n <= 200):
            bad.append("projected %d runs for a 0.10 ppm target, which does not"
                       " match the between-run SD" % n)
    return bad


def case_an_anchor_that_shares_the_offset_cannot_find_it(tmp):
    """The subtle case: a reference carrying the same shift agrees, and agreeing
    is not evidence of correctness."""
    f = tmp / "test.jsonl"
    a = tmp / "anchor_same.jsonl"
    plant(f)
    plant(a, seed=99)                    # same planted offset
    rc, out = run(f, "--anchor", a)
    bad = []
    if "anchor mean of run means" not in out:
        bad.append("the anchor was not summarised")
    if "cannot find it" not in out:
        bad.append("an anchor sharing the offset was not called out")
    if "in the mass computation or shared by both" not in out:
        bad.append("did not say where a shared offset would then live")
    # and an anchor that does NOT share it must locate the offset
    b = tmp / "anchor_clean.jsonl"
    plant(b, offset=0.0, seed=11)
    rc2, out2 = run(f, "--anchor", b)
    if "this pair locates it" not in out2:
        bad.append("an anchor without the offset did not locate it")
    return bad


def case_an_anchor_with_no_shared_run_is_two_cohorts(tmp):
    f = tmp / "t2.jsonl"
    a = tmp / "a2.jsonl"
    plant(f)
    rows = plant(a, offset=0.0, seed=5)
    # rename every anchor run so none is shared
    Path(a).write_text("".join(
        json.dumps({**r, "run": "OTHER_" + r["run"]}) + "\n" for r in rows))
    rc, out = run(f, "--anchor", a)
    bad = []
    if "NO RUN IS IN BOTH" not in out:
        bad.append("did not notice the anchor shares no run with the test set")
    if "mixes any cohort difference" not in out:
        bad.append("did not say what that costs")
    return bad


def case_the_sign_test_arithmetic_is_right(tmp):
    bad = []
    for n, k, want in ((24, 24, 1.1920928955078125e-07), (10, 10, 0.001953125),
                       (10, 5, 1.0), (4, 3, 0.625)):
        vals = [1.0] * k + [-1.0] * (n - k)
        s = O.sign_test(vals)
        if abs(s["p"] - want) > 1e-9:
            bad.append("sign test %d of %d gave p=%.6g, expected %.6g"
                       % (k, n, s["p"], want))
    if O.sign_test([0.0, 0.0]) is not None:
        bad.append("all-zero input produced a sign test")
    return bad


CASES = [case_a_planted_offset_is_found,
         case_no_offset_is_not_reported_as_one,
         case_the_spread_controls_are_bit_identical,
         case_the_ratio_control_is_non_monotone,
         case_computed_rows_are_refused,
         case_one_run_cannot_separate_offset_from_mean,
         case_more_samples_are_shown_to_make_it_worse,
         case_an_anchor_that_shares_the_offset_cannot_find_it,
         case_an_anchor_with_no_shared_run_is_two_cohorts,
         case_the_sign_test_arithmetic_is_right]


def main():
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for case in CASES:
            try:
                bad = case(tmp)
            except Exception as exc:
                bad = ["%s: %s" % (type(exc).__name__, exc)]
            failed += bool(bad)
            print("%s %s" % ("ok  " if not bad else "FAIL", case.__name__))
            for line in bad:
                print("     %s" % line)
    print("\n%d/%d passed" % (len(CASES) - failed, len(CASES)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
