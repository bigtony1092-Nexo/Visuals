#!/usr/bin/env python3
"""Prove the chain stops at a gate, and that "not reached" never reads as "passed".

The whole value of chaining is the ordering: a flat fragment reference makes a
later pattern comparison unreadable, so the chain must stop rather than print a
number from it. And a step that did not run must be visibly distinct from one
that passed — the five failures in docs/verification-run.md were all numbers
that looked like results.

    python3 tools/test_aria_check.py
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "aria_check.py"


def run(*args):
    p = subprocess.run([sys.executable, str(TOOL)] + list(args),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def flat_fragments(path):
    """A source with no proline enrichment: the chemistry gate must refuse it."""
    rows = []
    for seq in ("AAGSTVLKP", "SDEGTVAPR"):
        for kind in "by":
            for n in range(1, len(seq)):
                rows.append({"source": "flat", "sequence": seq, "charge": 2,
                             "fragment": "%s%d" % (kind, n), "intensity": 1.0})
    path.write_text("\n".join(json.dumps(r) for r in rows))


def good_fragments(path):
    rows = []
    for seq in ("AAGSTVLKP", "SDEGTVAPR", "GGSPTTAAK"):
        L = len(seq)
        for kind in "by":
            for n in range(1, L):
                k = n if kind == "b" else L - n
                inten = 4.0 if seq[k] == "P" else 1.0
                rows.append({"source": "real", "sequence": seq, "charge": 2,
                             "fragment": "%s%d" % (kind, n), "intensity": inten})
    path.write_text("\n".join(json.dumps(r) for r in rows))


def ions(path, n=6):
    rows = []
    for i in range(n):
        for j, jit in enumerate((-0.4, 0.0, 0.4)):
            rows.append({"run": "r%d" % j, "uid": "P%d|2" % i, "mz": 500.0 + 60 * i,
                         "z": 2, "drift": 40.0 + 12.0 * i + jit,
                         "rt": 10.0 + 3.0 * i + 0.05 * jit})
    path.write_text("\n".join(json.dumps(r) for r in rows))


def asymmetric_counts(path):
    path.write_text(json.dumps({
        "gates": [{"name": "5+", "real_total": 643, "real_again": 373,
                   "rev_total": 93, "rev_again": 4}],
        "ions": [{"uid": "i%d" % n, "real_ladder": 20, "rev_ladder": 18}
                 for n in range(30)]}))


def symmetric_counts(path):
    path.write_text(json.dumps({
        "gates": [{"name": "5+", "real_total": 643, "real_again": 373,
                   "rev_total": 93, "rev_again": 4}],
        "ions": [{"uid": "i%d" % n, "real_ladder": 20, "rev_ladder": 20}
                 for n in range(30)]}))


def a_flat_reference_stops_the_chain(tmp):
    f, i, cnt = tmp / "flat.jsonl", tmp / "ions.json", tmp / "counts.json"
    flat_fragments(f); ions(i); symmetric_counts(cnt)
    rc, out = run("--fragments", str(f), "--ions", str(i), "--counts", str(cnt),
                  "--rank-rt")
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1 when a gate failed" % rc)
    if "GATE FAILED" not in out:
        bad.append("did not mark the chemistry step as a failed gate")
    if "NOT REACHED" not in out:
        bad.append("later steps were not marked NOT REACHED")
    if "Not reached is not passed" not in out:
        bad.append("did not say that not reached is not passed")
    # The steps after the gate must not have produced numbers.
    if "uniquely resolved" in out:
        bad.append("ran the resolvability step past a failed gate")
    if "INTERVALS ON THE FOUND-AGAIN RATES" in out:
        bad.append("ran the stats step past a failed gate")
    return bad


def a_good_reference_lets_the_chain_through(tmp):
    f, i, cnt = tmp / "good.jsonl", tmp / "ions2.json", tmp / "sym.json"
    good_fragments(f); ions(i); symmetric_counts(cnt)
    # A stable elution order leaves rank with no scatter, so the tolerance is
    # chosen rather than measured: one rank position out of 100.
    rc, out = run("--fragments", str(f), "--ions", str(i), "--counts", str(cnt),
                  "--rank-rt", "--sigma-rt", "1.0")
    bad = []
    if "GATE FAILED" in out:
        bad.append("refused a source that shows the proline effect")
    if "NOT REACHED" in out:
        bad.append("stopped the chain on a passing gate")
    for want in ("uniquely resolved", "INTERVALS ON THE FOUND-AGAIN RATES", "SYMMETRIC"):
        if want not in out:
            bad.append("later step did not run: %r missing" % want)
    if rc != 0:
        bad.append("exit %d with every step passing" % rc)
    return bad


def asymmetric_ladders_are_a_gate_too(tmp):
    f, i, cnt = tmp / "good2.jsonl", tmp / "ions3.json", tmp / "asym.json"
    good_fragments(f); ions(i); asymmetric_counts(cnt)
    rc, out = run("--fragments", str(f), "--ions", str(i), "--counts", str(cnt),
                  "--rank-rt")
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1 on asymmetric ladders" % rc)
    if "GATE FAILED" not in out:
        bad.append("asymmetric ladders did not fail their gate")
    if "every rate above is a floor" not in out:
        bad.append("did not say why the intervals cannot be taken at face value")
    return bad


def a_skipped_step_is_not_a_passed_one(tmp):
    f = tmp / "good3.jsonl"
    good_fragments(f)
    rc, out = run("--fragments", str(f))
    bad = []
    if "skipped" not in out:
        bad.append("did not mark the steps with no inputs as skipped")
    if "A skipped step has not passed" not in out:
        bad.append("did not distinguish skipped from passed")
    if rc != 0:
        bad.append("exit %d: a skipped step should not fail the run" % rc)
    return bad


def the_shipped_examples_run_end_to_end(tmp):
    rc, out = run("--all-of", str(ROOT / "examples"))
    bad = []
    if "COULD NOT RUN" in out:
        bad.append("a tool could not be run on the shipped examples")
    if "1 EXPORT" not in out or "2 PANELS" not in out:
        bad.append("did not discover the example manifest and method")
    if "error(s)" in out and " 0 error(s)" not in out:
        bad.append("the shipped examples reported errors")
    return bad



def ppm_with_a_shared_offset(path, offset=-1.570, runs=24, per_run=53, seed=42):
    import random as _r
    rng = _r.Random(seed)
    rows = []
    for i in range(runs):
        mean = offset + rng.gauss(0, 0.643)
        for k in range(per_run):
            rows.append({"run": "CH%03d" % (163 + i), "uid": "ION_%03d_2" % k,
                         "ppm": mean + rng.gauss(0, 2.277)})
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))


def a_shared_offset_stops_the_chain(tmp):
    """Step 6 is a gate because nothing above it can see what it looks for."""
    f, p = tmp / "goodoff.jsonl", tmp / "ppm_off.jsonl"
    good_fragments(f)
    ppm_with_a_shared_offset(p)
    rc, out = run("--fragments", str(f), "--ppm", str(p))
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1 on a cohort with a shared offset" % rc)
    if "6 OFFSET" not in out:
        bad.append("step 6 did not appear")
    if "GATE FAILED" not in out:
        bad.append("a shared offset did not fail its gate")
    if "more confidently wrong" not in out:
        bad.append("did not say why a shared offset is a gate rather than a note")
    if "24 of 24 on one side" not in out:
        bad.append("the sign test did not run")
    return bad


def a_centred_cohort_passes_step_six(tmp):
    """The gate has to be able to pass, or it is not measuring anything."""
    f, p = tmp / "goodon.jsonl", tmp / "ppm_zero.jsonl"
    good_fragments(f)
    ppm_with_a_shared_offset(p, offset=0.0, seed=7)
    rc, out = run("--fragments", str(f), "--ppm", str(p))
    bad = []
    if rc != 0:
        bad.append("exit %d on a cohort centred at zero" % rc)
    if "GATE FAILED" in out:
        bad.append("failed the offset gate on a cohort planted without one")
    if "No shared offset detected" not in out:
        bad.append("did not report the absence of an offset")
    return bad


def step_six_is_skipped_not_assumed(tmp):
    f = tmp / "goodskip.jsonl"
    good_fragments(f)
    rc, out = run("--fragments", str(f))
    bad = []
    if "6 OFFSET" not in out:
        bad.append("step 6 vanished when it had no input")
    line = [l for l in out.splitlines() if "6 OFFSET" in l and "skipped" in l]
    if not line:
        bad.append("step 6 with no --ppm was not marked skipped: %r"
                   % [l.strip() for l in out.splitlines() if "6 OFFSET" in l])
    if "A skipped step has not passed" not in out:
        bad.append("did not distinguish a skipped offset check from a passed one")
    return bad


CASES = [a_flat_reference_stops_the_chain, a_good_reference_lets_the_chain_through,
         asymmetric_ladders_are_a_gate_too, a_skipped_step_is_not_a_passed_one,
         the_shipped_examples_run_end_to_end,
         a_shared_offset_stops_the_chain,
         a_centred_cohort_passes_step_six,
         step_six_is_skipped_not_assumed]


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
