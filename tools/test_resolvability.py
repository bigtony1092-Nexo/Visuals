#!/usr/bin/env python3
"""Prove the separation is measured in sigma, never in absolute units.

The point this tool exists to enforce: a distance means nothing without the
measurement scale beside it. The same 7-unit gap must read as resolved on a
precise instrument and unresolved on a coarse one, from identical coordinates.

    python3 tools/test_resolvability.py
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "resolvability.py"
sys.path.insert(0, str(ROOT / "tools"))
import resolvability as rv  # noqa: E402


def write(tmp, rows, name="ions.json"):
    f = tmp / name
    f.write_text("\n".join(json.dumps(r) for r in rows))
    return f


def run(path, *extra):
    p = subprocess.run([sys.executable, str(TOOL), str(path)] + list(extra),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def two_ions(gap_drift, gap_rt=10.0):
    return [{"run": "r1", "uid": "A|2", "mz": 600.0, "z": 2, "drift": 70.0, "rt": 20.0},
            {"run": "r1", "uid": "B|2", "mz": 700.0, "z": 2,
             "drift": 70.0 + gap_drift, "rt": 20.0 + gap_rt}]


def the_same_gap_flips_with_the_measurement_scale(tmp):
    """Identical coordinates, two instruments. This is the whole idea."""
    f = write(tmp, two_ions(7.4, 0.0))
    bad = []
    rc_fine, out_fine = run(f, "--sigma-drift", "0.87", "--sigma-rt", "0.15")
    if rc_fine != 0 or "uniquely resolved" not in out_fine or "(100.0%)" not in out_fine:
        bad.append("a 7.4-unit gap at sigma 0.87 (8.5 sigma) was not called resolved")
    rc_coarse, out_coarse = run(f, "--sigma-drift", "6.0", "--sigma-rt", "0.15")
    if rc_coarse != 1:
        bad.append("the same gap at sigma 6.0 (1.2 sigma) was not called unresolved")
    if "sharing with 1+ others" not in out_coarse or "(100.0%)" not in out_coarse:
        bad.append("coarse instrument did not report both ions as colliding")
    return bad


def separation_on_either_axis_resolves(tmp):
    """Both axes must be within k sigma to collide; one is enough to separate."""
    bad = []
    # Same drift, far apart in rt -> resolved.
    f = write(tmp, two_ions(0.0, 5.0), "rt_only.json")
    rc, out = run(f, "--sigma-drift", "0.87", "--sigma-rt", "0.15")
    if rc != 0:
        bad.append("ions separated in rt alone were called unresolved")
    # Same rt, far apart in drift -> resolved.
    f = write(tmp, two_ions(7.4, 0.0), "drift_only.json")
    rc, out = run(f, "--sigma-drift", "0.87", "--sigma-rt", "0.15")
    if rc != 0:
        bad.append("ions separated in drift alone were called unresolved")
    # Close on both -> collide.
    f = write(tmp, two_ions(0.5, 0.05), "both.json")
    rc, out = run(f, "--sigma-drift", "0.87", "--sigma-rt", "0.15")
    if rc != 1:
        bad.append("ions close on both axes were not called unresolved")
    return bad


def sigma_comes_from_within_uid_not_between(tmp):
    """The scale is one precursor's own scatter, not the spread of different ones.

    Planted: each precursor repeats with SD about 0.5, while different
    precursors are spread over 60 units. Using the between-ion spread would give
    a sigma two orders too large and call everything unresolved.
    """
    rows = []
    for i in range(6):
        base = 40.0 + 12.0 * i
        for j, jitter in enumerate((-0.5, 0.0, 0.5)):
            rows.append({"run": "r%d" % j, "uid": "P%d|2" % i, "mz": 500.0 + 50 * i,
                         "z": 2, "drift": base + jitter,
                         "rt": 10.0 + 3.0 * i + 0.1 * jitter})
    est, used = rv.estimate_sigma(rows, "drift")
    bad = []
    if used != 6:
        bad.append("used %d precursors for the estimate, expected 6" % used)
    if not 0.3 < est < 0.8:
        bad.append("sigma_drift estimated at %.3f, planted about 0.5" % est)
    rc, out = run(write(tmp, rows, "rep.json"), "--estimate-sigma")
    if "uniquely resolved" not in out or "(100.0%)" not in out:
        bad.append("with a within-uid sigma, the 6 precursors should all resolve")
    if "upper bound" not in out:
        bad.append("did not warn that a cross-run estimate overstates the peak width")
    return bad


def charge_and_mz_do_not_separate_fragments(tmp):
    """Fragments carry parent rt and drift, not parent m/z or charge.

    Two precursors far apart in m/z and different in charge, but sitting on the
    same rt and drift, must still be reported as colliding.
    """
    rows = [{"run": "r1", "uid": "A|2", "mz": 400.0, "z": 2, "drift": 70.0, "rt": 20.0},
            {"run": "r1", "uid": "B|5", "mz": 1200.0, "z": 5, "drift": 70.1, "rt": 20.02}]
    rc, out = run(write(tmp, rows, "mz.json"), "--sigma-drift", "0.87", "--sigma-rt", "0.15")
    bad = []
    if rc != 1:
        bad.append("800 m/z and 3 charges apart was treated as separating fragments")
    if "sharing with 1+ others" not in out or "(100.0%)" not in out:
        bad.append("did not report both as colliding")
    return bad


def a_reproducible_order_is_not_reported_as_a_defect(tmp):
    """--rank-rt plus a stable order gives sigma 0, which is the BEST case.

    The axis has no scatter to measure because the order never changed. Saying
    "a run-level constant or a coarsely rounded value does this" would blame the
    data for being good.
    """
    rows = []
    for i in range(6):
        for j, jit in enumerate((-0.4, 0.0, 0.4)):
            rows.append({"run": "r%d" % j, "uid": "P%d|2" % i, "mz": 500.0 + 60 * i,
                         "z": 2, "drift": 40.0 + 12.0 * i + jit,
                         "rt": 10.0 + 3.0 * i + 0.05 * jit})
    rc, out = run(write(tmp, rows, "stable.json"), "--estimate-sigma", "--rank-rt")
    bad = []
    if "best case" not in out:
        bad.append("blamed the data for a perfectly reproducible elution order")
    if "coarsely rounded" in out:
        bad.append("gave the wrong diagnosis for a rank axis with no scatter")
    if "choose the tolerance" not in out:
        bad.append("did not say the tolerance must be chosen rather than measured")
    # And with a chosen tolerance it must then work.
    rc2, out2 = run(write(tmp, rows, "stable2.json"), "--estimate-sigma",
                    "--rank-rt", "--sigma-rt", "1.0")
    if "uniquely resolved" not in out2:
        bad.append("did not run once a rank tolerance was supplied")
    return bad


def a_zero_sigma_is_refused_not_divided_by(tmp):
    """An axis that did not vary gives SD 0. That is no scale, not a perfect one.

    Dividing by it would make every pair infinitely separated and report a
    flawless run. Real data does this whenever an axis is recorded per run or
    rounded flat.
    """
    rows = []
    for i in range(5):
        for j, jit in enumerate((-0.5, 0.0, 0.5)):
            rows.append({"run": "r%d" % j, "uid": "P%d|2" % i, "mz": 500.0 + 50 * i,
                         "z": 2, "drift": 40.0 + 12.0 * i + jit,
                         "rt": 10.0 + 3.0 * i})          # identical across runs
    rc, out = run(write(tmp, rows, "zero.json"), "--estimate-sigma")
    bad = []
    if rc != 1:
        bad.append("exit %d: ran with a zero measurement scale" % rc)
    if "ZERO" not in out:
        bad.append("did not say the estimate came out at zero")
    if "did not" not in out or "vary" not in out:
        bad.append("did not explain that the axis never varied")
    if "Traceback" in out:
        bad.append("crashed instead of refusing: %r" % out[-200:])
    return bad


def a_missing_scale_is_refused(tmp):
    """No sigma, no answer. An absolute gap is not a result."""
    rc, out = run(write(tmp, two_ions(7.4), "nos.json"))
    bad = []
    if rc != 1:
        bad.append("ran without a measurement scale")
    if "Need a measurement scale" not in out:
        bad.append("did not say why it refused")
    rc, out = run(write(tmp, two_ions(7.4), "nos2.json"), "--estimate-sigma")
    if "could not estimate" not in out:
        bad.append("--estimate-sigma on single-observation data did not say it failed")
    return bad


def crowding_is_reported_not_hidden(tmp):
    """Ten precursors stacked in one cell: the count is the deliverable."""
    rows = [{"run": "r1", "uid": "X%d|2" % i, "mz": 600.0 + i, "z": 2,
             "drift": 70.0 + 0.05 * i, "rt": 20.0 + 0.005 * i} for i in range(10)]
    rc, out = run(write(tmp, rows, "crowd.json"), "--sigma-drift", "0.87", "--sigma-rt", "0.15")
    bad = []
    if rc != 1:
        bad.append("a fully stacked cell was not flagged")
    if "max 9" not in out:
        bad.append("did not report 9 co-residents for each of 10 stacked ions")
    if "fragment-physics prior" not in out:
        bad.append("did not say what has to resolve the collisions")
    return bad


def the_verdict_does_not_depend_on_line_order(tmp):
    """One coordinate per precursor, as a median, not whichever row came first.

    Keeping the first row seen compared coordinates from different runs, so
    shuffling the input could flip a pair between resolved and colliding.
    """
    rows = []
    for i, uid in enumerate(("A|2", "B|2")):
        for j, jit in enumerate((-3.0, 0.0, 3.0)):
            rows.append({"run": "r%d" % j, "uid": uid, "mz": 600.0 + 100 * i,
                         "z": 2, "drift": 70.0 + 7.4 * i + jit, "rt": 20.0 + 0.1 * jit})
    import itertools
    verdicts = set()
    for k, perm in enumerate(itertools.permutations(rows)):
        if k > 20:
            break
        rc, _ = run(write(tmp, list(perm), "ord%d.json" % k),
                    "--sigma-drift", "0.87", "--sigma-rt", "0.15")
        verdicts.add(rc)
    return (["verdict changed with input order: exit codes %s" % sorted(verdicts)]
            if len(verdicts) > 1 else [])


CASES = [the_same_gap_flips_with_the_measurement_scale,
         the_verdict_does_not_depend_on_line_order,
         a_reproducible_order_is_not_reported_as_a_defect,
         a_zero_sigma_is_refused_not_divided_by,
         separation_on_either_axis_resolves,
         sigma_comes_from_within_uid_not_between,
         charge_and_mz_do_not_separate_fragments,
         a_missing_scale_is_refused,
         crowding_is_reported_not_hidden]


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
