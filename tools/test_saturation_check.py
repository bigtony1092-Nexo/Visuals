#!/usr/bin/env python3
"""Can it tell a scan-level effect from a channel-level one, and refuse both?

The tool exists to separate droplet competition from detector dead time. The
only way to know it does that is to plant each one and check it says so --
and to plant the shapes that would let it answer by construction and check it
refuses instead.
"""

import os
import random
import subprocess
import sys
import tempfile

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "saturation_check.py")
HDR = "Function,ScanIndex,ScanTimeMin,Mz,MzFwhm,Drift,DriftFwhm,QuadMass,Counts,Intensity,IsSaturated,PeakGroupIndex,LinkCount"

SAT_THRESHOLD = 6.0e6


def write(path, rows):
    """rows: (scan, mz, counts, intensity, sat, group)"""
    with open(path, "w") as fh:
        fh.write(HDR + "\n")
        for scan, mz, counts, inten, sat, grp in rows:
            fh.write("1,%d,%.3f,%.6f,0.0329,45.0,3.19,0.0,%d,%.1f,%d,%s,3\n"
                     % (scan, scan * 0.01, mz, counts, inten, 1 if sat else 0, grp))


def peak(apex_inten, nscans=7, centre=None, width=2.0):
    """A chromatographic profile: rises and falls, so an ion has a shape."""
    mid = (nscans - 1) / 2.0 if centre is None else centre
    return [apex_inten * math.exp(-((i - mid) ** 2) / (2.0 * width ** 2))
            for i in range(nscans)]


import math


def cohort(path, mechanism, seed=11, nblocks=60, per_block=8, nscans=13,
           suppress=0.70):
    """Plant a cohort where suppression is scan-level or channel-level.

    Both cohorts hold identical ions in identical scans with identical
    intensities before suppression. The ONLY difference is who gets depressed
    when a bright ion is present in the scan.

    Intensities are log-uniform and continuous, because a bimodal fixture
    would leave nothing near the flag threshold and so no honest overlap
    between flagged and unflagged -- which is itself a thing the tool refuses.
    """
    rng = random.Random(seed)
    rows = []
    gid = 0
    for block in range(nblocks):
        start_scan = 100 + block * (nscans + 3)
        hot = (block % 2 == 0)
        apexes = []
        for j in range(per_block):
            if hot and j < 2:
                apexes.append(math.exp(rng.uniform(math.log(8e6), math.log(4e7))))
            else:
                apexes.append(math.exp(rng.uniform(math.log(2e5), math.log(5e6))))
        # staggered apexes: ions do not all co-elute with the bright one
        centres = [(nscans - 1) / 2.0 + rng.uniform(-2.0, 2.0)
                   for _ in apexes]
        centres[0] = (nscans - 1) / 2.0
        if len(centres) > 1:
            centres[1] = (nscans - 1) / 2.0 + rng.uniform(-1.0, 1.0)
        profiles = [peak(a, nscans, c) for a, c in zip(apexes, centres)]

        # Which scans does the block's bright ion push over the limit? Those
        # are the scans where droplet charge is spent.
        spent = set()
        if hot:
            for i, v in enumerate(profiles[0]):
                if v > SAT_THRESHOLD:
                    spent.add(start_scan + i)

        for j, prof in enumerate(profiles):
            for i, base in enumerate(prof):
                scan = start_scan + i
                counts = max(5.0, base / 2.0e4)
                inten = base * (1.0 + rng.gauss(0, 1.0 / math.sqrt(counts)))
                if inten <= 0:
                    continue
                if scan in spent:
                    if mechanism == "source":
                        inten *= suppress       # everything in the scan
                    elif mechanism == "detector":
                        if base > SAT_THRESHOLD:
                            inten *= suppress   # only the bright channel
                # A measured flag, not a threshold: the margin varies, so two
                # rows at one intensity can disagree.
                sat = inten > SAT_THRESHOLD * rng.uniform(0.6, 1.6)
                rows.append((scan, 500.0 + gid + j, counts, inten, sat, gid + j))
        gid += per_block
    write(path, rows)
    return rows


def mean_factors(out):
    """The (factor, p) pair for each test, read off the MEAN line only.

    A looser regex also matches the median line that follows it, which then
    zips four factors against two p-values and reports a number belonging to
    neither test.
    """
    import re
    pairs = []
    lines = out.split("\n")
    for i, l in enumerate(lines):
        m = re.search(r"mean marked-minus-unmarked departure [-+][0-9.]+ "
                      r"-> factor ([0-9.]+)", l)
        if not m:
            continue
        fac = float(m.group(1))
        pv = None
        for j in range(i + 1, min(i + 6, len(lines))):
            q = re.search(r"t = [-+][0-9.]+ on \d+ df,  p = ([0-9.e+-]+)", lines[j])
            if q:
                pv = float(q.group(1))
                break
        if pv is not None:
            pairs.append((fac, pv))
    return pairs


def verdict_of(out):
    """The verdict keyword, matched precisely. "SOURCE-LEVEL" is a substring
    of "NOT SOURCE-LEVEL", so a naive `in` test reads one as the other."""
    if "VERDICT" not in out:
        return "none"
    tail = out[out.index("VERDICT"):]
    for key in ("NOT SOURCE-LEVEL", "SOURCE-LEVEL", "DETECTOR-LEVEL",
                "BOTH.", "UNDECIDED", "ABOVE their own"):
        if key in tail:
            return key
    return "unrecognised: " + " ".join(tail.split()[:12])


def run(*args):
    p = subprocess.run([sys.executable, "-I", TOOL] + list(args),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


# ---- the two mechanisms -----------------------------------------------------

def source_level_is_identified(tmp):
    f = os.path.join(tmp, "source.csv")
    cohort(f, "source", seed=11)
    rc, out = run(f)
    bad = []
    v = verdict_of(out)
    if v != "SOURCE-LEVEL":
        bad.append("planted scan-level suppression read as %r" % v)
    if rc != 1:
        bad.append("exit %d on a planted source-level cohort, expected 1 "
                   "(a finding)" % rc)
    return bad


def detector_level_is_identified(tmp):
    f = os.path.join(tmp, "detector.csv")
    cohort(f, "detector", seed=11)
    rc, out = run(f)
    bad = []
    v = verdict_of(out)
    if v != "DETECTOR-LEVEL":
        bad.append("planted dead-time suppression read as %r -- calling it "
                   "source-level would license a per-scan correction on ions "
                   "that were never suppressed" % v)
    if rc != 1:
        bad.append("exit %d on a planted dead-time cohort, expected 1 "
                   "(undercounted bright ions are a finding)" % rc)
    return bad


def a_clean_cohort_passes(tmp):
    """The tool has to be able to say nothing is wrong, or it measures nothing."""
    f = os.path.join(tmp, "clean.csv")
    cohort(f, "none", seed=23)
    rc, out = run(f)
    bad = []
    v = verdict_of(out)
    if v != "NOT SOURCE-LEVEL":
        bad.append("a cohort with no suppression planted read as %r" % v)
    if rc != 0:
        bad.append("exit %d on a clean cohort, expected 0" % rc)
    return bad


def the_null_stays_at_one(tmp):
    """The guard on the whole result.

    Four separate biases in this measurement all push the same way, toward
    the ion appearing ABOVE its own profile, which reads as interference in
    data containing none. An earlier version of this tool measured a null of
    1.089 and called it significant at p = 0.002. If the null drifts off 1.0
    again, every verdict above becomes a statement about the estimator.
    """
    bad = []
    for seed in (23, 31, 47, 59):
        f = os.path.join(tmp, "null%d.csv" % seed)
        cohort(f, "none", seed=seed)
        rc, out = run(f)
        got = mean_factors(out)
        if not got:
            bad.append("seed %d: no factor reported" % seed)
            continue
        for fac, pv in got:
            # The guard that matters is that nothing is DECLARED. A noisy
            # point estimate with a wide interval is not a false positive;
            # calling one significant is.
            if pv < 0.025:
                bad.append("seed %d: factor %.4f called significant at "
                           "p = %.3g in data with nothing suppressed"
                           % (seed, fac, pv))
            if abs(math.log(fac)) > 0.20:
                bad.append("seed %d: null factor %.4f is far enough off 1.0 "
                           "that the estimator, not the data, is moving"
                           % (seed, fac))
        if rc != 0:
            bad.append("seed %d: exit %d on a clean cohort" % (seed, rc))
    return bad


def a_known_factor_is_recovered(tmp):
    """Unbiased is not enough; it has to get the SIZE roughly right."""
    bad = []
    for planted in (0.70, 0.85):
        f = os.path.join(tmp, "amp%d.csv" % int(planted * 100))
        cohort(f, "source", seed=11, suppress=planted)
        rc, out = run(f)
        got = mean_factors(out)
        if not got:
            bad.append("planted %.2f: test B reported no factor" % planted)
            continue
        fac = got[0][0]
        # Mismarking attenuates toward 1, so the recovered factor must sit
        # between the planted value and 1 -- never below it, which would mean
        # the estimator is inventing magnitude.
        if not (planted - 0.05 <= fac <= 1.0):
            bad.append("planted %.2f recovered as %.3f, outside [%.2f, 1.0]"
                       % (planted, fac, planted - 0.05))
    return bad


def the_two_are_distinguished(tmp):
    """The whole point. Same ions, same scans, different mechanism."""
    a, b = os.path.join(tmp, "a.csv"), os.path.join(tmp, "b.csv")
    cohort(a, "source", seed=11)
    cohort(b, "detector", seed=11)
    _, oa = run(a)
    _, ob = run(b)
    bad = []
    va, vb = verdict_of(oa), verdict_of(ob)
    if va == vb:
        bad.append("identical verdict %r for both mechanisms: a tool that "
                   "cannot separate them adds nothing" % va)
    return bad


# ---- the refusals, each in isolation ---------------------------------------

def refuses_constant_counts(tmp):
    f = os.path.join(tmp, "flatcounts.csv")
    rows = cohort(os.path.join(tmp, "_src.csv"), "source", seed=3)
    write(f, [(s, m, 119, i, t, g) for (s, m, c, i, t, g) in rows])
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on constant Counts, expected 2" % rc)
    if "not a count" not in out:
        bad.append("did not say why a constant Counts column is unusable")
    return bad


def refuses_flat_profiles(tmp):
    """No peak shape: every scan of an ion holds one intensity."""
    f = os.path.join(tmp, "flatprof.csv")
    rows = cohort(os.path.join(tmp, "_src2.csv"), "source", seed=5)
    first = {}
    for (s, m, c, i, t, g) in rows:
        first.setdefault(g, i)
    write(f, [(s, m, c, first[g], t, g) for (s, m, c, i, t, g) in rows])
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on flat elution profiles, expected 2" % rc)
    if "no shape" not in out and "ONE intensity" not in out:
        bad.append("did not name the missing peak shape")
    return bad


def refuses_threshold_flag(tmp):
    """IsSaturated reproduced exactly by an intensity cut."""
    f = os.path.join(tmp, "thresh.csv")
    rows = cohort(os.path.join(tmp, "_src3.csv"), "source", seed=7)
    write(f, [(s, m, c, i, i > SAT_THRESHOLD, g) for (s, m, c, i, t, g) in rows])
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on a flag that is a pure intensity cut, expected 2" % rc)
    if "threshold on intensity" not in out and "no" not in out.lower():
        bad.append("did not name the circularity")
    return bad


def refuses_no_flags(tmp):
    f = os.path.join(tmp, "noflag.csv")
    rows = cohort(os.path.join(tmp, "_src4.csv"), "source", seed=9)
    write(f, [(s, m, c, i, False, g) for (s, m, c, i, t, g) in rows])
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d with no flags set, expected 2" % rc)
    if "not a pass" not in out:
        bad.append("an unpopulated flag column must not read as clean")
    return bad


def a_zero_intensity_does_not_abort_the_run(tmp):
    """math.log(0) raised ValueError and took the whole tool down.

    Worse than a crash: the traceback exits 1, which is the same code the tool
    uses for "ran and found something", so extract_profiles' self-check printed
    "reads it" and aria_check recorded a finding.
    """
    f = os.path.join(tmp, "zero.csv")
    rows = cohort(os.path.join(tmp, "_z.csv"), "source", seed=19)
    rows = [(s_, m, c, i, t, g) for (s_, m, c, i, t, g) in rows]
    rows[len(rows) // 2] = (rows[len(rows) // 2][0], rows[len(rows) // 2][1],
                            rows[len(rows) // 2][2], 0.0,
                            rows[len(rows) // 2][4], rows[len(rows) // 2][5])
    write(f, rows)
    rc, out = run(f)
    bad = []
    if "Traceback" in out or "math domain error" in out:
        bad.append("one zero-intensity row aborted the run")
    if rc not in (0, 1, 2):
        bad.append("exit %d" % rc)
    if rc != 2 and "zero or less" not in out:
        bad.append("the ion carrying it was dropped without saying so")
    return bad


def an_export_with_the_pax8_defects_is_refused(tmp):
    """The shipped PAX8 export, reproduced from what was measured in it.

    This used to glob a path under /tmp. That directory existed in the session
    that wrote the test and nowhere else, so on any other machine the glob
    matched nothing, the case returned no complaints, and it reported ok
    WITHOUT RUNNING THE TOOL -- a test that cannot fail, guarding the rule that
    a tool must be able to fail.

    The defects are planted instead, exactly as they were counted in the real
    file: Counts the single value 119 while Intensity spans 733x, one intensity
    repeated across every scan of a peak group, and IsSaturated separating
    intensity with no overlap.
    """
    f = os.path.join(tmp, "pax8_shape.csv")
    rng = random.Random(31)
    rows = []
    for g in range(96):
        inten = math.exp(rng.uniform(math.log(9.0e4), math.log(6.6e7)))
        sat = inten > 6.48e6                     # a pure threshold, as measured
        for k in range(3):                       # three scans, one intensity
            rows.append((110 + g * 10 + k, 400.0 + g, 119, inten, sat, g))
    write(f, rows)
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on an export with the PAX8 defects, expected 2" % rc)
    if "CANNOT RUN" not in out:
        bad.append("did not refuse it")
    if "not a count" not in out:
        bad.append("did not name the constant Counts column")
    return bad


CASES = [source_level_is_identified,
         detector_level_is_identified,
         a_clean_cohort_passes,
         the_null_stays_at_one,
         a_known_factor_is_recovered,
         the_two_are_distinguished,
         refuses_constant_counts,
         refuses_flat_profiles,
         refuses_threshold_flag,
         refuses_no_flags,
         a_zero_intensity_does_not_abort_the_run,
         an_export_with_the_pax8_defects_is_refused]


def main():
    ok = 0
    with tempfile.TemporaryDirectory() as tmp:
        for c in CASES:
            try:
                bad = c(tmp)
            except Exception as e:
                bad = ["raised %s: %s" % (type(e).__name__, e)]
            if bad:
                print("FAIL %s" % c.__name__)
                for b in bad:
                    print("     %s" % b)
            else:
                ok += 1
                print("ok   %s" % c.__name__)
    print("\n%d/%d passed" % (ok, len(CASES)))
    return 0 if ok == len(CASES) else 1


if __name__ == "__main__":
    sys.exit(main())
