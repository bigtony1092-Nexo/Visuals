#!/usr/bin/env python3
"""Can it tell counting noise from multiplicative noise, and refuse neither?

The tool's whole claim is that the form of a tolerance band is a measurable
question rather than a choice. So it has to get both answers right on planted
data, find the crossing when both are present, and refuse input whose
scan-to-scan variation has already been removed -- while still accepting real
instrument data, since a refusal set that rejects everything measures nothing.
"""

import json
import math
import os
import random
import subprocess
import sys
import tempfile

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "counting_noise.py")
HELA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "examples", "hela-scans.example.jsonl")


def run(*args):
    p = subprocess.run([sys.executable, "-I", TOOL] + list(args),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def pois(lam, rng):
    if lam <= 0:
        return 0.0
    if lam < 30:
        L = math.exp(-lam); k = 0; p = 1.0
        while True:
            p *= rng.random()
            if p <= L:
                return float(k)
            k += 1
            if k > 400:
                return float(k)
    return max(0.0, lam + math.sqrt(lam) * rng.gauss(0, 1))


def trend(n, lo, hi):
    """A chromatographic shape, so the estimator has curvature to cope with."""
    out = []
    for i in range(n):
        x = (i / float(n - 1)) * 4 * math.pi
        out.append(lo + (hi - lo) * (0.5 + 0.5 * math.sin(x - math.pi / 2)) ** 2)
    return out


def write(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")


def counting_series(path, n=900, lo=80, hi=4000, seed=3):
    rng = random.Random(seed)
    t = trend(n, lo, hi)
    write(path, [{"run": "p", "scan": i, "sig": pois(v, rng)}
                 for i, v in enumerate(t)])


def multiplicative_series(path, n=900, lo=80, hi=4000, cv=0.25, seed=3):
    rng = random.Random(seed)
    t = trend(n, lo, hi)
    rows = []
    for i, v in enumerate(t):
        rows.append({"run": "m", "scan": i,
                     "sig": max(0.0, v * (1 + rng.gauss(0, cv)))})
    write(path, rows)


def mixed_series(path, n=1400, seed=5):
    """Poisson at the bottom, multiplicative at the top, so a crossing exists."""
    rng = random.Random(seed)
    t = trend(n, 60, 60000)
    rows = []
    for i, v in enumerate(t):
        if v < 600:
            s = pois(v, rng)
        else:
            s = max(0.0, v * (1 + rng.gauss(0, 0.25)))
        rows.append({"run": "x", "scan": i, "sig": s})
    write(path, rows)


def two_channel(path, mode, n=1200, seed=7):
    """A part inside a whole, with a planted relationship between them."""
    rng = random.Random(seed)
    wt = trend(n, 2.0e6, 2.0e7)
    bt = trend(n, 2.0e3, 6.0e4)
    rows = []
    for i in range(n):
        shared = rng.gauss(0, 0.05)        # scan-to-scan delivery
        own = rng.gauss(0, 0.20)
        whole = wt[i] * (1 + shared)
        if mode == "common":
            part = bt[i] * (1 + shared + own)
        elif mode == "competition":
            part = bt[i] * (1 - shared + own)
        else:
            part = bt[i] * (1 + own)
        part = max(1.0, min(part, whole * 0.4))
        rows.append({"run": "c", "scan": i, "bpi": part, "tic": max(whole, part)})
    write(path, rows)


# ---- the two noise forms ---------------------------------------------------

def pure_counting_reads_as_counting(tmp):
    f = os.path.join(tmp, "pois.jsonl")
    counting_series(f)
    rc, out = run(f)
    bad = []
    if "counting describes it at every level" not in out:
        bad.append("Poisson data was not read as counting noise")
    if rc != 0:
        bad.append("exit %d on Poisson data, expected 0" % rc)
    if "1/sqrt(N) is" not in out:
        bad.append("did not say which band form follows")
    return bad


def multiplicative_reads_as_flat(tmp):
    f = os.path.join(tmp, "mult.jsonl")
    multiplicative_series(f)
    rc, out = run(f)
    bad = []
    if "nowhere does counting describe it" not in out:
        bad.append("constant-CV data was not read as multiplicative")
    if "FLAT" not in out:
        bad.append("did not say a flat percentage is the right band")
    if rc != 1:
        bad.append("exit %d on multiplicative data, expected 1" % rc)
    return bad


def the_null_lands_at_one_on_counting(tmp):
    """The guard on everything else.

    The Fano factor is read off second differences, which carry whatever trend
    curvature the data has. If the Poisson null does not come back near 1 on
    data that IS Poisson, the estimator is contributing and no ratio above it
    means anything.
    """
    import re
    f = os.path.join(tmp, "nullcheck.jsonl")
    counting_series(f, seed=11)
    _, out = run(f)
    nulls = [float(m) for m in
             re.findall(r"^\s+\S+\s+sig\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+%",
                        out, re.M)]
    bad = []
    if not nulls:
        bad.append("could not read the null column from the output")
    for v in nulls:
        if not (0.5 <= v <= 2.0):
            bad.append("Poisson null came back %.2f, outside [0.5, 2.0]" % v)
    return bad


def per_ion_profiles_pool_correctly(tmp):
    """One ion's elution is 20-40 scans, not a whole run.

    Differencing has to stay inside one ion -- a step between ions is not
    noise -- while the level bins fill from every ion at once, which is what
    puts a Fano factor on an intensity axis. Planted Poisson must read as
    counting noise at every level.
    """
    f = os.path.join(tmp, "ions.jsonl")
    rng = random.Random(17)
    rows = []
    for ion in range(40):
        apex = math.exp(rng.uniform(math.log(40), math.log(9000)))
        c = 14 + rng.random() * 6
        for k in range(28):
            lam = apex * math.exp(-((k - c) ** 2) / (2 * 4.0 ** 2))
            if lam < 0.5:
                continue
            rows.append({"run": "ion%02d" % ion, "scan": 100 + k,
                         "counts": pois(lam, rng)})
    write(f, rows)
    rc, out = run(f, "--pool")
    bad = []
    if "pooled" not in out:
        bad.append("did not report pooling")
    if "counting describes it at every level" not in out:
        bad.append("planted Poisson per-ion profiles were not read as counting "
                   "noise; got: " + " ".join(
                       l.strip() for l in out.split("\n") if "counting" in l)[:150])
    if rc != 0:
        bad.append("exit %d on Poisson per-ion profiles" % rc)
    return bad


def a_gap_is_not_differenced_across(tmp):
    """The bug this guards.

    An extraction with a count floor drops the scans where an ion fell below
    it, so a profile arrives with holes. Differencing across one spans a jump
    in retention time, where a step is elution and not noise. Left in, planted
    Poisson read as 7x to 260x counting noise.
    """
    rng = random.Random(23)
    def build(path, drop):
        rows = []
        for ion in range(40):
            apex = math.exp(rng.uniform(math.log(200), math.log(9000)))
            for k in range(30):
                lam = apex * math.exp(-((k - 15) ** 2) / (2 * 4.5 ** 2))
                if lam < 0.5:
                    continue
                # drop a block from the middle of the rise, as a floor would
                if drop and 8 <= k <= 11:
                    continue
                rows.append({"run": "ion%02d" % ion, "scan": 100 + k,
                             "counts": pois(lam, rng)})
        write(path, rows)
    whole = os.path.join(tmp, "whole.jsonl")
    holed = os.path.join(tmp, "holed.jsonl")
    build(whole, False); build(holed, True)
    _, a = run(whole, "--pool")
    _, b = run(holed, "--pool")
    bad = []
    for lab, o in (("ungapped", a), ("gapped", b)):
        if "counting describes it at every level" not in o:
            bad.append("%s Poisson profiles were not read as counting noise -- "
                       "a gap is being differenced across" % lab)
    return bad


def the_crossover_is_found(tmp):
    f = os.path.join(tmp, "mixed.jsonl")
    mixed_series(f)
    rc, out = run(f)
    bad = []
    if "counting describes it at and below" not in out and \
       "crossing is not clean" not in out:
        bad.append("a planted crossover was not reported either way")
    if "every level" in out:
        bad.append("called a mixed series pure counting")
    if "nowhere does counting" in out:
        bad.append("called a mixed series pure multiplicative")
    return bad


# ---- the shared term, both directions and neither --------------------------

def common_mode_is_identified(tmp):
    f = os.path.join(tmp, "common.jsonl")
    two_channel(f, "common")
    rc, out = run(f, "--part", "bpi", "--whole", "tic")
    bad = []
    if "COMMON-MODE" not in out:
        bad.append("a planted shared term was not called common-mode")
    if "COMPETITION" in out:
        bad.append("also called it competition")
    return bad


def competition_is_identified(tmp):
    f = os.path.join(tmp, "compete.jsonl")
    two_channel(f, "competition")
    rc, out = run(f, "--part", "bpi", "--whole", "tic")
    bad = []
    if "COMPETITION" not in out:
        bad.append("a planted anti-correlation was not called competition -- "
                   "missing it would hide the one mechanism this is for")
    if "COMMON-MODE" in out:
        bad.append("also called it common-mode")
    return bad


def no_shared_term_when_none_planted(tmp):
    f = os.path.join(tmp, "none.jsonl")
    two_channel(f, "none")
    rc, out = run(f, "--part", "bpi", "--whole", "tic")
    bad = []
    if "NO SHARED TERM FOUND" not in out:
        bad.append("found a shared term in data planted without one")
    return bad


# ---- refusals --------------------------------------------------------------

def refuses_a_smoothed_series(tmp):
    """Scan-to-scan variation already removed: the factor would measure the
    smoothing, not the measurement."""
    f = os.path.join(tmp, "smooth.jsonl")
    rng = random.Random(2)
    t = trend(900, 80, 4000)
    y = [pois(v, rng) for v in t]
    for i in range(1, len(y)):          # forward-fill every other scan
        if i % 2:
            y[i] = y[i - 1]
    write(f, [{"run": "s", "scan": i, "sig": v} for i, v in enumerate(y)])
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on a forward-filled series, expected 2" % rc)
    if "smoothed" not in out and "forward-filled" not in out:
        bad.append("did not name the missing scan-to-scan variation")
    return bad


def refuses_a_constant_channel(tmp):
    f = os.path.join(tmp, "const.jsonl")
    write(f, [{"run": "k", "scan": i, "sig": 119} for i in range(300)])
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on a constant channel, expected 2" % rc)
    if "single value" not in out:
        bad.append("did not name the constant channel")
    return bad


def refuses_too_few_scans(tmp):
    f = os.path.join(tmp, "few.jsonl")
    rng = random.Random(1)
    write(f, [{"run": "f", "scan": i, "sig": pois(500, rng)} for i in range(30)])
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on 30 scans, expected 2" % rc)
    if "below the" not in out:
        bad.append("did not say how many scans are needed")
    return bad


def refuses_negative_values(tmp):
    f = os.path.join(tmp, "neg.jsonl")
    rng = random.Random(1)
    rows = [{"run": "n", "scan": i, "sig": pois(500, rng) - 600}
            for i in range(300)]
    write(f, rows)
    rc, out = run(f)
    bad = []
    if rc != 2:
        bad.append("exit %d on negative values, expected 2" % rc)
    if "cannot be counts" not in out:
        bad.append("did not say why negatives have no Poisson expectation")
    return bad


def refuses_part_exceeding_whole(tmp):
    f = os.path.join(tmp, "over.jsonl")
    rng = random.Random(1)
    write(f, [{"run": "o", "scan": i,
               "bpi": pois(5000, rng), "tic": pois(900, rng)}
              for i in range(300)])
    rc, out = run(f, "--part", "bpi", "--whole", "tic")
    bad = []
    if rc != 2:
        bad.append("exit %d when the part exceeds the whole, expected 2" % rc)
    if "not contained in it" not in out:
        bad.append("did not say the containment claim is false")
    return bad


def refuses_part_that_is_most_of_whole(tmp):
    """Remove 80% of the whole and the remainder cannot speak."""
    f = os.path.join(tmp, "most.jsonl")
    rng = random.Random(1)
    rows = []
    for i in range(400):
        w = pois(100000, rng)
        rows.append({"run": "m", "scan": i, "bpi": w * 0.8, "tic": w})
    write(f, rows)
    rc, out = run(f, "--part", "bpi", "--whole", "tic")
    bad = []
    if rc != 2:
        bad.append("exit %d when the part is 80%% of the whole, expected 2" % rc)
    if "too little independent signal" not in out:
        bad.append("did not say why subtraction leaves nothing")
    return bad


def refuses_part_without_whole(tmp):
    f = os.path.join(tmp, "halfpair.jsonl")
    counting_series(f)
    rc, out = run(f, "--part", "sig")
    bad = []
    if rc != 2:
        bad.append("exit %d with --part and no --whole, expected 2" % rc)
    if "go together" not in out:
        bad.append("did not say the two flags are a pair")
    return bad


def the_real_hela_runs_are_usable(tmp):
    """The mirror of every refusal above.

    A refusal set that rejects real instrument data measures nothing. These
    are three Waters acquisitions read from the .raw, and the tool has to get
    all the way to a verdict on them.
    """
    if not os.path.exists(HELA):
        return []
    rc, out = run(HELA, "--part", "bpi", "--whole", "tic")
    bad = []
    if rc == 2:
        bad.append("refused real Waters data: %s"
                   % " ".join(out.split("CANNOT RUN")[-1].split())[:160])
    if "A  IS IT COUNTING NOISE?" not in out:
        bad.append("did not reach the counting test on real data")
    if "COMMON-MODE" not in out and "COMPETITION" not in out and \
       "NO SHARED TERM" not in out:
        bad.append("did not reach a shared-term verdict on real data")
    return bad


CASES = [pure_counting_reads_as_counting,
         multiplicative_reads_as_flat,
         the_null_lands_at_one_on_counting,
         the_crossover_is_found,
         per_ion_profiles_pool_correctly,
         a_gap_is_not_differenced_across,
         common_mode_is_identified,
         competition_is_identified,
         no_shared_term_when_none_planted,
         refuses_a_smoothed_series,
         refuses_a_constant_channel,
         refuses_too_few_scans,
         refuses_negative_values,
         refuses_part_exceeding_whole,
         refuses_part_that_is_most_of_whole,
         refuses_part_without_whole,
         the_real_hela_runs_are_usable]


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
