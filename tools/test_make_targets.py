#!/usr/bin/env python3
"""Does the target list avoid choosing ions by the thing under test?

A list ranked by brightness selects on the quantity the extraction then
measures. The tool's job is to make that hard to do by accident and impossible
to do silently.
"""

import csv
import math
import os
import random
import subprocess
import sys
import tempfile

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "make_targets.py")


def run(*a):
    p = subprocess.run([sys.executable, "-I", TOOL] + list(a),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def make_run(tmp, name, n=400, seed=1, shift_mz=0.0):
    """A matched psm.tsv and feature table describing one run."""
    rng = random.Random(seed)
    psm = os.path.join(tmp, name + "_psm.tsv")
    apex = os.path.join(tmp, name + "_apex.csv")
    ions = []
    for i in range(n):
        mz = rng.uniform(400, 1400)
        ions.append({"pep": "PEP%04d" % i, "z": rng.choice([2, 2, 3]),
                     "mz": mz, "rt": rng.uniform(8, 50),
                     "drift": 40 + (mz - 400) * 0.05 + rng.gauss(0, 1.5),
                     "inten": math.exp(rng.uniform(math.log(5e3), math.log(6e7)))})
    with open(psm, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Spectrum", "Peptide", "Charge", "Retention",
                    "Calculated M/Z", "Qvalue", "Is Decoy"])
        for i in ions:
            for _ in range(rng.randint(1, 3)):
                w.writerow(["s", i["pep"], i["z"], "%.3f" % i["rt"],
                            "%.4f" % i["mz"], "0.001", "false"])
        for k in range(30):                       # decoys that must be dropped
            w.writerow(["s", "DECOY%d" % k, 2, "20.0", "700.0", "0.44", "true"])
    with open(apex, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Function", "ScanIndex", "ScanTimeMin", "Mz", "Drift",
                    "Counts", "Intensity", "IsSaturated", "PeakGroupIndex"])
        for n_, i in enumerate(ions):
            w.writerow([1, 100 + n_, "%.3f" % i["rt"],
                        "%.4f" % (i["mz"] + shift_mz), "%.2f" % i["drift"],
                        100, "%.0f" % i["inten"], 0, n_])
    return psm, apex, ions


def read_targets(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def builds_a_list_stratified_by_mz(tmp):
    psm, apex, _ = make_run(tmp, "a")
    out_csv = os.path.join(tmp, "t_mz.csv")
    rc, out = run("--psm", psm, "--apex", apex, "--out", out_csv,
                  "--bins", "8", "--per-bin", "20")
    bad = []
    if rc != 0:
        bad.append("exit %d: %s" % (rc, out.strip().splitlines()[-1:]))
        return bad
    t = read_targets(out_csv)
    if len(t) < 20:
        bad.append("only %d target(s)" % len(t))
    mz = sorted(float(r["mz"]) for r in t)
    if mz[-1] - mz[0] < 500:
        bad.append("m/z span %.0f is too narrow for a stratified set"
                   % (mz[-1] - mz[0]))
    for col in ("drift", "rt_lo", "rt_hi", "mz_ppm", "drift_bins"):
        if not all(r.get(col) for r in t):
            bad.append("column %s is not populated" % col)
    return bad


def the_set_is_not_ranked_by_brightness(tmp):
    """The point of the tool.

    An m/z-stratified set must not come out as the brightest ions. If its
    median intensity sits near the top of the run's own distribution, the list
    has selected on the quantity the extraction will measure.
    """
    psm, apex, ions = make_run(tmp, "b", seed=4)
    out_csv = os.path.join(tmp, "t_b.csv")
    rc, _ = run("--psm", psm, "--apex", apex, "--out", out_csv,
                "--bins", "8", "--per-bin", "20")
    if rc != 0:
        return ["could not build the list"]
    t = read_targets(out_csv)
    picked = sorted(float(r["median_intensity"]) for r in t)
    allit = sorted(i["inten"] for i in ions)
    pm = picked[len(picked) // 2]
    rank = sum(1 for v in allit if v < pm) / float(len(allit))
    bad = []
    if not (0.3 <= rank <= 0.7):
        bad.append("the set's median intensity sits at the %.0fth percentile "
                   "of the run; an unbiased stratification should land near "
                   "the middle" % (100 * rank))
    return bad


def refuses_intensity_stratification_without_another_run(tmp):
    psm, apex, _ = make_run(tmp, "c")
    rc, out = run("--psm", psm, "--apex", apex, "--out", os.path.join(tmp, "x.csv"),
                  "--stratify", "intensity")
    bad = []
    if rc != 2:
        bad.append("exit %d, expected 2" % rc)
    if "selects on the quantity under test" not in out:
        bad.append("did not say why ranking by brightness is the failure")
    return bad


def refuses_selecting_from_the_run_being_extracted(tmp):
    psm, apex, _ = make_run(tmp, "d")
    rc, out = run("--psm", psm, "--apex", apex, "--out", os.path.join(tmp, "y.csv"),
                  "--stratify", "intensity", "--from-run", apex)
    bad = []
    if rc != 2:
        bad.append("exit %d when the selection run IS the target run, expected 2" % rc)
    if "share their data" not in out:
        bad.append("did not name the circularity")
    return bad


def accepts_intensity_from_a_different_run(tmp):
    psm, apex, _ = make_run(tmp, "e", seed=7)
    _, other, _ = make_run(tmp, "e_other", seed=7)   # same ions, separate file
    out_csv = os.path.join(tmp, "t_e.csv")
    rc, out = run("--psm", psm, "--apex", apex, "--out", out_csv,
                  "--stratify", "intensity", "--from-run", other,
                  "--bins", "6", "--per-bin", "20")
    bad = []
    if rc != 0:
        bad.append("exit %d with a separate selection run: %s"
                   % (rc, out.strip().splitlines()[-2:]))
    elif len(read_targets(out_csv)) < 20:
        bad.append("too few targets")
    return bad


def seconds_and_minutes_are_reconciled(tmp):
    """The bug this guards, found on a real-shaped file pair.

    FragPipe writes retention in seconds; an ApexRT table writes minutes. Left
    alone the join matched nothing and the tool then said the two files
    described different runs -- wrong, and the most misleading thing it could
    have said, since the m/z values agreed to 0 ppm.
    """
    psm, apex, _ = make_run(tmp, "sec", seed=13)
    secs = os.path.join(tmp, "sec_seconds_psm.tsv")
    with open(psm, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    hdr, body = rows[0], rows[1:]
    rt = hdr.index("Retention")
    with open(secs, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        for r in body:
            r = list(r)
            r[rt] = "%.1f" % (float(r[rt]) * 60.0)     # minutes -> seconds
            w.writerow(r)
    out_csv = os.path.join(tmp, "t_sec.csv")
    rc, out = run("--psm", secs, "--apex", apex, "--out", out_csv,
                  "--bins", "6", "--per-bin", "20")
    bad = []
    if "SECONDS" not in out:
        bad.append("did not detect that the PSM table was in seconds")
    if rc != 0:
        bad.append("exit %d on a seconds/minutes pair that describes one run" % rc)
    else:
        n = len(read_targets(out_csv))
        if n < 20:
            bad.append("only %d target(s) after unit conversion" % n)
    # and the override has to work too
    rc2, out2 = run("--psm", secs, "--apex", apex, "--out",
                    os.path.join(tmp, "t_sec2.csv"), "--psm-rt-units", "sec",
                    "--bins", "6", "--per-bin", "20")
    if rc2 != 0:
        bad.append("exit %d with --psm-rt-units sec given explicitly" % rc2)
    return bad


def a_real_unit_mismatch_still_names_units(tmp):
    """When m/z agrees and retention does not, say so rather than blaming the
    files for being different runs."""
    psm, apex, _ = make_run(tmp, "mm", seed=21)
    bad_rt = os.path.join(tmp, "mm_badrt_psm.tsv")
    with open(psm, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    hdr, body = rows[0], rows[1:]
    rt = hdr.index("Retention")
    with open(bad_rt, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        for r in body:
            r = list(r)
            r[rt] = "%.3f" % (float(r[rt]) * 7.0)   # not 60: no clean unit story
            w.writerow(r)
    rc, out = run("--psm", bad_rt, "--apex", apex, "--out",
                  os.path.join(tmp, "t_mm.csv"))
    bad = []
    if rc != 2:
        bad.append("exit %d on a retention mismatch, expected 2" % rc)
    if "check the units" not in out:
        bad.append("did not point at units when m/z matched and retention did not")
    return bad


def refuses_mismatched_files(tmp):
    """A search and a feature table from different runs would point the
    extraction at coordinates nothing was measured at."""
    psm, _, _ = make_run(tmp, "f", seed=2)
    _, apex, _ = make_run(tmp, "g", seed=99, shift_mz=3.0)
    rc, out = run("--psm", psm, "--apex", apex, "--out", os.path.join(tmp, "z.csv"))
    bad = []
    if rc != 2:
        bad.append("exit %d on a search and features from different runs" % rc)
    if "may not describe the same run" not in out:
        bad.append("did not say the two files disagree")
    return bad


def reports_a_narrow_intensity_span(tmp):
    """A set that holds no bright ion answers nothing about a ceiling."""
    rng = random.Random(3)
    psm = os.path.join(tmp, "n_psm.tsv")
    apex = os.path.join(tmp, "n_apex.csv")
    rows = [{"pep": "P%04d" % i, "z": 2, "mz": rng.uniform(400, 1400),
             "rt": rng.uniform(8, 50)} for i in range(300)]
    with open(psm, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Peptide", "Charge", "Retention", "Calculated M/Z", "Qvalue"])
        for r in rows:
            w.writerow([r["pep"], 2, "%.3f" % r["rt"], "%.4f" % r["mz"], "0.001"])
    with open(apex, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Function", "ScanIndex", "ScanTimeMin", "Mz", "Drift",
                    "Counts", "Intensity", "IsSaturated", "PeakGroupIndex"])
        for n_, r in enumerate(rows):
            w.writerow([1, n_, "%.3f" % r["rt"], "%.4f" % r["mz"], 60.0, 100,
                        "%.0f" % rng.uniform(9e5, 1.1e6), 0, n_])
    rc, out = run("--psm", psm, "--apex", apex, "--out", os.path.join(tmp, "n.csv"),
                  "--bins", "6", "--per-bin", "20")
    bad = []
    if "span is narrow" not in out:
        bad.append("did not warn that a narrow intensity span answers nothing "
                   "about a ceiling")
    return bad


def the_output_feeds_the_extractor(tmp):
    """The two tools have to actually fit together."""
    psm, apex, _ = make_run(tmp, "h", seed=11)
    out_csv = os.path.join(tmp, "t_h.csv")
    rc, _ = run("--psm", psm, "--apex", apex, "--out", out_csv,
                "--bins", "6", "--per-bin", "12")
    if rc != 0:
        return ["could not build the list"]
    ex = os.path.join(os.path.dirname(TOOL), "extract_profiles.py")
    p = subprocess.run([sys.executable, "-I", ex, "demo", "--targets", out_csv,
                        "--out", os.path.join(tmp, "prof"), "--name", "T"],
                       capture_output=True, text=True)
    body = p.stdout + p.stderr
    bad = []
    if "targets" not in body:
        bad.append("extract_profiles.py could not read the list: %s"
                   % body.strip().splitlines()[-2:])
    if "wrote" not in body:
        bad.append("no profile was written from the generated list")
    return bad


CASES = [builds_a_list_stratified_by_mz,
         the_set_is_not_ranked_by_brightness,
         refuses_intensity_stratification_without_another_run,
         refuses_selecting_from_the_run_being_extracted,
         accepts_intensity_from_a_different_run,
         seconds_and_minutes_are_reconciled,
         a_real_unit_mismatch_still_names_units,
         refuses_mismatched_files,
         reports_a_narrow_intensity_span,
         the_output_feeds_the_extractor]


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
