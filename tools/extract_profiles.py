#!/usr/bin/env python3
"""Pull one ion's intensity scan by scan, in the shape the two tools already read.

Everything in this repository that would answer whether noise falls at
saturation, or whether suppression is scan-level or channel-level, needs the
same thing: a single ion followed across the scans it elutes over, in raw
counts. docs/per-scan-profiles.md says why neither a total ion current nor a
base peak can stand in for it.

This writes that file. The vendor-specific part -- opening a .raw and handing
back the points in a scan -- is ONE interface, so the reader that already works
on the acquisition machine can be used unchanged rather than rewritten here:

    class Reader:
        def scan_count(self): ...                 # how many scans
        def scan_rt(self, i): ...                 # retention time, minutes
        def scan_points(self, i): ...             # (mz, drift_bin, counts)
                                                  # or (mz, drift_bin, counts, saturated)

Point `--reader` at a module exposing `open_run(path) -> Reader`. The repository
ships `--reader demo`, which synthesises a run so the writer and its checks can
be exercised without an instrument.

    python3 tools/extract_profiles.py RUN.raw --targets examples/pax8-targets.csv \
        --reader my_masslynx_reader --out profiles/

Counts are summed, never scaled: the Poisson expectation is the whole point of
the exercise and a scaled intensity has none. Before writing, the output is put
through the refusals of saturation_check.py and counting_noise.py, so a file
that those tools would reject is reported here rather than after it travels.
"""

import argparse
import csv
import importlib
import json
import math
import os
import subprocess
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))

# A scan whose summed counts are below this is background for that ion rather
# than a measurement of it, and carrying it would put the elution's tails into
# every variance.
MIN_COUNTS = 3

# An ion needs this many scans for the leave-one-out reference the tests build,
# and at least counting_noise.MIN_SCANS_POOLED so the file does not carry ions
# its own consumer will skip.
MIN_SCANS = 8


class DemoReader:
    """A synthetic run, so the writer and its checks can be tested with no
    instrument. Counts are Poisson about a Gaussian elution profile, which is
    what the tools expect to find and what the shipped export did not have."""

    def __init__(self, targets, seed=5, scans=900, span=60.0):
        import random
        self.rng = random.Random(seed)
        self.n = scans
        self.span = span
        self.targets = targets

    def scan_count(self):
        return self.n

    def scan_rt(self, i):
        return self.span * i / float(self.n - 1)

    def scan_points(self, i):
        rt = self.scan_rt(i)
        pts = []
        for t in self.targets:
            centre = (t["rt_lo"] + t["rt_hi"]) / 2.0
            width = max(0.05, (t["rt_hi"] - t["rt_lo"]) / 6.0)
            lam = t["_apex"] * math.exp(-((rt - centre) ** 2) / (2 * width ** 2))
            if lam < 0.2:
                continue
            c = self._pois(lam)
            if c > 0:
                pts.append((t["mz"], t["drift"], float(c)))
        return pts

    def _pois(self, lam):
        if lam < 30:
            L = math.exp(-lam)
            k, p = 0, 1.0
            while True:
                p *= self.rng.random()
                if p <= L:
                    return k
                k += 1
                if k > 500:
                    return k
        return max(0, int(lam + math.sqrt(lam) * self.rng.gauss(0, 1)))


def read_targets(path):
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                out.append({
                    "label": r["label"],
                    "mz": float(r["mz"]),
                    "mz_ppm": float(r.get("mz_ppm") or 10.0),
                    "drift": float(r["drift"]),
                    "drift_bins": float(r.get("drift_bins") or 3.0),
                    "rt_lo": float(r["rt_lo"]),
                    "rt_hi": float(r["rt_hi"]),
                    "_apex": float(r.get("median_intensity") or 1e5) / 2.0e3,
                })
            except (KeyError, ValueError) as e:
                raise SystemExit("targets file: %s in row %r" % (e, r))
    if not out:
        raise SystemExit("targets file holds no usable rows")
    return out


def extract(reader, targets, out):
    """Sum raw counts inside each target's window, one row per ion per scan."""
    rows = []
    n = reader.scan_count()
    out("reading %d scan(s)" % n)
    for i in range(n):
        rt = reader.scan_rt(i)
        live = [t for t in targets if t["rt_lo"] <= rt <= t["rt_hi"]]
        if not live:
            continue
        pts = list(reader.scan_points(i))
        if not pts:
            continue
        for t in live:
            tol = t["mz"] * t["mz_ppm"] / 1e6
            c = 0.0
            sat = None
            for pt in pts:
                mz, drift, counts = pt[0], pt[1], pt[2]
                if abs(mz - t["mz"]) <= tol and abs(drift - t["drift"]) <= t["drift_bins"]:
                    c += counts
                    if len(pt) > 3:
                        # The reader knows; carry its answer rather than
                        # inventing one. A flag derived from the counts here
                        # would be an intensity threshold wearing the vendor's
                        # name, which is exactly what saturation_check refuses.
                        sat = bool(sat) or bool(pt[3])
            if c >= MIN_COUNTS:
                rows.append({"label": t["label"], "scan": i, "rt": rt,
                             "mz": t["mz"], "drift": t["drift"], "counts": c,
                             "sat": sat})
    return rows


def write_outputs(rows, run, outdir, out):
    """Two files: the ApexRT shape saturation_check.py reads, and a per-ion
    series counting_noise.py reads with the ion as the series."""
    by = {}
    for r in rows:
        by.setdefault(r["label"], []).append(r)
    keep = {k: v for k, v in by.items() if len(v) >= MIN_SCANS}
    dropped = len(by) - len(keep)
    if dropped:
        out("%d ion(s) had fewer than %d scans and were left out" % (dropped, MIN_SCANS))
    if not keep:
        raise SystemExit("no ion reached %d scans; nothing to write" % MIN_SCANS)

    sat = os.path.join(outdir, run, "Summary")
    os.makedirs(sat, exist_ok=True)
    csv_path = os.path.join(sat, "%s_profiles.csv" % run)
    idx = {k: n + 1 for n, k in enumerate(sorted(keep))}
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Function", "ScanIndex", "ScanTimeMin", "Mz", "MzFwhm", "Drift",
                    "DriftFwhm", "QuadMass", "Counts", "Intensity", "IsSaturated",
                    "PeakGroupIndex", "LinkCount", "Label"])
        for k in sorted(keep):
            for r in sorted(keep[k], key=lambda x: x["scan"]):
                # Intensity is the counts here. Nothing is scaled, because a
                # scaled intensity has no Poisson expectation to be held to.
                # An unpopulated flag must not read as "not saturated": left
                # empty, saturation_check says it cannot run, which is true.
                flag = "" if r["sat"] is None else ("1" if r["sat"] else "0")
                w.writerow([1, r["scan"], "%.4f" % r["rt"], "%.4f" % r["mz"], "",
                            "%.3f" % r["drift"], "", 0.0, "%.0f" % r["counts"],
                            "%.0f" % r["counts"], flag, idx[k], len(keep[k]), k])
    jsonl_path = os.path.join(outdir, "%s_profiles.jsonl" % run)
    with open(jsonl_path, "w") as fh:
        for k in sorted(keep):
            for r in sorted(keep[k], key=lambda x: x["scan"]):
                # The ION is the series here, so the Fano factor comes out per
                # ion rather than per run -- which is the measurement that a
                # total ion current cannot give.
                fh.write(json.dumps({"run": k, "scan": r["scan"],
                                     "counts": r["counts"]}) + "\n")
    out("wrote %s  (%d ion(s), %d row(s))"
        % (csv_path, len(keep), sum(len(v) for v in keep.values())))
    out("wrote %s" % jsonl_path)
    if all(r["sat"] is None for r in rows):
        out("")
        out("NOTE  the reader gave no saturation flag, so that column is empty.")
        out("      counting_noise.py does not need it and will still answer")
        out("      whether the Fano factor falls at the apex of a bright peak.")
        out("      saturation_check.py does need it and will say so.")
    return csv_path, jsonl_path, keep


def self_check(csv_path, jsonl_path, out):
    """Run both consumers' refusals here, so a file they would reject is known
    about before it travels."""
    out("")
    out("SELF-CHECK  the two tools that will read this")
    ok = True
    for name, argv in (("counting_noise.py", [jsonl_path, "--pool"]),
                       ("saturation_check.py", [csv_path])):
        p = subprocess.run([sys.executable, os.path.join(TOOLS, name)] + argv,
                           capture_output=True, text=True)
        body = (p.stdout or "") + (p.stderr or "")
        if p.returncode == 2:
            ok = False
            why = body.split("CANNOT RUN", 1)[-1].strip().splitlines()
            out("  %-22s WOULD REFUSE" % name)
            for line in why[:4]:
                if line.strip():
                    out("      %s" % line.strip())
        else:
            out("  %-22s reads it (exit %d)" % (name, p.returncode))
    out("")
    if ok:
        out("  Both tools accept this file. It is worth sending.")
    else:
        out("  Fix the export before sending it. A refusal is not a result,")
        out("  and neither is a number computed from a file like this.")
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run", help="the acquisition to read (a .raw, or 'demo')")
    p.add_argument("--targets", required=True, help="target list CSV")
    p.add_argument("--reader", default="demo",
                   help="module exposing open_run(path) -> Reader, or 'demo'")
    p.add_argument("--out", default="profiles", help="output directory")
    p.add_argument("--name", help="run name for the output path (default: the file's)")
    args = p.parse_args(argv)

    lines = []
    out = lambda s="": lines.append(s)
    targets = read_targets(args.targets)
    out("targets   %d ion(s) from %s" % (len(targets), args.targets))

    if args.reader == "demo":
        out("reader    demo (synthetic; no instrument involved)")
        reader = DemoReader(targets)
    else:
        try:
            mod = importlib.import_module(args.reader)
        except ImportError as e:
            raise SystemExit("cannot import reader %r: %s" % (args.reader, e))
        if not hasattr(mod, "open_run"):
            raise SystemExit("reader %r has no open_run(path)" % args.reader)
        out("reader    %s" % args.reader)
        reader = mod.open_run(args.run)

    name = args.name or os.path.splitext(os.path.basename(args.run))[0]
    rows = extract(reader, targets, out)
    out("extracted %d row(s) over %d ion(s)"
        % (len(rows), len(set(r["label"] for r in rows))))
    csv_path, jsonl_path, keep = write_outputs(rows, name, args.out, out)
    ok = self_check(csv_path, jsonl_path, out)
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
