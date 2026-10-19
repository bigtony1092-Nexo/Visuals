#!/usr/bin/env python3
"""Build an extraction target list from a search, without selecting on brightness.

tools/extract_profiles.py needs coordinates to extract at. Forty-five of them
came from a protein whose positions were already measured; there is no reason
to stop there, because reading a .raw costs the same whether the pass pulls 45
ions or 500, and the tests that follow were underpowered at 45.

The one thing that cannot be done casually is choosing which ions. Ranking by
intensity selects on the quantity under test: a Fano factor against level, or a
suppression measured as a departure from an ion's own profile, both move when
the set is picked for being bright. On planted data that selection alone
contributed +1.7% of a false result.

So the default stratifies by m/z, which is orthogonal to what is being
measured, and reports the intensity span it happened to get. Where the span is
too narrow to carry the question, intensity may be used instead -- but only
from a DIFFERENT run than the one being extracted, so the selection data and
the test data are separate. Stratifying on intensity from the same run is
refused rather than warned about.

    python3 tools/make_targets.py --psm psm.tsv --apex RUN_ApexRT_temp1.csv \
        --out targets.csv --bins 10 --per-bin 50

Identity comes from the search; the coordinate comes from the feature table,
because a search reports a calculated m/z and this needs the observed one.
"""

import argparse
import collections
import csv
import math
import os
import sys

JOIN_PPM = 15.0        # search m/z to feature m/z
JOIN_RT_MIN = 0.3      # search retention to feature apex
MIN_JOIN_SHARE = 0.25  # below this the two files do not describe one run
RT_PAD = 1.0           # widen each target's window either side of its apex
MZ_PPM = 10.0
DRIFT_BINS = 3.0


class Unusable(Exception):
    pass


def read_psms(path, qcut=0.01):
    rows, kept = [], 0
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (rd.fieldnames or [])}
        need = ("peptide", "charge", "retention", "calculated m/z")
        miss = [n for n in need if n not in cols]
        if miss:
            raise Unusable("%s lacks column(s): %s" % (path, ", ".join(miss)))
        qc = cols.get("qvalue") or cols.get("q-value")
        dec = cols.get("is decoy")
        for r in rd:
            if dec and str(r[dec]).strip().lower() in ("true", "1"):
                continue
            if qc:
                try:
                    if float(r[qc]) > qcut:
                        continue
                except (TypeError, ValueError):
                    pass
            try:
                rows.append({"pep": r[cols["peptide"]].strip(),
                             "z": int(float(r[cols["charge"]])),
                             "rt": float(r[cols["retention"]]),
                             "mz": float(r[cols["calculated m/z"]])})
                kept += 1
            except (TypeError, ValueError):
                continue
    if not rows:
        raise Unusable("%s yielded no usable PSMs at q <= %g" % (path, qcut))
    return rows


def read_features(path):
    feats = []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        need = ("Mz", "Drift")
        miss = [n for n in need if n not in (rd.fieldnames or [])]
        if miss:
            raise Unusable("%s lacks column(s): %s" % (path, ", ".join(miss)))
        rt = "ScanTimeMin" if "ScanTimeMin" in rd.fieldnames else None
        it = "Intensity" if "Intensity" in rd.fieldnames else None
        for r in rd:
            try:
                feats.append({"mz": float(r["Mz"]), "drift": float(r["Drift"]),
                              "rt": float(r[rt]) if rt else None,
                              "inten": float(r[it]) if it and r[it] else 0.0})
            except (TypeError, ValueError):
                continue
    if not feats:
        raise Unusable("%s yielded no usable features" % path)
    return feats


def rt_units(psms, feats, out, forced="auto"):
    """Search engines and feature finders disagree about retention units.

    FragPipe's psm.tsv reports seconds; an ApexRT table reports minutes. The
    two differ by a factor of 60, which is far enough from 1 to detect without
    guessing -- and worth detecting, because comparing the two directly joins
    nothing and then reports that the files describe different runs, which is
    both wrong and the most misleading thing it could say.
    """
    if forced == "min":
        return 1.0
    if forced == "sec":
        out("retention  psm.tsv read as SECONDS (told to)")
        return 1.0 / 60.0
    frt = [f["rt"] for f in feats if f["rt"] is not None]
    if not frt:
        return 1.0
    pm, fm = med([p["rt"] for p in psms]), med(frt)
    if fm <= 0 or pm <= 0:
        return 1.0
    ratio = pm / fm
    if 20 <= ratio <= 200:
        out("retention  psm.tsv is in SECONDS (median %.0f against the feature "
            "table's %.2f min); converting" % (pm, fm))
        return 1.0 / 60.0
    if 0.005 <= ratio <= 0.05:
        out("retention  psm.tsv is in HOURS (median %.3f against the feature "
            "table's %.2f min); converting" % (pm, fm))
        return 60.0
    return 1.0


def join(psms, feats, out):
    """Identity from the search, coordinate from the feature table."""
    feats = sorted(feats, key=lambda f: f["mz"])
    mzs = [f["mz"] for f in feats]
    import bisect
    hits, by = 0, {}
    for p in psms:
        tol = p["mz"] * JOIN_PPM / 1e6
        lo = bisect.bisect_left(mzs, p["mz"] - tol)
        hi = bisect.bisect_right(mzs, p["mz"] + tol)
        best = None
        for f in feats[lo:hi]:
            if f["rt"] is not None and abs(f["rt"] - p["rt"]) > JOIN_RT_MIN:
                continue
            d = abs(f["mz"] - p["mz"])
            if best is None or d < best[0]:
                best = (d, f)
        if best is None:
            continue
        hits += 1
        p["_joined"] = True
        f = best[1]
        key = (p["pep"], p["z"])
        e = by.setdefault(key, {"pep": p["pep"], "z": p["z"], "mz": [], "rt": [],
                                "drift": [], "inten": []})
        e["mz"].append(f["mz"]); e["drift"].append(f["drift"])
        e["rt"].append(f["rt"] if f["rt"] is not None else p["rt"])
        e["inten"].append(f["inten"])
    share = float(hits) / len(psms)
    out("joined %d of %d PSM(s) to a feature (%.0f%%) -> %d distinct ion(s)"
        % (hits, len(psms), 100 * share, len(by)))
    if share < MIN_JOIN_SHARE:
        # Reuse the sorted index the join itself uses. Scanning every feature
        # per PSM was quadratic and looked like a hang on real inputs -- on the
        # very path that exists to explain why nothing joined. Count only PSMs
        # that did NOT join, or the number overstates what retention rejected.
        near = 0
        for p_ in psms:
            if p_.get("_joined"):
                continue
            tol2 = p_["mz"] * JOIN_PPM / 1e6
            a = bisect.bisect_left(mzs, p_["mz"] - tol2)
            b = bisect.bisect_right(mzs, p_["mz"] + tol2)
            if b > a:
                near += 1
        extra = ""
        if near:
            extra = (". %d unjoined PSM(s) DID match on m/z and were rejected "
                     "on retention alone, so check the units: --psm-rt-units "
                     "sec or min" % near)
        raise Unusable(
            "only %.0f%% of PSMs found a feature within %g ppm and %g min. The "
            "two files may not describe the same run, and a target list built "
            "from them would point the extraction at coordinates nothing was "
            "measured at%s" % (100 * share, JOIN_PPM, JOIN_RT_MIN, extra))
    return list(by.values())


def med(v):
    s = sorted(v)
    return s[len(s) // 2]


def stratify(ions, key, bins, per_bin, out):
    """Take per_bin ions from each of `bins` equal-count strata of `key`."""
    ions = sorted(ions, key=lambda i: i[key])
    n = len(ions)
    picked, size = [], max(1, n // bins)
    for b in range(bins):
        chunk = ions[b * size: (b + 1) * size if b < bins - 1 else n]
        if not chunk:
            continue
        step = max(1, len(chunk) // per_bin)
        picked.extend(chunk[::step][:per_bin])
    out("stratified by %s into %d bin(s), up to %d each -> %d target(s)"
        % (key, bins, per_bin, len(picked)))
    return picked


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--psm", required=True, help="search PSM table (psm.tsv)")
    p.add_argument("--apex", required=True, help="feature table for the SAME run")
    p.add_argument("--out", required=True, help="target list CSV to write")
    p.add_argument("--bins", type=int, default=10)
    p.add_argument("--per-bin", type=int, default=50)
    p.add_argument("--q-cut", type=float, default=0.01)
    p.add_argument("--psm-rt-units", choices=("auto", "min", "sec"), default="auto",
                   help="retention units in the PSM table (default: detect)")
    p.add_argument("--stratify", choices=("mz", "intensity"), default="mz",
                   help="what to spread the set across. mz is orthogonal to "
                        "everything being measured; intensity needs --from-run")
    p.add_argument("--from-run",
                   help="with --stratify intensity: the feature table of a "
                        "DIFFERENT run, so selection and test are separate data")
    p.add_argument("--extract-run",
                   help="the run that will be extracted, to check --from-run "
                        "is not the same file")
    args = p.parse_args(argv)

    lines = []
    out = lambda s="": lines.append(s)
    try:
        if args.stratify == "intensity":
            if not args.from_run:
                raise Unusable(
                    "--stratify intensity needs --from-run. Ranking ions by "
                    "their brightness in the run being extracted selects on the "
                    "quantity under test, which is the failure this list exists "
                    "to avoid")
            same = args.extract_run and \
                os.path.abspath(args.from_run) == os.path.abspath(args.extract_run)
            if same or os.path.abspath(args.from_run) == os.path.abspath(args.apex):
                raise Unusable(
                    "--from-run is the same file the targets are being built "
                    "from, so selection and test would share their data. Point "
                    "it at a different run")

        psms = read_psms(args.psm, args.q_cut)
        out("PSMs      %d at q <= %g from %s" % (len(psms), args.q_cut, args.psm))
        feats = read_features(args.apex)
        out("features  %d from %s" % (len(feats), args.apex))
        scale = rt_units(psms, feats, out, args.psm_rt_units)
        if scale != 1.0:
            for p_ in psms:
                p_["rt"] *= scale
        ions = join(psms, feats, out)

        for i in ions:
            i["_mz"] = med(i["mz"]); i["_rt"] = med(i["rt"])
            i["_drift"] = med(i["drift"]); i["_inten"] = med(i["inten"])
        key = "_mz"
        if args.stratify == "intensity":
            other = read_features(args.from_run)
            omz = sorted(other, key=lambda f: f["mz"])
            import bisect
            mzs = [f["mz"] for f in omz]
            for i in ions:
                tol = i["_mz"] * JOIN_PPM / 1e6
                lo = bisect.bisect_left(mzs, i["_mz"] - tol)
                hi = bisect.bisect_right(mzs, i["_mz"] + tol)
                i["_sel"] = max([f["inten"] for f in omz[lo:hi]] or [0.0])
            before = len(ions)
            ions = [i for i in ions if i["_sel"] > 0]
            out("matched %d of %d ion(s) in the selection run" % (len(ions), before))
            key = "_sel"

        picked = stratify(ions, key, args.bins, args.per_bin, out)
        if len(picked) < 20:
            raise Unusable("only %d target(s) survived; too few to measure "
                           "anything against" % len(picked))

        it = sorted(i["_inten"] for i in picked if i["_inten"] > 0)
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["label", "peptide", "charge", "mz", "mz_ppm", "drift",
                        "drift_bins", "rt_lo", "rt_hi", "seen_in_runs",
                        "median_intensity"])
            for i in sorted(picked, key=lambda a: a["_mz"]):
                w.writerow(["%s_z%d" % (i["pep"], i["z"]), i["pep"], i["z"],
                            "%.4f" % i["_mz"], "%.1f" % MZ_PPM,
                            "%.2f" % i["_drift"], "%.1f" % DRIFT_BINS,
                            "%.3f" % (min(i["rt"]) - RT_PAD),
                            "%.3f" % (max(i["rt"]) + RT_PAD),
                            len(i["mz"]), "%.0f" % i["_inten"]])
        out("wrote %s" % args.out)
        out("")
        if it:
            span = it[-1] / it[0] if it[0] > 0 else float("inf")
            out("intensity span of the set: %.0f to %.0f (%.0fx)"
                % (it[0], it[-1], span))
            if span < 50:
                out("")
                out("That span is narrow. The saturation question needs ions")
                out("both well below and into the non-linear range; a set that")
                out("spans less than about 50x may hold no saturated ion at all")
                out("and will answer nothing about a ceiling.")
        out("m/z span: %.1f to %.1f"
            % (min(i["_mz"] for i in picked), max(i["_mz"] for i in picked)))
    except Unusable as e:
        print("\n".join(lines))
        print("")
        print("CANNOT BUILD")
        print("  %s." % str(e).rstrip("."))
        return 2
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
