#!/usr/bin/env python3
"""Merge single-run exports into one cohort export.

Every cross-sample panel needs this, and so does the atlas_consensus tier: one
run cannot observe agreement across runs. The spine writes one export per run,
so without a merge the cohort view does not exist.

    python3 tools/merge_exports.py export_a/ export_b/ export_c/ -o cohort/

What it will not do, because each would be a quiet lie:

  Assign atlas_consensus. Agreement across samples is a decision with a rule
  behind it -- how many samples, how close -- and this tool has no such rule.
  It reports which ions recur and leaves the tier alone.

  Merge two exports that disagree about one sample_id. The same id with a
  different vendor, acquisition method or drift_kind is two different samples.

  Average anything. Nothing here combines intensities across runs, because the
  runs are not normalised to each other and the merge is not the place to decide
  how they should be.

  Pool proximity controls measured through different windows. The bare_proximity
  false-positive rate is a property of a window; summing probe counts from runs
  that used different tolerances, or compared different axes, would report a rate
  for a window that never existed. Where the windows agree the counts pool, which
  is what independent Bernoulli trials on one window allow; where they differ the
  merged manifest carries no control and says which runs disagreed.

  Declare a table it did not write. tables{} comes from the files produced.

Exit 0 on a clean merge, 1 on a finding, 2 if nothing could be merged.
"""

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ladder_stats import wilson, katz_ratio  # noqa: E402

TABLES = ("ions", "spectra", "traces", "proteins", "fragment_hits")
SAMPLE_KEYS = ("vendor", "acquisition_method", "drift_kind", "has_ims",
               "search_engine", "batch", "sample_class")


# The window a proximity control measured. Two runs whose controls differ on any
# of these measured different things, and their counts do not add.
WINDOW_KEYS = ("precursor_ppm", "rt_tol", "drift_tol", "min_scans")


def _label_findings(have):
    """One finding per text, naming the runs that raised it when not all did."""
    seen = {}
    for name, pc in have:
        for f in (pc.get("findings") or []):
            seen.setdefault(f, []).append(name)
    out = []
    for f, runs in sorted(seen.items()):
        out.append(f if len(runs) == len(have)
                   else "%s [%s]" % (f, ", ".join(sorted(runs))))
    return out


def pool_proximity_control(manifests):
    """Sum the control arms across runs, or refuse and say why.

    Pooling k and n across runs is legitimate only when every run queried the
    SAME window: the arms are counts of independent Bernoulli trials against one
    window definition. Different tolerances, or a run where retention time was
    compared and another where it was not, make the sum a rate for a window
    nobody ran. Returns (pooled_or_None, note).
    """
    have = [(Path(d).name, (m.get("proximity_control") or None)) for d, m in manifests]
    missing = [n for n, pc in have if not pc]
    if missing:
        return None, ("%d of %d run(s) carry no measured bare_proximity"
                      " false-positive rate (%s), so the cohort has none. A"
                      " cohort rate cannot be assembled from runs that did not"
                      " measure one."
                      % (len(missing), len(have), ", ".join(sorted(missing))))
    windows = {}
    for n, pc in have:
        key = (tuple(pc.get(k) for k in WINDOW_KEYS),
               tuple(pc.get("axes_matched") or []))
        windows.setdefault(key, []).append(n)
    if len(windows) > 1:
        detail = "; ".join(
            "%s on %s" % (", ".join("%s=%s" % (k, v) for k, v in zip(WINDOW_KEYS, key[0]))
                          + " over " + "/".join(key[1]), ", ".join(sorted(runs)))
            for key, runs in windows.items())
        return None, ("the runs measured their bare_proximity false-positive rate"
                      " through DIFFERENT windows, so the counts do not add and the"
                      " cohort carries no control: %s" % detail)

    (window_key, axes), _ = next(iter(windows.items())), None
    arms = {}
    for _, pc in have:
        for name, a in (pc.get("arms") or {}).items():
            if not a.get("queried"):
                continue
            acc = arms.setdefault(name, {"arm": name, "queried": 0, "matched": 0,
                                         "distinct_peaks_matched": 0})
            acc["queried"] += a["queried"]
            acc["matched"] += a["matched"]
            acc["distinct_peaks_matched"] += a.get("distinct_peaks_matched") or 0
    if not arms:
        return None, ("no run queried a single probe, so there is nothing to pool."
                      " That is not a cohort rate of zero.")
    for a in arms.values():
        lo, hi = wilson(a["matched"], a["queried"])
        a["rate"] = a["matched"] / a["queried"]
        a["ci"] = [lo, hi]
        a["axes_matched"] = list(axes)
    tgt = arms.get("target")
    for name, a in arms.items():
        a["ratio_vs_target"] = None if (name == "target" or not tgt) else katz_ratio(
            tgt["matched"], tgt["queried"], a["matched"], a["queried"])
    pooled = {"pooled_from": sorted(n for n, _ in have),
              "axes_matched": list(axes), "arms": arms,
              "cohort_false_rate": (arms["entrapment"]["rate"]
                                    if "entrapment" in arms else None),
              # Carried from the runs rather than recomputed: the strata are
              # per-run m/z quartiles and do not align across runs, so a pooled
              # stratum would be a different partition from any ion's own.
              "per_run": {n: {"arms": {k: {"queried": v.get("queried"),
                                           "matched": v.get("matched")}
                                       for k, v in (pc.get("arms") or {}).items()},
                              "calibration": pc.get("calibration")}
                          for n, pc in have},
              # Labelled by run, the same way this tool collapses warnings: a
              # finding that says "on this run" is misleading in a cohort
              # manifest unless it names which run raised it.
              "findings": _label_findings(have)}
    for k, v in zip(WINDOW_KEYS, window_key):
        pooled[k] = v
    return pooled, None


def read_jsonl(path):
    if not path.exists():
        return []
    return [json.loads(l) for l in path.read_text().splitlines() if l.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("exports", nargs="+", help="export directories to merge")
    ap.add_argument("-o", "--out", required=True, help="cohort export directory")
    ap.add_argument("--export-id", help="id for the merged export")
    args = ap.parse_args(argv)

    manifests, findings = [], []
    for d in args.exports:
        p = Path(d)
        mp = p if p.is_file() else p / "manifest.json"
        if not mp.exists():
            findings.append("no manifest at %s" % mp)
            continue
        manifests.append((mp.parent, json.loads(mp.read_text())))
    if not manifests:
        for f in findings:
            print("  %s" % f, file=sys.stderr)
        print("nothing to merge", file=sys.stderr)
        return 2

    samples, rows, uid_from = {}, {k: [] for k in TABLES}, {}
    total_before = 0
    for d, man in manifests:
        for s in man.get("samples", []):
            sid = s["sample_id"]
            if sid in samples:
                clash = [k for k in SAMPLE_KEYS if samples[sid].get(k) != s.get(k)]
                if clash:
                    findings.append(
                        "sample_id %s appears in two exports with different %s."
                        " The same id describing two different acquisitions is not"
                        " one sample, and merging them would put both under one"
                        " label." % (sid, ", ".join(clash)))
                    continue
            samples[sid] = s
        total_before += (man.get("counts") or {}).get("ions_total", 0)
        for table in TABLES:
            for r in read_jsonl(d / ("%s.jsonl" % table)):
                if table == "ions":
                    uid = r["ion_uid"]
                    if uid in uid_from and uid_from[uid] != str(d):
                        # The uid is content-addressed over sample_id and the
                        # source file, so a collision across runs means two
                        # different detections hashed the same.
                        findings.append(
                            "ion_uid %s is in both %s and %s. The address covers"
                            " sample_id and source_file, so this is a hash"
                            " collision, not the same ion seen twice."
                            % (uid, uid_from[uid], d))
                        continue
                    uid_from[uid] = str(d)
                rows[table].append(r)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tables = {}
    for table in TABLES:
        if rows[table]:
            with (out / ("%s.jsonl" % table)).open("w") as fh:
                for r in rows[table]:
                    fh.write(json.dumps(r) + "\n")
        tables[table] = bool(rows[table])
    kinds = sorted({r.get("kind") for r in rows["spectra"]
                    if r.get("kind") in ("ms2_experimental", "ms2_library",
                                         "ms1_isotope")})

    tiers = {}
    for r in rows["ions"]:
        t = r["confidence"]["tier"]
        tiers[t] = tiers.get(t, 0) + 1

    # Which identities recur across samples. Reported, never promoted: the tier
    # above needs a rule for how many samples and how close, and that rule is
    # not here.
    seen = {}
    for r in rows["ions"]:
        ann = r.get("annotation") or {}
        key = (ann.get("peptide_sequence"), r["coord"].get("charge"))
        if key[0]:
            seen.setdefault(key, set()).add(r["sample_id"])
    recurring = {k: v for k, v in seen.items() if len(v) > 1}

    # The same warning from every run is one fact about the cohort, not N facts.
    # Repeating it buries the warnings that apply to only one run.
    seen_warn = {}
    for d, man in manifests:
        for w in (man.get("warnings") or []):
            seen_warn.setdefault(w, []).append(Path(d).name)
    warnings = []
    for w, where in seen_warn.items():
        warnings.append(w if len(where) == len(manifests)
                        else "%s [%s]" % (w, ", ".join(where)))
    kinds_of_drift = {s["drift_kind"] for s in samples.values()}
    if len(kinds_of_drift) > 1:
        warnings.append(
            "the cohort mixes drift kinds (%s). A drift axis carries one kind, so"
            " the drift panels will exclude the others."
            % ", ".join(sorted(kinds_of_drift)))
    if recurring:
        warnings.append(
            "%d identity(ies) appear in more than one sample. That is the"
            " observation atlas_consensus is about, but this merge does not assign"
            " the tier: how many samples and how close is a rule, and it is not"
            " in this tool." % len(recurring))
    pooled_control, control_note = pool_proximity_control(manifests)
    if control_note:
        warnings.append(control_note)
    unassigned = sum(1 for r in rows["ions"] if not (r.get("annotation") or {}))
    if unassigned == len(rows["ions"]) and rows["ions"]:
        warnings.append(
            "every ion in the cohort is unassigned, so no identity can recur and"
            " no cross-sample agreement can be observed at all.")

    manifest = {
        "export_id": args.export_id or ("cohort-%d-runs" % len(manifests)),
        "created_at": max(m.get("created_at", "") for _, m in manifests),
        "method_id": "merge_exports",
        "method_version": "1.0.0",
        "sources": {"raw_dir": None,
                    "consensus_db": next((m.get("sources", {}).get("consensus_db")
                                          for _, m in manifests
                                          if (m.get("sources") or {}).get("consensus_db")),
                                         None)},
        "pipeline": {"merge_exports": "1.0.0",
                     "merged": ", ".join(str(d) for d, _ in manifests)},
        "samples": [samples[k] for k in sorted(samples)],
        "tables": dict(tables, spectrum_kinds=kinds),
        "counts": {"ions_total": total_before or len(rows["ions"]),
                   "ions_exported": len(rows["ions"]), "by_tier": tiers},
        "proximity_control": pooled_control,
        "decimation": {"strategy": "none", "kept": len(rows["ions"]),
                       "of": total_before or len(rows["ions"]), "seed": 0},
        "warnings": warnings,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    resolve = []
    for r in rows["ions"]:
        ann = r.get("annotation") or {}
        if ann.get("peptide_sequence") and r["coord"].get("drift") is not None \
                and r["coord"].get("charge") is not None:
            resolve.append({"run": r["sample_id"],
                            "uid": "%s_%d" % (ann["peptide_sequence"],
                                              r["coord"]["charge"]),
                            "mz": r["coord"]["mz"], "z": r["coord"]["charge"],
                            "drift": r["coord"]["drift"], "rt": r["coord"]["rt_min"]})
    with (out / "resolve_ions.json").open("w") as fh:
        for r in resolve:
            fh.write(json.dumps(r) + "\n")

    print("merged %d export(s) into %s" % (len(manifests), out))
    print("  %d sample(s): %s" % (len(samples), ", ".join(sorted(samples))))
    for table in TABLES:
        if rows[table]:
            print("  %-14s %d row(s)" % (table, len(rows[table])))
    for tier, n in sorted(tiers.items()):
        print("    %-18s %d" % (tier, n))
    print("  %d identity(ies) recur across samples" % len(recurring))
    print("  %d resolvability observation(s)" % len(resolve))
    for w in warnings:
        print("  warning: %s" % w)
    for f in findings:
        print("  FINDING: %s" % f)
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
