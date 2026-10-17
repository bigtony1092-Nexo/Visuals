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

  Declare a table it did not write. tables{} comes from the files produced.

Exit 0 on a clean merge, 1 on a finding, 2 if nothing could be merged.
"""

import argparse
import json
import sys
from pathlib import Path

TABLES = ("ions", "spectra", "traces", "proteins", "fragment_hits")
SAMPLE_KEYS = ("vendor", "acquisition_method", "drift_kind", "has_ims",
               "search_engine", "batch", "sample_class")


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
