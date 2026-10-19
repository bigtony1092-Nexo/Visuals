#!/usr/bin/env python3
"""Prove the merge does not create agreement that was not observed.

A cohort export is where atlas_consensus would be tempting: the ions are finally
side by side. Promoting a tier because two rows look alike would be the whole
project's failure mode in one line, so these cases check that the merge reports
recurrence and assigns nothing.

    python3 tools/test_merge_exports.py
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "merge_exports.py"
sys.path.insert(0, str(ROOT / "tools"))
from ladder_stats import wilson  # noqa: E402


def run(out, *dirs):
    p = subprocess.run([sys.executable, str(TOOL)] + [str(d) for d in dirs]
                       + ["-o", str(out)], capture_output=True, text=True)
    man = {}
    if (Path(out) / "manifest.json").exists():
        man = json.loads((Path(out) / "manifest.json").read_text())
    return p.returncode, p.stdout + p.stderr, man


def ion(uid, sid, pep="LSIDSQSSSSGPR", z=2, tier="pass2_confirmed",
        drift=0.80, kind="ion_mobility"):
    return {"ion_uid": uid, "sample_id": sid,
            "coord": {"mz": 650.3123, "rt_min": 20.0, "charge": z, "drift": drift,
                      "drift_kind": kind},
            "intensity": {"matched_intensity": 1e6, "summed_precursor_intensity": None,
                          "cluster_rank": 1, "n_consecutive_scans": 7},
            "annotation": None if pep is None else
                          {"peptide_sequence": pep, "prot_accession": "Q06710",
                           "protein_name": None, "seq_start": 11, "mods": None,
                           "frag_type": None},
            "confidence": {"tier": tier, "decoy_controlled": tier != "bare_proximity",
                           "ppm_error": -1.5, "score": 0.001, "score_kind": "q_value"},
            "provenance": {"source_file": sid + ".raw", "vendor": "Waters",
                           "acquisition_method": "HDMSE",
                           "has_ims": kind == "ion_mobility", "search_engine": "ARiA",
                           "batch": "A", "sample_class": "HeLa"},
            "lineage": {"parent_ion_uids": [], "cluster_uid": None, "atlas_uid": None},
            "events": [{"stage": "detected", "at": "2026-10-17T00:00:00Z",
                        "tool": "t", "tool_version": None, "params_hash": "a",
                        "note": None, "coord_snapshot": {"mz": 650.3123,
                                                         "rt_min": 20.0,
                                                         "drift": drift}}]}


def control(ppm=10.0, rt_tol=0.5, axes=("mz", "charge", "rt"),
            ent=(0, 40), tgt=(38, 40)):
    """A proximity_control block in the shape aria_run writes."""
    def arm(name, k_n):
        k, n = k_n
        lo, hi = wilson(k, n)
        return {"arm": name, "queried": n, "matched": k, "rate": k / n,
                "ci": [lo, hi], "distinct_peaks_matched": k,
                "axes_matched": list(axes), "excluded": {}, "strata": []}
    return {"axes_matched": list(axes), "precursor_ppm": ppm, "rt_tol": rt_tol,
            "drift_tol": None, "min_scans": 2, "n_peaks_queryable": 100,
            "mz_range": [400.0, 1200.0], "charges_seen": [2],
            "calibration": {"offset_ppm": -0.8, "n_targets_used": 40,
                            "spread_ppm": 2.1},
            "arms": {"entrapment": arm("entrapment", ent),
                     "target": arm("target", tgt)},
            "cohort_false_rate": ent[0] / ent[1], "findings": []}


def export(tmp, name, ions, kind="ion_mobility", method="HDMSE", warnings=(),
           traces=None, proximity_control=None):
    d = Path(tmp) / name
    d.mkdir(parents=True, exist_ok=True)
    sids = sorted({i["sample_id"] for i in ions})
    tiers = {}
    for i in ions:
        tiers[i["confidence"]["tier"]] = tiers.get(i["confidence"]["tier"], 0) + 1
    man = {"export_id": name, "created_at": "2026-10-17T00:00:00Z", "method_id": "m",
           "samples": [{"sample_id": s, "vendor": "Waters",
                        "acquisition_method": method, "drift_kind": kind,
                        "has_ims": kind == "ion_mobility", "search_engine": "ARiA",
                        "batch": "A", "sample_class": "HeLa"} for s in sids],
           "tables": {"ions": True, "proteins": False, "spectra": False,
                      "traces": bool(traces), "spectrum_kinds": [],
                      "fragment_hits": False},
           "counts": {"ions_total": len(ions), "ions_exported": len(ions),
                      "by_tier": tiers},
           "warnings": list(warnings)}
    if proximity_control is not None:
        man["proximity_control"] = proximity_control
    (d / "manifest.json").write_text(json.dumps(man))
    (d / "ions.jsonl").write_text("".join(json.dumps(i) + "\n" for i in ions))
    if traces:
        (d / "traces.jsonl").write_text("".join(json.dumps(t) + "\n" for t in traces))
    return d


def case_recurrence_is_reported_not_promoted(tmp):
    """The same identity in two samples is an observation, not a tier."""
    a = export(tmp, "m_a", [ion("a" * 16, "CH215")])
    b = export(tmp, "m_b", [ion("b" * 16, "CH216")])
    rc, out, man = run(Path(tmp) / "m_out", a, b)
    bad = []
    if rc != 0:
        bad.append("exit %d on a clean merge" % rc)
    if man.get("counts", {}).get("by_tier", {}).get("atlas_consensus"):
        bad.append("assigned atlas_consensus during a merge")
    for row in json.loads("[" + ",".join(
            (Path(tmp) / "m_out" / "ions.jsonl").read_text().splitlines()) + "]"):
        if row["confidence"]["tier"] == "atlas_consensus":
            bad.append("an ion came out of the merge at atlas_consensus")
        if row["lineage"].get("atlas_uid"):
            bad.append("the merge wrote an atlas_uid")
    if "1 identity(ies) recur" not in out:
        bad.append("did not report the recurrence it did observe")
    if "does not assign the tier" not in out:
        bad.append("did not say it leaves the tier alone")
    return bad


def case_one_sample_id_two_acquisitions_is_refused(tmp):
    a = export(tmp, "c_a", [ion("c" * 16, "CH215")], kind="ion_mobility",
               method="HDMSE")
    b = export(tmp, "c_b", [ion("d" * 16, "CH215", drift=120.0, kind="quad_window")],
               kind="quad_window", method="SONAR")
    rc, out, man = run(Path(tmp) / "c_out", a, b)
    bad = []
    if rc == 0:
        bad.append("merged one sample_id describing two different acquisitions")
    if "not one sample" not in out:
        bad.append("did not say why the two cannot share an id")
    for key in ("acquisition_method", "drift_kind"):
        if key not in out:
            bad.append("did not name %s as a thing they disagree on" % key)
    return bad


def case_a_uid_collision_is_a_finding(tmp):
    a = export(tmp, "e_a", [ion("e" * 16, "CH215")])
    b = export(tmp, "e_b", [ion("e" * 16, "CH216")])
    rc, out, man = run(Path(tmp) / "e_out", a, b)
    bad = []
    if rc == 0:
        bad.append("accepted the same ion_uid from two exports")
    if "hash collision" not in out:
        bad.append("did not identify the repeat as a collision")
    if "not the same ion seen twice" not in out:
        bad.append("did not rule out the benign reading")
    n = len((Path(tmp) / "e_out" / "ions.jsonl").read_text().splitlines())
    if n != 1:
        bad.append("%d ions written; the colliding one must not be merged in" % n)
    return bad


def case_mixed_drift_kinds_are_declared(tmp):
    a = export(tmp, "f_a", [ion("f" * 16, "CH215")])
    b = export(tmp, "f_b", [ion("g" * 16, "CO021", drift=120.0, kind="quad_window")],
               kind="quad_window", method="SONAR")
    rc, out, man = run(Path(tmp) / "f_out", a, b)
    bad = []
    if not any("mixes drift kinds" in w for w in man.get("warnings", [])):
        bad.append("a cohort mixing mobility and quad_window did not say so")
    if not any("will exclude" in w for w in man.get("warnings", [])):
        bad.append("did not say the drift panels will exclude one kind")
    return bad


def case_tables_come_from_the_files(tmp):
    tr = [{"ion_uid": "h" * 16, "sample_id": "CH215", "kind": "xic",
           "points": [[20.0, 5.0], [20.1, 9.0]], "apex_rt_min": 20.1}]
    a = export(tmp, "g_a", [ion("h" * 16, "CH215")], traces=tr)
    b = export(tmp, "g_b", [ion("i" * 16, "CH216")])
    rc, out, man = run(Path(tmp) / "g_out", a, b)
    bad = []
    if not man["tables"]["traces"]:
        bad.append("traces came in from one export but tables says false")
    if man["tables"]["proteins"] or man["tables"]["spectra"]:
        bad.append("declared a table nothing wrote")
    if (Path(tmp) / "g_out" / "proteins.jsonl").exists():
        bad.append("wrote an empty proteins.jsonl")
    return bad


def case_identical_warnings_collapse(tmp):
    w = "every ion is bare_proximity"
    a = export(tmp, "h_a", [ion("j" * 16, "CH215")], warnings=[w, "only in A"])
    b = export(tmp, "h_b", [ion("k" * 16, "CH216")], warnings=[w])
    rc, out, man = run(Path(tmp) / "h_out", a, b)
    bad = []
    hits = [x for x in man["warnings"] if x.startswith(w)]
    if len(hits) != 1:
        bad.append("%d copies of a warning every run made" % len(hits))
    one = [x for x in man["warnings"] if x.startswith("only in A")]
    if not one:
        bad.append("a warning from one run only was dropped")
    elif "h_a" not in one[0]:
        bad.append("a warning from one run only does not say which: %r" % one[0])
    return bad


def case_an_all_unassigned_cohort_says_agreement_is_impossible(tmp):
    a = export(tmp, "i_a", [ion("l" * 16, "CH215", pep=None)])
    b = export(tmp, "i_b", [ion("m" * 16, "CH216", pep=None)])
    rc, out, man = run(Path(tmp) / "i_out", a, b)
    bad = []
    if not any("no identity can recur" in w for w in man.get("warnings", [])):
        bad.append("an all-unassigned cohort did not say agreement is unobservable")
    if "0 identity(ies) recur" not in out:
        bad.append("did not report zero recurrence")
    return bad


def case_nothing_to_merge_exits_two(tmp):
    rc, out, man = run(Path(tmp) / "j_out", Path(tmp) / "no_such_dir")
    bad = []
    if rc != 2:
        bad.append("exit %d when no export could be read" % rc)
    if "nothing to merge" not in out:
        bad.append("did not say nothing could be merged")
    return bad


def case_resolvability_needs_an_identity_and_a_charge(tmp):
    """A chargeless or unassigned ion has no identity to repeat."""
    a = export(tmp, "k_a", [ion("n" * 16, "CH215"),
                            ion("o" * 16, "CH215", pep=None),
                            ion("p" * 16, "CH215", z=None, tier="bare_proximity")])
    b = export(tmp, "k_b", [ion("q" * 16, "CH216")])
    rc, out, man = run(Path(tmp) / "k_out", a, b)
    rows = [json.loads(l) for l in
            (Path(tmp) / "k_out" / "resolve_ions.json").read_text().splitlines()]
    bad = []
    if len(rows) != 2:
        bad.append("%d resolvability row(s); only the two identified, charged ions"
                   " can have a within-ion scale" % len(rows))
    for r in rows:
        if "None" in r["uid"]:
            bad.append("wrote a uid built from a missing charge: %r" % r["uid"])
    return bad


def case_controls_pool_only_when_the_window_is_the_same(tmp):
    """The bare_proximity false-positive rate is a property of a WINDOW.

    Summing probe counts across runs is legitimate when every run queried the
    same window -- they are independent Bernoulli trials against one definition.
    Runs with different tolerances, or one that compared retention time and one
    that did not, measured different things, and their sum would be a rate for a
    window nobody ran.
    """
    bad = []
    a = export(tmp, "pc_a", [ion("a" * 16, "S1")], proximity_control=control())
    b = export(tmp, "pc_b", [ion("b" * 16, "S2")], proximity_control=control())
    rc, out, man = run(Path(tmp) / "pc_same", a, b)
    pc = man.get("proximity_control")
    if not pc:
        bad.append("two runs with the same window produced no pooled control")
    else:
        e = pc["arms"]["entrapment"]
        if (e["queried"], e["matched"]) != (80, 0):
            bad.append("pooled entrapment is %d/%d, expected 0/80"
                       % (e["matched"], e["queried"]))
        if pc["arms"]["target"]["queried"] != 80:
            bad.append("the target comparator did not pool with it")
        if sorted(pc.get("pooled_from") or []) != ["pc_a", "pc_b"]:
            bad.append("the pooled block does not name the runs it came from")

    # Different window: must refuse, and say which runs disagreed.
    c = export(tmp, "pc_c", [ion("c" * 16, "S3")],
               proximity_control=control(ppm=20.0))
    rc2, out2, man2 = run(Path(tmp) / "pc_diff", a, c)
    if man2.get("proximity_control") is not None:
        bad.append("runs with different precursor_ppm were pooled anyway")
    if "DIFFERENT windows" not in out2:
        bad.append("pooling was refused without saying the windows differed")

    # Different axes is the same problem, and the easier one to miss.
    d = export(tmp, "pc_d", [ion("d" * 16, "S4")],
               proximity_control=control(axes=("mz", "charge")))
    rc3, out3, man3 = run(Path(tmp) / "pc_axes", a, d)
    if man3.get("proximity_control") is not None:
        bad.append("a run that never compared rt was pooled with one that did")
    return bad


def case_a_run_with_no_control_leaves_the_cohort_without_one(tmp):
    """A cohort rate cannot be assembled from runs that did not measure one.

    Inventing one from the runs that did would report the measured runs' window
    as the cohort's, which is how the borrowed 33-70% got into this project in
    the first place.
    """
    bad = []
    a = export(tmp, "nc_a", [ion("a" * 16, "S1")], proximity_control=control())
    b = export(tmp, "nc_b", [ion("b" * 16, "S2")])          # never measured
    rc, out, man = run(Path(tmp) / "nc_out", a, b)
    if man.get("proximity_control") is not None:
        bad.append("a cohort control was assembled although one run measured none")
    if "carry no measured" not in out:
        bad.append("the merge did not say which run carried no measured rate")
    if "nc_b" not in out:
        bad.append("the warning does not name the run that was missing it")
    return bad


CASES = [case_controls_pool_only_when_the_window_is_the_same,
         case_a_run_with_no_control_leaves_the_cohort_without_one,
         case_recurrence_is_reported_not_promoted,
         case_one_sample_id_two_acquisitions_is_refused,
         case_a_uid_collision_is_a_finding,
         case_mixed_drift_kinds_are_declared,
         case_tables_come_from_the_files,
         case_identical_warnings_collapse,
         case_an_all_unassigned_cohort_says_agreement_is_impossible,
         case_nothing_to_merge_exits_two,
         case_resolvability_needs_an_identity_and_a_charge]


def main():
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
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
