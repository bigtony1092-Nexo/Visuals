#!/usr/bin/env python3
"""Prove every gate in the spine can fail, and that the ion list is honest.

This repository's recurring defect is a check that cannot fail. So each case
here feeds the input that SHOULD stop the run and asserts that it stopped, with
the reason named. A case that passes because the gate is toothless is worse than
no case at all.

    python3 tools/test_aria_run.py
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "aria_run.py"
sys.path.insert(0, str(ROOT / "tools"))
import aria_run  # noqa: E402

PEP = "LSIDSQSSSSGPR"
TARGET = {"peptide_sequence": PEP, "charge": 2, "mz": 650.3133, "rt": 20.0,
          "drift": 0.80, "prot_accession": "Q06710", "seq_start": 11}


def run(*args):
    p = subprocess.run([sys.executable, str(TOOL)] + [str(a) for a in args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def waters(tmp, name, method="HDMSE"):
    d = Path(tmp) / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "_FUNC001.DAT").write_bytes(b"\0" * 8)
    (d / "_FUNC002.DAT").write_bytes(b"\0" * 8)
    if method:
        (d / "_extern.inf").write_text("Acquisition: %s\n" % method)
    return d


def jsonl(path, rows):
    Path(path).write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    return Path(path)


def good_peak(**over):
    r = {"mz": 650.3123, "rt": 20.02, "drift": 0.801, "charge": 2,
         "intensity": 1.2e6, "n_scans": 7, "score_kind": "q_value", "q_value": 0.001}
    r.update(over)
    return r


def case_no_extractor_refuses(tmp):
    """The adapter boundary. Nothing downstream may run on invented peaks."""
    raw = waters(tmp, "a.raw")
    rc, out = run("--raw", raw, "--out", Path(tmp) / "out_a")
    bad = []
    if rc == 0:
        bad.append("exit 0 with no peak source at all")
    if "will not invent peaks" not in out:
        bad.append("did not say it refuses to invent peaks")
    summary = out.split("SUMMARY")[-1]
    for stage in ("3 ANNOTATE", "4 IONLIST", "5 EXPORT", "6 CHECK", "7 PANELS"):
        line = next((l for l in summary.splitlines() if stage in l), None)
        if line is None:
            bad.append("%s is missing from the summary entirely" % stage)
        elif "NOT REACHED" not in line:
            bad.append("%s reads %r, not NOT REACHED"
                       % (stage, line.split(stage)[-1].strip()))
    if (Path(tmp) / "out_a" / "ions.jsonl").exists():
        bad.append("wrote an ion list despite having no peaks")
    return bad


def case_unknown_drift_kind_is_a_gate(tmp):
    """A tolerance on an unknown drift column is how SONAR becomes 'mobility'."""
    raw = waters(tmp, "b.raw", method=None)
    peaks = jsonl(Path(tmp) / "pk_b.jsonl", [good_peak()])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", Path(tmp) / "out_b")
    bad = []
    if rc == 0:
        bad.append("exit 0 with drift_kind unknown")
    if "drift_kind is unknown" not in out:
        bad.append("did not name drift_kind as the blocker")
    if "2 EXTRACT" not in out or "NOT REACHED" not in out:
        bad.append("extract was not marked NOT REACHED")
    return bad


def case_declared_drift_kind_proceeds_and_is_recorded(tmp):
    raw = waters(tmp, "c.raw", method=None)
    peaks = jsonl(Path(tmp) / "pk_c.jsonl", [good_peak()])
    out_dir = Path(tmp) / "out_c"
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--drift-kind", "ion_mobility", "--targets",
                  jsonl(Path(tmp) / "t_c.jsonl", [TARGET]), "--drift-tol", 0.02)
    bad = []
    if not (out_dir / "ions.jsonl").exists():
        bad.append("no ion list written after the drift kind was declared")
        return bad
    ion = json.loads((out_dir / "ions.jsonl").read_text().splitlines()[0])
    if ion["coord"]["drift_kind"] != "ion_mobility":
        bad.append("declared drift kind did not reach the ion record")
    run_rec = json.loads((out_dir / "run.json").read_text())
    vendor = run_rec.get("vendor") or {}
    if not vendor.get("drift_kind_declared"):
        bad.append("run.json does not mark the drift kind as declared, so a"
                   " declaration is indistinguishable from a reading of the files")
    if vendor.get("drift_kind_detected") != "unknown":
        bad.append("run.json records the detected drift kind as %r; the folder had"
                   " no _extern.inf, so it was unknown"
                   % vendor.get("drift_kind_detected"))
    if not any("nothing in the data confirms it" in w.lower()
               for w in run_rec["warnings"]):
        bad.append("the warnings do not say the declaration is unconfirmed: %r"
                   % run_rec["warnings"])
    return bad


def case_extractor_exit_zero_with_no_output(tmp):
    """Exit 0 and no file is a broken adapter, not an empty result."""
    raw = waters(tmp, "d.raw")
    rc, out = run("--raw", raw, "--extract", "true", "--out", Path(tmp) / "out_d")
    bad = []
    if rc != 2:
        bad.append("exit %d; a silent extractor should be COULD NOT RUN (2)" % rc)
    if "not an empty result" not in out:
        bad.append("did not distinguish exit-0-no-output from an empty result")
    return bad


def case_extractor_failure_stops_the_run(tmp):
    raw = waters(tmp, "e.raw")
    rc, out = run("--raw", raw, "--extract", "echo boom >&2; exit 3",
                  "--out", Path(tmp) / "out_e")
    bad = []
    if rc != 2:
        bad.append("exit %d for a failed extractor" % rc)
    if "exited 3" not in out:
        bad.append("did not report the extractor's exit code")
    if "boom" not in out:
        bad.append("swallowed the extractor's own error output")
    return bad


def case_an_extractor_that_works_feeds_the_whole_run(tmp):
    """The adapter's happy path -- the one a real machine uses.

    Both failure modes were covered and this was not, which is how a boundary
    that only ever gets exercised by its errors stops working unnoticed.
    """
    raw = waters(tmp, "ex.raw")
    out_dir = Path(tmp) / "out_ex"
    # Stand in for ApexRT: reads the raw folder, writes a peak list where asked.
    script = Path(tmp) / "fake_extract.py"
    script.write_text(
        "import json, os, sys\n"
        "raw, out = sys.argv[1], sys.argv[2]\n"
        "n = len([f for f in os.listdir(raw) if f.endswith('.DAT')])\n"
        "rows = [{'mz': 650.3123, 'rt': 20.0 + 0.3 * i, 'drift': 0.801,\n"
        "         'charge': 2, 'intensity': 1e6, 'n_scans': 7,\n"
        "         'score_kind': 'q_value', 'q_value': 0.002} for i in range(n)]\n"
        "open(out, 'w').write(''.join(json.dumps(r) + chr(10) for r in rows))\n"
        "print('wrote %d peaks' % len(rows))\n")
    cmd = "%s %s {raw} {peaks}" % (sys.executable, script)
    rc, out = run("--raw", raw, "--extract", cmd, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_ex.jsonl", [TARGET]),
                  "--rt-tol", 3.0, "--drift-tol", 0.02, "--sample-id", "CH215")
    bad = []
    if "wrote 2 peaks" not in out:
        bad.append("the extractor's own stdout was not shown to the operator")
    if "{raw}" in out or "{peaks}" in out:
        bad.append("a placeholder was left unsubstituted in the command")
    if not (out_dir / "peaks.jsonl").exists():
        bad.append("the extractor's peak list was not kept in the export, so the"
                   " run cannot be re-checked against what it actually read")
    if not (out_dir / "ions.jsonl").exists():
        return bad + ["no ion list from a working extractor"]
    ions = [json.loads(l) for l in (out_dir / "ions.jsonl").read_text().splitlines()]
    if len(ions) != 2:
        bad.append("%d ions from 2 extracted peaks" % len(ions))
    if ions and ions[0]["confidence"]["tier"] != "pass2_confirmed":
        bad.append("tier is %r; the extractor supplied a decoy-controlled q-value"
                   % ions[0]["confidence"]["tier"])
    man = json.loads((out_dir / "manifest.json").read_text())
    if script.name not in str(man["pipeline"].get("extract", "")):
        bad.append("the manifest does not record which extractor produced the"
                   " peaks: %r" % man["pipeline"].get("extract"))
    if not (out_dir / "panels.html").exists():
        bad.append("the run did not reach a drawn page")
    if rc != 0:
        bad.append("exit %d on a clean extract-to-page run" % rc)
    return bad


def case_empty_peak_list_is_a_gate(tmp):
    raw = waters(tmp, "f.raw")
    peaks = Path(tmp) / "pk_f.jsonl"
    peaks.write_text("")
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", Path(tmp) / "out_f")
    bad = []
    if rc == 0:
        bad.append("exit 0 on an empty peak list")
    if "no peak records" not in out:
        bad.append("did not say the peak list was empty")
    return bad


def case_peak_without_a_coordinate_cannot_be_identified(tmp):
    raw = waters(tmp, "g.raw")
    peaks = jsonl(Path(tmp) / "pk_g.jsonl", [{"mz": 650.31, "intensity": 1e5}])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", Path(tmp) / "out_g")
    bad = []
    if rc != 2:
        bad.append("exit %d for a peak with no rt or charge" % rc)
    if "has no rt" not in out and "has no charge" not in out:
        bad.append("did not name the missing field")
    return bad


def case_ion_uid_is_content_addressed(tmp):
    """Same detection, same uid. A changed coordinate, a different uid."""
    a = aria_run.ion_uid("CH215", "x.raw", 650.3123, 20.02, 0.801, 2)
    b = aria_run.ion_uid("CH215", "x.raw", 650.31234999, 20.0201, 0.8010499, 2)
    c = aria_run.ion_uid("CH215", "x.raw", 650.3133, 20.02, 0.801, 2)
    d = aria_run.ion_uid("CH216", "x.raw", 650.3123, 20.02, 0.801, 2)
    e = aria_run.ion_uid("CH215", "x.raw", 650.3123, 20.02, None, 2)
    bad = []
    if a != b:
        bad.append("rounding at the stated precision did not give the same uid")
    if a == c:
        bad.append("a 1 mDa shift gave the same uid")
    if a == d:
        bad.append("a different sample gave the same uid")
    if a == e:
        bad.append("null drift and a real drift gave the same uid")
    if len(a) != 16 or any(ch not in "0123456789abcdef" for ch in a):
        bad.append("uid %r does not match the schema pattern" % a)
    return bad


def case_sonar_drift_is_carried_but_never_matched(tmp):
    """SONAR's drift is a quadrupole position. Matching on it is a category error."""
    raw = waters(tmp, "h.raw", method="SONAR")
    out_dir = Path(tmp) / "out_h"
    # drift 400 units from the target's: if drift were a matching dimension this
    # peak would be rejected. It must not be.
    peaks = jsonl(Path(tmp) / "pk_h.jsonl", [good_peak(drift=400.0)])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_h.jsonl", [TARGET]),
                  "--drift-tol", 0.02)
    bad = []
    if not (out_dir / "ions.jsonl").exists():
        bad.append("no ion list written for a SONAR run")
        return bad
    lines = (out_dir / "ions.jsonl").read_text().splitlines()
    if not lines:
        bad.append("SONAR drift was used as a matching dimension and killed the match")
        return bad
    ion = json.loads(lines[0])
    if ion["coord"]["drift_kind"] != "quad_window":
        bad.append("drift_kind is %r for a SONAR run" % ion["coord"]["drift_kind"])
    if ion["coord"]["drift"] != 400.0:
        bad.append("the quadrupole position was not carried through")
    if ion["provenance"]["has_ims"]:
        bad.append("has_ims is true for a SONAR run")
    if "not mobility" not in out:
        bad.append("did not say the drift column is not mobility")
    if ion["coord"].get("one_over_k0") is not None:
        bad.append("wrote a 1/K0 for a non-mobility run")
    return bad


def case_dda_drift_is_null(tmp):
    raw = waters(tmp, "i.raw", method="MSE")   # MSE: no mobility
    out_dir = Path(tmp) / "out_i"
    peaks = jsonl(Path(tmp) / "pk_i.jsonl", [good_peak(drift=0.801)])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_i.jsonl", [TARGET]))
    bad = []
    if not (out_dir / "ions.jsonl").exists():
        return ["no ion list written"]
    ion = json.loads((out_dir / "ions.jsonl").read_text().splitlines()[0])
    if ion["coord"]["drift_kind"] != "absent":
        bad.append("drift_kind is %r on an acquisition with no mobility"
                   % ion["coord"]["drift_kind"])
    if ion["coord"]["drift"] is not None:
        bad.append("carried a drift value of %r on a run with no drift axis"
                   % ion["coord"]["drift"])
    return bad


def case_tiers_follow_the_evidence(tmp):
    raw = waters(tmp, "j.raw")
    out_dir = Path(tmp) / "out_j"
    peaks = jsonl(Path(tmp) / "pk_j.jsonl", [
        good_peak(),                                            # decoy-controlled, under cut
        good_peak(rt=20.5, score_kind=None, q_value=None),      # window only
        good_peak(rt=21.0, q_value=0.4),                        # over the cut
        good_peak(rt=21.5, score_kind="hyperscore", score=38.1, q_value=None),
    ])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_j.jsonl", [TARGET]),
                  "--rt-tol", 3.0, "--drift-tol", 0.02)
    bad = []
    ions = [json.loads(l) for l in (out_dir / "ions.jsonl").read_text().splitlines()]
    tiers = sorted(i["confidence"]["tier"] for i in ions)
    if tiers != ["bare_proximity", "bare_proximity", "pass2_confirmed"]:
        bad.append("tiers are %r; expected one pass2 and two bare" % tiers)
    for i in ions:
        if i["confidence"]["tier"] == "bare_proximity" and i["confidence"]["decoy_controlled"]:
            bad.append("a bare_proximity ion claims decoy control")
        if i["confidence"]["tier"] == "pass2_confirmed" and not i["confidence"]["decoy_controlled"]:
            bad.append("a pass2_confirmed ion has no decoy control")
    if any(i["confidence"]["tier"] == "atlas_consensus" for i in ions):
        bad.append("assigned atlas_consensus from a single run")
    rej = [json.loads(l) for l in (out_dir / "rejected.jsonl").read_text().splitlines()]
    if not any("above the 0.01 cut" in r["reason"] for r in rej):
        bad.append("the q_value above the cut was not rejected with its reason")
    # A hyperscore is a score, but it carries no decoy model.
    hs = [i for i in ions if i["confidence"]["score_kind"] == "hyperscore"]
    if hs and hs[0]["confidence"]["tier"] != "bare_proximity":
        bad.append("a hyperscore lifted an ion above bare_proximity")
    return bad


def case_no_scan_count_is_rejected(tmp):
    """bare_proximity's definition includes the consecutive scans. No count, no tier."""
    raw = waters(tmp, "k.raw")
    out_dir = Path(tmp) / "out_k"
    peaks = jsonl(Path(tmp) / "pk_k.jsonl", [good_peak(n_scans=None),
                                             good_peak(rt=20.3, n_scans=1)])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_k.jsonl", [TARGET]), "--rt-tol", 3.0)
    bad = []
    if rc == 0:
        bad.append("exit 0 with no ion surviving")
    if "no ion survived" not in out:
        bad.append("did not say the ion list came out empty")
    rej = [json.loads(l) for l in (out_dir / "rejected.jsonl").read_text().splitlines()]
    if len(rej) != 2:
        bad.append("%d rejected, expected 2" % len(rej))
    if not any("no scan count" in r["reason"] for r in rej):
        bad.append("a peak with no scan count was not rejected for that reason")
    if not any("below the 2" in r["reason"] for r in rej):
        bad.append("a single-scan peak was not rejected for that reason")
    if (out_dir / "manifest.json").exists():
        bad.append("wrote a manifest for an empty ion list")
    return bad


def case_all_bare_warns_about_the_decoy_rate(tmp):
    raw = waters(tmp, "l.raw")
    out_dir = Path(tmp) / "out_l"
    peaks = jsonl(Path(tmp) / "pk_l.jsonl",
                  [good_peak(score_kind=None, q_value=None)])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_l.jsonl", [TARGET]), "--drift-tol", 0.02)
    bad = []
    if "33-70%" not in out:
        bad.append("an all-bare_proximity export did not carry the decoy warning")
    man = json.loads((out_dir / "manifest.json").read_text())
    if not any("33-70%" in w for w in man["warnings"]):
        bad.append("the warning did not reach manifest.warnings")
    if man["counts"]["by_tier"]["pass2_confirmed"]:
        bad.append("counted a pass2 ion that does not exist")
    return bad


def case_counts_match_what_was_written(tmp):
    raw = waters(tmp, "m.raw")
    out_dir = Path(tmp) / "out_m"
    rows = [good_peak(rt=20.0 + 0.1 * i) for i in range(5)]
    rows.append(good_peak(rt=25.0, n_scans=1))
    peaks = jsonl(Path(tmp) / "pk_m.jsonl", rows)
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_m.jsonl", [TARGET]), "--rt-tol", 6.0)
    man = json.loads((out_dir / "manifest.json").read_text())
    n = len((out_dir / "ions.jsonl").read_text().splitlines())
    bad = []
    if man["counts"]["ions_exported"] != n:
        bad.append("manifest says %d exported, file holds %d"
                   % (man["counts"]["ions_exported"], n))
    if man["counts"]["ions_total"] != len(rows):
        bad.append("ions_total %d, peak list held %d"
                   % (man["counts"]["ions_total"], len(rows)))
    if sum(man["counts"]["by_tier"].values()) != n:
        bad.append("by_tier sums to %d, not %d"
                   % (sum(man["counts"]["by_tier"].values()), n))
    for name, present in man["tables"].items():
        if name == "spectrum_kinds":
            continue
        if present and not (out_dir / ("%s.jsonl" % name)).exists():
            bad.append("tables.%s is true but %s.jsonl is not there" % (name, name))
    return bad


def case_resolvability_is_not_fed_from_one_observation(tmp):
    """Within-ion sigma needs repeats. One run cannot supply them."""
    raw = waters(tmp, "n.raw")
    out_dir = Path(tmp) / "out_n"
    peaks = jsonl(Path(tmp) / "pk_n.jsonl", [good_peak()])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_n.jsonl", [TARGET]), "--drift-tol", 0.02)
    bad = []
    if "resolvability NOT fed" not in out:
        bad.append("fed the collision step a scale it could not estimate")
    if "measured 3 times" not in out:
        bad.append("did not say why the scale could not be estimated")
    if "4 RESOLVE" in out and "sigma_drift" in out:
        bad.append("the collision step reported sigmas it had no basis for")
    return bad


def case_mobility_without_a_tolerance_says_so(tmp):
    raw = waters(tmp, "o.raw", method="HDMSE")
    out_dir = Path(tmp) / "out_o"
    peaks = jsonl(Path(tmp) / "pk_o.jsonl", [good_peak()])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_o.jsonl", [TARGET]))
    bad = []
    if "sat unused" not in out:
        bad.append("a mobility run matched without --drift-tol and said nothing")
    if "drift unused" not in out:
        bad.append("did not state in the window line that drift was unused")
    return bad


def case_the_written_ions_satisfy_the_contract(tmp):
    """The spine's output must pass the repo's own validator, not just look right."""
    raw = waters(tmp, "p.raw")
    out_dir = Path(tmp) / "out_p"
    # Two identities, four detections each, with jitter: the shape resolvability
    # can actually measure a within-ion scale from. One identity cannot, because
    # then the axis range is that ion's own scatter.
    t2 = dict(TARGET, peptide_sequence="QHYPEAAGPAR", mz=700.3500, rt=30.0,
              drift=0.90, prot_accession="P04406")
    rows = []
    for i in range(4):
        rows.append(good_peak(rt=20.0 + 0.03 * i, drift=0.801 + 0.004 * i,
                              one_over_k0=0.801 + 0.004 * i))
        rows.append(good_peak(mz=700.3493, rt=30.0 + 0.03 * i,
                              drift=0.901 + 0.004 * i, one_over_k0=0.901 + 0.004 * i))
    peaks = jsonl(Path(tmp) / "pk_p.jsonl", rows)
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_p.jsonl", [TARGET, t2]), "--rt-tol", 3.0,
                  "--drift-tol", 0.05, "--sample-id", "CH215")
    bad = []
    v = subprocess.run([sys.executable, str(ROOT / "tools" / "validate_export.py"),
                        str(out_dir)], capture_output=True, text=True)
    body = v.stdout + v.stderr
    if v.returncode != 0:
        bad.append("validate_export exited %d on the spine's own export" % v.returncode)
        for line in [l for l in body.splitlines() if "ERROR" in l][:4]:
            bad.append("  %s" % line.strip())
    if "ERROR" in body:
        bad.append("the export carries contract errors")
    if rc != 0:
        bad.append("the run exited %d on input that should pass cleanly" % rc)
    return bad


def case_duplicate_detections_merge_on_their_uid(tmp):
    raw = waters(tmp, "q.raw")
    out_dir = Path(tmp) / "out_q"
    peaks = jsonl(Path(tmp) / "pk_q.jsonl", [good_peak(), good_peak()])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_q.jsonl", [TARGET]), "--drift-tol", 0.02)
    lines = (out_dir / "ions.jsonl").read_text().splitlines()
    bad = []
    if len(lines) != 1:
        bad.append("%d ions written for one detection listed twice" % len(lines))
    if "DROPPED, not combined" not in out:
        bad.append("dropped the duplicate silently, or claimed to have merged it")
    ion = json.loads(lines[0]) if lines else {}
    if lines and ion["intensity"]["matched_intensity"] != 1.2e6:
        bad.append("intensity is %r: the duplicate was summed in, which is a merge"
                   " the tool does not do" % ion["intensity"]["matched_intensity"])
    return bad


def case_events_are_ordered_and_append_only(tmp):
    raw = waters(tmp, "r.raw")
    out_dir = Path(tmp) / "out_r"
    peaks = jsonl(Path(tmp) / "pk_r.jsonl", [good_peak()])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_r.jsonl", [TARGET]), "--drift-tol", 0.02)
    ion = json.loads((out_dir / "ions.jsonl").read_text().splitlines()[0])
    bad = []
    stages = [e["stage"] for e in ion["events"]]
    if stages[0] != "detected":
        bad.append("the first event is %r, not detected" % stages[0])
    if "annotated" not in stages:
        bad.append("an assigned ion carries no annotated event")
    if stages != sorted(stages, key=lambda s: ["detected", "annotated"].index(s)
                        if s in ("detected", "annotated") else 9):
        bad.append("events are out of order: %r" % stages)
    first = ion["events"][0]
    if not first.get("coord_snapshot") or first["coord_snapshot"]["mz"] is None:
        bad.append("the detected event carries no coordinate snapshot")
    if not first.get("params_hash"):
        bad.append("no params_hash, so the run is not reproducible from the ion")
    return bad


def case_ambiguous_matches_are_declared(tmp):
    """Nearest-in-ppm is a tie-break. It must never read as evidence."""
    raw = waters(tmp, "s.raw")
    out_dir = Path(tmp) / "out_s"
    t2 = dict(TARGET, peptide_sequence="QHYPEAAGPAR", mz=650.3135,
              prot_accession="P04406")
    peaks = jsonl(Path(tmp) / "pk_s.jsonl", [good_peak()])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_s.jsonl", [TARGET, t2]),
                  "--drift-tol", 0.02)
    bad = []
    if "more than one candidate" not in out:
        bad.append("did not report that a peak had several candidates")
    man = json.loads((out_dir / "manifest.json").read_text())
    if not any("tie-break, not" in w for w in man["warnings"]):
        bad.append("the manifest does not record that a tie-break was used")
    return bad


def case_no_targets_means_no_identifications(tmp):
    raw = waters(tmp, "t.raw")
    out_dir = Path(tmp) / "out_t"
    peaks = jsonl(Path(tmp) / "pk_t.jsonl", [good_peak()])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir, "--drift-tol", 0.02)
    bad = []
    ion = json.loads((out_dir / "ions.jsonl").read_text().splitlines()[0])
    if ion["annotation"] is not None:
        bad.append("produced an annotation with no target list")
    if ion["confidence"]["tier"] != "bare_proximity":
        bad.append("tier is %r with nothing to match against"
                   % ion["confidence"]["tier"])
    if "nothing here is an identification" not in out:
        bad.append("did not say that no identification was made")
    return bad


def case_peg_and_z1_filters_are_recorded_not_silent(tmp):
    raw = waters(tmp, "u.raw")
    out_dir = Path(tmp) / "out_u"
    # Doubly charged PEG, so the z1 filter cannot consume it first and the PEG
    # filter is the only thing that can reject it.
    peg2 = (44.026215 * 22 + 18.010565 + 2 * 1.007276) / 2
    peaks = jsonl(Path(tmp) / "pk_u.jsonl", [
        good_peak(),
        good_peak(mz=peg2, charge=2, rt=12.0),
        good_peak(mz=700.4000, charge=1, rt=13.0)])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--drop-z1", "--drop-peg",
                  "--targets", jsonl(Path(tmp) / "t_u.jsonl", [TARGET]), "--drift-tol", 0.02)
    rej = [json.loads(l) for l in (out_dir / "rejected.jsonl").read_text().splitlines()]
    stages = sorted(r["events"][-1]["stage"] for r in rej)
    bad = []
    if stages != ["peg_filtered", "z1_filtered"]:
        bad.append("rejected stages are %r; expected one of each" % stages)
    for r in rej:
        if not r["events"][-1].get("note"):
            bad.append("a filtered ion carries no reason on its last event")
    if not aria_run.looks_like_peg(peg2, 2):
        bad.append("the PEG lattice test misses a doubly charged PEG ion")
    # The end group is the part that is easy to leave out, so test a lattice
    # position that is only right if the water is included.
    if aria_run.looks_like_peg(44.026215 * 8 + 1.007276, 1):
        bad.append("the PEG test matches the lattice WITHOUT its water end group,"
                   " which is a lattice no real PEG ion sits on")
    if aria_run.looks_like_peg(650.3123, 2):
        bad.append("the PEG lattice test flags an ordinary peptide")
    if aria_run.looks_like_peg(44.026215 * 2 + 18.010565 + 1.007276, 1):
        bad.append("flagged a 2-unit oligomer; the filter claims n >= 3")
    return bad


def case_csv_peaks_are_read(tmp):
    raw = waters(tmp, "v.raw")
    out_dir = Path(tmp) / "out_v"
    (Path(tmp) / "pk_v.csv").write_text(
        "m/z,RT,Charge,Intensity,n_scans\n650.3123,20.02,2,1200000,7\n")
    rc, out = run("--raw", raw, "--peaks", Path(tmp) / "pk_v.csv", "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_v.jsonl", [TARGET]))
    bad = []
    if rc not in (0, 1):
        bad.append("exit %d reading a CSV peak list" % rc)
    if not (out_dir / "ions.jsonl").exists():
        return bad + ["no ion written from a CSV peak list"]
    ion = json.loads((out_dir / "ions.jsonl").read_text().splitlines()[0])
    if abs(ion["coord"]["mz"] - 650.3123) > 1e-6:
        bad.append("CSV m/z came through as %r" % ion["coord"]["mz"])
    if ion["intensity"]["matched_intensity"] != 1200000.0:
        bad.append("CSV intensity came through as %r"
                   % ion["intensity"]["matched_intensity"])
    return bad


def case_a_score_may_not_certify_a_different_sequence(tmp):
    """The target list selects ions. It does not identify them.

    If the row's own search result says VTIMPKDIQLAR and the target list says
    LSIDSQSSSSGPR, exporting the target's sequence with the row's q-value writes
    a decoy-controlled identification for a peptide nothing scored.
    """
    raw = waters(tmp, "w.raw")
    out_dir = Path(tmp) / "out_w"
    peaks = jsonl(Path(tmp) / "pk_w.jsonl", [
        good_peak(peptide="VTIMPKDIQLAR", protein="P12345")])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_w.jsonl", [TARGET]), "--drift-tol", 0.02)
    bad = []
    ion = json.loads((out_dir / "ions.jsonl").read_text().splitlines()[0])
    ann = ion["annotation"] or {}
    if ann.get("peptide_sequence") != "VTIMPKDIQLAR":
        bad.append("exported %r, but the score was computed for VTIMPKDIQLAR"
                   % ann.get("peptide_sequence"))
    if ann.get("prot_accession") != "P12345":
        bad.append("exported protein %r, not the row's own P12345"
                   % ann.get("prot_accession"))
    if ion["confidence"]["tier"] == "pass2_confirmed" and \
            ann.get("peptide_sequence") != "VTIMPKDIQLAR":
        bad.append("promoted an ion to pass2_confirmed under a sequence its score"
                   " does not belong to")
    if "proposed" not in json.dumps(ion["events"]):
        bad.append("the disagreement left no trace on the ion's events")
    man = json.loads((out_dir / "manifest.json").read_text())
    if not any("own identification" in w for w in man["warnings"]):
        bad.append("the manifest does not record the conflict")
    return bad


def case_tier_definitions_cannot_be_loosened_by_flag(tmp):
    """A gate a flag can switch off is not a gate."""
    raw = waters(tmp, "x.raw")
    peaks = jsonl(Path(tmp) / "pk_x.jsonl", [good_peak(n_scans=1, q_value=0.95)])
    targets = jsonl(Path(tmp) / "t_x.jsonl", [TARGET])
    bad = []
    for flag, value in (("--min-scans", "0"), ("--min-scans", "1"),
                        ("--q-cut", "1.0"), ("--q-cut", "0.5"),
                        ("--precursor-ppm", "5000")):
        rc, out = run("--raw", raw, "--peaks", peaks, "--out", Path(tmp) / "out_x",
                      "--targets", targets, "--drift-tol", 0.02, flag, value)
        if rc == 0:
            bad.append("%s %s was accepted" % (flag, value))
        if "error" not in out.lower() and "usage" not in out.lower():
            bad.append("%s %s was refused without saying why" % (flag, value))
    return bad


def case_an_unknown_vendor_is_not_laundered_by_a_declaration(tmp):
    """--drift-kind declares the drift axis. It says nothing about the vendor."""
    mixed = Path(tmp) / "mixed_folder"
    mixed.mkdir(exist_ok=True)
    (mixed / "a.raw").write_bytes(b"\x01" * 8)
    (mixed / "b.d").mkdir(exist_ok=True)
    (mixed / "b.d" / "analysis.tdf").write_bytes(b"\0" * 8)
    peaks = jsonl(Path(tmp) / "pk_y.jsonl", [good_peak()])
    out_dir = Path(tmp) / "out_y"
    rc, out = run("--raw", mixed, "--peaks", peaks, "--out", out_dir,
                  "--drift-kind", "ion_mobility")
    bad = []
    if rc == 0:
        bad.append("exit 0 on a mixed-vendor folder with a declared drift kind")
    if "vendor could not be established" not in out and "more than one vendor" not in out:
        bad.append("did not name the vendor as the blocker")
    if (out_dir / "ions.jsonl").exists():
        bad.append("wrote ions for a folder whose vendor is unknown")
    # and a folder of runs: the method is per run and must not be assumed shared
    folder = Path(tmp) / "runs_folder"
    folder.mkdir(exist_ok=True)
    for n in ("r1.raw", "r2.raw"):
        (folder / n).mkdir(exist_ok=True)
        (folder / n / "_FUNC001.DAT").write_bytes(b"\0" * 8)
    rc2, out2 = run("--raw", folder, "--peaks", peaks, "--out", Path(tmp) / "out_y2",
                    "--drift-kind", "ion_mobility")
    if rc2 == 0:
        bad.append("exit 0 on a folder of runs with a declared drift kind")
    if "folder of runs" not in out2:
        bad.append("did not refuse a folder of runs")
    return bad


def case_nothing_matching_stops_the_run_without_crashing(tmp):
    raw = waters(tmp, "z.raw")
    out_dir = Path(tmp) / "out_z"
    peaks = jsonl(Path(tmp) / "pk_z.jsonl", [good_peak(mz=1200.5000)])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_z.jsonl", [TARGET]), "--drift-tol", 0.02)
    bad = []
    if "Traceback" in out:
        bad.append("crashed instead of reporting: %s"
                   % out.split("Traceback")[-1].strip().splitlines()[-1][:80])
    if "SUMMARY" not in out:
        bad.append("never printed a summary")
    if not (out_dir / "run.json").exists():
        bad.append("wrote no run.json, so the run left no record of why it stopped")
    if rc == 0:
        bad.append("exit 0 when no peak matched any target")
    if "no peak matched any target" not in out:
        bad.append("did not say that nothing matched")
    return bad


def case_targets_and_atlas_together_are_refused(tmp):
    raw = waters(tmp, "aa.raw")
    peaks = jsonl(Path(tmp) / "pk_aa.jsonl", [good_peak()])
    out_dir = Path(tmp) / "out_aa"
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_aa.jsonl", [TARGET]),
                  "--atlas", Path(tmp) / "nonexistent.db", "--drift-tol", 0.02)
    bad = []
    if rc == 0:
        bad.append("accepted two target lists at once")
    if "false audit trail" not in out:
        bad.append("did not say why it refuses: the manifest would name an unopened DB")
    if (out_dir / "manifest.json").exists():
        man = json.loads((out_dir / "manifest.json").read_text())
        if man["sources"].get("consensus_db"):
            bad.append("wrote a manifest naming a consensus DB it never opened")
    return bad




def case_an_undetermined_charge_is_tracked_not_invented(tmp):
    """Real case, from the HeLa wash landmark at m/z 352.340.

    An MS1 detection whose isotope envelope was not resolved has no determined
    charge. The contract tracks unassigned ions on purpose, so refusing the
    measurement is wrong -- and writing a guessed charge to get it in is worse.
    """
    raw = waters(tmp, "nz.raw")
    out_dir = Path(tmp) / "out_nz"
    rows = [{"mz": 352.340, "rt": 50.047, "intensity": 68782.0, "n_scans": 154},
            {"mz": 352.340, "rt": 50.528, "intensity": 61077.0, "n_scans": 164,
             "score_kind": "q_value", "q_value": 0.001}]
    peaks = jsonl(Path(tmp) / "pk_nz.jsonl", rows)
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--sample-id", "HeLa_60_r1",
                  "--targets", jsonl(Path(tmp) / "t_nz.jsonl", [TARGET]))
    bad = []
    if rc != 0:
        bad.append("exit %d: a real detection with no determined charge was"
                   " refused rather than tracked" % rc)
    if not (out_dir / "ions.jsonl").exists():
        return bad + ["no ion list written"]
    ions = [json.loads(l) for l in (out_dir / "ions.jsonl").read_text().splitlines()]
    if len(ions) != 2:
        bad.append("%d ion(s) written from 2 chargeless peaks" % len(ions))
    for i in ions:
        if i["coord"]["charge"] is not None:
            bad.append("charge came out as %r; nothing measured it"
                       % i["coord"]["charge"])
        if i["annotation"] is not None:
            bad.append("annotated an ion whose charge is unknown: m/z cannot be"
                       " turned into a peptide mass without z")
        if i["confidence"]["tier"] != "bare_proximity":
            bad.append("tier %r on an ion with no charge; a decoy-controlled score"
                       " certifies an identification and there is none without z"
                       % i["confidence"]["tier"])
    # the second peak carried a q-value under the cut and must NOT be promoted
    scored = [i for i in ions if i["confidence"]["score_kind"]
              and "q_value" in i["confidence"]["score_kind"]]
    if scored and scored[0]["confidence"]["decoy_controlled"]:
        bad.append("a chargeless ion claims decoy control")
    if "no determined charge" not in out:
        bad.append("said nothing about the missing charge")
    if "never be identified" not in out:
        bad.append("did not say what the missing charge costs")
    # distinct uids despite identical m/z, so the null charge is in the address
    if len({i["ion_uid"] for i in ions}) != 2:
        bad.append("the two detections collapsed to one uid")
    if aria_run.ion_uid("S", "f", 352.34, 50.047, None, None) == \
            aria_run.ion_uid("S", "f", 352.34, 50.047, None, 1):
        bad.append("a null charge and charge 1 give the same address; null is not 1")
    return bad


def case_the_extraction_window_is_never_the_matching_window(tmp):
    """--precursor-ppm is what the run matched in, not what the trace came through."""
    raw = waters(tmp, "tw.raw")
    out_dir = Path(tmp) / "out_tw"
    rows = [{"mz": 650.3123, "rt": 20.0, "charge": 2, "intensity": 1e6,
             "n_scans": 7, "points": [[19.9, 3.0], [20.0, 9.0], [20.1, 4.0]]}]
    peaks = jsonl(Path(tmp) / "pk_tw.jsonl", rows)
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--precursor-ppm", "10",
                  "--targets", jsonl(Path(tmp) / "t_tw.jsonl", [TARGET]))
    bad = []
    if not (out_dir / "traces.jsonl").exists():
        return ["no trace written from a peak carrying elution points"]
    t = json.loads((out_dir / "traces.jsonl").read_text().splitlines()[0])
    if t.get("extraction_window_ppm") is not None:
        bad.append("wrote extraction_window_ppm %r, which no input stated"
                   % t["extraction_window_ppm"])
    if len(t["points"]) != 3:
        bad.append("%d point(s) carried through" % len(t["points"]))
    if t.get("apex_rt_min") != 20.0:
        bad.append("apex came out at %r, not the most intense point"
                   % t.get("apex_rt_min"))
    man = json.loads((out_dir / "manifest.json").read_text())
    if not man["tables"]["traces"]:
        bad.append("wrote traces.jsonl but tables.traces is false")
    if "not the window the extractor" not in out:
        bad.append("did not distinguish the matching window from the extraction one")
    # and a stated window IS carried
    rows[0]["extraction_window_ppm"] = 7.1
    peaks = jsonl(Path(tmp) / "pk_tw2.jsonl", rows)
    rc2, out2 = run("--raw", raw, "--peaks", peaks, "--out", Path(tmp) / "out_tw2",
                    "--targets", jsonl(Path(tmp) / "t_tw.jsonl", [TARGET]))
    t2 = json.loads((Path(tmp) / "out_tw2" / "traces.jsonl").read_text().splitlines()[0])
    if t2.get("extraction_window_ppm") != 7.1:
        bad.append("a stated extraction window was not carried: %r"
                   % t2.get("extraction_window_ppm"))
    return bad


def case_a_one_point_trace_is_not_written(tmp):
    raw = waters(tmp, "op.raw")
    out_dir = Path(tmp) / "out_op"
    peaks = jsonl(Path(tmp) / "pk_op.jsonl", [
        {"mz": 650.3123, "rt": 20.0, "charge": 2, "intensity": 1e6, "n_scans": 7,
         "points": [[20.0, 9.0]]}])
    rc, out = run("--raw", raw, "--peaks", peaks, "--out", out_dir,
                  "--targets", jsonl(Path(tmp) / "t_op.jsonl", [TARGET]))
    bad = []
    if (out_dir / "traces.jsonl").exists():
        bad.append("wrote a trace from a single point")
    if "One point is not a peak" not in out:
        bad.append("did not say why the trace was not written")
    man = json.loads((out_dir / "manifest.json").read_text())
    if man["tables"]["traces"]:
        bad.append("declared traces with no traces.jsonl")
    return bad


CASES = [case_no_extractor_refuses, case_unknown_drift_kind_is_a_gate,
         case_declared_drift_kind_proceeds_and_is_recorded,
         case_extractor_exit_zero_with_no_output, case_extractor_failure_stops_the_run,
         case_an_extractor_that_works_feeds_the_whole_run,
         case_empty_peak_list_is_a_gate, case_peak_without_a_coordinate_cannot_be_identified,
         case_ion_uid_is_content_addressed, case_sonar_drift_is_carried_but_never_matched,
         case_dda_drift_is_null, case_tiers_follow_the_evidence,
         case_no_scan_count_is_rejected, case_all_bare_warns_about_the_decoy_rate,
         case_counts_match_what_was_written,
         case_resolvability_is_not_fed_from_one_observation,
         case_mobility_without_a_tolerance_says_so,
         case_the_written_ions_satisfy_the_contract,
         case_duplicate_detections_merge_on_their_uid,
         case_events_are_ordered_and_append_only, case_ambiguous_matches_are_declared,
         case_no_targets_means_no_identifications,
         case_peg_and_z1_filters_are_recorded_not_silent, case_csv_peaks_are_read,
         case_a_score_may_not_certify_a_different_sequence,
         case_tier_definitions_cannot_be_loosened_by_flag,
         case_an_unknown_vendor_is_not_laundered_by_a_declaration,
         case_nothing_matching_stops_the_run_without_crashing,
         case_targets_and_atlas_together_are_refused,
         case_an_undetermined_charge_is_tracked_not_invented,
         case_the_extraction_window_is_never_the_matching_window,
         case_a_one_point_trace_is_not_written]


def main():
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        for case in CASES:
            try:
                bad = case(tmp)
            except Exception as exc:            # a crashing case is a failing case
                bad = ["%s: %s" % (type(exc).__name__, exc)]
            failed += bool(bad)
            print("%s %s" % ("ok  " if not bad else "FAIL", case.__name__))
            for line in bad:
                print("     %s" % line)
    print("\n%d/%d passed" % (len(CASES) - failed, len(CASES)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
