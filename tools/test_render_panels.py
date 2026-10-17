#!/usr/bin/env python3
"""Prove the page does not draw what is not there, and does not mislabel what is.

A chart is the last place an error can hide, because by then nobody re-reads the
data. So these cases attack the drawing itself: the axis label that makes a
quadrupole position look like mobility, the log axis that swallows a zero
without saying so, the blank panel that gets quietly omitted so the page looks
complete.

    python3 tools/test_render_panels.py
"""

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "render_panels.py"
sys.path.insert(0, str(ROOT / "tools"))
import render_panels  # noqa: E402

PEP = "LSIDSQSSSSGPR"


def ion(uid, mz=650.3123, rt=20.0, drift=0.80, kind="ion_mobility",
        tier="pass2_confirmed", inten=1.2e6, pep=PEP, rank=1, events=None,
        seq_start=11, summed=None):
    return {
        "ion_uid": uid, "sample_id": "CH215",
        "coord": {"mz": mz, "rt_min": rt, "charge": 2, "drift": drift,
                  "drift_kind": kind},
        "intensity": {"matched_intensity": inten,
                      "summed_precursor_intensity": summed,
                      "cluster_rank": rank, "n_consecutive_scans": 7},
        "annotation": None if pep is None else
                      {"peptide_sequence": pep, "prot_accession": "Q06710",
                       "protein_name": None, "seq_start": seq_start,
                       "mods": None, "frag_type": None},
        "confidence": {"tier": tier, "decoy_controlled": tier != "bare_proximity",
                       "ppm_error": -1.5, "score": 0.001, "score_kind": "q_value"},
        "provenance": {"source_file": "CH215.raw", "vendor": "Waters",
                       "acquisition_method": "HDMSE",
                       "has_ims": kind == "ion_mobility", "search_engine": "ARiA",
                       "batch": "A", "sample_class": "HeLa"},
        "lineage": {"parent_ion_uids": [], "cluster_uid": None, "atlas_uid": None},
        "events": events or [{"stage": "detected", "at": "2026-10-17T00:00:00Z",
                              "tool": "test", "tool_version": None,
                              "params_hash": "abc", "note": None,
                              "coord_snapshot": {"mz": mz, "rt_min": rt,
                                                 "drift": drift}}],
    }


def export(tmp, name, ions, kind="ion_mobility", method_name="HDMSE",
           warnings=(), decimation=None, absolute=False, tiers=None,
           n_samples=1):
    d = Path(tmp) / name
    d.mkdir(parents=True, exist_ok=True)
    by_tier = tiers or {}
    if not by_tier:
        for i in ions:
            by_tier[i["confidence"]["tier"]] = by_tier.get(i["confidence"]["tier"], 0) + 1
    samples = [{"sample_id": "CH%d" % (215 + n), "vendor": "Waters",
                "acquisition_method": method_name, "drift_kind": kind,
                "has_ims": kind == "ion_mobility", "search_engine": "ARiA",
                "batch": "A", "sample_class": "HeLa", "role": "target"}
               for n in range(n_samples)]
    man = {"export_id": "test-" + name, "created_at": "2026-10-17T00:00:00Z",
           "method_id": "m", "method_version": "1.0.0", "samples": samples,
           "tables": {"ions": True, "proteins": False, "spectra": False,
                      "traces": False, "spectrum_kinds": [], "fragment_hits": False},
           "counts": {"ions_total": len(ions), "ions_exported": len(ions),
                      "by_tier": by_tier},
           "warnings": list(warnings)}
    if decimation:
        man["decimation"] = decimation
    if absolute:
        man["quant"] = {"metric": "top3avg", "absolute_standard_declared": True}
    (d / "manifest.json").write_text(json.dumps(man))
    (d / "ions.jsonl").write_text("".join(json.dumps(i) + "\n" for i in ions))
    return d


def render(d, *extra):
    p = subprocess.run([sys.executable, str(TOOL), str(d)] + [str(e) for e in extra],
                       capture_output=True, text=True)
    html = ""
    out = d if Path(d).is_dir() else Path(d).parent
    if (out / "panels.html").exists():
        html = (out / "panels.html").read_text()
    return p.returncode, p.stdout + p.stderr, html


def case_sonar_drift_is_never_labelled_mobility(tmp):
    """The axis label IS the claim. 'Mobility' over a quad position is a false figure."""
    d = export(tmp, "sonar", [ion("a" * 16, drift=120.0, kind="quad_window"),
                              ion("b" * 16, drift=140.0, kind="quad_window", rt=24.0)],
               kind="quad_window", method_name="SONAR")
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    # The plan blanks the drift panels outright on a pure-SONAR export, so the
    # guarantee is stronger than a caption: no drift panel is drawn at all.
    for title in ("Drift vs retention time", "m/z vs drift"):
        if title not in html:
            bad.append("%r is not on the page" % title)
            continue
        chunk = html.split(title)[1].split("</div>")[0]
        if "<svg" in chunk:
            bad.append("%r was DRAWN on a SONAR export, where drift is a"
                       " quadrupole position" % title)
        if "quad_window" not in chunk and "quadrupole" not in chunk:
            bad.append("%r blanks without naming the quadrupole position" % title)
    if "quadrupole" not in html:
        bad.append("never names what the drift column actually is")
    # and the label itself, where it does get used
    if render_panels.drift_label("quad_window", False) != \
            "quadrupole window position - NOT mobility":
        bad.append("the quad_window axis label is %r"
                   % render_panels.drift_label("quad_window", False))
    if "mobility" in render_panels.drift_label("quad_window", False) \
            .replace("NOT mobility", ""):
        bad.append("the quad_window label still says mobility")
    return bad


def case_mixed_drift_kinds_refuse_a_shared_axis(tmp):
    d = export(tmp, "mixed", [ion("c" * 16), ion("d" * 16, kind="quad_window", rt=24.0)],
               kind="ion_mobility", n_samples=2)
    man = json.loads((d / "manifest.json").read_text())
    man["samples"][1]["drift_kind"] = "quad_window"
    man["samples"][1]["acquisition_method"] = "SONAR"
    (d / "manifest.json").write_text(json.dumps(man))
    rc, out, html = render(d)
    bad = []
    if "share no axis" not in html:
        bad.append("the banner does not say the kinds share no axis")
    if "other drift kinds excluded" not in html:
        bad.append("the axis does not say that other kinds were excluded")
    if "1 ion(s) excluded" not in html:
        bad.append("did not count the excluded ions: %r"
                   % re.findall(r"\d+ ion\(s\) excluded", html))
    if "different physical quantity" not in html:
        bad.append("does not say why they cannot share the axis")
    # the exclusion must be real, not a caption: one ion drawn, not two
    chunk = html.split("Drift vs retention time")[1].split("</div>")[0]
    if "1 of 1 point(s) drawn" not in chunk:
        bad.append("the drift panel drew %r, so the quad_window ion was plotted"
                   % re.findall(r"\d+ of \d+ point\(s\) drawn", chunk)[:1])
    return bad


def case_a_blank_panel_is_drawn_not_omitted(tmp):
    """A page that renders nothing and a page that renders everything must differ."""
    d = export(tmp, "blank", [ion("e" * 16)])
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "MS/MS head-to-tail" not in html:
        bad.append("a blanked panel was left off the page entirely")
    if "no spectra.jsonl" not in html:
        bad.append("the blank does not say why it is blank")
    if "not drawn" not in html:
        bad.append("the blank is not marked as not drawn")
    n_blank = html.count("blankbox")
    if n_blank < 5:
        bad.append("only %d blank(s) on an ions-only export" % n_blank)
    return bad


def case_points_a_log_axis_cannot_show_are_counted(tmp):
    """A zero on a log axis is invisible, and so is its absence."""
    d = export(tmp, "zeros", [ion("f" * 16, inten=1.2e6),
                              ion("g" * 16, inten=0.0, rt=21.0),
                              ion("h" * 16, inten=None, rt=22.0),
                              ion("i" * 16, inten=-5.0, rt=23.0)])
    rc, out, html = render(d)
    bad = []
    if "1 of 4 point(s) drawn" not in html:
        bad.append("the drawn count does not say 1 of 4: %r"
                   % re.findall(r"\d+ of \d+ point\(s\) drawn", html))
    if "log axis cannot show" not in html:
        bad.append("did not say why the non-positive points are missing")
    return bad


def case_an_all_bare_page_carries_the_decoy_warning(tmp):
    d = export(tmp, "bare", [ion("j" * 16, tier="bare_proximity")])
    rc, out, html = render(d)
    bad = []
    if "33-70%" not in html:
        bad.append("an all-bare_proximity page does not carry the decoy rate")
    if "no decoy control" not in html:
        bad.append("does not say the tier is uncontrolled")
    return bad


def case_mixed_tiers_are_marked_and_announced(tmp):
    d = export(tmp, "tiers", [ion("k" * 16, tier="pass2_confirmed"),
                              ion("l" * 16, tier="bare_proximity", rt=21.0),
                              ion("m" * 16, tier="atlas_consensus", rt=22.0)])
    rc, out, html = render(d)
    bad = []
    if "Several confidence tiers" not in html:
        bad.append("drew three tiers together without saying so")
    colours = {c for _, (c, _) in render_panels.TIER_MARK.items()}
    for colour in colours:
        if colour not in html:
            bad.append("tier colour %s never appears, so a tier is unmarked" % colour)
    shapes = ["<rect", "<polygon", "<circle"]
    missing = [s for s in shapes if s not in html]
    if missing:
        bad.append("tiers are not distinguishable without colour: %r missing"
                   % missing)
    if "legend" not in html:
        bad.append("no legend for the tier marks")
    return bad


def case_quant_is_relative_unless_a_standard_was_declared(tmp):
    d = export(tmp, "rel", [ion("n" * 16)])
    rc, out, html = render(d)
    bad = []
    if "RELATIVE abundance, not concentration" not in html:
        bad.append("a quant axis is not labelled relative")
    if "concentration" in html.replace("not concentration", "") \
            .replace("as a concentration", ""):
        bad.append("the word concentration appears outside its denial")
    d2 = export(tmp, "abs", [ion("o" * 16)], absolute=True)
    rc2, out2, html2 = render(d2)
    if "absolute (standard declared)" not in html2:
        bad.append("a declared absolute standard is not reflected on the axis")
    if "RELATIVE abundance" in html2:
        bad.append("still says relative after a standard was declared")
    return bad


def case_decimation_is_on_the_page(tmp):
    d = export(tmp, "dec", [ion("p" * 16)],
               decimation={"strategy": "grid_density", "kept": 60000,
                           "of": 2140883, "seed": 7})
    rc, out, html = render(d)
    bad = []
    if "grid_density" not in html:
        bad.append("the decimation strategy is not shown")
    if "2.80%" not in html:
        bad.append("the kept fraction is not shown: %r"
                   % re.findall(r"\d+\.\d\d%", html))
    if "not the data" not in html:
        bad.append("does not say the panels show the fraction, not the data")
    return bad


def case_cluster_rank_is_required_not_assumed(tmp):
    """Pooling clusters corrupts the pattern, so a missing rank is excluded."""
    d = export(tmp, "rank", [ion("q" * 16, rank=None, summed=5e5),
                             ion("r" * 16, rank=2, rt=21.0, summed=5e5)])
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "no cluster_rank" not in html:
        bad.append("excluded ions without a cluster_rank silently")
    if "primary-cluster only" not in html:
        bad.append("does not state the panel's own rule")
    if "SeqStart" in html and "not drawn" not in html.split("SeqStart")[1][:400]:
        bad.append("drew a SeqStart panel from non-primary clusters")
    return bad


def case_an_empty_export_refuses_rather_than_drawing_a_blank_page(tmp):
    d = export(tmp, "empty", [])
    (d / "ions.jsonl").write_text("")
    rc, out, html = render(d)
    bad = []
    if rc != 2:
        bad.append("exit %d on an export with no ions" % rc)
    if "nothing to draw" not in out:
        bad.append("did not say there is nothing to draw")
    if (d / "panels.html").exists():
        bad.append("wrote a page for an export with no ions")
    return bad


def case_the_page_needs_no_network(tmp):
    """A viewer that needs a CDN is not a deliverable."""
    d = export(tmp, "offline", [ion("s" * 16)])
    rc, out, html = render(d)
    bad = []
    for ref in re.findall(r'(?:src|href)\s*=\s*"([^"]+)"', html):
        if not ref.startswith("#"):
            bad.append("external reference: %s" % ref)
    for word in ("cdn.", "https://", "http://", "@import", "fonts.googleapis"):
        if word in html:
            bad.append("page contains %r" % word)
    if "<script" in html:
        bad.append("page contains a script tag")
    return bad


def case_untrusted_text_is_escaped(tmp):
    """Peptide names and warnings come from data. Data is not markup."""
    evil = '<script>alert(1)</script>'
    d = export(tmp, "escape", [ion("t" * 16)], warnings=[evil])
    man = json.loads((d / "manifest.json").read_text())
    man["samples"][0]["sample_id"] = evil
    (d / "manifest.json").write_text(json.dumps(man))
    rc, out, html = render(d)
    bad = []
    if "<script>alert" in html:
        bad.append("a warning rendered as live markup")
    if "&lt;script&gt;" not in html:
        bad.append("the text was dropped instead of escaped")
    return bad


def case_the_plan_is_the_single_source_of_truth(tmp):
    """A panel the plan blanks must be blank here, whatever the renderer thinks."""
    d = export(tmp, "truth", [ion("u" * 16, kind="quad_window", drift=120.0)],
               kind="quad_window", method_name="SONAR")
    method = {"method_id": "m", "name": "drift only", "version": "1.0.0",
              "tier_policy": {"min_tier": "pass2_confirmed", "mix_tiers": False},
              "panels": [{"id": "drift_vs_time"}, {"id": "mz_vs_drift"}]}
    mp = Path(tmp) / "m_truth.json"
    mp.write_text(json.dumps(method))
    rc, out, html = render(d, "--method", mp)
    catalog = {p["id"]: p for p in json.loads(
        (ROOT / "schema" / "panels.json").read_text())["panels"]}
    import plan_panels
    bad = []
    for pid in ("drift_vs_time", "mz_vs_drift"):
        status, _ = plan_panels.assess(catalog[pid], method,
                                       json.loads((d / "manifest.json").read_text()))
        title = catalog[pid]["title"]
        if title not in html:
            bad.append("%s is not on the page at all" % pid)
            continue
        chunk = html.split(title)[1][:600]
        drawn = "<svg" in chunk
        if status == plan_panels.BLANK and drawn:
            bad.append("%s is BLANK in the plan but drawn on the page" % pid)
        if status == plan_panels.RENDER and not drawn:
            bad.append("%s is RENDER in the plan but blank on the page" % pid)
    return bad


def case_ion_flight_needs_a_path(tmp):
    one = ion("v" * 16)
    two = ion("w" * 16, rt=21.0)
    two["events"] = two["events"] + [
        {"stage": "merged", "at": "2026-10-17T00:01:00Z", "tool": "test",
         "tool_version": None, "params_hash": "abc", "note": "recalibrated",
         "coord_snapshot": {"mz": 650.3140, "rt_min": 21.4, "drift": 0.802}}]
    d = export(tmp, "flight1", [one])
    rc, out, html = render(d)
    bad = []
    if "no path to fly" not in html:
        bad.append("drew a flight panel from a single coordinate snapshot")
    d2 = export(tmp, "flight2", [one, two])
    rc2, out2, html2 = render(d2)
    chunk = html2.split("Ion flight")[1][:900] if "Ion flight" in html2 else ""
    if "1 of 2 ion(s) moved" not in html2:
        bad.append("does not say how many ions actually moved: %r"
                   % re.findall(r"\d+ of \d+ ion\(s\) moved", html2))
    if "<svg" not in chunk:
        bad.append("an ion with two snapshots did not produce a flight panel")
    return bad


def case_no_method_offers_every_drawable_panel(tmp):
    d = export(tmp, "nomethod", [ion("x" * 16)])
    rc, out, html = render(d)
    bad = []
    if "no --method" not in out:
        bad.append("did not say it was offering all drawable panels")
    for title in ("Intensity vs retention time", "Drift vs retention time",
                  "m/z vs drift", "Forensic track"):
        if title not in html:
            bad.append("%r missing from the no-method page" % title)
    return bad


CASES = [case_sonar_drift_is_never_labelled_mobility,
         case_mixed_drift_kinds_refuse_a_shared_axis,
         case_a_blank_panel_is_drawn_not_omitted,
         case_points_a_log_axis_cannot_show_are_counted,
         case_an_all_bare_page_carries_the_decoy_warning,
         case_mixed_tiers_are_marked_and_announced,
         case_quant_is_relative_unless_a_standard_was_declared,
         case_decimation_is_on_the_page,
         case_cluster_rank_is_required_not_assumed,
         case_an_empty_export_refuses_rather_than_drawing_a_blank_page,
         case_the_page_needs_no_network, case_untrusted_text_is_escaped,
         case_the_plan_is_the_single_source_of_truth,
         case_ion_flight_needs_a_path, case_no_method_offers_every_drawable_panel]


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
