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
    page = ""
    out = d if Path(d).is_dir() else Path(d).parent
    if (out / "panels.html").exists():
        # Entities are how the page stores text; a test that matches on prose
        # should read the prose, not the encoding.
        page = __import__("html").unescape((out / "panels.html").read_text())
    return p.returncode, p.stdout + p.stderr, page


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


def case_an_all_bare_page_says_whether_the_rate_was_measured(tmp):
    """The caption used to quote 33-70% from elsewhere in this project, which
    described retired code. It now reads the rate out of the export, so a page
    built from an export that never measured it has to say so rather than
    borrow a number from a different run."""
    d = export(tmp, "bare", [ion("j" * 16, tier="bare_proximity")])
    rc, out, html = render(d)
    bad = []
    if "was NOT measured on this" not in html:
        bad.append("the export carries no proximity_control and the page did not"
                   " say the window's false-positive rate was unmeasured")
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
    if "Nothing here to draw" not in out:
        bad.append("did not say there is nothing to draw")
    if (d / "panels.html").exists():
        bad.append("wrote a page for an export with no ions")

    # A directory holding quant but no ion list is a different refusal: panels
    # COULD be drawn from it, and the contract is the reason they are not.
    d2 = export(tmp, "quantonly", [])
    (d2 / "ions.jsonl").write_text("")
    write(d2, "proteins.jsonl", [
        {"sample_id": "CH215", "prot_accession": "Q99836", "protein_name": "M",
         "top3avg": 7.209, "n_peptides_used": 3, "n_peptides_matched": 5,
         "peptides_used": ["X"], "rank": None, "percentile": 94.0,
         "normalized": None, "tier": "bare_proximity", "decoy_controlled": False}])
    declare(d2, proteins=True)
    rc2, out2, _ = render(d2)
    if rc2 != 2:
        bad.append("exit %d on a quant-only directory" % rc2)
    if "proteins.jsonl" not in out2:
        bad.append("did not say what the directory does hold")
    if "const true" not in out2:
        bad.append("did not name the contract rule that makes this a refusal")
    if "no evidence trail" not in out2:
        bad.append("did not say why quant without the ion list is not an export")
    if "Nothing here to draw" in out2:
        bad.append("claimed there was nothing to draw while holding a quant table")
    if (d2 / "panels.html").exists():
        bad.append("wrote a page from quant with no ion list behind it")
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
    """Peptide names and warnings come from data. Data is not markup.

    This one reads the file as written, not unescaped: render() unescapes so the
    other cases can match on prose, and a test for escaping that reads the
    unescaped page is a test that cannot pass.
    """
    evil = '<script>alert(1)</script>'
    d = export(tmp, "escape", [ion("t" * 16)], warnings=[evil])
    man = json.loads((d / "manifest.json").read_text())
    man["samples"][0]["sample_id"] = evil
    (d / "manifest.json").write_text(json.dumps(man))
    render(d)
    raw = (d / "panels.html").read_text()
    bad = []
    if "<script>alert" in raw:
        bad.append("a warning rendered as live markup")
    if "&lt;script&gt;" not in raw:
        bad.append("the text was dropped instead of escaped")
    # and the only <script in the file should be the escaped one
    if raw.count("<script") != 0:
        bad.append("%d live script tag(s) in the page" % raw.count("<script"))
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




# ---------------------------------------------------------------- new panels

def spec(uid, kind, peaks, match=None, charge=2):
    r = {"ion_uid": uid, "sample_id": "CH215", "kind": kind,
         "precursor_mz": 745.8821, "charge": charge,
         "scan_id": "CH215:%s:1" % kind, "peaks": peaks}
    if match:
        r["match"] = match
    return r


def pk(mz, inten, ion=None, ppm=None, idx=None):
    return {"mz": mz, "intensity": inten, "ion": ion, "idx": idx,
            "charge": 1, "ppm_error": ppm, "ar1": None}


def write(d, name, rows):
    (Path(d) / name).write_text("".join(json.dumps(r) + "\n" for r in rows))


def declare(d, **flags):
    man = json.loads((Path(d) / "manifest.json").read_text())
    man["tables"].update(flags)
    (Path(d) / "manifest.json").write_text(json.dumps(man))


def case_a_one_sided_mirror_is_not_a_comparison(tmp):
    """Drawing experimental alone against nothing is not a head-to-tail plot."""
    d = export(tmp, "mirror1", [ion("a1" + "0" * 14)])
    write(d, "spectra.jsonl", [spec("a1" + "0" * 14, "ms2_experimental",
                                    [pk(262.1, 18.4, "b2"), pk(449.2, 100.0, "y4")])])
    # The manifest claims both kinds. Only one is in the file -- the lie the plan
    # cannot catch, because the plan reads the manifest.
    declare(d, spectra=True, spectrum_kinds=["ms2_experimental", "ms2_library"])
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    chunk = html.split("MS/MS head-to-tail")[1].split("</div>")[0] \
        if "MS/MS head-to-tail" in html else ""
    if "<svg" in chunk:
        bad.append("drew a mirror with no library side")
    if "no row in spectra.jsonl carries them" not in html:
        bad.append("a declared spectrum kind with no rows was not reported")
    if "ms2_library" not in html:
        bad.append("did not name the missing kind")
    # and the renderer's own reason, for an export that declares honestly
    declare(d, spectra=True, spectrum_kinds=["ms2_experimental"])
    rc2, out2, html2 = render(d, "--method", ROOT / "examples" / "method.example.json")
    if "missing kind(s): ms2_library" not in html2:
        bad.append("an honest manifest missing the library kind was not blanked"
                   " with that reason")
    return bad


def case_the_mirror_reports_what_it_did_not_explain(tmp):
    uid = "b1" + "0" * 14
    d = export(tmp, "mirror2", [ion(uid)])
    exp = [pk(262.1292, 18.4, "b2", 1.9, 2), pk(449.2245, 100.0, "y4", 0.8, 4),
           pk(512.2399, 60.0, None, None, None)]
    lib = [pk(262.1287, 22.0, "b2", None, 2), pk(449.2241, 100.0, "y4", None, 4),
           pk(846.4705, 12.3, "y7", None, 7)]
    write(d, "spectra.jsonl", [
        spec(uid, "ms2_experimental", exp,
             {"dot_product": 0.914, "basis": "raw_intensity", "n_matched": 2,
              "n_library": 3, "library_source": "atlas-2026-10-09"}),
        spec(uid, "ms2_library", lib)])
    declare(d, spectra=True, spectrum_kinds=["ms2_experimental", "ms2_library"])
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "<svg" not in html.split("MS/MS head-to-tail")[1][:600]:
        return ["the mirror did not draw with both sides present"]
    if "34% of the intensity" not in html and "33% of the intensity" not in html:
        bad.append("the unexplained intensity share is wrong or missing: %r"
                   % re.findall(r"\d+% of the intensity", html))
    if "y7" not in html:
        bad.append("a library ion absent from the experimental spectrum is not named")
    if "dominated by the base peak" not in html:
        bad.append("a dot product NOT on sqrt intensities was not flagged, so a"
                   " number computed on a different basis reads as the panel's own")
    if "-2.2" in html:
        bad.append("reported a ppm range from peaks that are not in this spectrum")
    if "+0.8 to +1.9 ppm" not in html:
        bad.append("the matched-peak ppm range is missing or wrong: %r"
                   % re.findall(r"[-+][\d.]+ to [-+][\d.]+ ppm", html))
    # the colour convention practitioners read
    if '#1b6ca8' not in html or '#b3312c' not in html or '#9a9a92' not in html:
        bad.append("b / y / unassigned are not drawn in the conventional colours")
    return bad


def case_an_envelope_with_wrong_spacing_is_not_one_species(tmp):
    uid = "c1" + "0" * 14
    i = ion(uid)
    # 1.00335/2 = 0.5017 apart is one species at 2+. These are 1.0 apart.
    i["isotopes"] = [{"k": 0, "mz": 745.8821, "rel": 100.0, "observed": True},
                     {"k": 1, "mz": 746.8821, "rel": 79.6, "observed": True}]
    i["isotope_fit"] = {"expected_basis": "composition",
                        "residual": [{"k": 0, "delta": 0.0}, {"k": 1, "delta": -0.8}],
                        "rms": 0.79, "ratio_cv": [{"k": 1, "cv": 0.40}], "n_scans": 7}
    d = export(tmp, "iso1", [i])
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "0.5017" not in html:
        bad.append("did not state the spacing a 2+ envelope should have")
    if "not one species" not in html:
        bad.append("a 1.0 Th gap at charge 2 was not flagged: the envelope spans"
                   " two species and the panel called it one")
    if "varies 40%" not in html:
        bad.append("an isotope ratio drifting 40% across the peak was not flagged")
    return bad


def case_a_modelled_envelope_is_never_drawn_as_observed(tmp):
    uid = "d1" + "0" * 14
    i = ion(uid)
    i["isotopes"] = [{"k": 0, "mz": 745.8821, "rel": 100.0, "observed": True},
                     {"k": 1, "mz": 746.3838, "rel": 79.6, "observed": False}]
    d = export(tmp, "iso2", [i])
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "MODELLED" not in html:
        bad.append("a modelled (averagine) isotope is not marked as such")
    if "1 of 2 isotope(s) are MODELLED" not in html:
        bad.append("the modelled isotopes are not counted: %r"
                   % re.findall(r"\d+ of \d+ isotope\(s\) are MODELLED", html))
    if "must not be read as observation" not in html:
        bad.append("does not say a modelled envelope is not an observation")
    return bad


def case_a_single_point_is_not_a_peak(tmp):
    uid = "e1" + "0" * 14
    d = export(tmp, "xic1", [ion(uid)])
    write(d, "traces.jsonl", [{"ion_uid": uid, "sample_id": "CH215", "kind": "xic",
                               "extraction_window_ppm": 7.1,
                               "points": [[24.1, 19.5]], "apex_rt_min": 24.1,
                               "fwhm_min": None}])
    declare(d, traces=True)
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    chunk = html.split("Extracted ion chromatogram")[1].split("</div>")[0]
    if "<svg" in chunk:
        bad.append("drew an XIC from a single point")
    if "at least two consecutive scans" not in html:
        bad.append("did not say why one point is not a peak")
    return bad


def case_the_xic_states_its_extraction_window(tmp):
    uid = "f1" + "0" * 14
    d = export(tmp, "xic2", [ion(uid)])
    pts = [[23.9, 3.2], [24.0, 11.7], [24.1, 19.5], [24.2, 8.0]]
    write(d, "traces.jsonl", [
        {"ion_uid": uid, "sample_id": "CH215", "kind": "xic",
         "extraction_window_ppm": 7.1, "points": pts, "apex_rt_min": 24.1,
         "fwhm_min": 0.21},
        {"ion_uid": uid, "sample_id": "CH215", "kind": "xic",
         "extraction_window_ppm": 7.1, "points": pts, "apex_rt_min": 24.1,
         "fwhm_min": 0.21}])
    declare(d, traces=True)
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "7.1 ppm" not in html:
        bad.append("the extraction window is not on the panel")
    if "part of the peak's shape" not in html:
        bad.append("does not say the window shapes the peak")
    if html.count("FWHM 0.210 min") != 1:
        bad.append("repeated the same FWHM note per identical trace (%d times)"
                   % html.count("FWHM 0.210 min"))
    if "2 traces identical" not in html:
        bad.append("did not say the identical traces are repeats")
    if "repeats of one extraction, not a comparison" not in html:
        bad.append("does not say same-sample traces are not a cross-sample view")
    return bad


def case_a_cross_sample_panel_checks_its_own_table(tmp):
    """The manifest's sample list is what the export DESCRIBES. The table is what
    it HOLDS, and a cross-sample figure drawn from one sample is false."""
    uid = "g1" + "0" * 14
    d = export(tmp, "prot1", [ion(uid)], n_samples=5)
    write(d, "proteins.jsonl", [
        {"sample_id": "CH215", "prot_accession": "Q99836",
         "protein_name": "MYD88_HUMAN", "top3avg": 7.209, "n_peptides_used": 3,
         "n_peptides_matched": 6, "peptides_used": ["FITVCDYTNPCTK"], "rank": 214,
         "percentile": 94.0, "normalized": {"median_scaled": 1.42},
         "tier": "pass2_confirmed", "decoy_controlled": True}])
    declare(d, proteins=True)
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    chunk = html.split("Top3avg percentile")[1].split("</div>")[0]
    if "<svg" in chunk:
        bad.append("drew a cross-sample panel from one sample's rows")
    if "covers 1 sample(s)" not in html:
        bad.append("did not say how many samples the table actually covers")
    if "manifest declares 5 samples" not in html:
        bad.append("did not contrast the table's coverage with the manifest's claim")
    if "is false" not in html:
        bad.append("did not say the figure would be false")
    return bad


def case_flux_needs_a_normalization_that_was_applied(tmp):
    uid = "h1" + "0" * 14
    rows = [{"sample_id": "CH%d" % (215 + n), "prot_accession": "Q99836",
             "protein_name": "M", "top3avg": 7.2 + n, "n_peptides_used": 3,
             "n_peptides_matched": 6, "peptides_used": ["X"], "rank": 2,
             "percentile": 94.0, "normalized": {}, "tier": "pass2_confirmed",
             "decoy_controlled": True} for n in range(2)]
    d = export(tmp, "flux1", [ion(uid)], n_samples=2)
    write(d, "proteins.jsonl", rows)
    declare(d, proteins=True)
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "declares no normalization" not in html:
        bad.append("drew or blanked flux without naming the missing normalization")
    if "mixes biology with loading" not in html:
        bad.append("did not say why raw Top3avg across runs will not do")
    # now declare one that no row carries
    man = json.loads((d / "manifest.json").read_text())
    man["quant"] = {"metric": "top3avg", "absolute_standard_declared": False,
                    "normalizations": ["median_scaled"]}
    (d / "manifest.json").write_text(json.dumps(man))
    rc2, out2, html2 = render(d, "--method", ROOT / "examples" / "method.example.json")
    if "was not actually applied" not in html2:
        bad.append("a declared normalization absent from every row was not caught")
    # and now actually apply it
    for n, r in enumerate(rows):
        r["normalized"] = {"median_scaled": 1.0 + 0.4 * n}
    write(d, "proteins.jsonl", rows)
    rc3, out3, html3 = render(d, "--method", ROOT / "examples" / "method.example.json")
    chunk = html3.split("Flux across runs")[1].split("</div>")[0]
    if "<svg" not in chunk:
        bad.append("did not draw flux once a declared normalization was present")
    if "normalization median_scaled, as declared" not in html3:
        bad.append("drew flux without naming the normalization it used")
    return bad


def case_stoichiometry_stays_inside_one_run(tmp):
    uid = "i1" + "0" * 14
    d = export(tmp, "stoich", [ion(uid)], n_samples=2)
    rows = [{"sample_id": "CH215", "prot_accession": "Q99836", "protein_name": "A",
             "top3avg": 8.0, "n_peptides_used": 3, "n_peptides_matched": 6,
             "peptides_used": ["X"], "rank": 1, "percentile": 99.0,
             "normalized": {}, "tier": "pass2_confirmed", "decoy_controlled": True},
            {"sample_id": "CH215", "prot_accession": "P04406", "protein_name": "B",
             "top3avg": 2.0, "n_peptides_used": 1, "n_peptides_matched": 1,
             "peptides_used": ["Y"], "rank": 9, "percentile": 70.0,
             "normalized": {}, "tier": "pass2_confirmed", "decoy_controlled": True},
            {"sample_id": "CH216", "prot_accession": "Q99836", "protein_name": "A",
             "top3avg": 9.0, "n_peptides_used": 3, "n_peptides_matched": 6,
             "peptides_used": ["X"], "rank": 1, "percentile": 99.0,
             "normalized": {}, "tier": "pass2_confirmed", "decoy_controlled": True}]
    write(d, "proteins.jsonl", rows)
    declare(d, proteins=True)
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    chunk = html.split("Protein stoichiometry")[1].split("</div>")[0]
    if "<svg" not in chunk:
        bad.append("did not draw a ratio from a run holding two proteins")
    if "1 ratio(s) from 1 run(s)" not in html:
        bad.append("ratio/run count wrong: CH216 holds one protein and must not"
                   " pair across runs -- got %r"
                   % re.findall(r"\d+ ratio\(s\) from \d+ run\(s\)", html))
    if "rest on a single peptide" not in html:
        bad.append("a ratio with one peptide on one side was not flagged")
    if "not the same measurement" not in html:
        bad.append("does not say a 1-peptide vs 3-peptide ratio differs in kind")
    return bad


def case_fragment_hits_show_their_controls(tmp):
    uid = "j1" + "0" * 14
    d = export(tmp, "frag", [ion(uid)])
    write(d, "fragment_hits.jsonl", [
        {"hit_id": "fh-1", "sample_id": "CH215",
         "peptide": {"sequence": "FITVCDYTNPCTK", "prot_accession": "Q99836",
                     "seq_start": 171, "boundary": "tryptic"},
         "queried": {"n_fragments": 8, "fragment_labels": ["b3"], "charges": [1],
                     "ppm_tolerance": 20, "selection": "both_termini"},
         "matched": [], "missed": ["b7"],
         "k_of_n": {"matched": 4, "queried": 8, "required": 3},
         "coelution": {"rt_min": 24.1, "rt_span_min": 0.012,
                       "consecutive_scans": 4, "drift": 0.941,
                       "drift_offset_bins_applied": -3.25},
         "precursor": {"found": True}, "evidence": "fragment_and_precursor",
         "control": {"reversed_matched": None, "entrapment_matched": None}}])
    declare(d, fragment_hits=True)
    method = {"method_id": "m", "name": "frag", "version": "1.0.0",
              "tier_policy": {}, "panels": [{"id": "fragment_hits"}]}
    mp = Path(tmp) / "m_frag.json"
    mp.write_text(json.dumps(method))
    rc, out, html = render(d, "--method", mp)
    bad = []
    if "4/8 (need 3)" not in html:
        bad.append("k-of-n is not shown as matched/queried with the requirement")
    if "not a false-discovery rate" not in html:
        bad.append("an uncontrolled k-of-n was not distinguished from an FDR")
    if "-3.25" not in html or "partly constructed" not in html:
        bad.append("a drift offset applied to the coelution test was not disclosed")
    return bad


def case_a_declared_table_that_is_missing_is_reported(tmp):
    d = export(tmp, "liar", [ion("k1" + "0" * 14)])
    declare(d, spectra=True, spectrum_kinds=["ms2_experimental", "ms2_library"],
            traces=True, proteins=True)
    rc, out, html = render(d, "--method", ROOT / "examples" / "method.example.json")
    bad = []
    if "does not describe itself" not in html:
        bad.append("a manifest declaring tables the directory lacks was not called"
                   " out on the page")
    for name in ("spectra", "traces", "proteins"):
        if "%s.jsonl" % name not in html:
            bad.append("%s is declared and missing but not named" % name)
    chunk = html.split("MS/MS head-to-tail")[1].split("</div>")[0]
    if "<svg" in chunk:
        bad.append("drew a panel from a table that is not there")
    return bad


CASES = [case_sonar_drift_is_never_labelled_mobility,
         case_mixed_drift_kinds_refuse_a_shared_axis,
         case_a_blank_panel_is_drawn_not_omitted,
         case_points_a_log_axis_cannot_show_are_counted,
         case_an_all_bare_page_says_whether_the_rate_was_measured,
         case_mixed_tiers_are_marked_and_announced,
         case_quant_is_relative_unless_a_standard_was_declared,
         case_decimation_is_on_the_page,
         case_cluster_rank_is_required_not_assumed,
         case_an_empty_export_refuses_rather_than_drawing_a_blank_page,
         case_the_page_needs_no_network, case_untrusted_text_is_escaped,
         case_the_plan_is_the_single_source_of_truth,
         case_ion_flight_needs_a_path, case_no_method_offers_every_drawable_panel,
         case_a_one_sided_mirror_is_not_a_comparison,
         case_the_mirror_reports_what_it_did_not_explain,
         case_an_envelope_with_wrong_spacing_is_not_one_species,
         case_a_modelled_envelope_is_never_drawn_as_observed,
         case_a_single_point_is_not_a_peak,
         case_the_xic_states_its_extraction_window,
         case_a_cross_sample_panel_checks_its_own_table,
         case_flux_needs_a_normalization_that_was_applied,
         case_stoichiometry_stays_inside_one_run,
         case_fragment_hits_show_their_controls,
         case_a_declared_table_that_is_missing_is_reported]


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
