#!/usr/bin/env python3
"""Prove the PAX8 join reads the real file shapes, and refuses what it cannot measure.

Fixtures are written in the formats Scott-PC reported: ApexRT Summary temp1 with
its thirteen columns, the iDTs plainFP spectrum table, and a FragPipe psm.tsv
whose Intensity column is zero throughout, as the real ones are.

    python3 tools/test_pax8_load_response.py
"""

import json
import math
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "pax8_load_response.py"
LOADS = ROOT / "examples" / "pax8_recorded_loads.json"
sys.path.insert(0, str(ROOT / "tools"))
import pax8_load_response as P  # noqa: E402

PEPS = [("KMDDSDQDSCRPGLK", 2), ("MDDSDQDSCR", 2), ("VLTEIADQFGVNK", 2),
        ("SGPYTNFLQEAK", 2), ("ATLTPSNTPLGR", 2), ("LSIDSQSSSSGPR", 2),
        ("QHYPEAAGPAR", 2), ("FITVCDYTNPCTK", 2), ("GALQNIIPASTGAAK", 3),
        ("LISWYDNEFGYSNR", 2), ("VPTANVSVVDLTCR", 2), ("AAVLDLTCR", 2)]
RES = dict(zip("ACDEFGHIKLMNPQRSTVWY",
               [71.03711, 103.00919 + 57.02146, 115.02694, 129.04259, 147.06841,
                57.02146, 137.05891, 113.08406, 128.09496, 113.08406, 131.04049,
                114.04293, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
                99.06841, 186.07931, 163.06333]))
PROTON, WATER = 1.007276, 18.010565


def mz_of(seq, z):
    return (sum(RES[c] for c in seq) + WATER + z * PROTON) / z


def record():
    doc = json.loads(LOADS.read_text())
    return doc, {r["run"]: r for r in doc["runs"]}


def plant(root, no_psm=(), alt_table=(), load_scale=None, vol_slope=0.0,
          saturate_above=None, detect_floor=None, seed=4,
          saturate_grad=None, max_peps=None, decoy_group=False, psm_repeats=1):
    """Write the three real file shapes for every recorded run.

    load_scale   intensity multiplier per ng, so the planted load response is known
    vol_slope    log10 intensity change per decade of injection volume at fixed load
    saturate_above  set IsSaturated on detections above this intensity
    saturate_grad   {gradient: threshold}, so one LC method can flag and the
                    other not -- which is how a cross-gradient comparison gets
                    caught
    max_peps        use only the first N peptides, to starve a series of ions
    decoy_group     add a BRIGHTER feature further away inside every window, so
                    a join that takes the brightest instead of the nearest
                    picks the wrong one
    psm_repeats     write N PSM rows per peptide, as a real search does when an
                    ion is sequenced more than once
    """
    rng = random.Random(seed)
    doc, rec = record()
    apex = Path(root) / "apex"
    search = Path(root) / "search"
    amps = {(p, z): 10 ** rng.uniform(3.5, 5.5) for p, z in PEPS}
    for run, r in rec.items():
        a = apex / run
        (a / "Summary").mkdir(parents=True, exist_ok=True)
        (a / "idts_out").mkdir(parents=True, exist_ok=True)
        rows, table = [], []
        peps_here = PEPS[:max_peps] if max_peps else PEPS
        for n, (pep, z) in enumerate(peps_here, start=1):
            mz = mz_of(pep, z)
            rt = 5.0 + 3.1 * n
            drift = 40.0 + 2.5 * n
            scale = (load_scale or (lambda ng: ng))(r["load_ng"])
            trend = 10 ** (vol_slope * math.log10(r["recorded_uL"]))
            inten = amps[(pep, z)] * scale * trend * math.exp(rng.gauss(0, 0.05))
            if detect_floor and inten < detect_floor:
                # below the floor the feature is simply absent from the ApexRT
                # table, which is how a low-load run loses its weakest ions
                continue
            sat = 1 if (saturate_above and inten > saturate_above) else 0
            if saturate_grad:
                thr = saturate_grad.get(r["gradient"])
                if thr is not None and inten > thr:
                    sat = 1
            # three detections in one peak group, as the real tables carry
            for k in range(3):
                rows.append([1, 100 + n * 10 + k, round(rt + 0.01 * k, 6),
                             round(mz + 0.0001 * k, 6), 0.0329, round(drift, 6),
                             3.19, 0.0, 119, int(inten / 3), sat, n, 3])
            if decoy_group:
                # brighter, and further from the PSM inside the same window:
                # a join that ranks by intensity takes this one
                # FLAT in load, so a join that takes it reports a ratio near 1
                # instead of the planted 5. A decoy that merely scaled with load
                # would cancel in the ratio and catch nothing.
                flat = amps[(pep, z)] * 400.0
                for k in range(3):
                    rows.append([1, 100 + n * 10 + k, round(rt + 0.09 + 0.01 * k, 6),
                                 round(mz + 0.0000045 * mz, 6), 0.0329,
                                 round(drift + 1.6, 6), 3.19, 0.0, 119,
                                 int(flat / 3), sat, 500 + n, 3])
            table.append([n, "CPPIS: %d PseudoID: Z%dMZ" % (n, n),
                          round(rt, 3), round(drift, 3), round(mz * z - z * PROTON, 4),
                          "%d+" % z, 42])
        with (a / "Summary" / ("%s_ApexRT_temp1.csv" % run)).open("w") as fh:
            fh.write(",".join(P.APEX_COLS) + "\n")
            for row in rows:
                fh.write(",".join(str(v) for v in row) + "\n")
        name = ("%s.plain_table.csv" if run in alt_table
                else "%s_plainFP_table.csv") % run
        with (a / "idts_out" / name).open("w") as fh:
            fh.write("idx,title,rt_min,drift,pepmass,charge,n_peaks\n")
            for row in table:
                fh.write(",".join(str(v) for v in row) + "\n")
        if run in no_psm:
            (search / ("%s_apexrt_idts" % run) / ("%s_apexrt_idts" % run)).mkdir(
                parents=True, exist_ok=True)
            continue
        d = search / ("%s_apexrt_idts" % run) / ("%s_apexrt_idts" % run)
        d.mkdir(parents=True, exist_ok=True)
        cols = ["Spectrum", "Spectrum File", "Peptide", "Charge", "Retention",
                "Calculated M/Z", "Expectation", "Hyperscore", "Qvalue",
                "Intensity", "Is Decoy", "Protein"]
        with (d / "psm.tsv").open("w") as fh:
            fh.write("\t".join(cols) + "\n")
            peps_here = PEPS[:max_peps] if max_peps else PEPS
            for n, (pep, z) in enumerate(peps_here, start=1):
              rt = 5.0 + 3.1 * n
              for rep_i in range(psm_repeats):
                fh.write("\t".join(str(v) for v in [
                    "%s_plainFP.%05d.%05d.%d" % (run, n, n, z), run + ".mgf",
                    pep, z, round(rt * 60.0, 2),      # SECONDS, as the real file
                    round(mz_of(pep, z), 4), 3.0e-6, 28.1, 0.001,
                    0,                                 # Intensity: zero, as it is
                    "false", "sp|TEST|TEST"]) + "\n")
            # a decoy and a poor-q row that must both be dropped
            fh.write("\t".join(str(v) for v in [
                "%s_plainFP.00099.00099.2" % run, run + ".mgf", "DECOYPEPTIDEK", 2,
                600.0, 700.1234, 0.5, 9.0, 0.001, 0, "true", "rev_sp|X|X"]) + "\n")
            fh.write("\t".join(str(v) for v in [
                "%s_plainFP.00098.00098.2" % run, run + ".mgf", "BADQPEPTIDEK", 2,
                600.0, 700.1234, 0.5, 9.0, 0.44, 0, "false", "sp|Y|Y"]) + "\n")
    return str(apex), str(search)


def run_tool(*args):
    p = subprocess.run([sys.executable, str(TOOL)] + [str(a) for a in args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def case_the_join_recovers_a_planted_load_ratio(tmp):
    """Intensity proportional to load: the 500/100 ratio must come out near 5."""
    d = tmp / "ratio"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "LOAD RESPONSE" not in out:
        return ["no load response section"]
    # Anchor on the load-response section and on the per-gradient marker, not on
    # a bare "60 min" -- that now appears in the saturation table above it.
    section = out.split("LOAD RESPONSE", 1)[1]
    for grad in (60, 120):
        marker = "%d min: 100 ng" % grad
        if marker not in section:
            bad.append("no %d min load-response block" % grad)
            continue
        chunk = section.split(marker, 1)[1].split(" min: 100 ng")[0]
        m = [l for l in chunk.splitlines() if "median I(" in l]
        if not m:
            bad.append("no ratio reported for %d min" % grad)
            continue
        val = float(m[0].split("=")[1].strip())
        if abs(val - 5.0) > 0.4:
            bad.append("%d min ratio %.3f, planted 5.0" % (grad, val))
    if "SINGLE run" not in out:
        bad.append("did not flag the unreplicated 100 ng level")
    if "NO curvature" not in out:
        bad.append("did not refuse a curvature verdict on two levels")
    return bad


def case_gradients_are_reported_separately(tmp):
    d = tmp / "grad"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if out.count("median I(") < 2:
        bad.append("the two gradients were not reported separately")
    if "60 min: 100 ng x1" not in out and "60 min: " not in out:
        bad.append("no per-gradient load breakdown")
    return bad


def case_fragpipe_intensity_is_never_used(tmp):
    """It is 0 in every row of every searched run."""
    d = tmp / "zero"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "every psm.tsv Intensity in this run is 0" not in out:
        bad.append("did not notice that the FragPipe intensities are all zero")
    if "ApexRT's" not in out:
        bad.append("did not say where the intensity it used came from")
    # and a bare --psm must refuse outright
    rc2, out2 = run_tool("--psm", str(Path(search)))
    if rc2 != 2:
        bad.append("exit %d for --psm alone; it cannot measure anything" % rc2)
    if "REFUSED" not in out2 or "zeros as" not in out2:
        bad.append("did not refuse --psm with the reason")
    return bad


def case_runs_without_a_psm_file_are_named(tmp):
    """Three real runs have a search folder and no psm.tsv."""
    missing = ("PAX8_JP210_500ng_60_1uL_20260924_1",
               "PAX8_JP212_500ng_60_2uL_20260924_3",
               "PAX8_JP215_500ng_60_5uL_20260924_3")
    d = tmp / "nopsm"
    apex, search = plant(d, no_psm=missing)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    for run_name in missing:
        if run_name not in out:
            bad.append("%s was dropped without being named" % run_name)
    if "SKIPPED: no psm.tsv" not in out:
        bad.append("did not give the reason those runs contributed nothing")
    if "can lose its only run this way" not in out:
        bad.append("did not warn that a load level can lose its only run")
    if rc == 0:
        bad.append("exit 0 while three runs contributed nothing")
    return bad


def case_the_alternate_table_name_is_found(tmp):
    """Three runs name it <run>.plain_table.csv instead."""
    alt = ("PAX8_JP210_500ng_60_1uL_20260924_1",
           "PAX8_JP210_500ng_60_1uL_20260924_2",
           "PAX8_JP210_500ng_60_1uL_20260924_3")
    d = tmp / "alt"
    apex, search = plant(d, alt_table=alt)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "plain_table.csv, not" not in out:
        bad.append("the alternate spectrum-table name was not reported")
    for run_name in alt:
        chunk = [l for l in out.splitlines() if run_name in l and "ion(s) of" in l]
        if not chunk:
            bad.append("%s produced no ions despite its table being present"
                       % run_name)
    return bad


def case_saturation_comes_from_the_flag(tmp):
    """ApexRT sets IsSaturated, so saturation is counted rather than inferred."""
    d = tmp / "sat"
    # saturate the bright ions at 500 ng only
    apex, search = plant(d, saturate_above=2e7)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "SATURATION" not in out:
        bad.append("no saturation section")
    if "IsSaturated" not in out:
        bad.append("did not name the flag it used")
    if "evidence, not ground truth" not in out:
        bad.append("did not say the flag is the instrument's, not ground truth")
    if "can clip without setting it" not in out:
        bad.append("did not say how the flag can be wrong in either direction")
    if "FLOOR on the true ratio" not in out:
        bad.append("did not say what saturation does to the load ratio")
    if "per run, over ALL function-1 features" not in out:
        bad.append("did not report the flagged share per run over all features")
    if "intensity bracket the flag starts in" not in out:
        bad.append("did not report where the flag turns on")
    if "lowest flagged intensity" not in out:
        bad.append("did not give the onset bracket's lower edge")
    if "stay OUT of any within-run normalization set" not in out:
        bad.append("did not say flagged ions are excluded from normalization")
    # and with nothing saturated it must say so rather than stay silent
    d2 = tmp / "nosat"
    apex2, search2 = plant(d2)
    rc2, out2 = run_tool("--apex", apex2, "--search", search2)
    if "does not rise with load" not in out2:
        bad.append("did not report the absence of load-dependent saturation")
    return bad


def case_a_constant_load_trend_is_found_and_not_called_load(tmp):
    d = tmp / "voltrend"
    apex, search = plant(d, vol_slope=0.30)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "CONSTANT-LOAD SERIES" not in out:
        bad.append("no constant-load section")
    m = [l for l in out.splitlines() if "median" in l and "quartiles" in l
         and "+0.3" in l or "+0.2" in l]
    if "per-ion slope of log I against log volume" not in out:
        bad.append("did not report the per-ion slope against volume")
    if "expected 0.000" not in out:
        bad.append("did not state the expectation at constant load")
    if "not the load" not in out:
        bad.append("did not rule out reading the volume trend as load")
    if "vial-to-vial spread" not in out:
        bad.append("did not report the replicate noise the trend must clear")
    return bad


def case_a_flat_series_clears_no_verdict(tmp):
    d = tmp / "volflat"
    apex, search = plant(d, vol_slope=0.0)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "CONSTANT-LOAD SERIES" not in out:
        return ["no constant-load section"]
    chunk = out.split("CONSTANT-LOAD SERIES", 1)[1]
    if "clears the replicate noise" in chunk and "does not clear" not in chunk:
        bad.append("a flat series was reported as a real trend")
    return bad


def case_the_disputed_run_is_named_before_any_number(tmp):
    d = tmp / "disp"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "DISPUTED" not in out:
        bad.append("the disputed run was not named")
    if "inherits the dispute" not in out:
        bad.append("did not say what else the dispute affects")
    if out.index("DISPUTED") > out.index("LOAD RESPONSE"):
        bad.append("the dispute is reported after the numbers that use it")
    return bad


def case_the_name_record_check_never_infers_load(tmp):
    d = tmp / "names"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "NAME VERSUS RECORD" not in out:
        bad.append("no name-versus-record check")
    if "3 of 24 run name(s) disagree" not in out:
        bad.append("wrong count of disagreeing names: %r"
                   % [l for l in out.splitlines() if "disagree" in l])
    if "a name can never" not in out:
        bad.append("did not say a name cannot establish load")
    return bad


def case_the_published_ratio_is_rechecked(tmp):
    d = tmp / "recheck"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "RECHECK" not in out:
        bad.append("no recheck of the published ratio")
    if "with the recorded loads" not in out:
        bad.append("did not report the correctly grouped figure")
    if "with the named-100ng runs grouped" not in out:
        bad.append("did not report the mis-grouped figure for comparison")
    if "Consistent with, not proof" not in out and "does not explain it" not in out:
        bad.append("drew a conclusion stronger than the comparison supports")
    return bad


def case_units_are_stated(tmp):
    d = tmp / "units"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search, "--schema-only")
    bad = []
    if "BINS" not in out:
        bad.append("did not state that drift is in TWIMS bins")
    if "1/K0" not in out:
        bad.append("did not rule out 1/K0")
    if "min" not in out:
        bad.append("did not state the retention unit")
    return bad


def case_schema_only_measures_nothing(tmp):
    d = tmp / "schema"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search, "--schema-only")
    bad = []
    if rc != 0:
        bad.append("exit %d for a schema-only run" % rc)
    for banned in ("LOAD RESPONSE", "CONSTANT-LOAD SERIES", "RECHECK"):
        if banned in out:
            bad.append("--schema-only computed %r" % banned)
    if "ion(s) of" not in out:
        bad.append("did not report what it found per run")
    return bad


def case_the_load_record_is_internally_consistent(tmp):
    doc, rec = record()
    bad = []
    if len(rec) != 24:
        bad.append("%d runs in the record, the cohort has 24" % len(rec))
    for run, r in rec.items():
        want = round(r["recorded_uL"] * r["conc_ug_per_uL"] * 1000, 3)
        if abs(r["load_ng"] - want) > 1e-6:
            bad.append("%s: %g uL x %g should be %g ng, record says %g"
                       % (run, r["recorded_uL"], r["conc_ug_per_uL"], want,
                          r["load_ng"]))
        if r["gradient"] not in (60, 120):
            bad.append("%s has gradient %r" % (run, r["gradient"]))
    # the constant-load design: every 500 ng run at any volume is 500 ng
    at500 = {r["recorded_uL"] for r in rec.values() if r["load_ng"] == 500.0}
    if at500 != {1.0, 2.0, 5.0}:
        bad.append("500 ng volumes are %r, the design is 1, 2 and 5 uL" % at500)
    # exactly one 100 ng run per gradient
    for g in (60, 120):
        n = sum(1 for r in rec.values() if r["gradient"] == g
                and r["load_ng"] == 100.0)
        if n != 1:
            bad.append("%d min has %d runs at 100 ng, the record says 1" % (g, n))
    return bad




def case_load_and_volume_are_separated(tmp):
    """The design confounds them and also separates them.

    100 ng is 0.2 uL and 500 ng is 1.0 uL of the SAME 0.5 ug/uL vial, so the 5x
    in load is delivered entirely as 5x in volume. The constant-load series
    spans the same 5x at a fixed load, so its slope can be subtracted. Planted
    volume slopes of 0, 0.3 and 0.6 must all give back the planted 5x load.
    """
    bad = []
    seen = []
    for n, vs in enumerate((0.0, 0.30, 0.60)):
        d = tmp / ("sep%d" % n)
        apex, search = plant(d, vol_slope=vs)
        rc, out = run_tool("--apex", apex, "--search", search)
        if "THE LOAD CONTRAST IS A VOLUME CONTRAST" not in out:
            bad.append("slope %g: did not declare the confound" % vs)
        i = out.find("60 min: 100 ng")
        chunk = out[i:i + 1600]
        obs = [l for l in chunk.splitlines() if "median I(" in l]
        adj = [l for l in chunk.splitlines() if "LOAD ONLY" in l]
        if not obs or not adj:
            bad.append("slope %g: no observed or corrected ratio" % vs)
            continue
        o = float(obs[0].split("=")[1].strip())
        a = float(adj[0].split(":")[1].split(",")[0].strip())
        seen.append((vs, o, a))
        if abs(a - 5.0) > 0.5:
            bad.append("slope %g: load-only %.3f, planted 5.0 (observed %.3f)"
                       % (vs, a, o))
    # the observed ratio must actually move with the planted slope, or the
    # fixture is not exercising the correction at all
    if len(seen) == 3:
        if seen[2][1] - seen[0][1] < 2.0:
            bad.append("the observed ratio barely moved across planted slopes %r,"
                       " so the correction was not tested" % [s[1] for s in seen])
        spread = max(s[2] for s in seen) - min(s[2] for s in seen)
        if spread > 0.1:
            bad.append("the corrected ratio moved %.3f across planted volume"
                       " slopes; it should not move at all" % spread)
    return bad


def case_the_contrast_is_restricted_to_one_concentration(tmp):
    d = tmp / "matched"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "restricted to the 0.5 ug/uL vial" not in out:
        bad.append("did not restrict the high side to the matched concentration")
    if "only the VOLUME differs" not in out:
        bad.append("did not say why the restriction matters")
    if "mix the injection into the ratio" not in out:
        bad.append("did not say what the excluded runs would do to the ratio")
    return bad


def case_the_extrapolation_is_declared(tmp):
    d = tmp / "extrap"
    apex, search = plant(d, vol_slope=0.30)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "extrapolates a slope measured between" not in out:
        bad.append("did not say the correction extrapolates")
    if "where it was not measured" not in out:
        bad.append("did not say the low volume is outside the measured range")
    if "nothing here can tell" not in out:
        bad.append("did not admit the correction cannot be checked")
    return bad


def case_without_a_volume_slope_the_confound_stays_declared(tmp):
    """Strip the constant-load series: the ratio must say it is still confounded."""
    d = tmp / "noconst"
    apex, search = plant(d)
    # remove every 500 ng run that is not at the matched 1.0 uL, which leaves no
    # volume series to measure a slope from
    doc, rec = record()
    import shutil
    for run_name, r in rec.items():
        if r["load_ng"] == 500.0 and r["recorded_uL"] != 1.0:
            shutil.rmtree(Path(apex) / run_name, ignore_errors=True)
            shutil.rmtree(Path(search) / ("%s_apexrt_idts" % run_name),
                          ignore_errors=True)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "LOAD ONLY" in out:
        bad.append("subtracted a volume effect with no series to measure one from")
    if "stay confounded" not in out and "NOT subtracted" not in out:
        bad.append("did not say load and injection remain confounded")
    return bad




def case_the_sanitised_inner_folder_is_found(tmp):
    """FragPipe rewrites 0-5uL to 0_5uL in the INNER folder only.

    The outer folder and the psm Spectrum prefix keep the hyphen. A path built
    from the run stem misses it, and the disputed run would be reported as
    unsearched -- the one run whose load is already in question.
    """
    d = tmp / "sanitised"
    apex, search = plant(d)
    run = "PAX8_JP213_100ng_60_0-5uL_20260924_4"
    inner = run.replace("0-5uL", "0_5uL")
    outer = Path(search) / ("%s_apexrt_idts" % run)
    (outer / ("%s_apexrt_idts" % inner)).mkdir(parents=True, exist_ok=True)
    (outer / ("%s_apexrt_idts" % run) / "psm.tsv").rename(
        outer / ("%s_apexrt_idts" % inner) / "psm.tsv")
    (outer / ("%s_apexrt_idts" % run)).rmdir()
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "SKIPPED: no psm.tsv" in out and run in out.split("SKIPPED: no psm.tsv")[0][-200:]:
        bad.append("the sanitised run was reported as unsearched")
    hit = [l for l in out.splitlines() if run in l and "ion(s) of" in l]
    if not hit:
        bad.append("the sanitised run produced no ions")
    if "sanitised the name" not in out:
        bad.append("found the folder but did not say the name was rewritten")
    if "would have called this run" not in out:
        bad.append("did not say what building the path from the stem would do")
    return bad



def case_the_ratio_carries_a_bootstrap_interval(tmp):
    """~30 ions per run is thin, and a bare median hides that."""
    d = tmp / "boot"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "bootstrap 95% interval" not in out:
        bad.append("no interval on the median ratio")
    if "over IONS only" not in out:
        bad.append("did not say what the interval covers")
    if "unreplicated" not in out or "not in it" not in out:
        bad.append("did not say the single low-load run's error is outside it")
    line = [l for l in out.splitlines() if "bootstrap 95% interval" in l]
    if line:
        import re as _re
        nums = _re.findall(r"\d+\.\d+", line[0])
        if len(nums) >= 2 and float(nums[0]) >= float(nums[1]):
            bad.append("interval is inverted: %r" % line[0].strip())
    return bad




def case_a_yield_gap_between_levels_is_a_selection_bias(tmp):
    """Fewer matches at the low load is not just a smaller n.

    With a detection floor the 100 ng side loses its weakest ions, so the ones
    it keeps are its strongest. That inflates the low side and DEPRESSES the
    ratio, and the tool must say so and say which way the true value lies.
    """
    d = tmp / "bias"
    # the floor sits between the 100 ng and 500 ng intensity ranges, so it bites
    # only on the low-load runs
    apex, search = plant(d, detect_floor=2e6)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "JOIN YIELD" not in out:
        bad.append("no join-yield table")
    yl = [l for l in out.splitlines() if "join yield" in l and "ng" in l]
    if not yl:
        return bad + ["no per-level yield beside the ratio"]
    if "SELECTION BIAS" not in out:
        bad.append("a yield gap between levels was not called a selection bias: %r"
                   % [l.strip() for l in yl])
        return bad
    if "not merely a smaller n" not in out:
        bad.append("did not rule out the benign reading")
    if "BRIGHTEST ions' ratio" not in out:
        bad.append("did not say what the ratio is now measured over")
    if "DIRECTION is not knowable" not in out:
        bad.append("claimed a direction the gap alone does not give")
    if "paired" not in out:
        bad.append("did not say the ratio is paired, which is why it is re-aimed"
                   " rather than depressed")
    for banned in ("DEPRESSES the ratio", "true ratio is higher",
                   "true ratio is lower"):
        if banned in out:
            bad.append("asserts a direction this ratio does not have: %r" % banned)
    # The measurement behind that: a 58-point yield gap must NOT move the paired
    # median much. If it does, the claim in the output is wrong.
    section = out.split("LOAD RESPONSE", 1)[1]
    m = [l for l in section.splitlines() if "median I(" in l]
    n = [l for l in section.splitlines() if "ion(s) at both levels" in l]
    if m:
        val = float(m[0].split("=")[1].strip())
        if abs(val - 5.0) > 0.4:
            bad.append("paired ratio %.3f under a 58-point yield gap; the output"
                       " says such a gap re-aims the ratio rather than shifting"
                       " it" % val)
    if n:
        count = int(n[0].strip().split()[0])
        if count >= 12:
            bad.append("%d ion(s) survived, so the fixture's floor did not create"
                       " a yield gap at all" % count)
    return bad


def case_a_low_overall_yield_points_at_the_units(tmp):
    """The join's two unit traps are where a collapsed yield would come from."""
    d = tmp / "tight"
    apex, search = plant(d)
    # a drift window far below the planted spacing kills most matches
    rc, out = run_tool("--apex", apex, "--search", search, "--drift-tol", "0.0001",
                       "--rt-tol", "0.0001")
    bad = []
    if "BELOW HALF" not in out:
        bad.append("a collapsed yield was not flagged")
    if "SECONDS" not in out or "MINUTES" not in out:
        bad.append("did not name the retention unit trap")
    if "TWIMS BINS" not in out or "1/K0" not in out:
        bad.append("did not name the drift unit trap")
    if "does not improve gradually" not in out:
        bad.append("did not say how to tell a unit error from a tight window")
    return bad


def case_the_cross_check_band_is_reported_not_enforced(tmp):
    d = tmp / "band"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "cross-check" not in out:
        bad.append("the independent precursor count is not used as a cross-check")
    if "26-42" not in out:
        bad.append("did not quote the cross-check band")
    # the fixture plants 12 peptides, so every run falls outside the band and
    # must still be used -- a cross-check is not a filter
    if "LOAD RESPONSE" not in out or "median I(" not in out:
        bad.append("runs outside the cross-check band were excluded; it is a"
                   " cross-check, not a gate")
    return bad


def case_the_owner_facing_ranking_is_printed(tmp):
    d = tmp / "verdict"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "WHAT TO TRUST IN THIS COHORT" not in out:
        bad.append("no owner-facing ranking of the outputs")
    tail = out.split("WHAT TO TRUST IN THIS COHORT", 1)[1]
    if "MOST: the constant-load series" not in tail:
        bad.append("did not name the constant-load series as the most trustworthy")
    if "LEAST: the load ratio" not in tail:
        bad.append("did not name the load ratio as the least trustworthy")
    if "opposite of how the work was framed" not in tail:
        bad.append("did not flag that this inverts the brief")
    if "partial-loop" not in tail or "0.2 uL" not in tail:
        bad.append("did not give the reason the load ratio is weakest")
    if "detector linearity is not" not in tail:
        bad.append("did not say which conclusion the cohort cannot support")
    return bad




def case_a_gap_plus_a_rising_flag_is_the_case_to_watch(tmp):
    """The conjunction, not either alone.

    A yield gap re-aims the ratio at the brighter ions. That only changes the
    NUMBER where response varies with brightness, and a flagged share rising
    with load is the evidence that it does. The two together is the case where
    the ratio should not be read; either alone is not.
    """
    bad = []
    # the saturation threshold must sit between the two load levels' intensity
    # ranges, or the surviving bright low-load ions are flagged too and the
    # pattern comes out backwards
    d = tmp / "conj"
    apex, search = plant(d, detect_floor=2e6, saturate_above=4e7)
    rc, out = run_tool("--apex", apex, "--search", search)
    if "SELECTION BIAS" not in out:
        return ["the fixture produced no yield gap"]
    if "AND IT DOES HERE" not in out:
        bad.append("a yield gap with a rising flagged share was not joined")
    if "where saturation bends" not in out:
        bad.append("did not connect the re-aimed ion set to saturation")
    if "Do not read the ratio" not in out:
        bad.append("did not say what to do about the conjunction")

    # a gap with NO rising flag must get the other branch, not this one
    d2 = tmp / "gaponly"
    apex2, search2 = plant(d2, detect_floor=2e6)
    rc2, out2 = run_tool("--apex", apex2, "--search", search2)
    if "SELECTION BIAS" not in out2:
        bad.append("the gap-only fixture produced no yield gap")
    if "AND IT DOES HERE" in out2:
        bad.append("claimed the conjunction on a fixture with no rising flag")
    if "does NOT rise with load" not in out2:
        bad.append("did not report the absence of a rising flagged share")
    if "without visibly" not in out2:
        bad.append("did not say the re-aiming changed the ion set but not"
                   " visibly what it measures")
    if "Do not read the ratio" in out2:
        bad.append("told the reader to discard a ratio on a gap alone")
    return bad


def case_a_volume_slope_corrects_only_its_own_gradient(tmp):
    """The 60 min constant-load slope must not reach the 120 min ratio.

    slope_out was a list and the caller took slope_out[0] for every gradient,
    so a slope measured over the 60 min volume series was subtracted from the
    120 min ratio -- where every run is 1 uL and there is no volume series to
    measure anything from.
    """
    d = tmp / "gradslope"
    apex, search = plant(d, vol_slope=0.25)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    n = out.count("constant-load series gives")
    if n != 1:
        bad.append("a volume correction was applied %d time(s); only the 60 min"
                   " gradient has a volume series, so it must be applied once" % n)
    if "constant-load volume series of its" not in out:
        bad.append("the gradient without a volume series is not told that its"
                   " injection was never subtracted")
    sec = out.split("LOAD RESPONSE", 1)[1]
    after = sec.split("120 min: 100 ng", 1)
    if len(after) > 1 and "constant-load series gives" in after[1]:
        bad.append("the 120 min block quotes a volume correction it has no"
                   " series for")
    return bad


def case_the_flagged_share_is_compared_within_a_gradient(tmp):
    """sorted() on (gradient, ng) made levels[0]=(60,100), levels[-1]=(120,500).

    Flagging nothing at either 60 min level and half at both 120 min levels is
    a pure METHOD difference with no load effect at all. Compared across
    gradients it reads as "the higher load flags more" and counts a finding.
    """
    d = tmp / "satgrad"
    apex, search = plant(d, saturate_grad={120: 0.0})
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "120 min: the flagged share does not rise with load" not in out:
        bad.append("the 120 min gradient is not compared against itself")
    if "60 min: the flagged share does not rise with load" not in out:
        bad.append("the 60 min gradient is not compared against itself")
    for line in out.splitlines():
        if "the higher load flags more" in line:
            bad.append("a rise in flagged share is reported where every rise is"
                       " between gradients: %r" % line.strip())
    return bad


def case_too_few_volume_ions_do_not_crash(tmp):
    """q = [None]*3 and then "%+.3f" % None took the whole run down."""
    d = tmp / "fewions"
    apex, search = plant(d, max_peps=2, vol_slope=0.2)
    rc, out = run_tool("--apex", apex, "--search", search)
    bad = []
    if "Traceback" in out or "TypeError" in out:
        bad.append("crashed on a volume series with too few ions: %s"
                   % [l for l in out.splitlines() if "Error" in l][:2])
    if rc == 2 and "Traceback" in out:
        bad.append("exit 2 from an unhandled exception")
    if "too few" not in out and "no quartiles" not in out:
        bad.append("did not say why the quartiles are missing")
    return bad


def case_the_yield_covers_the_runs_the_ratio_uses(tmp):
    """level_yield pooled every run at a load level; the ratio uses six of them.

    The verdict on selection bias was computed over runs the ratio never
    touches, so it answered about the wrong set.
    """
    d = tmp / "yieldset"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    sec = out.split("LOAD RESPONSE", 1)[1]
    chunk = sec.split("60 min: 100 ng", 1)[1].split("120 min: 100 ng")[0]
    bad = []
    restrict = [l for l in chunk.splitlines() if "restricted to the" in l]
    yline = [l for l in chunk.splitlines() if "join yield" in l]
    if not restrict or not yline:
        return ["could not find the restriction and yield lines"]
    import re as _re
    kept = int(_re.search(r"vial\(s\), (\d+) of", restrict[0]).group(1))
    both = _re.findall(r"\((\d+) of (\d+)\)", yline[0])
    if len(both) != 2:
        return ["could not read both yield denominators from %r" % yline[0].strip()]
    lo_total, hi_total = int(both[0][1]), int(both[1][1])
    if "SINGLE run" not in chunk:
        return ["this fixture is meant to have one run on the low side"]
    # One run on the low side, so its denominator IS the per-run precursor
    # count. The restricted high side must be exactly that many times the
    # number of runs the ratio kept -- no more.
    if hi_total != kept * lo_total:
        bad.append("the 500 ng yield denominator is %d; the ratio uses %d run(s)"
                   " of %d precursors, so it should be %d. The yield is being"
                   " judged over runs the ratio never touches."
                   % (hi_total, kept, lo_total, kept * lo_total))
    return bad


def case_the_disputed_run_is_not_on_both_sides(tmp):
    """The recheck put the suspect run in the high group AND moved it to the
    low group, so it sat on both sides of the same ratio."""
    d = tmp / "recheck2"
    apex, search = plant(d)
    rc, out = run_tool("--apex", apex, "--search", search)
    if "RECHECK OF THE PUBLISHED" not in out:
        return ["no recheck section"]
    sec = out.split("RECHECK OF THE PUBLISHED", 1)[1]
    bad = []
    import re as _re
    m = _re.search(r"low group is (\d+) run\(s\) and the high group (\d+)", sec)
    if not m:
        return ["could not read the group sizes"]
    lo_n, hi_n = int(m.group(1)), int(m.group(2))
    sus = _re.search(r"(\d+) run\(s\) named 100ng at 60 min record a different",
                     sec)
    n_sus = int(sus.group(1)) if sus else 0
    if "restricted to the" not in sec:
        bad.append("the high side is not restricted to the low side's vial"
                   " concentration, so the contrast mixes volume in")
    if hi_n > 6:
        bad.append("the high group holds %d run(s); the load ratio uses 6, and a"
                   " wider group puts other volumes into the comparison" % hi_n)
    # the suspect must not be inside the high group it is also moved out of
    if n_sus and hi_n + n_sus <= 6:
        bad.append("the suspect still appears to be counted in the high group")
    return bad


def case_the_join_takes_the_nearest_not_the_brightest(tmp):
    """A brighter feature further away in the same window must not win.

    The join kept whichever candidate had the larger intensity, which chooses
    by the quantity the whole study measures. Nothing stopped two ions claiming
    one feature either.
    """
    clean = tmp / "join_clean"
    decoy = tmp / "join_decoy"
    a1, s1 = plant(clean)
    a2, s2 = plant(decoy, decoy_group=True)
    _, out1 = run_tool("--apex", a1, "--search", s1)
    _, out2 = run_tool("--apex", a2, "--search", s2)

    def ratio(o, grad=60):
        sec = o.split("LOAD RESPONSE", 1)[1]
        chunk = sec.split("%d min: 100 ng" % grad, 1)[1]
        line = [l for l in chunk.splitlines() if "median I(" in l]
        return float(line[0].split("=")[1].strip()) if line else None

    bad = []
    r1, r2 = ratio(out1), ratio(out2)
    if r1 is None or r2 is None:
        return ["could not read the ratio from one of the runs"]
    # The decoy is 4x brighter at every load, so a brightest-wins join still
    # reports a ratio -- but it is measuring the decoy. The nearest-wins join
    # must give the same answer as the cohort without decoys.
    if abs(r1 - r2) > 0.25:
        bad.append("adding a brighter, further feature inside the window moved"
                   " the ratio from %.3f to %.3f, so the join is still taking"
                   " the brightest" % (r1, r2))
    return bad


def case_several_psms_of_one_ion_share_a_feature(tmp):
    """One feature per ION, not per PSM.

    The one-to-one rule contended between PSMs of the SAME peptide and charge,
    which `ions` collapses anyway. A search that sequences an ion four times
    then had three of its four PSMs rejected: 12 ion(s) of 48 PSM(s) with 36
    unmatched, where one PSM each matched cleanly.
    """
    import re as _re
    one = tmp / "psm1"
    many = tmp / "psm4"
    a1, s1 = plant(one, psm_repeats=1)
    a4, s4 = plant(many, psm_repeats=4)
    _, o1 = run_tool("--apex", a1, "--search", s1)
    _, o4 = run_tool("--apex", a4, "--search", s4)

    def counts(o):
        m = _re.search(r"(\d+) ion\(s\) of\s+(\d+) PSM\(s\); (\d+) unmatched", o)
        return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else None

    c1, c4 = counts(o1), counts(o4)
    bad = []
    if not c1 or not c4:
        return ["could not read the join counts"]
    if c4[2] != 0:
        bad.append("%d of %d PSM(s) went unmatched when the same ion was"
                   " sequenced four times; repeated PSMs of one ion must share"
                   " its feature" % (c4[2], c4[1]))
    if c4[0] != c1[0]:
        bad.append("repeating the PSMs changed the ion count from %d to %d"
                   % (c1[0], c4[0]))
    return bad


CASES = [case_the_join_recovers_a_planted_load_ratio,
         case_gradients_are_reported_separately,
         case_fragpipe_intensity_is_never_used,
         case_runs_without_a_psm_file_are_named,
         case_the_alternate_table_name_is_found,
         case_saturation_comes_from_the_flag,
         case_a_constant_load_trend_is_found_and_not_called_load,
         case_a_flat_series_clears_no_verdict,
         case_the_disputed_run_is_named_before_any_number,
         case_the_name_record_check_never_infers_load,
         case_the_published_ratio_is_rechecked,
         case_units_are_stated,
         case_schema_only_measures_nothing,
         case_the_load_record_is_internally_consistent,
         case_load_and_volume_are_separated,
         case_the_contrast_is_restricted_to_one_concentration,
         case_the_extrapolation_is_declared,
         case_without_a_volume_slope_the_confound_stays_declared,
         case_the_sanitised_inner_folder_is_found,
         case_the_ratio_carries_a_bootstrap_interval,
         case_a_yield_gap_between_levels_is_a_selection_bias,
         case_a_low_overall_yield_points_at_the_units,
         case_the_cross_check_band_is_reported_not_enforced,
         case_the_owner_facing_ranking_is_printed,
         case_a_gap_plus_a_rising_flag_is_the_case_to_watch,
         case_a_volume_slope_corrects_only_its_own_gradient,
         case_the_flagged_share_is_compared_within_a_gradient,
         case_too_few_volume_ions_do_not_crash,
         case_the_yield_covers_the_runs_the_ratio_uses,
         case_the_disputed_run_is_not_on_both_sides,
         case_the_join_takes_the_nearest_not_the_brightest,
         case_several_psms_of_one_ion_share_a_feature]


def main():
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
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
