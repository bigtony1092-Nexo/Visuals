#!/usr/bin/env python3
"""Prove the proximity control can report a bad window, and can report itself.

A control on a confidence tier is the easiest thing in this repository to write
so that it cannot fail: every way of breaking it -- an uncorrected calibration
offset, a charge no peak carries, a probe outside the acquired range, a thin
peak table -- drives the measured false-positive rate toward zero, and a rule
that passes a low rate rewards all of them. So the cases here plant each of
those breakages and assert the tool SAYS SO, and plant a window that does not
discriminate and assert the tool refuses it.

    python3 tools/test_proximity_control.py
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "proximity_control.py"
sys.path.insert(0, str(ROOT / "tools"))
import proximity_control as pc  # noqa: E402
from validate_export import peptide_mass  # noqa: E402

PROTON = 1.007276

# Real human tryptic peptides, used as the sample's targets.
HUMAN = ["LSIDSQSSSSGPR", "KMDDSDQDSCR", "SAMPLEDGTVK", "LLGATVLDALR",
         "AVFVDLEPTVIDEVR", "TIGGGDDSFNTFFSETGAGK", "YMVIQGEPGAVIR",
         "NLDIERPTYTNLNR", "SLDIEALFR", "VGINYQPPTVVPGGDLAK",
         "EIVHIQAGQCGNQIGAK", "DVNAAIATIK", "FDLMYAK", "QLFHPEQLITGK",
         "TLKFVTVDDGK", "ALGIDVTEVR", "WVDPNDPTK", "LMTEFDAHGR",
         "IVSNASCTTNCLAPLAK", "GYSFTTTAER", "ASGNYATVISHNPETK",
         "VQALEEANNDLENK", "LDGTEELLVK", "TPVITGAPYEYR", "SGPFGQIFR",
         "NVPLYKHIIR", "AQLVDMLGHR", "DSYVGDEAQSK", "ETLQQLQALR",
         "ILTERGYSFTTTAER", "EAQELLDHYGR", "TTGIVMDSGDGVTHTVPIYEGYALPHAIL",
         "MLSDFQVLR", "GILTLKYPIEHGIVTNWDDMEK", "VAPEEHPVLLTEAPLNPK"]

# E. coli tryptic peptides, used as entrapment: they cannot be in a human
# sample, and they are real sequences with real peptide masses rather than
# invented strings, because a mass off the peptide lattice would under-hit for
# a chemistry reason and the rate would read low for free.
ECOLI = ["AAVLTSHGK", "LTDAEVNSVLR", "GVDTSLIER", "VLNEPTAAALAYGLDK",
         "SIDLQTLTR", "AIFQGEVIR", "TLSDYNIQK", "ESTLHLVLR", "LGEHNIDVLEGNEQFINAAK",
         "ADIVLSHVK", "NLQDAMLVR", "VTDAEVNSLLR", "GIEPTLVER", "SLADELALVDVLEDK",
         "TIDLQSLTR", "AVFEGQVLR", "TLTDYNIQK", "ESMLHLVLR", "LGDHNIEVLEGNEQFLNAAK",
         "ADLVLSHIK", "NLQEAMLIR", "VSDAEVNTLLR", "GVEPSLIER", "SLSDELALVEVLEDK",
         "TVDLQTLSR", "AVFQGEVLR", "TLSEYNIQK", "ESTLHLVIR", "LGEHNIDVLEGNDQFINAAK",
         "ADIVLTHVK", "NLQDAMLIR", "VTDAEVNTLLR", "GIEPTLIER", "SLADELSLVDVLEDK",
         "TIDLQTLTR", "AVFEGEVIR", "TLTEYNIQK", "ESMLHLVIR", "LGDHNIDVLEGNEQFINAAK",
         "ADLVLSHVK", "NLQEAMLVR", "VSDAEVNSLLR"]


def mz(seq, z, ppm_off=0.0):
    m = peptide_mass(seq)
    return (m + z * PROTON) / z * (1.0 + ppm_off / 1e6)


def targets_file(tmp, name, seqs=HUMAN, z=2, rt_step=1.0, with_rt=True,
                 with_seq=True, ppm_off=0.0):
    """A target list in the shape make_targets.py emits: OBSERVED m/z."""
    lines = ["peptide,charge,mz,rt,drift"]
    for i, s in enumerate(seqs):
        lines.append("%s,%d,%.6f,%s,%.3f" % (s if with_seq else "", z,
                                             mz(s, z, ppm_off),
                                             ("%.3f" % (5.0 + i * rt_step)) if with_rt else "",
                                             50.0 + i))
    p = Path(tmp) / name
    p.write_text("\n".join(lines) + "\n")
    return p


def peaks_file(tmp, name, rows):
    import json
    p = Path(tmp) / name
    p.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    return p


def peaks_for(seqs, z=2, rt_step=1.0, ppm_off=0.0, n_scans=7, extra=()):
    """One peak per sequence, at the target's own coordinates."""
    rows = []
    for i, s in enumerate(seqs):
        rows.append({"mz": mz(s, z, ppm_off), "charge": z, "rt": 5.0 + i * rt_step,
                     "drift": 50.0 + i, "intensity": 1e6, "n_scans": n_scans})
    rows.extend(extra)
    return rows


def entrap_file(tmp, name, seqs=ECOLI):
    p = Path(tmp) / name
    p.write_text("\n".join(seqs) + "\n")
    return p


def run(*args):
    p = subprocess.run([sys.executable, str(TOOL)] + [str(a) for a in args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def control(tmp, peaks, targets, seqs=ECOLI, ppm_tol=10.0, rt_tol=None,
            drift_tol=None, min_scans=2):
    """Call the library directly, so a case can read the numbers."""
    from aria_run import read_peaks, load_targets
    rep, err = pc.run_control(read_peaks(str(peaks)), load_targets(str(targets)),
                              seqs, ppm_tol, rt_tol, drift_tol,
                              drift_tol is not None, min_scans)
    return rep, err


# ------------------------------------------------------------------ the cases

def case_a_planted_rate_is_recovered(tmp):
    """Plant a peak on a known share of entrapment coordinates; read it back."""
    bad = []
    tg = targets_file(tmp, "t1.csv")
    # Entrapment probes are built at the E. coli m/z on a human target's rt.
    # Plant a peak at the first 10 entrapment coordinates so the rate is known.
    rows = peaks_for(HUMAN)
    planted = 0
    for i, s in enumerate(ECOLI):
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
        planted += 1
        if planted >= 10:
            break
    pk = peaks_file(tmp, "p1.jsonl", rows)
    rep, err = control(tmp, pk, tg, rt_tol=0.2)
    if rep is None:
        return ["control did not run: %s" % err]
    e = rep["arms"]["entrapment"]
    want = 10
    if e["matched"] != want:
        bad.append("planted %d entrapment hits, tool reported %d of %d probes"
                   % (want, e["matched"], e["queried"]))
    if not (e["ci"][0] <= want / e["queried"] <= e["ci"][1]):
        bad.append("the 95%% interval %.3f-%.3f does not cover the planted rate %.3f"
                   % (e["ci"][0], e["ci"][1], want / e["queried"]))
    return bad


def case_a_window_that_does_not_discriminate_is_refused(tmp):
    """Entrapment matching as often as targets must not pass."""
    bad = []
    tg = targets_file(tmp, "t2.csv")
    rows = peaks_for(HUMAN)
    for i, s in enumerate(ECOLI):               # a peak on EVERY probe
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
    pk = peaks_file(tmp, "p2.jsonl", rows)
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    ok, why = pc.discriminated(rep)
    if ok:
        bad.append("every entrapment probe found a peak and the control still"
                   " reported that the window discriminated: %s" % why)
    rc, out = run(pk, tg, "--entrapment", entrap_file(tmp, "e2.txt"),
                  "--rt-tol", 0.2)
    if rc == 0:
        bad.append("exit 0 on a window that separated nothing")
    return bad


def case_the_same_window_is_used_as_the_match(tmp):
    """Widen --precursor-ppm and both rates must move.

    A control that reimplemented the window would not track the matcher, and the
    rate would describe code that never ran.
    """
    bad = []
    # Peaks sit 30 ppm off their target, so a 10 ppm window misses and 50 hits.
    tg = targets_file(tmp, "t3.csv")
    rows = peaks_for(HUMAN, ppm_off=30.0)
    for i, s in enumerate(ECOLI):
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2, 30.0), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
    pk = peaks_file(tmp, "p3.jsonl", rows)
    narrow, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=0.2)
    wide, _ = control(tmp, pk, tg, ppm_tol=50.0, rt_tol=0.2)
    for arm in ("target", "entrapment"):
        n, w = narrow["arms"][arm], wide["arms"][arm]
        if not (w["matched"] > n["matched"]):
            bad.append("%s arm did not move when the window widened from 10 to 50"
                       " ppm (%d then %d): the control is not measuring the"
                       " matcher's window" % (arm, n["matched"], w["matched"]))
    return bad


def case_axes_never_compared_are_not_claimed(tmp):
    """No --rt-tol, or no rt column, means the rate is global m/z crowding."""
    bad = []
    tg = targets_file(tmp, "t4.csv")
    pk = peaks_file(tmp, "p4.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=None)          # tolerance never passed
    if "rt" in rep["axes_matched"]:
        bad.append("rt listed as compared with no --rt-tol")
    if not any("GLOBAL" in f for f in rep["findings"]):
        bad.append("no --rt-tol and the report did not say the rate is a global"
                   " m/z-crowding rate")
    tg2 = targets_file(tmp, "t4b.csv", with_rt=False)   # tolerance passed, column null
    rep2, _ = control(tmp, pk, tg2, rt_tol=0.2)
    if "rt" in rep2["axes_matched"]:
        bad.append("rt listed as compared although no target carried one")
    return bad


def case_an_unmeasured_calibration_offset_is_declared(tmp):
    """Targets with no sequence: the offset cannot be measured, so say so."""
    bad = []
    tg = targets_file(tmp, "t5.csv", with_seq=False)
    pk = peaks_file(tmp, "p5.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    if rep["calibration"]["offset_ppm"] is not None:
        bad.append("an offset was reported from targets that carry no sequence")
    if not any("UNCORRECTED" in f for f in rep["findings"]):
        bad.append("the entrapment m/z was uncorrected and the report did not"
                   " say the rate is a floor")
    return bad


def case_a_real_calibration_offset_is_measured_and_applied(tmp):
    """Observed m/z 8 ppm off calculated: the probes must be moved with it.

    This is the bias that would otherwise make the control unfalsifiable. The
    planted offset is 18 ppm against a 10 ppm window, and it has to be: at 8 ppm
    an UNCORRECTED probe still lands inside a 10 ppm window, so the fixture
    could not tell a tool that applies the correction from one that ignores it.
    That is how this case first passed against code that never applied it.
    """
    bad = []
    tg = targets_file(tmp, "t6.csv", ppm_off=18.0)
    rows = peaks_for(HUMAN, ppm_off=18.0)
    for i, s in enumerate(ECOLI):               # peaks also sit at +18 ppm
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2, 18.0), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
    pk = peaks_file(tmp, "p6.jsonl", rows)
    rep, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=0.2)
    off = rep["calibration"]["offset_ppm"]
    if off is None or abs(off - 18.0) > 0.5:
        bad.append("planted a +18 ppm observed-minus-calculated offset, tool"
                   " measured %s" % ("none" if off is None else "%.2f" % off))
    if rep["arms"]["entrapment"]["matched"] == 0:
        bad.append("every entrapment probe missed although a peak sits on each"
                   " one: the measured offset was not applied to the probes")
    return bad


def case_a_calculated_coordinate_penalty_is_reported(tmp):
    """When calculated m/z costs the targets matches, the probes pay it too."""
    bad = []
    # Observed m/z 9 ppm off calculated, but the offset is NOT uniform: half the
    # targets are shifted, so the median correction leaves a residual and the
    # calc_target arm loses matches the observed-m/z arm keeps.
    lines = ["peptide,charge,mz,rt,drift"]
    for i, s in enumerate(HUMAN):
        off = 0.0 if i % 2 else 40.0
        lines.append("%s,2,%.6f,%.3f,%.3f" % (s, mz(s, 2, off), 5.0 + i, 50.0 + i))
    tg = Path(tmp) / "t7.csv"
    tg.write_text("\n".join(lines) + "\n")
    rows = []
    for i, s in enumerate(HUMAN):
        off = 0.0 if i % 2 else 40.0
        rows.append({"mz": mz(s, 2, off), "charge": 2, "rt": 5.0 + i,
                     "drift": 50.0 + i, "intensity": 1e6, "n_scans": 7})
    pk = peaks_file(tmp, "p7.jsonl", rows)
    rep, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=0.2)
    t, c = rep["arms"]["target"], rep["arms"]["calc_target"]
    if not (c["rate"] < t["rate"]):
        bad.append("calculated m/z cost the targets nothing (%.2f vs %.2f) in data"
                   " built so that it must" % (c["rate"], t["rate"]))
    elif not any("floor" in f for f in rep["findings"]):
        bad.append("calc_target matched at %.1f%% against target %.1f%% and the"
                   " report did not call the entrapment rate a floor"
                   % (100 * c["rate"], 100 * t["rate"]))
    return bad


def case_impossible_probes_do_not_pad_the_denominator(tmp):
    """A probe outside the acquired range, or at a charge no peak carries,
    contributes 0/1 and flatters the rate. Exclude and count it."""
    bad = []
    tg = targets_file(tmp, "t8.csv", z=4)   # probes inherit z=4
    pk = peaks_file(tmp, "p8.jsonl", peaks_for(HUMAN, z=2))  # no peak at z=4
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    if e["queried"]:
        bad.append("%d probe(s) at a charge no peak carries were counted as"
                   " trials" % e["queried"])
    if not e["excluded"].get("unseen_charge"):
        bad.append("probes at an unseen charge were not counted as excluded")
    return bad


def case_a_thin_peak_table_is_not_a_good_window(tmp):
    """No target matched either, so a zero probe rate says nothing."""
    bad = []
    tg = targets_file(tmp, "t9.csv")
    pk = peaks_file(tmp, "p9.jsonl", [{"mz": 999.9999, "charge": 2, "rt": 99.0,
                                       "drift": 1.0, "intensity": 1e4, "n_scans": 7}])
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    ok, why = pc.discriminated(rep)
    if ok:
        bad.append("nothing matched at all and the control reported that the"
                   " window discriminated: %s" % why)
    if not any("thin" in f for f in rep["findings"]):
        bad.append("no target matched and the report did not say a low probe"
                   " rate means a thin peak table")
    return bad


def case_too_few_probes_cannot_certify(tmp):
    """n below the point where the Wilson bound supports a claim."""
    bad = []
    tg = targets_file(tmp, "t10.csv")
    pk = peaks_file(tmp, "p10.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, seqs=ECOLI[:6], rt_tol=0.2)
    ok, why = pc.discriminated(rep)
    if ok:
        bad.append("6 probes was enough to certify the window: %s" % why)
    if not any("too wide" in f for f in rep["findings"]):
        bad.append("6 probes and the report did not say the interval is too wide")
    return bad


def case_no_entrapment_list_refuses(tmp):
    """The measurement is the entrapment arm. Without it, exit 2, never 0."""
    bad = []
    tg = targets_file(tmp, "t11.csv")
    pk = peaks_file(tmp, "p11.jsonl", peaks_for(HUMAN))
    rc, out = run(pk, tg, "--rt-tol", 0.2)
    if rc != 2:
        bad.append("no --entrapment gave exit %d, not 2" % rc)
    if "not measured" not in out:
        bad.append("no --entrapment and the output did not say the rate was not"
                   " measured")
    rc2, out2 = run(pk, tg, "--entrapment", "/dev/null", "--rt-tol", 0.2)
    if rc2 != 2:
        bad.append("an empty entrapment file gave exit %d, not 2" % rc2)
    return bad


def case_swap_is_not_reported_when_no_coordinate_is_compared(tmp):
    """The swap arm's whole question is whether rt/drift add anything. With
    neither compared there is no question, and a number would be a fiction."""
    bad = []
    tg = targets_file(tmp, "t12.csv")
    pk = peaks_file(tmp, "p12.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=None)
    s = rep["arms"]["swap"]
    if s.get("queried"):
        bad.append("the swap arm reported %d probes with no axis to swap on"
                   % s["queried"])
    if "no coordinate to swap" not in (s.get("note") or ""):
        bad.append("the swap arm was empty without saying why")
    rep2, _ = control(tmp, pk, tg, rt_tol=0.2)
    if not rep2["arms"]["swap"].get("queried"):
        bad.append("rt was compared and the swap arm still did not run")
    return bad


def case_a_ratio_against_zero_is_null_not_perfect(tmp):
    """katz_ratio is undefined against zero events. Reporting that as a huge
    separation would turn the best possible news into the strongest claim."""
    bad = []
    tg = targets_file(tmp, "t13.csv")
    pk = peaks_file(tmp, "p13.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    if e["matched"] != 0:
        bad.append("fixture no longer has a zero entrapment rate (%d matched)"
                   % e["matched"])
    if e["ratio_vs_target"] is not None:
        bad.append("a ratio was reported against zero entrapment events: %s"
                   % (e["ratio_vs_target"],))
    out = pc.render(rep)
    if "null rather than as perfect separation" not in out:
        bad.append("the rendered report does not explain what a null ratio means")
    return bad


def case_a_probe_that_is_its_own_target_is_excluded(tmp):
    """An entrapment peptide mass-degenerate with a listed target is testing
    that target, and counting it as a false positive is circular."""
    bad = []
    # Give the entrapment list a sequence that IS one of the targets.
    tg = targets_file(tmp, "t14.csv")
    pk = peaks_file(tmp, "p14.jsonl", peaks_for(HUMAN))
    seqs = ECOLI[:20] + [HUMAN[0]]
    rep, _ = control(tmp, pk, tg, seqs=seqs, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    n_excl = sum(e["excluded"].get(k, 0) for k in
                 ("self_match", "degenerate_with_a_target"))
    if not n_excl:
        bad.append("an entrapment sequence identical to a target was queried as"
                   " a decoy")
    return bad


def case_clustered_trials_are_declared(tmp):
    """Many probes matching one peak are not independent trials.

    The clustering is built from PERMUTED sequences, which is the cleanest way
    to get several distinct probes at one address: permuting a peptide keeps its
    residue composition and therefore its monoisotopic mass exactly, so these
    four probes share an m/z to the last decimal and one peak answers all of
    them. (It is also why reversed-tryptic decoys cannot test an m/z window at
    all -- reversal is a permutation.)

    An earlier version of this case used a 200000 ppm window instead. Every
    probe then fell inside its own paired target's window, was excluded as a
    self-match, and the case passed on zero probes against a tool that had no
    clustering check at all. Hence the explicit assertion that the fixture
    actually clustered.
    """
    bad = []
    perms = ["AAVLTSHGK", "VAALTSHGK", "LAAVTSHGK", "TAAVLSHGK"]
    if len({round(peptide_mass(s), 6) for s in perms}) != 1:
        return ["the fixture's permuted sequences no longer share one mass"]
    tg = targets_file(tmp, "t15.csv", seqs=HUMAN[:1])
    rows = [{"mz": mz(HUMAN[0], 2), "charge": 2, "rt": 5.0, "drift": 50.0,
             "intensity": 1e6, "n_scans": 7},
            {"mz": mz(perms[0], 2), "charge": 2, "rt": 5.0, "drift": 50.0,
             "intensity": 1e5, "n_scans": 7}]
    pk = peaks_file(tmp, "p15.jsonl", rows)
    rep, _ = control(tmp, pk, tg, seqs=perms, ppm_tol=10.0, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    if e["matched"] <= max(1, e.get("distinct_peaks_matched") or 1):
        return ["the fixture did not cluster: %d probe(s) matched %s distinct"
                " peak(s), so this case cannot test anything"
                % (e["matched"], e.get("distinct_peaks_matched"))]
    if not any("clustered" in f for f in rep["findings"]):
        bad.append("%d probes matched %d distinct peak(s) and the report did not"
                   " say the trials are clustered"
                   % (e["matched"], e["distinct_peaks_matched"]))
    return bad


def case_the_scan_floor_applies_to_probes_too(tmp):
    """A probe hit that would have been rejected as a real ion is not a false
    positive, so the same n_scans floor has to apply on both sides."""
    bad = []
    tg = targets_file(tmp, "t16.csv")
    rows = peaks_for(HUMAN)
    for i, s in enumerate(ECOLI):           # probe peaks below the floor
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 1})
    pk = peaks_file(tmp, "p16.jsonl", rows)
    rep, _ = control(tmp, pk, tg, rt_tol=0.2, min_scans=2)
    if rep["arms"]["entrapment"]["matched"]:
        bad.append("%d probe(s) matched a peak with 1 scan, below the floor that"
                   " bare_proximity requires"
                   % rep["arms"]["entrapment"]["matched"])
    rep2, _ = control(tmp, pk, tg, rt_tol=0.2, min_scans=1)
    if not rep2["arms"]["entrapment"]["matched"]:
        bad.append("the same peaks matched nothing with the floor at 1, so the"
                   " floor is not what excluded them")
    return bad


def case_a_sparse_table_says_clearing_the_bar_means_little(tmp):
    """The expectation is derived here, not borrowed.

    The 33-70% this project recorded came from reversed-tryptic decoys, and
    reversal permutes a peptide's residues and so preserves its mass exactly --
    a reversed decoy matches an m/z window 100% by construction, so that figure
    cannot describe this window and disagreeing with it would prove nothing. In
    its place the tool computes what chance alone would give from this run's own
    peak-window geometry, and when that is vanishingly small it says so: a 0%
    probe rate against a 0.001% expectation is not an achievement.

    It is reported, not enforced. See expected_collision_rate on why no gate
    hangs on it.
    """
    bad = []
    tg = targets_file(tmp, "t17.csv")
    pk = peaks_file(tmp, "p17.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=0.2)
    exp = rep.get("expected_by_chance")
    if exp is None:
        return ["no expectation was computed at all"]
    if exp >= 0.01:
        return ["the fixture is not sparse enough to test this: chance alone is"
                " %.4f" % exp]
    if not any("says very little" in f for f in rep["findings"]):
        bad.append("chance alone was %.4g%% and the report did not say a low"
                   " measured rate means very little on a table this sparse"
                   % (100 * exp))
    if "Chance alone" not in pc.render(rep):
        bad.append("the expectation is computed but never shown")
    return bad


def case_a_planted_rate_is_recovered(tmp):
    """Plant a peak on a known share of entrapment coordinates; read it back."""
    bad = []
    tg = targets_file(tmp, "t1.csv")
    # Entrapment probes are built at the E. coli m/z on a human target's rt.
    # Plant a peak at the first 10 entrapment coordinates so the rate is known.
    rows = peaks_for(HUMAN)
    planted = 0
    for i, s in enumerate(ECOLI):
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
        planted += 1
        if planted >= 10:
            break
    pk = peaks_file(tmp, "p1.jsonl", rows)
    rep, err = control(tmp, pk, tg, rt_tol=0.2)
    if rep is None:
        return ["control did not run: %s" % err]
    e = rep["arms"]["entrapment"]
    want = 10
    if e["matched"] != want:
        bad.append("planted %d entrapment hits, tool reported %d of %d probes"
                   % (want, e["matched"], e["queried"]))
    if not (e["ci"][0] <= want / e["queried"] <= e["ci"][1]):
        bad.append("the 95%% interval %.3f-%.3f does not cover the planted rate %.3f"
                   % (e["ci"][0], e["ci"][1], want / e["queried"]))
    return bad


def case_a_window_that_does_not_discriminate_is_refused(tmp):
    """Entrapment matching as often as targets must not pass."""
    bad = []
    tg = targets_file(tmp, "t2.csv")
    rows = peaks_for(HUMAN)
    for i, s in enumerate(ECOLI):               # a peak on EVERY probe
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
    pk = peaks_file(tmp, "p2.jsonl", rows)
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    ok, why = pc.discriminated(rep)
    if ok:
        bad.append("every entrapment probe found a peak and the control still"
                   " reported that the window discriminated: %s" % why)
    rc, out = run(pk, tg, "--entrapment", entrap_file(tmp, "e2.txt"),
                  "--rt-tol", 0.2)
    if rc == 0:
        bad.append("exit 0 on a window that separated nothing")
    return bad


def case_the_same_window_is_used_as_the_match(tmp):
    """Widen --precursor-ppm and both rates must move.

    A control that reimplemented the window would not track the matcher, and the
    rate would describe code that never ran.
    """
    bad = []
    # Peaks sit 30 ppm off their target, so a 10 ppm window misses and 50 hits.
    tg = targets_file(tmp, "t3.csv")
    rows = peaks_for(HUMAN, ppm_off=30.0)
    for i, s in enumerate(ECOLI):
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2, 30.0), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
    pk = peaks_file(tmp, "p3.jsonl", rows)
    narrow, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=0.2)
    wide, _ = control(tmp, pk, tg, ppm_tol=50.0, rt_tol=0.2)
    for arm in ("target", "entrapment"):
        n, w = narrow["arms"][arm], wide["arms"][arm]
        if not (w["matched"] > n["matched"]):
            bad.append("%s arm did not move when the window widened from 10 to 50"
                       " ppm (%d then %d): the control is not measuring the"
                       " matcher's window" % (arm, n["matched"], w["matched"]))
    return bad


def case_axes_never_compared_are_not_claimed(tmp):
    """No --rt-tol, or no rt column, means the rate is global m/z crowding."""
    bad = []
    tg = targets_file(tmp, "t4.csv")
    pk = peaks_file(tmp, "p4.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=None)          # tolerance never passed
    if "rt" in rep["axes_matched"]:
        bad.append("rt listed as compared with no --rt-tol")
    if not any("GLOBAL" in f for f in rep["findings"]):
        bad.append("no --rt-tol and the report did not say the rate is a global"
                   " m/z-crowding rate")
    tg2 = targets_file(tmp, "t4b.csv", with_rt=False)   # tolerance passed, column null
    rep2, _ = control(tmp, pk, tg2, rt_tol=0.2)
    if "rt" in rep2["axes_matched"]:
        bad.append("rt listed as compared although no target carried one")
    return bad


def case_an_unmeasured_calibration_offset_is_declared(tmp):
    """Targets with no sequence: the offset cannot be measured, so say so."""
    bad = []
    tg = targets_file(tmp, "t5.csv", with_seq=False)
    pk = peaks_file(tmp, "p5.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    if rep["calibration"]["offset_ppm"] is not None:
        bad.append("an offset was reported from targets that carry no sequence")
    if not any("UNCORRECTED" in f for f in rep["findings"]):
        bad.append("the entrapment m/z was uncorrected and the report did not"
                   " say the rate is a floor")
    return bad


def case_a_real_calibration_offset_is_measured_and_applied(tmp):
    """Observed m/z 8 ppm off calculated: the probes must be moved with it.

    This is the bias that would otherwise make the control unfalsifiable. The
    planted offset is 18 ppm against a 10 ppm window, and it has to be: at 8 ppm
    an UNCORRECTED probe still lands inside a 10 ppm window, so the fixture
    could not tell a tool that applies the correction from one that ignores it.
    That is how this case first passed against code that never applied it.
    """
    bad = []
    tg = targets_file(tmp, "t6.csv", ppm_off=18.0)
    rows = peaks_for(HUMAN, ppm_off=18.0)
    for i, s in enumerate(ECOLI):               # peaks also sit at +18 ppm
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2, 18.0), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 7})
    pk = peaks_file(tmp, "p6.jsonl", rows)
    rep, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=0.2)
    off = rep["calibration"]["offset_ppm"]
    if off is None or abs(off - 18.0) > 0.5:
        bad.append("planted a +18 ppm observed-minus-calculated offset, tool"
                   " measured %s" % ("none" if off is None else "%.2f" % off))
    if rep["arms"]["entrapment"]["matched"] == 0:
        bad.append("every entrapment probe missed although a peak sits on each"
                   " one: the measured offset was not applied to the probes")
    return bad


def case_a_calculated_coordinate_penalty_is_reported(tmp):
    """When calculated m/z costs the targets matches, the probes pay it too."""
    bad = []
    # Observed m/z 9 ppm off calculated, but the offset is NOT uniform: half the
    # targets are shifted, so the median correction leaves a residual and the
    # calc_target arm loses matches the observed-m/z arm keeps.
    lines = ["peptide,charge,mz,rt,drift"]
    for i, s in enumerate(HUMAN):
        off = 0.0 if i % 2 else 40.0
        lines.append("%s,2,%.6f,%.3f,%.3f" % (s, mz(s, 2, off), 5.0 + i, 50.0 + i))
    tg = Path(tmp) / "t7.csv"
    tg.write_text("\n".join(lines) + "\n")
    rows = []
    for i, s in enumerate(HUMAN):
        off = 0.0 if i % 2 else 40.0
        rows.append({"mz": mz(s, 2, off), "charge": 2, "rt": 5.0 + i,
                     "drift": 50.0 + i, "intensity": 1e6, "n_scans": 7})
    pk = peaks_file(tmp, "p7.jsonl", rows)
    rep, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=0.2)
    t, c = rep["arms"]["target"], rep["arms"]["calc_target"]
    if not (c["rate"] < t["rate"]):
        bad.append("calculated m/z cost the targets nothing (%.2f vs %.2f) in data"
                   " built so that it must" % (c["rate"], t["rate"]))
    elif not any("floor" in f for f in rep["findings"]):
        bad.append("calc_target matched at %.1f%% against target %.1f%% and the"
                   " report did not call the entrapment rate a floor"
                   % (100 * c["rate"], 100 * t["rate"]))
    return bad


def case_impossible_probes_do_not_pad_the_denominator(tmp):
    """A probe outside the acquired range, or at a charge no peak carries,
    contributes 0/1 and flatters the rate. Exclude and count it."""
    bad = []
    tg = targets_file(tmp, "t8.csv", z=4)   # probes inherit z=4
    pk = peaks_file(tmp, "p8.jsonl", peaks_for(HUMAN, z=2))  # no peak at z=4
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    if e["queried"]:
        bad.append("%d probe(s) at a charge no peak carries were counted as"
                   " trials" % e["queried"])
    if not e["excluded"].get("unseen_charge"):
        bad.append("probes at an unseen charge were not counted as excluded")
    return bad


def case_a_thin_peak_table_is_not_a_good_window(tmp):
    """No target matched either, so a zero probe rate says nothing."""
    bad = []
    tg = targets_file(tmp, "t9.csv")
    pk = peaks_file(tmp, "p9.jsonl", [{"mz": 999.9999, "charge": 2, "rt": 99.0,
                                       "drift": 1.0, "intensity": 1e4, "n_scans": 7}])
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    ok, why = pc.discriminated(rep)
    if ok:
        bad.append("nothing matched at all and the control reported that the"
                   " window discriminated: %s" % why)
    if not any("thin" in f for f in rep["findings"]):
        bad.append("no target matched and the report did not say a low probe"
                   " rate means a thin peak table")
    return bad


def case_too_few_probes_cannot_certify(tmp):
    """n below the point where the Wilson bound supports a claim."""
    bad = []
    tg = targets_file(tmp, "t10.csv")
    pk = peaks_file(tmp, "p10.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, seqs=ECOLI[:6], rt_tol=0.2)
    ok, why = pc.discriminated(rep)
    if ok:
        bad.append("6 probes was enough to certify the window: %s" % why)
    if not any("too wide" in f for f in rep["findings"]):
        bad.append("6 probes and the report did not say the interval is too wide")
    return bad


def case_no_entrapment_list_refuses(tmp):
    """The measurement is the entrapment arm. Without it, exit 2, never 0."""
    bad = []
    tg = targets_file(tmp, "t11.csv")
    pk = peaks_file(tmp, "p11.jsonl", peaks_for(HUMAN))
    rc, out = run(pk, tg, "--rt-tol", 0.2)
    if rc != 2:
        bad.append("no --entrapment gave exit %d, not 2" % rc)
    if "not measured" not in out:
        bad.append("no --entrapment and the output did not say the rate was not"
                   " measured")
    rc2, out2 = run(pk, tg, "--entrapment", "/dev/null", "--rt-tol", 0.2)
    if rc2 != 2:
        bad.append("an empty entrapment file gave exit %d, not 2" % rc2)
    return bad


def case_swap_is_not_reported_when_no_coordinate_is_compared(tmp):
    """The swap arm's whole question is whether rt/drift add anything. With
    neither compared there is no question, and a number would be a fiction."""
    bad = []
    tg = targets_file(tmp, "t12.csv")
    pk = peaks_file(tmp, "p12.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=None)
    s = rep["arms"]["swap"]
    if s.get("queried"):
        bad.append("the swap arm reported %d probes with no axis to swap on"
                   % s["queried"])
    if "no coordinate to swap" not in (s.get("note") or ""):
        bad.append("the swap arm was empty without saying why")
    rep2, _ = control(tmp, pk, tg, rt_tol=0.2)
    if not rep2["arms"]["swap"].get("queried"):
        bad.append("rt was compared and the swap arm still did not run")
    return bad


def case_a_ratio_against_zero_is_null_not_perfect(tmp):
    """katz_ratio is undefined against zero events. Reporting that as a huge
    separation would turn the best possible news into the strongest claim."""
    bad = []
    tg = targets_file(tmp, "t13.csv")
    pk = peaks_file(tmp, "p13.jsonl", peaks_for(HUMAN))
    rep, _ = control(tmp, pk, tg, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    if e["matched"] != 0:
        bad.append("fixture no longer has a zero entrapment rate (%d matched)"
                   % e["matched"])
    if e["ratio_vs_target"] is not None:
        bad.append("a ratio was reported against zero entrapment events: %s"
                   % (e["ratio_vs_target"],))
    out = pc.render(rep)
    if "null rather than as perfect separation" not in out:
        bad.append("the rendered report does not explain what a null ratio means")
    return bad


def case_a_probe_that_is_its_own_target_is_excluded(tmp):
    """An entrapment peptide mass-degenerate with a listed target is testing
    that target, and counting it as a false positive is circular."""
    bad = []
    # Give the entrapment list a sequence that IS one of the targets.
    tg = targets_file(tmp, "t14.csv")
    pk = peaks_file(tmp, "p14.jsonl", peaks_for(HUMAN))
    seqs = ECOLI[:20] + [HUMAN[0]]
    rep, _ = control(tmp, pk, tg, seqs=seqs, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    n_excl = sum(e["excluded"].get(k, 0) for k in
                 ("self_match", "degenerate_with_a_target"))
    if not n_excl:
        bad.append("an entrapment sequence identical to a target was queried as"
                   " a decoy")
    return bad


def case_clustered_trials_are_declared(tmp):
    """Many probes matching one peak are not independent trials.

    The clustering is built from PERMUTED sequences, which is the cleanest way
    to get several distinct probes at one address: permuting a peptide keeps its
    residue composition and therefore its monoisotopic mass exactly, so these
    four probes share an m/z to the last decimal and one peak answers all of
    them. (It is also why reversed-tryptic decoys cannot test an m/z window at
    all -- reversal is a permutation.)

    An earlier version of this case used a 200000 ppm window instead. Every
    probe then fell inside its own paired target's window, was excluded as a
    self-match, and the case passed on zero probes against a tool that had no
    clustering check at all. Hence the explicit assertion that the fixture
    actually clustered.
    """
    bad = []
    perms = ["AAVLTSHGK", "VAALTSHGK", "LAAVTSHGK", "TAAVLSHGK"]
    if len({round(peptide_mass(s), 6) for s in perms}) != 1:
        return ["the fixture's permuted sequences no longer share one mass"]
    tg = targets_file(tmp, "t15.csv", seqs=HUMAN[:1])
    rows = [{"mz": mz(HUMAN[0], 2), "charge": 2, "rt": 5.0, "drift": 50.0,
             "intensity": 1e6, "n_scans": 7},
            {"mz": mz(perms[0], 2), "charge": 2, "rt": 5.0, "drift": 50.0,
             "intensity": 1e5, "n_scans": 7}]
    pk = peaks_file(tmp, "p15.jsonl", rows)
    rep, _ = control(tmp, pk, tg, seqs=perms, ppm_tol=10.0, rt_tol=0.2)
    e = rep["arms"]["entrapment"]
    if e["matched"] <= max(1, e.get("distinct_peaks_matched") or 1):
        return ["the fixture did not cluster: %d probe(s) matched %s distinct"
                " peak(s), so this case cannot test anything"
                % (e["matched"], e.get("distinct_peaks_matched"))]
    if not any("clustered" in f for f in rep["findings"]):
        bad.append("%d probes matched %d distinct peak(s) and the report did not"
                   " say the trials are clustered"
                   % (e["matched"], e["distinct_peaks_matched"]))
    return bad


def case_the_scan_floor_applies_to_probes_too(tmp):
    """A probe hit that would have been rejected as a real ion is not a false
    positive, so the same n_scans floor has to apply on both sides."""
    bad = []
    tg = targets_file(tmp, "t16.csv")
    rows = peaks_for(HUMAN)
    for i, s in enumerate(ECOLI):           # probe peaks below the floor
        t_i = i % len(HUMAN)
        rows.append({"mz": mz(s, 2), "charge": 2, "rt": 5.0 + t_i * 1.0,
                     "drift": 50.0 + t_i, "intensity": 1e5, "n_scans": 1})
    pk = peaks_file(tmp, "p16.jsonl", rows)
    rep, _ = control(tmp, pk, tg, rt_tol=0.2, min_scans=2)
    if rep["arms"]["entrapment"]["matched"]:
        bad.append("%d probe(s) matched a peak with 1 scan, below the floor that"
                   " bare_proximity requires"
                   % rep["arms"]["entrapment"]["matched"])
    rep2, _ = control(tmp, pk, tg, rt_tol=0.2, min_scans=1)
    if not rep2["arms"]["entrapment"]["matched"]:
        bad.append("the same peaks matched nothing with the floor at 1, so the"
                   " floor is not what excluded them")
    return bad


def case_a_rate_below_chance_is_declared(tmp):
    """A probe set suppressed for a reason other than the window.

    The expectation is derived from this run's own peak-window geometry, not
    borrowed: the 33-70% this project recorded came from reversed-tryptic
    decoys, and reversal preserves a peptide's mass exactly, so that figure
    cannot describe an m/z window. Here the peak table is made dense enough that
    chance alone should put most probes on a peak, and the probes are then given
    a charge no peak carries -- so the measured rate collapses for a reason that
    has nothing to do with the window, which is what the tool must say.
    """
    bad = []
    tg = targets_file(tmp, "t17.csv")
    # A dense peak table: many peaks at z=2 across a narrow m/z span, so the
    # 10 ppm windows cover a large fraction of the acquired range.
    rows = []
    base = 600.0
    for i in range(400):
        rows.append({"mz": base + i * 0.0005, "charge": 2, "rt": 10.0,
                     "drift": 50.0, "intensity": 1e5, "n_scans": 7})
    pk = peaks_file(tmp, "p17.jsonl", rows)
    rep, _ = control(tmp, pk, tg, ppm_tol=10.0, rt_tol=None)
    exp = rep.get("expected_by_chance")
    if not exp or exp < 0.2:
        return ["the fixture is not dense enough to test this: chance alone is"
                " %s, so a measured zero is unremarkable" % exp]
    e = rep["arms"]["entrapment"]
    if e.get("queried") and e["ci"][1] < 0.25 * exp:
        if not any("Chance alone should have put" in f for f in rep["findings"]):
            bad.append("chance alone was %.1f%% and the measured upper bound %.1f%%,"
                       " and the report did not say the probes are suppressed"
                       % (100 * exp, 100 * e["ci"][1]))
    else:
        bad.append("the fixture did not produce a rate below chance: expected"
                   " %.3f, measured upper bound %s"
                   % (exp, e["ci"][1] if e.get("queried") else "no probes"))
    return bad


CASES = [case_a_sparse_table_says_clearing_the_bar_means_little,
         case_a_planted_rate_is_recovered,
         case_a_window_that_does_not_discriminate_is_refused,
         case_the_same_window_is_used_as_the_match,
         case_axes_never_compared_are_not_claimed,
         case_an_unmeasured_calibration_offset_is_declared,
         case_a_real_calibration_offset_is_measured_and_applied,
         case_a_calculated_coordinate_penalty_is_reported,
         case_impossible_probes_do_not_pad_the_denominator,
         case_a_thin_peak_table_is_not_a_good_window,
         case_too_few_probes_cannot_certify,
         case_no_entrapment_list_refuses,
         case_swap_is_not_reported_when_no_coordinate_is_compared,
         case_a_ratio_against_zero_is_null_not_perfect,
         case_a_probe_that_is_its_own_target_is_excluded,
         case_clustered_trials_are_declared,
         case_the_scan_floor_applies_to_probes_too]


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
            for line in bad or []:
                print("     %s" % line)
    print("\n%d/%d passed" % (len(CASES) - failed, len(CASES)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
