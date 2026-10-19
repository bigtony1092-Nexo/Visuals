#!/usr/bin/env python3
"""Prove the proximity-control tests can fail, by breaking the code on purpose.

This repository's recurring defect is a check that cannot fail, and a test suite
is not exempt. Four of the cases written for the proximity control passed on the
first run against code that did not implement what they claimed to test: an 8
ppm planted calibration offset still lands inside a 10 ppm window; a clustering
case used a 200000 ppm window, which excluded every probe as a self-match and
passed on zero trials; and twice a fixture protein digested to nothing because
trypsin does not cut before proline. A green suite proved none of them.

So each entry below breaks one specific behaviour and names the case that must
notice. The suite is run once per mutation, which is why this is not wired into
the per-file suites: it takes a few minutes. Run it before committing a change
to any gate it covers.

    python3 tools/mutation_check.py
    python3 tools/mutation_check.py --only index
    python3 tools/mutation_check.py --list

Three properties worth stating, because each stops a way this tool could have
reported success without checking anything:

A MUTATION THAT NO LONGER APPLIES IS A FAILURE, not a skip. Its `old` text is
copied from the source, so rewording a line detaches it -- that happened twice
while this was being written, once on a rename and once on a capitalisation fix.
Detached, it silently proves nothing, which is the exact shape of defect this
tool exists to catch. It is reported as a failure and names the file.

THE BASELINE MUST BE GREEN FIRST. If a case already fails, every mutation
"catches" it and the whole run is meaningless. Each suite is run unmutated
before anything is touched, and a red baseline stops the run.

THE SOURCE IS RESTORED AND THE RESTORATION IS VERIFIED. Every file is restored
from content held in memory, in a finally block, and compared byte for byte
afterwards. If a restore ever fails, the tool says which file and tells you to
check it out of git rather than carrying on.

Not every property can be proved this way. A mutation whose effect a second,
independent guard masks will not be caught, and inventing a case to catch it
would be testing the guard rather than the behaviour. Two were dropped for that
reason rather than kept as evidence of something they do not show: the
sort-before-sample in make_entrapment.py guards against a dict ordering CPython
does not produce, and the header-comment skip in load_sequences is masked by the
mass filter behind it. Both are marked as defensive in the code instead.

Exit 0 when every mutation was caught by the case written for it, 1 otherwise.
"""

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "tools"

# (name, file, old, new, suite, cases that must fail)
# `old` is copied verbatim from the source and must appear exactly once.
MUTATIONS = [
    ('offset not applied to probes',
     'proximity_control.py',
     '        if offset_ppm is not None:\n            mz *= (1.0 + offset_ppm / 1e6)\n',
     '',
     'test_proximity_control',
     ['case_a_real_calibration_offset_is_measured_and_applied']),
    ('discriminated always says yes',
     'proximity_control.py',
     '    e = rep["arms"]["entrapment"]\n    t = rep["arms"]["target"]\n    if not e.get("queried")',
     '    return True, "x"\n    e = rep["arms"]["entrapment"]\n    t = rep["arms"]["target"]\n    if not e.get("queried")',
     'test_proximity_control',
     ['case_a_window_that_does_not_discriminate_is_refused', 'case_a_thin_peak_table_is_not_a_good_window', 'case_too_few_probes_cannot_certify']),
    ('scan floor ignored for probes',
     'proximity_control.py',
     '            if pk["n_scans"] is not None and pk["n_scans"] < min_scans:\n                continue          # would have been rejected as a real ion too\n',
     '',
     'test_proximity_control',
     ['case_the_scan_floor_applies_to_probes_too']),
    ('out-of-range and unseen-charge probes kept',
     'proximity_control.py',
     '        if not (mz_lo <= mz <= mz_hi):\n            excluded["out_of_range"] += 1\n            continue\n        if z not in charges:\n            excluded["unseen_charge"] += 1\n            continue\n',
     '',
     'test_proximity_control',
     ['case_impossible_probes_do_not_pad_the_denominator']),
    ('self/degenerate probes kept',
     'proximity_control.py',
     '        if ppm(mz, t["mz"]) <= ppm_tol * SELF_PPM_FACTOR:\n            excluded["self_match"] += 1\n            continue\n        if any(zz == z and ppm(mz, m) <= ppm_tol for m, zz in tmz):\n            excluded["degenerate_with_a_target"] += 1\n            continue\n',
     '',
     'test_proximity_control',
     ['case_a_probe_that_is_its_own_target_is_excluded']),
    ('MIN_CONTROL_PROBES bound removed',
     'proximity_control.py',
     '    if e["queried"] < MIN_CONTROL_PROBES:\n        return False, "%d probes is too few for the bound to mean anything" % e["queried"]\n',
     '',
     'test_proximity_control',
     ['case_too_few_probes_cannot_certify']),
    ('exit 0 with no entrapment list',
     'proximity_control.py',
     '        print(render(rep))\n        return 2\n',
     '        print(render(rep))\n        return 0\n',
     'test_proximity_control',
     ['case_no_entrapment_list_refuses']),
    ('the global-crowding finding dropped',
     'proximity_control.py',
     '    if "rt" not in rep["axes_matched"] and "drift" not in rep["axes_matched"]:',
     '    if False:',
     'test_proximity_control',
     ['case_axes_never_compared_are_not_claimed']),
    ('the uncorrected-calibration finding dropped',
     'proximity_control.py',
     '    if rep["calibration"]["offset_ppm"] is None:\n        out.append("No target carried a sequence',
     '    if False:\n        out.append("No target carried a sequence',
     'test_proximity_control',
     ['case_an_unmeasured_calibration_offset_is_declared']),
    ('the calculated-coordinate floor finding dropped',
     'proximity_control.py',
     '        if c["rate"] < 0.75 * t["rate"]:',
     '        if False:',
     'test_proximity_control',
     ['case_a_calculated_coordinate_penalty_is_reported']),
    ('the thin-table finding dropped',
     'proximity_control.py',
     '    if t.get("queried") and not t.get("matched"):',
     '    if False:',
     'test_proximity_control',
     ['case_a_thin_peak_table_is_not_a_good_window']),
    ('the clustered-trials finding dropped',
     'proximity_control.py',
     '    if e.get("queried") and e.get("distinct_peaks_matched") is not None \\\n            and e["matched"] > max(1, e["distinct_peaks_matched"]):',
     '    if False:',
     'test_proximity_control',
     ['case_clustered_trials_are_declared']),
    ('swap arm runs with no axis to swap on',
     'proximity_control.py',
     '    if "rt" in axes or "drift" in axes:',
     '    if True:',
     'test_proximity_control',
     ['case_swap_is_not_reported_when_no_coordinate_is_compared']),
    ('a ratio is invented against zero events',
     'proximity_control.py',
     '        a["ratio_vs_target"] = katz_ratio(t["matched"], t["queried"],\n                                          a["matched"], a["queried"])',
     '        a["ratio_vs_target"] = katz_ratio(t["matched"], t["queried"],\n                                          a["matched"], a["queried"]) or (99.0, 99.0)',
     'test_proximity_control',
     ['case_a_ratio_against_zero_is_null_not_perfect']),
    ('probe matches undercounted by one',
     'proximity_control.py',
     '        if best is not None:\n            k += 1',
     '        if best is not None:\n            k += 0',
     'test_proximity_control',
     ['case_a_planted_rate_is_recovered']),
    ('index bound too tight by one ulp of ppm (window silently narrowed)',
     'proximity_control.py',
     '    half = ppm_tol * probe["mz"] / 1e6',
     '    half = ppm_tol * probe["mz"] / 1e6 * 0.999',
     'test_proximity_control',
     ['case_the_candidate_index_changes_no_answer']),
    ('index ignores probes with no charge',
     'proximity_control.py',
     '        buckets = list(index.values())',
     '        buckets = []',
     'test_proximity_control',
     ['case_the_candidate_index_changes_no_answer']),
    ('index takes the first candidate, not the nearest',
     'proximity_control.py',
     '            if best is None or d < best[0]:',
     '            if best is None:',
     'test_proximity_control',
     ['case_the_candidate_index_changes_no_answer']),
    ('the sparse-table caveat dropped',
     'proximity_control.py',
     '    if e.get("queried") and exp is not None and exp < 0.01:',
     '    if False:',
     'test_proximity_control',
     ['case_a_sparse_table_says_clearing_the_bar_means_little']),
    ('the expectation never computed',
     'proximity_control.py',
     '    rep["expected_by_chance"] = expected_collision_rate(\n        prows, ent_probes, ppm_tol, rt_tol, axes)',
     '    rep["expected_by_chance"] = None',
     'test_proximity_control',
     ['case_a_sparse_table_says_clearing_the_bar_means_little']),
    ('window hardcoded instead of taken from the matcher',
     'proximity_control.py',
     '                         ppm_tol, rt_tol, drift_tol, use_drift)',
     '                         10.0, rt_tol, drift_tol, use_drift)',
     'test_proximity_control',
     ['case_the_same_window_is_used_as_the_match']),
    ('per-ion control requirement disabled',
     'validate_export.py',
     '    if tier == "bare_proximity" and ion.get("annotation"):\n        ctrl = conf.get("control")',
     '    if False:\n        ctrl = conf.get("control")',
     'test_validate_export',
     ['an_annotated_bare_ion_with_no_control', 'a_control_whose_rate_disagrees_with_its_counts', 'a_control_that_names_no_axes']),
    ('interval-overlap gate removed',
     'validate_export.py',
     '        if e_ci[1] >= t_ci[0]:',
     '        if False:',
     'test_validate_export',
     ['a_window_that_did_not_discriminate']),
    ('too-few-probes gate removed',
     'validate_export.py',
     '    if n < MIN_CONTROL_PROBES:',
     '    if False:',
     'test_validate_export',
     ['a_control_with_too_few_probes']),
    ('target-comparator gate removed',
     'validate_export.py',
     '    if not tgt.get("queried"):',
     '    if False and not tgt.get("queried"):',
     'test_validate_export',
     ['a_control_with_no_target_comparator']),
    ('thin-table gate removed',
     'validate_export.py',
     '    if not tgt.get("matched"):',
     '    if False:',
     'test_validate_export',
     ['a_control_whose_targets_matched_nothing']),
    ('missing-manifest-control gate removed',
     'validate_export.py',
     '    if not isinstance(pc, dict):',
     '    if False:',
     'test_validate_export',
     ['a_manifest_with_no_proximity_control']),
    ('proximity_control never reaches the manifest',
     'aria_run.py',
     '        "proximity_control": control_rep,',
     '        "proximity_control": None,',
     'test_aria_run',
     ['case_all_bare_says_whether_the_rate_was_measured']),
    ('the run always claims the rate was measured',
     'aria_run.py',
     '        if control_rep is None:\n            rate = ("That window\'s false-positive rate was NOT measured on this run,"\n                    " so there is no number behind the tier")',
     '        if False:\n            rate = ("x")',
     'test_aria_run',
     ['case_all_bare_says_whether_the_rate_was_measured']),
    ('the panel caption drops the rate entirely',
     'render_panels.py',
     '            tail = ("That window\'s false-positive rate was NOT measured on this"\n                    " run, so nothing here says how often it is wrong.")',
     '            tail = ("")',
     'test_render_panels',
     ['case_an_all_bare_page_says_whether_the_rate_was_measured']),
    ('controls pooled regardless of window',
     'merge_exports.py',
     '    if len(windows) > 1:',
     '    if False:',
     'test_merge_exports',
     ['case_controls_pool_only_when_the_window_is_the_same']),
    ('a missing control silently ignored',
     'merge_exports.py',
     '    if missing:',
     '    if False:',
     'test_merge_exports',
     ['case_a_run_with_no_control_leaves_the_cohort_without_one']),
    ('pooled control never reaches the manifest',
     'merge_exports.py',
     '        "proximity_control": pooled_control,',
     '        "proximity_control": None,',
     'test_merge_exports',
     ['case_controls_pool_only_when_the_window_is_the_same']),
    ('counts not actually summed',
     'merge_exports.py',
     '            acc["queried"] += a["queried"]\n            acc["matched"] += a["matched"]',
     '            acc["queried"] = a["queried"]\n            acc["matched"] = a["matched"]',
     'test_merge_exports',
     ['case_controls_pool_only_when_the_window_is_the_same']),
    ('the shared-sequence screen does nothing',
     'make_entrapment.py',
     '        peptides = [s for s in peptides if s not in shared]',
     '        peptides = list(peptides)',
     'test_make_entrapment',
     ['case_a_peptide_that_is_in_the_sample_is_not_entrapment']),
    ('unmassable peptides kept',
     'make_entrapment.py',
     '    return [s for s in seqs if peptide_mass(s) is not None]',
     '    return list(seqs)',
     'test_make_entrapment',
     ['case_unmassable_residues_are_dropped']),
    ('a too-small list reported as success',
     'make_entrapment.py',
     '    if len(keep) < MIN_CONTROL_PROBES:',
     '    if False:',
     'test_make_entrapment',
     ['case_a_list_too_small_to_certify_says_so']),
    ('provenance hash omitted',
     'make_entrapment.py',
     '        "# source sha256    %s" % sha256(args.fasta),',
     '        "# source sha256    (omitted)",',
     'test_make_entrapment',
     ['case_the_list_records_what_it_came_from']),
    ('an empty list written instead of refusing',
     'make_entrapment.py',
     '        if not peptides:\n            print("every candidate peptide also occurs in %s, so none of them is"\n                  " entrapment." % args.exclude_fasta, file=sys.stderr)\n            return 2',
     '        pass',
     'test_make_entrapment',
     ['case_everything_shared_is_refused_not_emptied']),
    ('non-FASTA input accepted',
     'make_entrapment.py',
     '    if not text.lstrip().startswith(">"):',
     '    if False:',
     'test_make_entrapment',
     ['case_a_file_that_is_not_a_fasta_is_refused']),
    ('the unrecorded-modification exclusion disabled',
     'proximity_control.py',
     '        if abs(err) > MAX_CALIBRATION_PPM:',
     '        if False:',
     'test_proximity_control',
     ['case_an_unrecorded_modification_is_kept_out_of_the_calibration']),
    ('calc_target keeps targets whose mass contradicts their sequence',
     'proximity_control.py',
     '        if (seq, t["charge"]) in skip:',
     '        if False:',
     'test_proximity_control',
     ['case_an_unrecorded_modification_is_kept_out_of_the_calibration']),
    ('a modification is named past the tolerance that distinguishes it',
     'proximity_control.py',
     '        if abs(delta - mass) <= TOL:',
     '        if abs(delta - mass) <= 8.0:',
     'test_proximity_control',
     ['case_a_difference_that_matches_nothing_is_not_named']),
    ('the careless-target-list finding dropped',
     'proximity_control.py',
     '    if mx and n_seq and len(mx) / n_seq > 0.10:',
     '    if False:',
     'test_proximity_control',
     ['case_a_list_careless_about_its_modifications_is_declared']),
]


def run_suite(suite):
    """(cases that failed, whether the suite ran at all)."""
    p = subprocess.run([sys.executable, str(TOOLS / (suite + ".py"))],
                       capture_output=True, text=True)
    out = p.stdout + p.stderr
    if "passed" not in out:
        return set(), out
    return set(re.findall(r"FAIL (\S+)", out)), None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="run only mutations whose name or file contains this")
    ap.add_argument("--list", action="store_true", help="list the mutations and exit")
    args = ap.parse_args(argv)

    muts = MUTATIONS
    if args.only:
        q = args.only.lower()
        muts = [m for m in muts if q in m[0].lower() or q in m[1].lower()]
        if not muts:
            print("no mutation matches %r" % args.only, file=sys.stderr)
            return 1
    if args.list:
        for name, f, _, _, suite, cases in muts:
            print("%-26s %-22s %s" % (f, suite, name))
        print("\n%d mutation(s)" % len(muts))
        return 0

    # A red baseline makes every mutation look caught, so check it first.
    print("baseline, before anything is touched:")
    for suite in sorted({m[4] for m in muts}):
        failing, broke = run_suite(suite)
        if broke is not None:
            print("  %-24s DID NOT RUN" % suite)
            print("     %s" % broke.strip().splitlines()[-1] if broke.strip() else "")
            return 1
        if failing:
            print("  %-24s ALREADY FAILING: %s" % (suite, ", ".join(sorted(failing))))
            print("\nFix the suite first. Against a red baseline every mutation"
                  " below would be reported as caught.")
            return 1
        print("  %-24s green" % suite)
    print()

    originals = {}
    problems = []
    t0 = time.time()
    try:
        for name, fname, old, new, suite, cases in muts:
            path = TOOLS / fname
            if fname not in originals:
                originals[fname] = path.read_text()
            src = originals[fname]
            n = src.count(old)
            if n != 1:
                problems.append(
                    "DETACHED  %-46s its `old` text appears %d times in %s, not once."
                    " The code was reworded and this mutation now proves nothing."
                    % (name, n, fname))
                print("DETACHED %s" % name)
                continue
            path.write_text(src.replace(old, new, 1))
            try:
                failing, broke = run_suite(suite)
            finally:
                path.write_text(src)
            if broke is not None:
                problems.append("CRASHED   %-46s %s could not run under this mutation"
                                % (name, suite))
                print("CRASHED  %s" % name)
                continue
            caught = set(cases) & failing
            if caught:
                print("caught   %-46s by %s" % (name, ", ".join(sorted(caught))))
            else:
                problems.append(
                    "NOT CAUGHT %-45s expected %s to fail; actually failing: %s"
                    % (name, ", ".join(cases), ", ".join(sorted(failing)) or "nothing"))
                print("NOT CAUGHT %s" % name)
    finally:
        bad_restore = []
        for fname, src in originals.items():
            path = TOOLS / fname
            if path.read_text() != src:
                path.write_text(src)
                if path.read_text() != src:
                    bad_restore.append(fname)
        if bad_restore:
            print("\nCOULD NOT RESTORE: %s" % ", ".join(bad_restore), file=sys.stderr)
            print("Check those files out of git before doing anything else.",
                  file=sys.stderr)
            return 1

    print("\n%d mutation(s) in %.0fs" % (len(muts), time.time() - t0))
    if problems:
        print()
        for p in problems:
            print("  " + p)
        print("\n%d of %d not accounted for." % (len(problems), len(muts)))
        return 1
    print("every mutation was caught by the case written for it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
