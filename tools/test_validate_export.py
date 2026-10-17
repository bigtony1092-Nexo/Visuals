#!/usr/bin/env python3
"""Prove `validate_export.py` catches what it claims to.

A validator that only ever says "ok" is worse than none, so each case here takes the
shipped examples, injects exactly one contract violation, and asserts the validator
names it and exits with the right status. The last case mutates nothing: the examples
must come back clean, or every other case is measuring noise.

    python3 tools/test_validate_export.py
    python3 tools/test_validate_export.py -v     # print each finding

No pytest, no dependencies. Exit status is 1 if any case fails.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"
VALIDATOR = ROOT / "tools" / "validate_export.py"
TABLES = ("ions", "proteins", "spectra", "traces", "fragment_hits")


def build(directory, mutate):
    """Write the examples into an export layout, with one mutation applied."""
    manifest = json.loads((EXAMPLES / "manifest.example.json").read_text())
    tables = {
        name: [
            json.loads(line)
            for line in (EXAMPLES / f"{name}.example.jsonl").read_text().splitlines()
            if line.strip()
        ]
        for name in TABLES
    }
    mutate(manifest, tables)
    (directory / "manifest.json").write_text(json.dumps(manifest))
    for name, rows in tables.items():
        (directory / f"{name}.jsonl").write_text(
            "".join(json.dumps(row) + "\n" for row in rows)
        )


def run(directory):
    proc = subprocess.run(
        [sys.executable, str(VALIDATOR), str(directory), "--excerpt"],
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout + proc.stderr


# --- one violation each ------------------------------------------------------------
# Section numbers refer to docs/export-contract.md.


def quad_window_carries_mobility(manifest, tables):
    """Section 3: a SONAR quadrupole position presented as a mobility."""
    ion = tables["ions"][0]
    ion["coord"]["drift_kind"] = "quad_window"
    ion["coord"]["one_over_k0"] = 0.94
    for sample in manifest["samples"]:
        if sample["sample_id"] == ion["sample_id"]:
            sample["drift_kind"] = "quad_window"


def bare_proximity_claims_decoy_control(manifest, tables):
    """Section 4: the tier that is not decoy-controlled saying it is."""
    tables["ions"][0]["confidence"] = {"tier": "bare_proximity", "decoy_controlled": True}


def validated_tier_disclaims_decoy_control(manifest, tables):
    """Section 4: the inverse — a validated tier that was not decoy-controlled."""
    tables["ions"][0]["confidence"].update(
        {"tier": "atlas_consensus", "decoy_controlled": False}
    )


def monoisotope_off_coord(manifest, tables):
    """Section 2: k=0 is the monoisotopic peak and sits at coord.mz."""
    tables["ions"][0]["isotopes"][0]["mz"] = 999.1234


def envelope_on_a_different_basis(manifest, tables):
    """Section 2: the most intense isotope is scaled to 100."""
    for isotope in tables["ions"][0]["isotopes"]:
        isotope["rel"] = isotope["rel"] / 100.0


def reused_ion_uid(manifest, tables):
    """Section 2: ion_uid is content-addressed and never reassigned."""
    tables["ions"][1]["ion_uid"] = tables["ions"][0]["ion_uid"]


def ledger_out_of_order(manifest, tables):
    """Section 2: events[] is append-only and ordered oldest first."""
    events = tables["ions"][0]["events"]
    events[0], events[1] = events[1], events[0]


def empty_ledger(manifest, tables):
    """Section 2: the ledger is also what the animation replays."""
    tables["ions"][0]["events"] = []


def ion_names_an_absent_sample(manifest, tables):
    tables["ions"][0]["sample_id"] = "NOPE-999"


def drift_kind_disagrees_with_sample(manifest, tables):
    """Section 3: drift_kind is keyed to vendor and method, so it cannot differ."""
    coord = tables["ions"][0]["coord"]
    coord.update({"drift_kind": "absent", "drift": None, "one_over_k0": None, "ccs": None})


def ar1_without_a_basis(manifest, tables):
    """Section 6: the two denominators give different numbers."""
    tables["spectra"][0].pop("ar1_basis")


def ar1_normalized_within_spectrum(manifest, tables):
    """Section 6: dividing by the largest cancels the precursor exactly."""
    peaks = [p for p in tables["spectra"][0]["peaks"] if p.get("ar1") is not None]
    top = max(p["ar1"] for p in peaks)
    for peak in peaks:
        peak["ar1"] = peak["ar1"] / top
    tables["spectra"][0]["ar1_sum"] = sum(p["ar1"] for p in peaks)


def spectrum_for_no_ion(manifest, tables):
    tables["spectra"][0]["ion_uid"] = "deadbeefdeadbeef"


def two_traces_for_one_isotope(manifest, tables):
    """Section 11: one trace per isotope, so the ratio can be followed scan by scan."""
    tables["traces"][2]["isotope_k"] = tables["traces"][1]["isotope_k"]


def apex_outside_its_own_trace(manifest, tables):
    tables["traces"][0]["apex_rt_min"] = 999.0


def top3avg_padded(manifest, tables):
    """Section 5: padding to 3 is forbidden."""
    tables["proteins"][0].update({"n_peptides_matched": 1, "n_peptides_used": 3})


def peptide_count_disagrees_with_list(manifest, tables):
    tables["proteins"][0]["n_peptides_used"] = 2


def declared_spectrum_kind_absent(manifest, tables):
    """A panel gated on this kind would be selected, then come back empty."""
    tables["spectra"] = [s for s in tables["spectra"] if s["kind"] != "ms1_isotope"]


def protein_twice_in_one_sample(manifest, tables):
    tables["proteins"][1] = dict(tables["proteins"][0])


def impossible_decimation(manifest, tables):
    """Section 10: kept cannot exceed of."""
    manifest["decimation"]["kept"] = manifest["decimation"]["of"] + 1


def table_present_but_undeclared(manifest, tables):
    """The data is there and the panel blanks anyway."""
    manifest["tables"]["traces"] = False


# --- fragment-first, section 7 -----------------------------------------------------


def _hit(tables):
    return tables["fragment_hits"][0]


def control_clears_the_gate(manifest, tables):
    _hit(tables)["control"]["reversed_matched"] = 3


def entrapment_clears_the_gate(manifest, tables):
    _hit(tables)["control"]["entrapment_matched"] = 5


def gate_demands_every_fragment(manifest, tables):
    _hit(tables)["k_of_n"]["required"] = 8


def matched_count_disagrees_with_list(manifest, tables):
    _hit(tables)["k_of_n"]["matched"] = 3


def label_both_matched_and_missed(manifest, tables):
    _hit(tables)["missed"].append("b3")
    _hit(tables)["k_of_n"]["queried"] = 9
    _hit(tables)["queried"]["n_fragments"] = 9


def evidence_disagrees_with_precursor(manifest, tables):
    _hit(tables)["precursor"]["found"] = False


def delta_mass_not_derivable(manifest, tables):
    _hit(tables)["delta_mass"]["value"] = 42.0


def candidate_needs_a_residue_the_peptide_lacks(manifest, tables):
    _hit(tables)["delta_mass"]["candidates"][0]["residues"] = "W"


def fragment_cites_an_absent_ion(manifest, tables):
    _hit(tables)["matched"][0]["ion_uid"] = "0000000000000bad"


def coelution_over_one_scan(manifest, tables):
    _hit(tables)["coelution"]["consecutive_scans"] = 1


def unmutated(manifest, tables):
    """The control. The shipped examples must validate clean."""


def a_sample_with_no_id(manifest, tables):
    """Must be reported, not raise before any finding is printed."""
    manifest["samples"][0].pop("sample_id", None)


def isotopes_repeat_k_with_one_missing(manifest, tables):
    """A repeated k plus an entry with no k used to raise TypeError in sorted()."""
    for ion in tables["ions"]:
        iso = ion.get("isotopes")
        if isinstance(iso, list) and len(iso) >= 2:
            iso[1] = dict(iso[1])
            iso[1]["k"] = iso[0].get("k")          # duplicate
            iso.append({"rel": 0.1})               # no k at all
            return


CASES = [
    (a_sample_with_no_id, "has no sample_id", 1),
    (isotopes_repeat_k_with_one_missing, "isotopes[] repeats k", 1),
    (quad_window_carries_mobility, "quad_window but one_over_k0 is set", 1),
    (bare_proximity_claims_decoy_control, "bare_proximity with decoy_controlled=true", 1),
    (validated_tier_disclaims_decoy_control, "atlas_consensus with decoy_controlled=false", 1),
    (monoisotope_off_coord, "k=0 is at m/z", 1),
    (envelope_on_a_different_basis, "rel peaks at", 1),
    (reused_ion_uid, "already used on line", 1),
    (ledger_out_of_order, "not ordered oldest first", 1),
    (empty_ledger, "non-empty", 1),
    (ion_names_an_absent_sample, "is not in manifest.samples", 1),
    (drift_kind_disagrees_with_sample, "but sample", 1),
    (ar1_without_a_basis, "ar1 but ar1_basis is unset", 1),
    (ar1_normalized_within_spectrum, "largest ar1 is exactly 1.0", 0),
    (spectrum_for_no_ion, "is not an ion in", 1),
    (two_traces_for_one_isotope, "a second trace for ion", 1),
    (apex_outside_its_own_trace, "falls outside the trace's own range", 1),
    (top3avg_padded, "padding to 3 is forbidden", 1),
    (peptide_count_disagrees_with_list, "peptides_used lists", 1),
    (declared_spectrum_kind_absent, "tables.spectrum_kinds declares", 1),
    (protein_twice_in_one_sample, "already appears on line", 1),
    (impossible_decimation, "exceeds decimation.of", 1),
    (table_present_but_undeclared, "does not declare it", 0),
    (control_clears_the_gate, "reversed sequence matched 3 fragments", 1),
    (entrapment_clears_the_gate, "entrapment peptide matched 5 fragments", 1),
    (gate_demands_every_fragment, "Never n of n", 1),
    (matched_count_disagrees_with_list, "k_of_n.matched is 3", 1),
    (label_both_matched_and_missed, "listed as both matched and missed", 1),
    (evidence_disagrees_with_precursor, "no precursor was found but evidence is", 1),
    (delta_mass_not_derivable, "delta_mass.value is 42.00000", 1),
    (candidate_needs_a_residue_the_peptide_lacks, "carries none of them", 1),
    (fragment_cites_an_absent_ion, "is not an ion in ions.jsonl", 1),
    (coelution_over_one_scan, "A single scan is a blip", 1),
    (unmutated, "", 0),
]


def residue_parsing_separates_positions_from_residues():
    """"N-term" is a position, not residue N.

    Scanned character by character it became residue N, so a legitimate
    N-terminal modification read as impossible on any peptide lacking one.
    Checked directly rather than through a mutation, because changing a
    fixture's sequence also breaks its delta-mass arithmetic and the resulting
    error hides the one under test.
    """
    sys.path.insert(0, str(ROOT / "tools"))
    import validate_export as ve
    bad = []
    for residues, want_sites, want_pos in (
            ("STY", ["S", "T", "Y"], []),
            ("S", ["S"], []),
            ("N-term", [], ["N-term"]),
            ("C-term", [], ["C-term"]),
            ("protein N-term", [], ["protein", "N-term"]),
            ("K, N-term", ["K"], ["N-term"]),
            ("", [], []),
            (None, [], []),
    ):
        sites, pos = ve.required_residues(residues)
        if sorted(sites) != sorted(want_sites):
            bad.append("required_residues(%r) sites %r, expected %r"
                       % (residues, sites, want_sites))
        if bool(pos) != bool(want_pos):
            bad.append("required_residues(%r) positional %r, expected %r"
                       % (residues, pos, want_pos))
    return bad


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-v", "--verbose", action="store_true", help="print the findings")
    args = ap.parse_args()

    failed = 0
    unit = residue_parsing_separates_positions_from_residues()
    failed += bool(unit)
    print("%s residue_parsing_separates_positions_from_residues"
          % ("ok  " if not unit else "FAIL"))
    for line in unit:
        print("     %s" % line)
    for mutate, expect, want_rc in CASES:
        with tempfile.TemporaryDirectory() as tmp:
            build(Path(tmp), mutate)
            rc, out = run(Path(tmp))
        named = expect in out if expect else True
        ok = named and rc == want_rc
        failed += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {mutate.__name__}")
        if not ok:
            if not named:
                print(f"     expected a finding containing {expect!r}")
            if rc != want_rc:
                print(f"     exit {rc}, expected {want_rc}")
        if args.verbose or not ok:
            for line in out.splitlines():
                print(f"     {line}")

    total = len(CASES) + 1
    print(f"\n{total - failed}/{total} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
