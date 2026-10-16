#!/usr/bin/env python3
"""Check a real export against the contract and say, field by field, what is wrong.

`plan_panels.py` asks whether a panel *can* render. This asks whether the export is
honest. Schema validation is the easy half and is done here when `jsonschema` is
installed; the half that matters is the rules the schemas cannot express — a SONAR
quadrupole position sitting in a mobility field, an isotope minted as its own ion,
`bare_proximity` claiming to be decoy-controlled, a 1-peptide Top3avg padded to 3,
`ar1` normalized inside its own spectrum until it carries no information.

Each finding is one of:

    ERROR  the export violates the contract; a panel built on it would mislead
    WARN   legal, but lossy or suspicious — read the reason before trusting the panel
    NOTE   worth knowing, no action implied

    python3 tools/validate_export.py export/
    python3 tools/validate_export.py export/manifest.json
    python3 tools/validate_export.py export/ --limit 100000   # spot-check a huge export
    python3 tools/validate_export.py export/ --strict         # warnings fail too
    python3 tools/validate_export.py examples/manifest.example.json --excerpt

Exit status is 1 if any ERROR was found, or any WARN under --strict.

Records stream one line at a time, but the referential checks hold one `ion_uid` per
exported ion in memory. For an export in the hundreds of millions, use --limit.
"""

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

ERROR, WARN, NOTE = "ERROR", "WARN", "NOTE"

# Tiers that the contract defines as decoy-controlled (section 4, confidence tier).
DECOY_CONTROLLED = {"bare_proximity": False, "pass2_confirmed": True, "atlas_consensus": True}

# The m/z grid the ion_uid hash rounds to (section 2, the forensic track), and so the finest
# distinction two m/z values in this export are allowed to claim.
MZ_TOL = 1e-4

TABLES = {
    "ions": "ion_track",
    "proteins": "protein_quant",
    "spectra": "spectrum",
    "traces": "trace",
    "fragment_hits": "fragment_hit",
}

# Monoisotopic residue masses, for recomputing what a record claims about mass
# (section 7, fragment-first). Carbamidomethyl is fixed on Cys, so the unmodified mass of a
# peptide is computable and a stated delta can be checked rather than trusted.
RESIDUE = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
WATER = 18.010565
CARBAMIDOMETHYL = 57.02146
DELTA_TOL = 0.01


def peptide_mass(sequence):
    """Unmodified neutral mass, carbamidomethyl fixed on every Cys."""
    if not sequence or any(c not in RESIDUE for c in sequence):
        return None
    return (WATER + sum(RESIDUE[c] for c in sequence)
            + CARBAMIDOMETHYL * sequence.count("C"))


class Report:
    """Collected findings, grouped by the file they were found in."""

    def __init__(self):
        self.findings = []
        self.counts = Counter()

    def add(self, level, where, message):
        self.findings.append((level, where, message))
        self.counts[level] += 1

    def error(self, where, message):
        self.add(ERROR, where, message)

    def warn(self, where, message):
        self.add(WARN, where, message)

    def note(self, where, message):
        self.add(NOTE, where, message)

    def print(self):
        if not self.findings:
            print("no findings — the export satisfies every rule this tool checks")
            return
        order = {ERROR: 0, WARN: 1, NOTE: 2}
        by_file = defaultdict(list)
        for level, where, message in self.findings:
            by_file[where].append((order[level], level, message))
        for where in sorted(by_file):
            print(where)
            for _, level, message in sorted(by_file[where]):
                lead = f"  {level:<6}"
                lines = message.split("\n")
                print(f"{lead} {lines[0]}")
                for extra in lines[1:]:
                    print(f"{'':<8} {extra}")
            print()


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def read_jsonl(path, limit=None):
    """Yield (line_number, record). A line that will not parse is reported, not fatal."""
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if limit is not None and n > limit:
                return
            line = line.strip()
            if not line:
                continue
            try:
                yield n, json.loads(line)
            except json.JSONDecodeError as exc:
                yield n, exc


def schema_validators(report):
    """One Draft 2020-12 validator per table, or None if jsonschema is absent."""
    try:
        import jsonschema
    except ImportError:
        report.note(
            "(setup)",
            "jsonschema not installed — schema validation skipped, contract rules still ran.\n"
            "pip install jsonschema  for the full check",
        )
        return None
    out = {}
    for key, schema_name in TABLES.items():
        schema = load_json(ROOT / "schema" / f"{schema_name}.schema.json")
        out[key] = jsonschema.Draft202012Validator(schema)
    out["manifest"] = jsonschema.Draft202012Validator(
        load_json(ROOT / "schema" / "manifest.schema.json")
    )
    return out


def check_schema(report, validator, where, record, line=None):
    """Report every schema error on one record. Returns True if it validated."""
    if validator is None:
        return True
    errors = sorted(validator.iter_errors(record), key=lambda e: list(e.path))
    at = f"line {line}: " if line else ""
    for err in errors:
        path = "/".join(str(p) for p in err.path) or "(root)"
        report.error(where, f"{at}{path}: {err.message}")
    return not errors


def parse_instant(value):
    """ISO-8601 to datetime, or None if it will not parse."""
    if not isinstance(value, str):
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None


# --------------------------------------------------------------------------- manifest


def check_manifest(report, manifest, export_dir, validators, names, where):
    check_schema(report, (validators or {}).get("manifest"), where, manifest)

    samples = {s["sample_id"]: s for s in manifest.get("samples", [])}
    if not samples:
        report.error(where, "samples[] is empty — no ion can name a sample that exists")

    tables = manifest.get("tables") or {}
    if not tables:
        report.warn(
            where,
            "no tables{} block. The viewer reads the manifest first to decide which panels\n"
            "can render; without it every optional panel has to guess or fail on a missing file",
        )

    # A declared table must be on disk, and a table on disk should be declared.
    declared = {}
    for key in TABLES:
        filename = names[key]
        on_disk = (export_dir / filename).exists()
        # `ions` is implicit: the contract requires it whether or not it is declared.
        is_declared = bool(tables.get(key)) if key != "ions" else tables.get("ions", True)
        declared[key] = bool(is_declared) and on_disk
        if is_declared and not on_disk:
            report.error(where, f"tables.{key} is declared but {filename} is not in the export")
        if on_disk and not is_declared:
            report.warn(
                where,
                f"{filename} is present but tables.{key} does not declare it — panels that\n"
                f"need {key} will blank although the data is right there",
            )
    if not (export_dir / names["ions"]).exists():
        report.error(where, f"{names['ions']} is missing; it is the one table the contract requires")

    # Decimation has to be disclosed and has to be arithmetically possible (section 11,
    # decimation is disclosed).
    dec = manifest.get("decimation")
    if dec and dec.get("strategy") != "none":
        kept, of = dec.get("kept"), dec.get("of")
        if isinstance(kept, int) and isinstance(of, int):
            if kept > of:
                report.error(where, f"decimation.kept ({kept:,}) exceeds decimation.of ({of:,})")
            elif of:
                pct = 100.0 * kept / of
                report.note(
                    where,
                    f"decimated: {dec['strategy']}, {kept:,} of {of:,} ({pct:.2f}%) — every"
                    " affected panel must show this ratio",
                )
        if dec.get("seed") is None and dec.get("strategy") == "uniform_random":
            report.warn(where, "uniform_random decimation without a seed cannot be reproduced")

    # Top3avg is relative unless a spike-in is declared (section 5, Top3avg).
    quant = manifest.get("quant") or {}
    if quant and not quant.get("absolute_standard_declared"):
        report.note(
            where,
            "no absolute standard declared — Top3avg is RELATIVE abundance and no panel"
            " may be labelled as concentration",
        )
    return samples, declared


# ------------------------------------------------------------------------------- ions


def check_ion(report, where, line, ion, samples, seen_uids):
    uid = ion.get("ion_uid")
    sample_id = ion.get("sample_id")
    at = f"line {line}"

    # ion_uid is assigned once and never reassigned (section 2, the forensic track).
    if uid in seen_uids:
        report.error(
            where,
            f"{at}: ion_uid {uid} already used on line {seen_uids[uid]}. The uid is"
            " content-addressed, so a repeat means two different ions hashed the same"
            " or one ion was exported twice",
        )
    elif uid is not None:
        seen_uids[uid] = line

    sample = samples.get(sample_id) if isinstance(sample_id, str) else None
    if sample_id is not None and sample is None:
        report.error(where, f"{at}: sample_id {sample_id!r} is not in manifest.samples")

    coord = ion.get("coord")
    coord = coord if isinstance(coord, dict) else {}
    kind = coord.get("drift_kind")
    drift = coord.get("drift")
    k0 = coord.get("one_over_k0")

    # Section 3. The whole reason drift_kind is a required field.
    if sample is not None and kind is not None and sample.get("drift_kind") != kind:
        report.error(
            where,
            f"{at}: coord.drift_kind is {kind!r} but sample {sample_id} declares"
            f" {sample.get('drift_kind')!r}",
        )
    if kind == "quad_window" and k0 is not None:
        report.error(
            where,
            f"{at}: drift_kind is quad_window but one_over_k0 is set ({k0}). A SONAR"
            " quadrupole window position is not a mobility; on a mobility axis it draws a"
            " clean, convincing, meaningless chart",
        )
    if kind == "quad_window" and coord.get("ccs") is not None:
        report.error(
            where,
            f"{at}: drift_kind is quad_window but a CCS is reported. There is no"
            " collision cross-section to derive from a quadrupole position",
        )
    if kind == "absent" and drift is not None:
        report.error(where, f"{at}: drift_kind is absent but drift is {drift}, not null")
    if kind == "ion_mobility" and drift is not None and k0 is None:
        report.warn(
            where,
            f"{at}: real mobility present but one_over_k0 is null — the mobility panels"
            " plot 1/K0, so this ion will be missing from them",
        )

    # Section 4. The tier decides whether decoys were controlled; an ion cannot claim
    # otherwise.
    conf = ion.get("confidence")
    conf = conf if isinstance(conf, dict) else {}
    tier = conf.get("tier")
    if tier in DECOY_CONTROLLED:
        want = DECOY_CONTROLLED[tier]
        got = conf.get("decoy_controlled")
        if got is not None and got != want:
            if tier == "bare_proximity":
                report.error(
                    where,
                    f"{at}: tier bare_proximity with decoy_controlled=true. That tier is a"
                    " ppm window plus consecutive scans with no MS2; the same mechanism"
                    " passed reversed-tryptic decoys at 33-70% in this project",
                )
            else:
                report.error(
                    where,
                    f"{at}: tier {tier} with decoy_controlled=false — the tier is defined"
                    " by being decoy-controlled",
                )

    # Section 2. Isotopes ride on the precursor; k=0 is the monoisotopic peak at coord.mz
    # and the most intense peak is scaled to 100.
    isotopes = ion.get("isotopes")
    if isinstance(isotopes, list) and isotopes:
        isotopes = [iso for iso in isotopes if isinstance(iso, dict)]
        ks = [iso.get("k") for iso in isotopes]
        if len(set(ks)) != len(ks):
            dupes = sorted({k for k in ks if ks.count(k) > 1})
            report.error(where, f"{at}: isotopes[] repeats k={dupes}")
        mono = next((iso for iso in isotopes if iso.get("k") == 0), None)
        if mono is None:
            report.warn(
                where,
                f"{at}: isotopes[] has no k=0 entry, so the envelope has no monoisotopic"
                " anchor at coord.mz",
            )
        elif isinstance(coord.get("mz"), (int, float)) and isinstance(mono.get("mz"), (int, float)):
            if abs(mono["mz"] - coord["mz"]) > MZ_TOL:
                report.error(
                    where,
                    f"{at}: isotope k=0 is at m/z {mono['mz']} but coord.mz is"
                    f" {coord['mz']}. k=0 is the monoisotopic peak and sits at coord.mz",
                )
        rels = [iso.get("rel") for iso in isotopes if isinstance(iso.get("rel"), (int, float))]
        if rels and abs(max(rels) - 100.0) > 0.5:
            report.error(
                where,
                f"{at}: isotopes[].rel peaks at {max(rels):g}; the most intense isotope is"
                " scaled to 100, so a different maximum means a different basis",
            )
        if any(iso.get("observed") is False for iso in isotopes):
            n = sum(1 for iso in isotopes if iso.get("observed") is False)
            report.note(
                where,
                f"{at}: {n} isotope peak(s) are modelled, not observed — a panel drawing"
                " them must say so",
            )

    # Section 11. expected_basis=averagine on an ion whose sequence is known throws away
    # the exact distribution that identity buys.
    fit = ion.get("isotope_fit")
    if isinstance(fit, dict) and fit.get("expected_basis") == "averagine":
        annotation = ion.get("annotation")
        if isinstance(annotation, dict) and annotation.get("peptide_sequence"):
            report.warn(
                where,
                f"{at}: isotope_fit uses averagine although the ion is annotated"
                f" ({annotation['peptide_sequence']}). With the sequence known the"
                " composition gives the exact distribution; averagine is for unidentified ions",
            )

    # Section 11. null and flag:none are different claims.
    inter = ion.get("interference")
    if isinstance(inter, dict):
        start, end = inter.get("rt_start"), inter.get("rt_end")
        if isinstance(start, (int, float)) and isinstance(end, (int, float)) and start > end:
            report.error(where, f"{at}: interference.rt_start {start} is after rt_end {end}")
        if inter.get("flag") == "none" and (start is not None or end is not None):
            report.warn(
                where,
                f"{at}: interference.flag is none but an rt window is given — if an"
                " interferer was found the flag should name how",
            )

    # Section 2. events[] is the audit trail and the animation; order carries meaning.
    events = ion.get("events")
    events = [e for e in events if isinstance(e, dict)] if isinstance(events, list) else []
    if not events:
        report.error(
            where,
            f"{at}: events[] is empty. It is both the audit trail and what the flight"
            " animation replays",
        )
    else:
        if events[0].get("stage") != "detected":
            report.warn(
                where,
                f"{at}: events[] starts at {events[0].get('stage')!r}, not 'detected' —"
                " the ion's first appearance is missing from the ledger",
            )
        instants = [parse_instant(e.get("at")) for e in events]
        known = [t for t in instants if t is not None]
        if len(known) > 1 and any(b < a for a, b in zip(known, known[1:])):
            report.error(
                where,
                f"{at}: events[] is not ordered oldest first. A correction is a new event"
                " appended, never an edit, so the order is the history",
            )
        for i, e in enumerate(events):
            if e.get("stage") == "depleted" and e.get("explained_by") is None:
                report.warn(
                    where,
                    f"{at}: events[{i}] is stage=depleted with no explained_by — depletion"
                    " means some other ion accounted for it, and that is the whole record of"
                    " which",
                )
    return tier


def check_ions(report, path, samples, validators, limit):
    where = path.name
    seen_uids = {}
    tiers = Counter()
    annotated_peptides = set()
    n = 0
    for line, ion in read_jsonl(path, limit):
        if isinstance(ion, json.JSONDecodeError):
            report.error(where, f"line {line}: not valid JSON — {ion.msg}")
            continue
        n += 1
        if not isinstance(ion, dict):
            report.error(where, f"line {line}: not a JSON object")
            continue
        # Schema errors are reported but do not stop the contract checks, and above all
        # do not stop the uid being registered: an unregistered uid turns one bad ion
        # into an orphan error on every spectrum and trace that refers to it.
        check_schema(report, (validators or {}).get("ions"), where, ion, line)
        tier = check_ion(report, where, line, ion, samples, seen_uids)
        if tier:
            tiers[tier] += 1
        pep = (ion.get("annotation") or {}).get("peptide_sequence")
        if isinstance(pep, str) and pep:
            annotated_peptides.add(pep)
    return n, seen_uids, tiers, annotated_peptides


# ------------------------------------------------------------------ spectra and traces


def check_spectra(report, path, uids, samples, validators, limit):
    where = path.name
    kinds = Counter()
    for line, spec in read_jsonl(path, limit):
        if isinstance(spec, json.JSONDecodeError):
            report.error(where, f"line {line}: not valid JSON — {spec.msg}")
            continue
        if not check_schema(report, (validators or {}).get("spectra"), where, spec, line):
            continue
        at = f"line {line}"
        uid = spec.get("ion_uid")
        if uids and uid not in uids:
            report.error(where, f"{at}: ion_uid {uid} is not an ion in ions.jsonl")
        if spec.get("sample_id") and spec["sample_id"] not in samples:
            report.error(where, f"{at}: sample_id {spec['sample_id']!r} is not in manifest.samples")
        kinds[spec.get("kind")] += 1

        peaks = spec.get("peaks") or []
        ar1s = [p.get("ar1") for p in peaks if p.get("ar1") is not None]

        # Section 6. ar1 is a ratio to the precursor, and which precursor number was used
        # changes the answer.
        if ar1s and not spec.get("ar1_basis"):
            report.error(
                where,
                f"{at}: peaks carry ar1 but ar1_basis is unset."
                " summed_precursor_intensity and matched_intensity are different"
                " denominators and give different numbers",
            )
        if ar1s and len(ar1s) > 1 and abs(max(ar1s) - 1.0) < 1e-9:
            report.warn(
                where,
                f"{at}: the largest ar1 is exactly 1.0, the signature of ar1 having been"
                " normalized inside the spectrum. Dividing by the largest cancels the"
                " precursor exactly and leaves base-peak relative intensity, which"
                " peaks[].intensity already holds",
            )
        if spec.get("ar1_sum") is not None and ar1s:
            total = sum(ar1s)
            if total and abs(spec["ar1_sum"] - total) / max(total, 1e-12) > 0.01:
                report.warn(
                    where,
                    f"{at}: ar1_sum is {spec['ar1_sum']:.4g} but the assigned peaks sum to"
                    f" {total:.4g}. ar1_sum is the share of the precursor that converted to"
                    " identified fragments, so it should be that sum",
                )
        if spec.get("kind") == "ms2_library" and ar1s:
            report.warn(
                where,
                f"{at}: a library spectrum carries ar1. Library sources record no precursor"
                " intensity, so their ar1 is null by construction (section 6, AR1)",
            )
    return kinds


def check_traces(report, path, uids, validators, limit):
    where = path.name
    seen = {}
    for line, tr in read_jsonl(path, limit):
        if isinstance(tr, json.JSONDecodeError):
            report.error(where, f"line {line}: not valid JSON — {tr.msg}")
            continue
        if not check_schema(report, (validators or {}).get("traces"), where, tr, line):
            continue
        at = f"line {line}"
        uid = tr.get("ion_uid")
        if uids and uid not in uids:
            report.error(where, f"{at}: ion_uid {uid} is not an ion in ions.jsonl")
        key = (uid, tr.get("isotope_k"))
        if key in seen:
            report.error(
                where,
                f"{at}: a second trace for ion {uid} isotope_k={tr.get('isotope_k')}"
                f" (first on line {seen[key]}). One trace per isotope per ion",
            )
        else:
            seen[key] = line
        points = tr.get("points") or []
        if not points:
            report.warn(where, f"{at}: points[] is empty — nothing to draw")
        apex = tr.get("apex_rt_min")
        if apex is not None and points:
            rts = [p[0] if isinstance(p, list) else p.get("rt_min") for p in points]
            rts = [r for r in rts if isinstance(r, (int, float))]
            if rts and not (min(rts) - 1e-9 <= apex <= max(rts) + 1e-9):
                report.error(
                    where,
                    f"{at}: apex_rt_min {apex} falls outside the trace's own range"
                    f" {min(rts)}-{max(rts)}",
                )


def check_fragment_hits(report, path, uids, samples, validators, limit):
    """Section 7. A fragment-first hit is a claim built from fragments whose
    masses are shared across peptides, so almost every rule here is about
    whether the claim can be told apart from one that cannot be true."""
    where = path.name
    for line, hit in read_jsonl(path, limit):
        if isinstance(hit, json.JSONDecodeError):
            report.error(where, f"line {line}: not valid JSON — {hit.msg}")
            continue
        if not check_schema(report, (validators or {}).get("fragment_hits"), where, hit, line):
            continue
        at = f"line {line}"

        sample = hit.get("sample_id")
        if samples and sample not in samples:
            report.error(where, f"{at}: sample_id {sample!r} is not a sample in the manifest")

        matched = hit.get("matched") or []
        missed = hit.get("missed") or []
        kn = hit.get("k_of_n") or {}
        need = kn.get("required")
        got = kn.get("matched")
        queried = kn.get("queried")

        # --- the control gate: a bar a decoy clears is not a bar -----------------
        control = hit.get("control") or {}
        for field, what in (("reversed_matched", "the reversed sequence"),
                            ("entrapment_matched", "an entrapment peptide")):
            n = control.get(field)
            if isinstance(n, int) and isinstance(need, int) and n >= need:
                report.error(
                    where,
                    f"{at}: {what} matched {n} fragments, reaching the same bar of {need}."
                    " A control that clears the gate means the gate did not discriminate;"
                    " this hit does not stand",
                )

        # --- the gate itself -----------------------------------------------------
        if isinstance(need, int) and isinstance(queried, int) and need >= queried:
            report.error(
                where,
                f"{at}: k_of_n requires {need} of {queried}. Never n of n — demanding every"
                " queried fragment fails exactly when the modification is real, because"
                " some queried fragments carry the modified residue",
            )
        if isinstance(need, int) and isinstance(got, int) and got < need:
            report.error(where, f"{at}: {got} fragments matched but {need} were required —"
                                " this is not a hit")

        # --- internal arithmetic -------------------------------------------------
        if isinstance(got, int) and got != len(matched):
            report.error(where, f"{at}: k_of_n.matched is {got} but matched[] holds"
                                f" {len(matched)} fragment(s)")
        if isinstance(queried, int) and queried != len(matched) + len(missed):
            report.error(where, f"{at}: k_of_n.queried is {queried} but matched[] and missed[]"
                                f" hold {len(matched)} + {len(missed)} together")
        n_declared = (hit.get("queried") or {}).get("n_fragments")
        if isinstance(n_declared, int) and isinstance(queried, int) and n_declared != queried:
            report.error(where, f"{at}: queried.n_fragments is {n_declared} but k_of_n.queried"
                                f" is {queried}")
        both = sorted({m.get("label") for m in matched} & set(missed))
        if both:
            report.error(where, f"{at}: {', '.join(both)} listed as both matched and missed")

        # --- fragments must be detections the ledger holds ----------------------
        for m in matched:
            uid = m.get("ion_uid")
            if uid and uids and uid not in uids:
                report.error(where, f"{at}: fragment {m.get('label')} cites ion_uid {uid},"
                                    " which is not an ion in ions.jsonl")

        co = hit.get("coelution") or {}
        scans = co.get("consecutive_scans")
        if isinstance(scans, int) and scans < 2:
            report.error(where, f"{at}: coelution over {scans} scan(s). A single scan is a blip,"
                                " not a peak")

        # --- evidence must match what was actually found -------------------------
        prec = hit.get("precursor") or {}
        found = bool(prec.get("found"))
        evidence = hit.get("evidence")
        if found and evidence != "fragment_and_precursor":
            report.error(where, f"{at}: a precursor was found but evidence is {evidence!r}")
        if not found and evidence != "fragment_only":
            report.error(where, f"{at}: no precursor was found but evidence is {evidence!r}")

        delta = hit.get("delta_mass")
        if delta is not None and not found:
            report.error(where, f"{at}: delta_mass is given with no precursor. The delta is"
                                " observed mass minus unmodified; without a precursor there is"
                                " nothing to subtract from")

        # --- the delta is derivable, so it is checkable --------------------------
        sequence = (hit.get("peptide") or {}).get("sequence")
        base = peptide_mass(sequence)
        if delta is not None and found and base is not None:
            observed = prec.get("neutral_mass")
            stated = delta.get("value")
            if isinstance(observed, (int, float)) and isinstance(stated, (int, float)):
                expect = observed - base
                if abs(expect - stated) > DELTA_TOL:
                    report.error(
                        where,
                        f"{at}: delta_mass.value is {stated:.5f} but the precursor's neutral"
                        f" mass {observed:.4f} less the unmodified peptide {base:.4f} is"
                        f" {expect:.5f}",
                    )
            for cand in (delta.get("candidates") or []):
                residues = cand.get("residues") or ""
                sites = [r for r in residues if r in RESIDUE]
                if sites and sequence and cand.get("plausible") is True:
                    if not any(r in sequence for r in sites):
                        report.error(
                            where,
                            f"{at}: {cand.get('name')} is marked plausible but needs one of"
                            f" {residues}, and {sequence} carries none of them",
                        )

        # --- localization only means something when something shifted -----------
        if hit.get("site_localization") and not missed:
            report.warn(where, f"{at}: a site is localized although no queried fragment is"
                               " missing. Nothing shifted, so nothing bounds the site")


def check_proteins(report, path, samples, peptides, validators, limit, excerpt):
    where = path.name
    seen = {}
    for line, rec in read_jsonl(path, limit):
        if isinstance(rec, json.JSONDecodeError):
            report.error(where, f"line {line}: not valid JSON — {rec.msg}")
            continue
        if not check_schema(report, (validators or {}).get("proteins"), where, rec, line):
            continue
        at = f"line {line}"
        key = (rec.get("sample_id"), rec.get("prot_accession"))
        if key in seen:
            report.error(
                where,
                f"{at}: {key[1]} in {key[0]} already appears on line {seen[key]}."
                " One record per (protein, sample)",
            )
        else:
            seen[key] = line
        if rec.get("sample_id") not in samples:
            report.error(where, f"{at}: sample_id {rec.get('sample_id')!r} is not in manifest.samples")

        # Section 5. Padding to 3 is forbidden, and a Top3avg cannot average more
        # peptides than it matched.
        used = rec.get("n_peptides_used")
        matched = rec.get("n_peptides_matched")
        if isinstance(used, int) and not 1 <= used <= 3:
            report.error(where, f"{at}: n_peptides_used is {used}; Top3avg averages 1, 2 or 3")
        if isinstance(used, int) and isinstance(matched, int):
            if used > matched:
                report.error(
                    where,
                    f"{at}: n_peptides_used {used} exceeds n_peptides_matched {matched} —"
                    " padding to 3 is forbidden; a 1-peptide Top3avg is a different"
                    " measurement, not a worse one",
                )
            elif used != min(3, matched):
                report.warn(
                    where,
                    f"{at}: {matched} peptides matched but {used} averaged; Top3avg takes"
                    f" the top {min(3, matched)}",
                )
        listed = rec.get("peptides_used")
        if listed is not None and isinstance(used, int) and len(listed) != used:
            report.error(
                where,
                f"{at}: peptides_used lists {len(listed)} peptides but n_peptides_used is {used}",
            )
        if listed and peptides and not excerpt:
            orphans = [p for p in listed if p not in peptides]
            if orphans:
                report.warn(
                    where,
                    f"{at}: peptides_used names {len(orphans)} peptide(s) with no annotated"
                    f" ion in this export, e.g. {orphans[0]}",
                )
        tier = rec.get("tier")
        if tier in DECOY_CONTROLLED and rec.get("decoy_controlled") is not None:
            if rec["decoy_controlled"] != DECOY_CONTROLLED[tier]:
                report.error(
                    where,
                    f"{at}: tier {tier} with decoy_controlled={rec['decoy_controlled']}",
                )


# ------------------------------------------------------------------------ cross-checks


def check_totals(report, manifest, n_ions, tiers, kinds, limited, excerpt, where):
    counts = manifest.get("counts") or {}
    exported = counts.get("ions_exported")
    total = counts.get("ions_total")

    if limited or excerpt:
        why = "--limit was given" if limited else "--excerpt was given"
        report.note(where, f"{why}, so record counts were not compared with counts{{}}")
        limited = True
    elif isinstance(exported, int) and exported != n_ions:
        report.error(
            where,
            f"counts.ions_exported is {exported:,} but ions.jsonl holds {n_ions:,} records",
        )
    if isinstance(total, int) and isinstance(exported, int) and exported > total:
        report.error(where, f"counts.ions_exported ({exported:,}) exceeds ions_total ({total:,})")

    by_tier = counts.get("by_tier") or {}
    if by_tier and not limited:
        for tier, declared in by_tier.items():
            actual = tiers.get(tier, 0)
            if declared != actual:
                report.error(
                    where,
                    f"counts.by_tier[{tier}] is {declared:,} but {actual:,} ions carry that tier",
                )
        undeclared = {t: c for t, c in tiers.items() if t not in by_tier}
        if undeclared:
            detail = ", ".join(f"{t}={c:,}" for t, c in sorted(undeclared.items()))
            report.warn(where, f"ions carry tiers absent from counts.by_tier: {detail}")

    dec = manifest.get("decimation")
    if dec and dec.get("strategy") != "none" and isinstance(dec.get("kept"), int):
        if isinstance(exported, int) and dec["kept"] != exported:
            report.warn(
                where,
                f"decimation.kept is {dec['kept']:,} but counts.ions_exported is"
                f" {exported:,}; the ratio the viewer must display comes from these two",
            )

    tables = manifest.get("tables") or {}
    declared_kinds = set(tables.get("spectrum_kinds") or [])
    present = {k for k in kinds if k}
    if declared_kinds or present:
        missing = sorted(declared_kinds - present)
        extra = sorted(present - declared_kinds)
        if missing:
            report.error(
                where,
                "tables.spectrum_kinds declares " + ", ".join(missing) +
                " but no such spectrum is in spectra.jsonl. A panel gated on it will be"
                " selected and then come back empty",
            )
        if extra and declared_kinds:
            report.warn(
                where,
                "spectra.jsonl holds " + ", ".join(extra) +
                " which tables.spectrum_kinds does not declare",
            )


# ------------------------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("export", help="an export directory, or its manifest.json")
    ap.add_argument("--limit", type=int, metavar="N",
                    help="read only the first N records of each table")
    ap.add_argument("--excerpt", action="store_true",
                    help="the tables are a sample of a larger export, so skip the"
                         " whole-export cross-checks (counts, peptide back-references)")
    ap.add_argument("--strict", action="store_true", help="exit non-zero on warnings too")
    ap.add_argument("--quiet", action="store_true", help="print findings only, no summary")
    args = ap.parse_args()

    target = Path(args.export)
    if target.is_dir():
        export_dir, manifest_path = target, target / "manifest.json"
    else:
        export_dir, manifest_path = target.parent, target
    if not manifest_path.exists():
        print(f"no manifest at {manifest_path}", file=sys.stderr)
        return 2

    # "manifest.json" -> tables are "ions.jsonl"; "manifest.example.json" -> they are
    # "ions.example.jsonl". So the shipped examples validate in place.
    infix = manifest_path.name[len("manifest"):-len(".json")]
    names = {key: f"{key}{infix}.jsonl" for key in TABLES}

    report = Report()
    manifest = load_json(manifest_path)
    validators = schema_validators(report)

    if not args.quiet:
        print(f"export   {manifest.get('export_id')}  ·  {manifest_path}")
        print(f"method   {manifest.get('method_id')} v{manifest.get('method_version')}")
        print(f"samples  {len(manifest.get('samples', []))}")
        print()

    samples, declared = check_manifest(
        report, manifest, export_dir, validators, names, manifest_path.name
    )

    n_ions, uids, tiers, peptides = 0, {}, Counter(), set()
    if (export_dir / names["ions"]).exists():
        n_ions, uids, tiers, peptides = check_ions(
            report, export_dir / names["ions"], samples, validators, args.limit
        )

    kinds = Counter()
    if declared.get("spectra"):
        kinds = check_spectra(
            report, export_dir / names["spectra"], uids, samples, validators, args.limit
        )
    if declared.get("traces"):
        check_traces(report, export_dir / names["traces"], uids, validators, args.limit)
    if declared.get("fragment_hits"):
        check_fragment_hits(
            report, export_dir / names["fragment_hits"], uids, samples, validators, args.limit
        )
    if declared.get("proteins"):
        check_proteins(
            report, export_dir / names["proteins"], samples, peptides, validators,
            args.limit, args.excerpt
        )

    check_totals(report, manifest, n_ions, tiers, kinds,
                 limited=args.limit is not None, excerpt=args.excerpt,
                 where=manifest_path.name)

    report.print()

    if not args.quiet:
        print(f"{n_ions:,} ions read"
              + (f" (limit {args.limit:,})" if args.limit else "")
              + f"  ·  {report.counts[ERROR]} error(s), {report.counts[WARN]} warning(s),"
              f" {report.counts[NOTE]} note(s)")

    if report.counts[ERROR]:
        return 1
    if args.strict and report.counts[WARN]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
