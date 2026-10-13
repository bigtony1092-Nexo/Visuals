#!/usr/bin/env python3
"""Check a method against an export manifest and report, panel by panel, what will
actually render.

The point is to find out *before* a run that a panel will come back empty, and to say
why in terms of the data rather than in terms of a stack trace. Run it with a method
alone to validate the method; add a manifest to get the feasibility report.

    python3 tools/plan_panels.py examples/method.example.json
    python3 tools/plan_panels.py examples/method.example.json examples/manifest.example.json
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TIER_ORDER = {"bare_proximity": 0, "pass2_confirmed": 1, "atlas_consensus": 2}

RENDER, PARTIAL, BLANK = "RENDER", "PARTIAL", "BLANK"


def load(path):
    with open(path) as fh:
        return json.load(fh)


def validate(method):
    """Schema-validate the method. Returns a list of human-readable errors."""
    try:
        import jsonschema
    except ImportError:
        return ["jsonschema not installed — skipped schema validation"]

    schema = load(ROOT / "schema" / "method.schema.json")
    validator = jsonschema.Draft202012Validator(schema)
    errors = []
    for err in sorted(validator.iter_errors(method), key=lambda e: list(e.path)):
        where = "/".join(str(p) for p in err.path) or "(root)"
        errors.append(f"{where}: {err.message}")
    return errors


def assess(panel, method, manifest):
    """Decide whether one panel can render. Returns (status, [reasons])."""
    req = panel.get("requires", {})
    reasons = []
    status = RENDER

    samples = manifest["samples"] if manifest else []

    # --- drift reality -----------------------------------------------------
    want_drift = req.get("drift_kind")
    if want_drift and samples:
        ok = [s for s in samples if s["drift_kind"] in want_drift]
        bad = [s for s in samples if s["drift_kind"] not in want_drift]
        if not ok:
            status = BLANK
            reasons.append(
                "no sample has real ion mobility ("
                + ", ".join(f"{s['sample_id']}={s['drift_kind']}" for s in bad)
                + ")"
            )
        elif bad:
            status = PARTIAL if status == RENDER else status
            for s in bad:
                why = (
                    "SONAR drift is a quadrupole window position, not mobility"
                    if s["drift_kind"] == "quad_window"
                    else "no drift column in the raw file"
                )
                reasons.append(f"excludes {s['sample_id']} ({s['acquisition_method']}): {why}")

    # --- confidence tier ---------------------------------------------------
    need_tier = req.get("min_tier")
    if need_tier:
        policy = (method.get("tier_policy") or {}).get("min_tier", "pass2_confirmed")
        if TIER_ORDER[policy] < TIER_ORDER[need_tier]:
            mix = (method.get("tier_policy") or {}).get("mix_tiers", False)
            if mix:
                status = PARTIAL if status != BLANK else status
                reasons.append(
                    f"method admits {policy} but panel needs {need_tier}; mix_tiers=true, "
                    "so results blend a decoy-controlled tier with an untested one"
                )
            else:
                status = BLANK
                reasons.append(
                    f"method admits {policy}, panel needs {need_tier}, mix_tiers=false"
                )
        if manifest:
            by_tier = manifest.get("counts", {}).get("by_tier", {})
            have = sum(v for k, v in by_tier.items() if TIER_ORDER[k] >= TIER_ORDER[need_tier])
            if by_tier and have == 0:
                status = BLANK
                reasons.append(f"export contains 0 ions at {need_tier} or above")

    # --- sample count ------------------------------------------------------
    need_n = req.get("min_samples")
    if need_n and samples and len(samples) < need_n:
        status = BLANK
        reasons.append(f"needs {need_n} samples, export has {len(samples)}")

    # --- stoichiometry pairs ----------------------------------------------
    if req.get("pairs"):
        pairs = ((method.get("stoichiometry") or {}).get("pairs")) or []
        if not pairs:
            status = BLANK
            reasons.append("no stoichiometry.pairs declared in the method")

    # --- classifier --------------------------------------------------------
    if req.get("labels"):
        cls = method.get("classification") or {}
        if not cls.get("enabled"):
            status = BLANK
            reasons.append("classification.enabled is false")
        elif samples:
            labelled = [s for s in samples if s.get("sample_class")]
            classes = {s["sample_class"] for s in labelled}
            if len(classes) < req.get("min_classes", 2):
                status = BLANK
                reasons.append(
                    f"needs {req.get('min_classes', 2)} labelled classes, found {len(classes)}"
                )
            else:
                need_b = req.get("min_batches_per_class", 2)
                thin = {}
                for c in classes:
                    batches = {s.get("batch") for s in labelled if s["sample_class"] == c}
                    batches.discard(None)
                    if len(batches) < need_b:
                        thin[c] = len(batches)
                if thin:
                    status = BLANK
                    detail = ", ".join(f"{c}: {n} batch(es)" for c, n in sorted(thin.items()))
                    reasons.append(
                        f"batch holdout needs {need_b} batches per class — {detail}. "
                        "Any accuracy here would measure the instrument, not the cell line"
                    )
                need_s = cls.get("min_samples_per_class", 3)
                small = {
                    c: sum(1 for s in labelled if s["sample_class"] == c)
                    for c in classes
                }
                short = {c: n for c, n in small.items() if n < need_s}
                if short and status != BLANK:
                    status = PARTIAL
                    reasons.append(
                        "under min_samples_per_class="
                        f"{need_s}: " + ", ".join(f"{c}={n}" for c, n in sorted(short.items()))
                    )

    return status, reasons


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("method")
    ap.add_argument("manifest", nargs="?")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any selected panel would be blank")
    args = ap.parse_args()

    method = load(args.method)
    manifest = load(args.manifest) if args.manifest else None
    catalog = {p["id"]: p for p in load(ROOT / "schema" / "panels.json")["panels"]}

    print(f"method   {method.get('name')}  ({method.get('method_id')} v{method.get('version')})")
    if manifest:
        print(f"export   {manifest['export_id']}  ·  {len(manifest['samples'])} samples")
        c = manifest.get("counts", {})
        if c:
            print(f"ions     {c.get('ions_exported'):,} exported of {c.get('ions_total'):,}")
    print()

    errors = validate(method)
    if errors:
        print("method validation")
        for e in errors:
            print(f"  ! {e}")
        print()

    # Quant honesty check, independent of any panel.
    quant = method.get("quant", {})
    if not quant.get("absolute_standard"):
        print("quant    no absolute_standard declared — Top3avg is RELATIVE abundance.")
        print("         Panels must not be labelled as concentration.")
        print()

    print(f"{'panel':<24} {'status':<9} detail")
    print("-" * 78)
    blanks = 0
    for sel in method.get("panels", []):
        pid = sel["id"]
        panel = catalog.get(pid)
        if panel is None:
            print(f"{pid:<24} {'UNKNOWN':<9} not in schema/panels.json")
            blanks += 1
            continue
        status, reasons = assess(panel, method, manifest)
        if status == BLANK:
            blanks += 1
        first = reasons[0] if reasons else ""
        print(f"{pid:<24} {status:<9} {first}")
        for extra in reasons[1:]:
            print(f"{'':<24} {'':<9} {extra}")
    print()

    if manifest and manifest.get("warnings"):
        print("export warnings")
        for w in manifest["warnings"]:
            print(f"  · {w}")
        print()

    dec = manifest.get("decimation") if manifest else None
    if dec and dec.get("strategy") != "none":
        pct = 100.0 * dec["kept"] / dec["of"] if dec["of"] else 0.0
        print(f"decimation  {dec['strategy']}: {dec['kept']:,} of {dec['of']:,} ({pct:.2f}%)")
        print("            the viewer must show this ratio on every affected panel")
        print()

    if args.strict and blanks:
        print(f"{blanks} panel(s) would not render", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
