#!/usr/bin/env python3
"""Run the pipeline end to end: a raw folder in, annotated ion lists and panels out.

The stages, and what each one is allowed to assume:

  1 VENDOR    detect_vendor.py names the vendor and acquisition method from the
              folder's own files.
              GATE: an unknown acquisition method means drift_kind is unknown,
              and a tolerance applied to an unknown drift column is how a SONAR
              quadrupole position gets treated as a mobility measurement.
  2 EXTRACT   the adapter boundary. Vendor peak picking happens OUTSIDE this
              repository -- ApexRT, match_targets_spatial_in_cube, a vendor SDK.
              This stage either runs the command it is given or reads the peak
              list that command already wrote.
              GATE: with neither, the run stops. It does not invent peaks.
  3 ANNOTATE  match peaks to targets in m/z, and in rt and drift where both the
              target and the acquisition carry them.
              Tier is assigned from the evidence present, never from the fact
              that a match was found.
  4 IONLIST   write ions.jsonl: one ion_track record per ion, ion_uid content
              addressed at first detection, events append-only. This is the
              annotated ion list, and it is the thing worth keeping.
  5 EXPORT    write manifest.json describing what was actually produced,
              including the tables that are NOT there.
  6 CHECK     aria_check.py over the export: contract, panels, and the gates.
  7 PANELS    plan_panels.py says which panels can render, and
              render_panels.py draws them to a single self-contained page --
              blanks included, with the reason on each one.

    python3 tools/aria_run.py --raw RUN.raw --peaks peaks.jsonl \
        --targets targets.json --out export/ --sample-id CH215

    python3 tools/aria_run.py --raw RUN.d --extract 'apexrt --in {raw} --out {peaks}' \
        --atlas ARiA_Consensus.db --out export/

A peak record (JSON-lines, or CSV with a header) needs at least:

    mz, rt (minutes), charge, intensity

and may carry: drift, n_scans, score, score_kind, q_value, one_over_k0, ccs,
peptide, protein, seq_start, mods.

Tiers are not negotiable. An m/z window plus two consecutive scans is
bare_proximity and nothing more -- that bare mechanism passed reversed-tryptic
decoys at 33-70% elsewhere in this project. pass2_confirmed requires a
decoy-controlled score under a declared cut. atlas_consensus requires agreement
ACROSS samples, which a single-run spine cannot observe, so it never assigns it.

Exit status 0 when every stage that had inputs ran and passed, 1 when a gate
failed or a stage reported findings, 2 when a stage could not be run at all.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import detect_vendor  # noqa: E402

RULE = "=" * 78
TOOL_VERSION = "aria_run/1.0.0"

# From the production ingest contract. Changing these changes what the Atlas
# accepts, so they are named here rather than hidden in a default.
PRECURSOR_PPM = 10.0
PRODUCT_PPM = 15.0

PROTON = 1.007276
C13 = 1.003355

# bare_proximity's own definition: a window is not a detection without a peak.
MIN_CONSECUTIVE_SCANS = 2

# Score kinds that come with a decoy model behind them. A score from anywhere
# else is recorded but cannot lift an ion out of bare_proximity.
DECOY_CONTROLLED_SCORES = {"q_value", "fdr", "peptide_prophet", "percolator_q", "diann_qvalue"}
DEFAULT_Q_CUT = 0.01
# A q-value above this is not a control, so no flag may promote on it.
MAX_Q_CUT = 0.05

PEG_SPACING = 44.026215          # ethylene oxide repeat
PEG_TOL = 0.01


def now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def params_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True,
                                     default=str).encode()).hexdigest()[:16]


def ion_uid(sample_id, source_file, mz, rt, drift, charge):
    """Content addressed exactly as schema/ion_track.schema.json states.

    sha256(sample_id | source_file | mz@4dp | rt@3dp | drift@4dp-or-null | charge)

    The point of content addressing is that re-merging or re-clustering the same
    detection lands on the same uid. So the canonical string is fixed here and
    must not drift: a change to it silently orphans every earlier ion list.
    """
    parts = [str(sample_id), str(source_file), "%.4f" % float(mz), "%.3f" % float(rt),
             "null" if drift is None else "%.4f" % float(drift),
             "null" if charge is None else str(int(charge))]
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


class Run:
    """The stage log. A stage that did not run never looks like one that passed."""

    def __init__(self, verbose=False):
        self.stages = []
        self.stopped_by = None
        self.verbose = verbose
        self.warnings = []
        self.vendor_rec = None

    def banner(self, name):
        print("%s\n%s\n%s" % (RULE, name, RULE))

    def ok(self, name, detail=""):
        self.stages.append((name, "passed", detail))

    def findings(self, name, detail, stop=False):
        self.stages.append((name, "findings", detail))
        if stop:
            self.stopped_by = name

    def gate(self, name, why):
        self.stages.append((name, "GATE FAILED", why))
        self.stopped_by = name
        print("  GATE FAILED: %s" % why)

    def cannot(self, name, why):
        self.stages.append((name, "COULD NOT RUN", why))
        self.stopped_by = name
        print("  COULD NOT RUN: %s" % why)

    def skip(self, name, why):
        self.stages.append((name, "skipped", why))

    def unreached(self, name):
        why = ("stopped by %s" % self.stopped_by.strip()) if self.stopped_by \
              else "an earlier stage did not produce its input"
        self.stages.append((name, "NOT REACHED", why))

    def warn(self, text):
        self.warnings.append(text)

    def report(self):
        print("%s\nSUMMARY\n%s" % (RULE, RULE))
        width = max(len(n) for n, _, _ in self.stages)
        for name, state, detail in self.stages:
            print("  %-*s  %-14s %s" % (width, name, state, detail))
        print()
        for w in self.warnings:
            print("  warning: %s" % w)
        if self.warnings:
            print()
        bad = [n for n, s, _ in self.stages
               if s in ("GATE FAILED", "findings", "COULD NOT RUN")]
        unreached = [n for n, s, _ in self.stages if s == "NOT REACHED"]
        if unreached and not self.stopped_by:
            print("  %d stage(s) were NOT REACHED: %s."
                  % (len(unreached), ", ".join(u.strip() for u in unreached)))
            print("  Not reached is not passed.")
        if self.stopped_by:
            print("  The run stopped at %s. %d stage(s) NOT REACHED: %s."
                  % (self.stopped_by.strip(), len(unreached),
                     ", ".join(u.strip() for u in unreached) or "none"))
            print("  Not reached is not passed.")
        elif bad:
            print("  Every stage ran. %d reported findings: %s."
                  % (len(bad), ", ".join(b.strip() for b in bad)))
        else:
            ran = [n for n, s, _ in self.stages if s == "passed"]
            print("  %d stage(s) passed, %d skipped for want of inputs."
                  % (len(ran), len(self.stages) - len(ran)))
        if any(s == "COULD NOT RUN" for _, s, _ in self.stages):
            return 2
        return 1 if bad else 0

    def as_json(self):
        return {"tool": TOOL_VERSION, "finished_at": now(),
                "vendor": self.vendor_rec,
                "stages": [{"stage": n.strip(), "state": s, "detail": d}
                           for n, s, d in self.stages],
                "stopped_by": self.stopped_by.strip() if self.stopped_by else None,
                "warnings": list(self.warnings)}


# ---------------------------------------------------------------- 1 VENDOR

def stage_vendor(run, raw, override_drift_kind=None):
    run.banner("1 VENDOR    read the folder")
    rec = detect_vendor.detect(str(raw))
    if override_drift_kind:
        was = rec["drift_kind"]
        if was == override_drift_kind:
            run.warn("--drift-kind %s repeats what the folder's own files already"
                     " said, so it changed nothing." % override_drift_kind)
        elif was == "unknown":
            run.warn("drift_kind %s was declared with --drift-kind. The folder's"
                     " files did not say, so nothing in the data confirms it and"
                     " every drift tolerance below rests on the declaration."
                     % override_drift_kind)
        else:
            run.warn("--drift-kind %s CONTRADICTS the %s read from the folder's own"
                     " files. The declaration wins and both are in run.json, but one"
                     " of the two is wrong about what the drift column means."
                     % (override_drift_kind, was))
        rec["drift_kind"] = override_drift_kind
        rec["drift_kind_declared"] = True
        rec["drift_kind_detected"] = was
        rec.setdefault("evidence", []).append(
            "drift_kind set by --drift-kind, not read from files (files said: %s)" % was)
        if rec["acquisition_method"] == "unknown":
            rec["acquisition_method"] = "declared:" + override_drift_kind
    for line in rec.get("evidence", []):
        print("  evidence: %s" % line)
    for line in rec.get("why_unknown", []):
        print("  unknown:  %s" % line)
    print("  vendor              %s" % rec["vendor"])
    print("  acquisition_method  %s" % rec["acquisition_method"])
    print("  drift_kind          %s" % rec["drift_kind"])
    print()
    run.vendor_rec = rec
    if rec["vendor"] == "unknown":
        run.gate("1 VENDOR    read the folder",
                 "the vendor could not be established: %s. --drift-kind declares the"
                 " drift axis, not the vendor, so it cannot carry this."
                 % ("; ".join(rec.get("why_unknown", [])) or "no reason recorded"))
        return rec, False
    # A mixed-vendor folder and a folder-of-runs are refusals in their own right:
    # each run needs its own tolerances, and the method is per run. A declared
    # drift kind must not launder either of them.
    blocking = [w for w in rec.get("why_unknown", [])
                if "more than one vendor" in w or "folder of runs" in w]
    if blocking:
        run.gate("1 VENDOR    read the folder", " ".join(blocking))
        return rec, False
    if rec["drift_kind"] == "unknown":
        run.gate("1 VENDOR    read the folder",
                 "drift_kind is unknown, so no tolerance can be applied to the drift"
                 " column. Name it with --drift-kind if you know it from the method.")
        return rec, False
    if rec["drift_kind"] == "quad_window":
        run.warn("drift_kind is quad_window (SONAR): the drift column is a quadrupole"
                 " position, not mobility. Drift is carried but never used as a"
                 " matching dimension.")
    run.vendor_rec = rec
    run.ok("1 VENDOR    read the folder",
           "%s / %s / %s%s" % (rec["vendor"], rec["acquisition_method"],
                               rec["drift_kind"],
                               " (declared)" if rec.get("drift_kind_declared") else ""))
    return rec, True


# --------------------------------------------------------------- 2 EXTRACT

def read_peaks(path):
    """JSON-lines, a JSON array, or CSV with a header. Numbers stay numbers."""
    text = Path(path).read_text()
    stripped = text.lstrip()
    rows = []
    if stripped.startswith("["):
        rows = json.loads(text)
    elif stripped.startswith("{"):
        for n, line in enumerate(text.splitlines(), 1):
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except ValueError as exc:
                raise ValueError("%s line %d: %s" % (path, n, exc))
    else:
        for row in csv.DictReader(text.splitlines()):
            out = {}
            for k, v in row.items():
                if k is None:
                    continue
                k = k.strip()
                if v is None or v.strip() == "":
                    out[k] = None
                    continue
                v = v.strip()
                try:
                    out[k] = int(v) if re.fullmatch(r"-?\d+", v) else float(v)
                except ValueError:
                    out[k] = v
            rows.append(out)
    return rows


ALIASES = {
    "mz": ("mz", "m/z", "precursor_mz", "mz_observed", "observed_mz", "PrecursorMz"),
    "rt": ("rt", "rt_min", "retention_time", "rt_minutes", "apex_rt", "RT"),
    "drift": ("drift", "drift_time", "mobility", "ion_mobility", "quad_position", "scan_number_1k0"),
    "charge": ("charge", "z", "precursor_charge", "Charge"),
    "intensity": ("intensity", "matched_intensity", "apex_intensity", "summed_intensity", "Intensity"),
    "n_scans": ("n_scans", "n_consecutive_scans", "scan_count", "n_points"),
    "score": ("score", "hyperscore", "expect", "q_value", "qvalue"),
    "score_kind": ("score_kind", "score_type"),
    "q_value": ("q_value", "qvalue", "QValue", "fdr"),
    "one_over_k0": ("one_over_k0", "1/k0", "inv_k0"),
    "ccs": ("ccs", "CCS"),
    "peptide": ("peptide", "peptide_sequence", "sequence", "Peptide"),
    "protein": ("protein", "prot_accession", "accession", "Protein"),
    "seq_start": ("seq_start", "peptide_start", "start"),
    "mods": ("mods", "modifications", "assigned_modifications"),
    "source_file": ("source_file", "file", "raw_file", "spectrum_file"),
    "cluster_rank": ("cluster_rank", "rank"),
    "points": ("points", "trace", "xic", "elution"),
    "window_ppm": ("window_ppm", "extraction_window_ppm", "xic_window_ppm"),
    "apex_rt": ("apex_rt", "apex_rt_min"),
    "fwhm": ("fwhm", "fwhm_min"),
}


def field(row, name):
    for key in ALIASES[name]:
        if key in row and row[key] not in (None, ""):
            return row[key]
    return None


def stage_extract(run, raw, peaks_path, extract_cmd, workdir):
    run.banner("2 EXTRACT   vendor peak picking (adapter boundary)")
    if extract_cmd:
        out = Path(workdir) / "peaks.jsonl"
        cmd = extract_cmd.replace("{raw}", shlex.quote(str(raw))) \
                         .replace("{peaks}", shlex.quote(str(out))) \
                         .replace("{out}", shlex.quote(str(workdir)))
        print("  running: %s" % cmd)
        t0 = time.time()
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        body = ((proc.stdout or "") + (proc.stderr or "")).rstrip()
        if body:
            print("  ---- extractor output ----")
            for line in body.splitlines()[-40:]:
                print("  %s" % line)
            print("  --------------------------")
        if proc.returncode != 0:
            run.cannot("2 EXTRACT   vendor peaks",
                       "the extractor exited %d. No peak list was produced, so"
                       " nothing downstream has input." % proc.returncode)
            return None, False
        if not out.exists():
            run.cannot("2 EXTRACT   vendor peaks",
                       "the extractor exited 0 but wrote no %s. Exit 0 with no"
                       " output is not an empty result, it is a broken adapter." % out.name)
            return None, False
        peaks_path = out
        print("  extractor finished in %.1fs" % (time.time() - t0))
    if not peaks_path:
        run.gate("2 EXTRACT   vendor peaks",
                 "no --peaks and no --extract. Vendor peak picking happens outside"
                 " this repository; without one of the two there is nothing to read."
                 " This stage will not invent peaks.")
        return None, False
    try:
        rows = read_peaks(peaks_path)
    except (OSError, ValueError) as exc:
        run.cannot("2 EXTRACT   vendor peaks", "could not read %s: %s" % (peaks_path, exc))
        return None, False
    if not rows:
        run.gate("2 EXTRACT   vendor peaks",
                 "%s holds no peak records. An empty extraction is a result, but it"
                 " is not one anything downstream can be read from." % peaks_path)
        return None, False
    missing, chargeless = [], 0
    for i, row in enumerate(rows, 1):
        for need in ("mz", "rt"):
            if field(row, need) is None:
                missing.append("record %d has no %s" % (i, need))
        if field(row, "charge") is None:
            chargeless += 1
    if missing:
        for line in missing[:5]:
            print("  %s" % line)
        run.cannot("2 EXTRACT   vendor peaks",
                   "%d record(s) lack an m/z or a retention time. A peak without a"
                   " coordinate is not a detection." % len(missing))
        return None, False
    if chargeless:
        # Not a refusal. An MS1 detection whose envelope was not resolved has no
        # determined charge, and the contract tracks unassigned ions on purpose.
        # What it cannot do is carry an identity or a tier, and that is enforced
        # downstream rather than by discarding the measurement here.
        print("  %d record(s) carry no charge: recorded as not determined, which"
              " bars them from an annotation and from any tier above"
              " bare_proximity" % chargeless)
        run.warn("%d peak(s) have no determined charge. They are tracked with"
                 " coord.charge null and can never be identified from this run:"
                 " m/z cannot be turned into a peptide mass without z."
                 % chargeless)
    print("  read %d peak record(s) from %s" % (len(rows), peaks_path))
    print()
    run.ok("2 EXTRACT   vendor peaks", "%d peaks from %s" % (len(rows), Path(peaks_path).name))
    return rows, True


# -------------------------------------------------------------- 3 ANNOTATE

def load_targets(path):
    rows = read_peaks(path)
    out = []
    for r in rows:
        mz = field(r, "mz")
        if mz is None:
            continue
        out.append({"mz": float(mz),
                    "charge": int(field(r, "charge") or 0) or None,
                    "rt": None if field(r, "rt") is None else float(field(r, "rt")),
                    "drift": None if field(r, "drift") is None else float(field(r, "drift")),
                    "peptide": field(r, "peptide"),
                    "protein": field(r, "protein"),
                    "seq_start": field(r, "seq_start"),
                    "mods": field(r, "mods"),
                    "atlas_uid": r.get("uid") or r.get("atlas_uid") or r.get("ion_uid")})
    return out


def load_atlas_targets(run, db_path, limit):
    """Read targets out of a consensus DB without assuming its schema.

    The DB on the owner's machine has been renamed and re-cut more than once in
    this project, so this discovers the table and columns rather than hardcoding
    them, and says exactly what it used.
    """
    try:
        con = sqlite3.connect("file:%s?mode=ro" % db_path, uri=True)
    except sqlite3.Error as exc:
        return None, "could not open %s read-only: %s" % (db_path, exc)
    try:
        tables = [r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type IN ('table','view')")]
        best = None
        for t in tables:
            cols = {r[1].lower(): r[1] for r in con.execute('PRAGMA table_info("%s")' % t)}
            mz = next((cols[c] for c in ("mz", "precursor_mz", "mz_calculated", "calc_mz")
                       if c in cols), None)
            z = next((cols[c] for c in ("charge", "z", "precursor_charge") if c in cols), None)
            if not mz or not z:
                continue
            pep = next((cols[c] for c in ("peptide_sequence", "peptide", "sequence")
                        if c in cols), None)
            score = (pep is not None) * 2 + ("uid" in cols)
            if best is None or score > best[0]:
                best = (score, t, cols, mz, z, pep)
        if best is None:
            return None, ("no table in %s carries both an m/z and a charge column."
                          " Tables seen: %s" % (db_path, ", ".join(tables[:12]) or "none"))
        _, t, cols, mz, z, pep = best
        rt = next((cols[c] for c in ("rt", "rt_min", "retention_time", "apex_rt")
                   if c in cols), None)
        dr = next((cols[c] for c in ("drift", "drift_time", "mobility") if c in cols), None)
        pr = next((cols[c] for c in ("prot_accession", "protein", "accession") if c in cols), None)
        uid = cols.get("uid") or cols.get("ion_uid")
        picked = [mz, z] + [c for c in (pep, rt, dr, pr, uid) if c]
        sql = 'SELECT %s FROM "%s"' % (", ".join('"%s"' % c for c in picked), t)
        total = None
        if limit:
            # A bare LIMIT over an unordered table returns whatever the engine
            # hands back first, which can differ between runs of the same query.
            # Order it, and say how much was left behind.
            try:
                total = con.execute('SELECT COUNT(*) FROM "%s"' % t).fetchone()[0]
            except sqlite3.Error:
                total = None
            sql += ' ORDER BY "%s", "%s" LIMIT %d' % (mz, z, int(limit))
        print("  table %s, columns %s" % (t, ", ".join(picked)))
        out = []
        dropped = 0
        for row in con.execute(sql):
            rec = dict(zip(picked, row))
            try:
                out.append({"mz": float(rec[mz]), "charge": int(rec[z]),
                            "rt": None if not rt or rec.get(rt) is None else float(rec[rt]),
                            "drift": None if not dr or rec.get(dr) is None else float(rec[dr]),
                            "peptide": rec.get(pep), "protein": rec.get(pr),
                            "seq_start": None, "mods": None,
                            "atlas_uid": None if not uid else str(rec.get(uid))})
            except (TypeError, ValueError):
                dropped += 1
        if dropped:
            run.warn("%d row(s) in %s.%s carried an m/z or charge that is not a"
                     " number and were not loaded as targets. Peaks that should"
                     " have matched them will read as unassigned." % (dropped, t, mz))
        if total is not None and total > len(out) + dropped:
            run.warn("--atlas-limit %d took %d of %d rows in %s, ordered by m/z then"
                     " charge. The rest were not loaded, so an ion that would have"
                     " matched one of them reads as unassigned, not as missing."
                     % (limit, len(out) + dropped, total, t))
        return out, None
    finally:
        con.close()


def ppm(a, b):
    return abs(a - b) / b * 1e6


def tier_of(row, q_cut):
    """What the evidence on THIS row supports. Nothing more.

    atlas_consensus is deliberately unreachable here: it means the ion held
    across samples, and one run cannot see across samples.
    """
    kind = field(row, "score_kind")
    q = field(row, "q_value")
    if kind and str(kind).lower() in DECOY_CONTROLLED_SCORES and q is not None:
        if float(q) <= q_cut:
            return "pass2_confirmed", True, float(q), str(kind)
        return None, True, float(q), str(kind)
    score = field(row, "score")
    return "bare_proximity", False, (None if score is None else float(score)),\
        (str(kind) if kind else None)


def stage_annotate(run, peaks, targets, vendor_rec, args):
    run.banner("3 ANNOTATE  match peaks to targets")
    drift_kind = vendor_rec["drift_kind"]
    use_drift = (drift_kind == "ion_mobility" and args.drift_tol is not None)
    if drift_kind == "ion_mobility" and args.drift_tol is None:
        run.warn("the acquisition carries real ion mobility but no --drift-tol was"
                 " given, so drift was not used to separate co-isolated peaks."
                 " The dimension that does the separating sat unused.")
    if drift_kind == "quad_window":
        print("  drift is a quadrupole position here: carried, never matched on")
    print("  m/z window %.1f ppm%s%s" % (
        args.precursor_ppm,
        ", rt +/- %.2f min" % args.rt_tol if args.rt_tol is not None else ", rt unused",
        ", drift +/- %g" % args.drift_tol if use_drift else ", drift unused"))

    assigned = unassigned = 0
    ambiguous = 0
    no_charge = 0
    for row in peaks:
        mz = float(field(row, "mz"))
        zf = field(row, "charge")
        z = None if zf is None else int(zf)
        rt = float(field(row, "rt"))
        dr = field(row, "drift")
        dr = None if dr is None else float(dr)
        if z is None:
            # Without z there is no mass, so a target match would rest on m/z
            # coincidence alone. Left unassigned, which is a tracked state.
            row["_match"] = None
            no_charge += 1
            unassigned += 1
            continue
        hits = []
        for t in targets:
            if t["charge"] and t["charge"] != z:
                continue
            if ppm(mz, t["mz"]) > args.precursor_ppm:
                continue
            if args.rt_tol is not None and t["rt"] is not None and abs(rt - t["rt"]) > args.rt_tol:
                continue
            if use_drift and t["drift"] is not None and dr is not None \
                    and abs(dr - t["drift"]) > args.drift_tol:
                continue
            hits.append(t)
        if not hits:
            row["_match"] = None
            unassigned += 1
            continue
        if len(hits) > 1:
            ambiguous += 1
            hits.sort(key=lambda t: ppm(mz, t["mz"]))
        row["_match"] = hits[0]
        row["_n_candidates"] = len(hits)
        row["_ppm"] = (mz - hits[0]["mz"]) / hits[0]["mz"] * 1e6
        assigned += 1

    print("  %d assigned, %d unassigned, %d had more than one candidate"
          % (assigned, unassigned, ambiguous))
    if no_charge:
        print("  %d of the unassigned have no determined charge, so no target could"
              " be tested against them" % no_charge)
    print()
    if ambiguous:
        run.warn("%d peak(s) matched more than one target inside the window and were"
                 " assigned the nearest in m/z. Nearest-in-ppm is a tie-break, not"
                 " evidence: see docs/resolvability.md for what the window can and"
                 " cannot separate." % ambiguous)
    if assigned == 0:
        if no_charge == len(peaks):
            print("  no peak had a determined charge, so no target could be tested"
                  " against any of them. This is not a failed match.")
            run.skip("3 ANNOTATE  match to targets",
                     "no peak carried a charge, so nothing was testable")
            run.warn("the target list was never consulted: every peak lacks a"
                     " determined charge, so there was nothing to test against."
                     " The ions below are tracked and unassigned, which is a state,"
                     " not a result.")
            return True
        run.findings("3 ANNOTATE  match to targets",
                     "no peak matched any target. Of %d peak(s), %d were testable"
                     " and none matched, so either the targets are wrong for this"
                     " run or the window is."
                     % (len(peaks), len(peaks) - no_charge), stop=True)
        return False
    run.ok("3 ANNOTATE  match to targets",
           "%d/%d assigned" % (assigned, len(peaks)))
    return True


# --------------------------------------------------------------- 4 IONLIST

WATER = 18.010565


def looks_like_peg(mz, z):
    """Does this sit on the PEG series lattice?

    PEG is H(OCH2CH2)nOH: neutral mass n * 44.026215 + a water end group. The
    water matters -- test the lattice without it and every real PEG ion is 18 Da
    off it, so the filter names a contaminant it can never catch.

    Only the plain protonated series. Sodiated and ammoniated PEG sit on their
    own lattices and are not tested here, so a clean result is not a clean run.
    """
    if not z:
        return False
    neutral = mz * z - z * PROTON - WATER
    k = round(neutral / PEG_SPACING)
    return k >= 3 and abs(neutral - k * PEG_SPACING) < PEG_TOL * z


def stage_ionlist(run, peaks, vendor_rec, args, out_dir):
    run.banner("4 IONLIST   write the annotated ion list")
    drift_kind = vendor_rec["drift_kind"]
    sample_id = args.sample_id or Path(args.raw).name
    source_file = str(args.raw)
    ph = params_hash({"precursor_ppm": args.precursor_ppm, "rt_tol": args.rt_tol,
                      "drift_tol": args.drift_tol, "q_cut": args.q_cut,
                      "drop_z1": args.drop_z1, "drop_peg": args.drop_peg,
                      "min_scans": args.min_scans,
                      # The drift semantics and the target list change what the
                      # same peaks mean, so they belong in the hash. Without them
                      # two runs that disagree about the drift axis are
                      # indistinguishable from the ion record.
                      "drift_kind": drift_kind,
                      "acquisition_method": vendor_rec["acquisition_method"],
                      "drift_kind_declared": bool(vendor_rec.get("drift_kind_declared")),
                      "targets": str(args.targets or args.atlas or "none")})
    stamp = now()
    kept, kept_rows, rejected, conflicts = [], [], [], []
    tiers = {"bare_proximity": 0, "pass2_confirmed": 0, "atlas_consensus": 0}
    seen = {}
    dup = 0

    for row in peaks:
        mz = float(field(row, "mz"))
        zf = field(row, "charge")
        z = None if zf is None else int(zf)
        rt = float(field(row, "rt"))
        dr = field(row, "drift")
        dr = None if dr is None else float(dr)
        if drift_kind == "absent":
            dr = None
        src = field(row, "source_file") or source_file
        uid = ion_uid(sample_id, src, mz, rt, dr, z)

        events = [{"stage": "detected", "at": stamp, "tool": TOOL_VERSION,
                   "tool_version": None, "params_hash": ph,
                   "coord_snapshot": {"mz": mz, "rt_min": rt, "drift": dr},
                   "note": "peak read at the extract boundary"}]

        coord = {"mz": mz, "rt_min": rt, "charge": z, "drift": dr,
                 "drift_kind": drift_kind}
        if drift_kind == "ion_mobility":
            k0 = field(row, "one_over_k0")
            ccs = field(row, "ccs")
            if k0 is not None:
                coord["one_over_k0"] = float(k0)
            if ccs is not None:
                coord["ccs"] = float(ccs)

        def reject(stage, why):
            events.append({"stage": stage, "at": stamp, "tool": TOOL_VERSION,
                           "tool_version": None, "params_hash": ph,
                           "coord_snapshot": None, "note": why})
            rejected.append({"ion_uid": uid, "sample_id": sample_id, "coord": coord,
                             "reason": why, "events": events})

        if args.drop_z1 and z == 1:
            reject("z1_filtered", "charge 1: excluded by --drop-z1")
            continue
        if args.drop_peg and looks_like_peg(mz, z):
            reject("peg_filtered", "sits on the %g/z PEG lattice" % PEG_SPACING)
            continue

        n_scans = field(row, "n_scans")
        n_scans = None if n_scans is None else int(n_scans)
        if n_scans is None:
            reject("rejected", "no scan count, so the >=%d consecutive scan evidence"
                               " that defines bare_proximity is absent"
                   % args.min_scans)
            continue
        if n_scans < args.min_scans:
            reject("rejected", "%d consecutive scan(s), below the %d that"
                               " bare_proximity requires" % (n_scans, args.min_scans))
            continue

        match = row.get("_match")
        tier, decoy, score, kind = tier_of(row, args.q_cut)
        if z is None and tier != "bare_proximity":
            # A decoy-controlled score is a score for an identification, and
            # there is no identification without z. The score is kept on the
            # record; the tier it would buy is not.
            run.warn("an ion with no determined charge carried a %s score; it is"
                     " recorded but the tier stays bare_proximity, because a score"
                     " certifies an identification and there is none without z."
                     % kind)
            tier, decoy = "bare_proximity", False
        if tier is None:
            reject("rejected", "decoy-controlled %s = %s is above the %g cut"
                   % (kind, score, args.q_cut))
            continue
        if match is None:
            tier = "bare_proximity"
            decoy = False

        annotation = None
        if match is not None:
            # The row's own sequence is what its score was computed for, so it
            # wins. The target list SELECTS ions; it does not identify them. Take
            # the target's sequence only where the row has none of its own.
            row_pep = field(row, "peptide")
            pep = row_pep or match.get("peptide")
            if row_pep and match.get("peptide") and \
                    str(row_pep).upper() != str(match["peptide"]).upper():
                conflicts.append((str(row_pep), str(match["peptide"])))
                events.append({"stage": "annotated", "at": stamp, "tool": TOOL_VERSION,
                               "tool_version": None, "params_hash": ph,
                               "coord_snapshot": None,
                               "note": "the target list proposed %s; the row's own"
                                       " identification %s was kept, because the"
                                       " score belongs to it"
                                       % (match["peptide"], row_pep)})
            if pep and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", str(pep)):
                prot = (field(row, "protein") if row_pep else None) \
                    or match.get("protein") or field(row, "protein")
                annotation = {"peptide_sequence": str(pep),
                              "prot_accession": str(prot) if prot else None,
                              "protein_name": None,
                              "seq_start": (int(match["seq_start"])
                                            if match.get("seq_start") and not row_pep
                                            else None),
                              "mods": field(row, "mods") or match.get("mods"),
                              "frag_type": None}
                events.append({"stage": "annotated", "at": stamp, "tool": TOOL_VERSION,
                               "tool_version": None, "params_hash": ph,
                               "coord_snapshot": None,
                               "note": "matched %s within %.2f ppm%s"
                                       % (pep, abs(row.get("_ppm") or 0.0),
                                          "" if row.get("_n_candidates", 1) < 2 else
                                          "; %d targets were inside the window and"
                                          " this one was nearest in m/z, which is a"
                                          " tie-break and not evidence"
                                          % row["_n_candidates"])})
            elif pep:
                run.warn("a target carried the sequence %r, which is not a bare"
                         " residue string; it was not written as an annotation."
                         % str(pep)[:24])

        inten = field(row, "intensity")
        rank = field(row, "cluster_rank")
        rec = {"ion_uid": uid, "sample_id": sample_id, "coord": coord,
               "intensity": {"matched_intensity": None if inten is None else float(inten),
                             "summed_precursor_intensity": None,
                             "cluster_rank": None if rank is None else int(rank),
                             "n_consecutive_scans": n_scans},
               "annotation": annotation,
               "confidence": {"tier": tier, "decoy_controlled": bool(decoy),
                              "ppm_error": row.get("_ppm"), "score": score,
                              "score_kind": ((kind or "") +
                                             ("" if row.get("_n_candidates", 1) < 2
                                              else " (%d candidates in window)"
                                              % row["_n_candidates"])) or None},
               "provenance": {"source_file": src,
                              "vendor": vendor_rec["vendor"],
                              "acquisition_method": vendor_rec["acquisition_method"],
                              "has_ims": drift_kind == "ion_mobility",
                              "search_engine": args.search_engine,
                              "batch": args.batch,
                              "sample_class": args.sample_class},
               "lineage": {"parent_ion_uids": [], "cluster_uid": None,
                           "atlas_uid": (match or {}).get("atlas_uid")},
               "events": events}
        if uid in seen:
            dup += 1
            continue
        seen[uid] = True
        tiers[tier] += 1
        kept.append(rec)
        kept_rows.append(row)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ions_path = out_dir / "ions.jsonl"
    with ions_path.open("w") as fh:
        for rec in kept:
            fh.write(json.dumps(rec) + "\n")
    rej_path = out_dir / "rejected.jsonl"
    with rej_path.open("w") as fh:
        for rec in rejected:
            fh.write(json.dumps(rec) + "\n")

    # The resolvability input, written from the ions themselves. uid here is the
    # IDENTITY (sequence and charge), not ion_uid: ion_uid is content addressed
    # per detection, so each one would have exactly one observation and the
    # within-ion SD that resolvability needs could never be estimated.
    resolve = []
    for rec in kept:
        ann = rec["annotation"]
        if not ann or rec["coord"]["drift"] is None:
            continue
        resolve.append({"run": sample_id,
                        "uid": "%s_%d" % (ann["peptide_sequence"], rec["coord"]["charge"]),
                        "mz": rec["coord"]["mz"], "z": rec["coord"]["charge"],
                        "drift": rec["coord"]["drift"], "rt": rec["coord"]["rt_min"]})
    # Traces, where the extractor supplied the elution points. Two points is the
    # floor: one point is not a peak, and the XIC panel refuses it anyway.
    traces, short, no_window = [], 0, 0
    for rec, row in zip(kept, kept_rows):
        pts = field(row, "points")
        if not isinstance(pts, list) or not pts:
            continue
        clean = [[float(q[0]), float(q[1])] for q in pts
                 if isinstance(q, (list, tuple)) and len(q) >= 2]
        if len(clean) < 2:
            short += 1
            continue
        t = {"ion_uid": rec["ion_uid"], "sample_id": rec["sample_id"], "kind": "xic",
             "points": clean}
        w = field(row, "window_ppm")
        if w is not None and float(w) > 0:
            t["extraction_window_ppm"] = float(w)
        else:
            # --precursor-ppm is the window this run MATCHED in. The window the
            # extractor pulled the trace through is a different number, and
            # substituting one for the other would put a figure on the panel
            # that nothing measured.
            no_window += 1
        a = field(row, "apex_rt")
        t["apex_rt_min"] = float(a) if a is not None else \
            max(clean, key=lambda q: q[1])[0]
        f = field(row, "fwhm")
        if f is not None:
            t["fwhm_min"] = float(f)
        traces.append(t)
    if traces:
        with (Path(out_dir) / "traces.jsonl").open("w") as fh:
            for t in traces:
                fh.write(json.dumps(t) + "\n")
        print("  wrote %d trace(s) to %s"
              % (len(traces), Path(out_dir) / "traces.jsonl"))
    if no_window:
        run.warn("%d trace(s) state no extraction window. --precursor-ppm is the"
                 " window this run matched in, not the window the extractor pulled"
                 " the trace through, so none was written. Without it there is no"
                 " telling how much of the peak is the ion." % no_window)
    if short:
        run.warn("%d peak(s) carried fewer than 2 elution points, so no trace was"
                 " written for them. One point is not a peak." % short)

    resolve_path = Path(out_dir) / "resolve_ions.json"
    with resolve_path.open("w") as fh:
        for r in resolve:
            fh.write(json.dumps(r) + "\n")

    why = {}
    for rec in rejected:
        key = rec["events"][-1]["stage"]
        why[key] = why.get(key, 0) + 1
    print("  wrote %d ion(s) to %s" % (len(kept), ions_path))
    print("  wrote %d rejected to %s%s" % (len(rejected), rej_path,
          (" (" + ", ".join("%s %d" % kv for kv in sorted(why.items())) + ")") if why else ""))
    for tier, n in sorted(tiers.items()):
        print("    %-18s %d" % (tier, n))
    print("  wrote %d resolvability observation(s) to %s" % (len(resolve), resolve_path))
    if dup:
        # Nothing is combined: the first record wins and the later one is dropped.
        # Saying "merged" would imply intensities were summed. They are not.
        print("  %d later peak(s) had the same content address as an earlier one"
              " and were DROPPED, not combined" % dup)
    print()
    if not kept:
        run.gate("4 IONLIST   annotated ions",
                 "no ion survived. %d peak(s) were rejected; the ion list is empty,"
                 " so there is nothing to export, check or draw." % len(rejected))
        return None, tiers, False
    if conflicts:
        shown = ", ".join("%s vs %s" % c for c in conflicts[:3])
        run.warn("%d ion(s) carried their own identification AND matched a different"
                 " sequence in the target list (%s%s). The row's own sequence was"
                 " kept, because the score that promotes the ion was computed for"
                 " it. A target list that disagrees this often is the wrong list"
                 " for this run." % (len(conflicts), shown,
                                     ", ..." if len(conflicts) > 3 else ""))
    if drift_kind == "ion_mobility" and not any(
            r["coord"].get("one_over_k0") is not None for r in kept):
        run.warn("the acquisition carries real ion mobility but no peak brought a"
                 " 1/K0 value, so every mobility panel plots an empty axis. The"
                 " extractor has the number; it is not reaching the peak list.")
    if tiers["bare_proximity"] == len(kept):
        run.warn("every ion is bare_proximity: an m/z window plus %d consecutive"
                 " scans, with no decoy control. That mechanism passed"
                 " reversed-tryptic decoys at 33-70%% elsewhere in this project."
                 % args.min_scans)
    run.ok("4 IONLIST   annotated ions", "%d ions, %d rejected" % (len(kept), len(rejected)))
    return kept, tiers, True


def wrote_traces(out_dir):
    return (Path(out_dir) / "traces.jsonl").exists()


def resolve_is_usable(path, args):
    """Can step 4 of the chain actually answer, or would it only look like it?

    Within-ion sigma needs a precursor measured at least three times. One run
    gives one observation per identity, so a single-run export can only feed the
    chain when the scale is supplied from outside with --sigma-drift/--sigma-rt.
    """
    if args.sigma_drift is not None and args.sigma_rt is not None:
        return True, "scales given on the command line"
    counts = {}
    try:
        for line in Path(path).read_text().splitlines():
            if line.strip():
                uid = json.loads(line)["uid"]
                counts[uid] = counts.get(uid, 0) + 1
    except (OSError, ValueError, KeyError):
        return False, "could not be read back"
    repeated = sum(1 for n in counts.values() if n >= 3)
    if not repeated:
        return False, ("no identity was measured 3 times in this export, so the"
                       " within-ion scale cannot be estimated. Give --sigma-drift"
                       " and --sigma-rt, or merge several runs first.")
    if len(counts) < 2:
        return False, ("only one distinct identity is present. With one precursor"
                       " the axis range IS its own scatter, so an estimated scale"
                       " can only come back as a confound -- and there is nothing"
                       " for it to collide with anyway.")
    return True, "%d identity(ies) measured 3+ times, %d distinct" % (repeated, len(counts))


# ----------------------------------------------------------------- 5 EXPORT

def stage_export(run, ions, tiers, vendor_rec, args, out_dir, peaks_n):
    run.banner("5 EXPORT    manifest")
    sample_id = args.sample_id or Path(args.raw).name
    drift_kind = vendor_rec["drift_kind"]
    manifest = {
        "export_id": "run-%s-%s" % (sample_id, now().replace(":", "").replace("-", "")),
        "created_at": now(),
        "method_id": args.method_id or "aria-run-spine",
        "method_version": "1.0.0",
        # sources is a closed set in the contract; the target list is a tool
        # input, so it is recorded under pipeline, which takes any named string.
        "sources": {"raw_dir": str(args.raw), "consensus_db": args.atlas},
        "pipeline": {"aria_run": TOOL_VERSION,
                     "extract": args.extract or ("peaks:" + str(args.peaks)),
                     "targets": str(args.targets) if args.targets else "none",
                     "detect_vendor": "detect_vendor@" + TOOL_VERSION},
        "samples": [{"sample_id": sample_id,
                     "vendor": vendor_rec["vendor"],
                     "acquisition_method": vendor_rec["acquisition_method"],
                     "drift_kind": drift_kind,
                     "has_ims": drift_kind == "ion_mobility",
                     "search_engine": args.search_engine,
                     "batch": args.batch,
                     "sample_class": args.sample_class,
                     "role": args.role,
                     "n_proteins": len({i["annotation"]["prot_accession"]
                                        for i in ions
                                        if i["annotation"]
                                        and i["annotation"]["prot_accession"]}),
                     "total_signal": sum(i["intensity"]["matched_intensity"] or 0.0
                                         for i in ions)}],
        # Declared from what was written, not from what was hoped for. A table
        # claimed but absent is the contract violation validate_export looks for.
        "tables": {"ions": True, "proteins": False, "spectra": False,
                   "traces": wrote_traces(out_dir), "spectrum_kinds": [],
                   "fragment_hits": False},
        "counts": {"ions_total": peaks_n, "ions_exported": len(ions),
                   "by_tier": dict(tiers)},
        "decimation": {"strategy": "none", "kept": len(ions), "of": peaks_n, "seed": 0},
        "warnings": list(run.warnings),
    }
    # No `quant` block. The contract's only quant metric is top3avg, which is a
    # per-protein roll-up this spine does not compute -- tables.proteins is
    # false. Declaring top3avg here would name a metric nothing produced, and
    # naming matched_intensity would invent an enum value. So quant is absent,
    # which is what "no quant claim" looks like in this contract.
    path = Path(out_dir) / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    print("  wrote %s" % path)
    print("  ions_total %d, ions_exported %d" % (peaks_n, len(ions)))
    print()
    run.ok("5 EXPORT    manifest", str(path))
    return manifest


# ------------------------------------------------------------ 6 CHECK / 7 PANELS

def stage_check(run, out_dir, args):
    run.banner("6 CHECK     aria_check over the export")
    extra = []
    if args.excerpt:
        extra.append("--excerpt")
    if args.limit:
        extra += ["--limit", str(args.limit)]
    usable, why = resolve_is_usable(Path(out_dir) / "resolve_ions.json", args)
    if usable:
        extra += ["--ions", str(Path(out_dir) / "resolve_ions.json")]
        if args.sigma_drift is not None:
            extra += ["--sigma-drift", str(args.sigma_drift)]
        if args.sigma_rt is not None:
            extra += ["--sigma-rt", str(args.sigma_rt)]
        if args.rank_rt:
            extra.append("--rank-rt")
        print("  resolvability fed from this export: %s" % why)
    else:
        print("  resolvability NOT fed: %s" % why)
        run.warn("the chain's collision step was not fed: %s" % why)
    if args.fragments:
        extra += ["--fragments"] + list(args.fragments)
    if args.counts:
        extra += ["--counts", args.counts]
    argv = [sys.executable, str(TOOLS / "aria_check.py"), "--export", str(out_dir)] + extra
    proc = subprocess.run(argv, capture_output=True, text=True)
    body = ((proc.stdout or "") + (proc.stderr or "")).rstrip()
    print(body or "(no output)")
    print()
    if proc.returncode == 2:
        run.cannot("6 CHECK     contract + gates", "aria_check.py exited 2")
        return False
    if proc.returncode != 0:
        run.findings("6 CHECK     contract + gates",
                     "aria_check reported findings; read its summary above")
        return False
    run.ok("6 CHECK     contract + gates", "")
    return True


def stage_panels(run, method, out_dir):
    run.banner("7 PANELS    plan, then draw")
    plan_ok = True
    if method:
        argv = [sys.executable, str(TOOLS / "plan_panels.py"), str(method),
                str(Path(out_dir) / "manifest.json")]
        proc = subprocess.run(argv, capture_output=True, text=True)
        body = ((proc.stdout or "") + (proc.stderr or "")).rstrip()
        print(body or "(no output)")
        (Path(out_dir) / "panels.txt").write_text(body + "\n")
        print("\n  wrote %s" % (Path(out_dir) / "panels.txt"))
        if proc.returncode == 2:
            run.cannot("7 PANELS    plan + draw", "plan_panels.py exited 2")
            return False
        plan_ok = proc.returncode == 0
    else:
        print("  no --method given: every panel the renderer can draw is offered")

    argv = [sys.executable, str(TOOLS / "render_panels.py"), str(out_dir)]
    if method:
        argv += ["--method", str(method)]
    proc = subprocess.run(argv, capture_output=True, text=True)
    body = ((proc.stdout or "") + (proc.stderr or "")).rstrip()
    for line in body.splitlines():
        print("  %s" % line)
    print()
    if proc.returncode != 0:
        run.cannot("7 PANELS    plan + draw",
                   "render_panels.py exited %d, so the run produced no page"
                   % proc.returncode)
        return False
    if not plan_ok:
        run.findings("7 PANELS    plan + draw",
                     "the page was drawn, but some requested panels are blank on it")
        return False
    run.ok("7 PANELS    plan + draw", str(Path(out_dir) / "panels.html"))
    return True


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw", required=True, help="the raw run folder")
    p.add_argument("--out", required=True, help="export directory to write")
    p.add_argument("--peaks", help="peak list the extractor already wrote (JSONL/JSON/CSV)")
    p.add_argument("--extract", metavar="CMD",
                   help="command to produce peaks; {raw} {peaks} {out} are substituted")
    p.add_argument("--targets", help="target list (JSONL/JSON/CSV) to annotate against")
    p.add_argument("--atlas", help="consensus DB to read targets from")
    p.add_argument("--atlas-limit", type=int, default=200000,
                   help="cap on targets read from --atlas (default 200000)")
    p.add_argument("--method", help="method JSON for plan_panels.py")
    p.add_argument("--method-id")
    p.add_argument("--sample-id")
    p.add_argument("--batch")
    p.add_argument("--sample-class")
    p.add_argument("--search-engine")
    p.add_argument("--role", default="target")
    p.add_argument("--drift-kind", choices=("ion_mobility", "quad_window", "absent"),
                   help="declare the drift kind when the folder's files do not say")
    p.add_argument("--precursor-ppm", type=float, default=PRECURSOR_PPM)
    p.add_argument("--rt-tol", type=float, help="retention tolerance, minutes")
    p.add_argument("--drift-tol", type=float, help="drift tolerance, same units as drift")
    p.add_argument("--q-cut", type=float, default=DEFAULT_Q_CUT)
    p.add_argument("--min-scans", type=int, default=MIN_CONSECUTIVE_SCANS)
    p.add_argument("--drop-z1", action="store_true", help="drop singly charged peaks")
    p.add_argument("--drop-peg", action="store_true", help="drop the PEG lattice")
    p.add_argument("--sigma-drift", type=float, help="drift resolution scale for the chain")
    p.add_argument("--sigma-rt", type=float, help="retention resolution scale for the chain")
    p.add_argument("--rank-rt", action="store_true", help="rank retention instead of minutes")
    p.add_argument("--fragments", nargs="+", metavar="JSONL",
                   help="fragment tables for the chain's chemistry gate")
    p.add_argument("--counts", help="gates/ions JSON for the chain's interval step")
    p.add_argument("--excerpt", action="store_true", help="pass --excerpt to the checker")
    p.add_argument("--limit", type=int, help="pass --limit to the checker")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    # These two define what the tiers MEAN. bare_proximity is a window plus at
    # least two consecutive scans; pass2_confirmed is a decoy-controlled score
    # under a credible cut. A flag that loosens either does not widen the export,
    # it renames a weaker thing as a stronger one.
    if args.min_scans < MIN_CONSECUTIVE_SCANS:
        p.error("--min-scans %d: bare_proximity is defined as a window plus at least"
                " %d consecutive scans. A lower floor does not keep more ions, it"
                " relabels single-scan noise as a detection."
                % (args.min_scans, MIN_CONSECUTIVE_SCANS))
    if args.q_cut > MAX_Q_CUT:
        p.error("--q-cut %g: above %g a decoy-controlled score is not evidence of"
                " anything, so promoting on it would print pass2_confirmed over"
                " ions that have not been controlled. Lower the cut, or accept"
                " bare_proximity." % (args.q_cut, MAX_Q_CUT))
    if args.precursor_ppm > 50.0:
        p.error("--precursor-ppm %g: an isotope pick sits 358-1003 ppm from the"
                " calculated m/z, and a window this wide stops being a mass"
                " measurement. See docs/verification-run.md."
                % args.precursor_ppm)

    run = Run(args.verbose)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    names = ["1 VENDOR    read the folder", "2 EXTRACT   vendor peaks",
             "3 ANNOTATE  match to targets", "4 IONLIST   annotated ions",
             "5 EXPORT    manifest", "6 CHECK     contract + gates",
             "7 PANELS    plan + draw"]

    def finish(stage_index):
        for n in names[stage_index:]:
            if not any(s[0] == n for s in run.stages):
                run.unreached(n)
        (out_dir / "run.json").write_text(json.dumps(run.as_json(), indent=2) + "\n")
        return run.report()

    vendor_rec, ok = stage_vendor(run, args.raw, args.drift_kind)
    if not ok:
        return finish(1)

    peaks, ok = stage_extract(run, args.raw, args.peaks, args.extract, out_dir)
    if not ok:
        return finish(2)

    if args.targets and args.atlas:
        print("  both --targets and --atlas were given. They are two different"
              " target lists and this does not guess between them.")
        run.cannot("3 ANNOTATE  match to targets",
                   "--targets and --atlas were both given. Pick one: a manifest that"
                   " records a consensus DB it never opened is a false audit trail.")
        return finish(2)

    targets, why = [], None
    if args.targets:
        try:
            targets = load_targets(args.targets)
        except (OSError, ValueError) as exc:
            why = "could not read %s: %s" % (args.targets, exc)
    elif args.atlas:
        run.banner("3 ANNOTATE  targets from the Atlas")
        targets, why = load_atlas_targets(run, args.atlas, args.atlas_limit)
        targets = targets or []
    if why:
        run.cannot("3 ANNOTATE  match to targets", why)
        return finish(3)
    if not targets:
        run.skip("3 ANNOTATE  match to targets",
                 "no --targets and no --atlas: ions stay unassigned")
        run.warn("no target list was given, so every ion is unassigned. Unassigned"
                 " ions still track and still render, but nothing here is an"
                 " identification.")
        for row in peaks:
            row["_match"] = None
    else:
        print("  %d target(s) loaded" % len(targets))
        if not stage_annotate(run, peaks, targets, vendor_rec, args):
            return finish(4)

    ions, tiers, ok = stage_ionlist(run, peaks, vendor_rec, args, out_dir)
    if not ok:
        return finish(5)

    stage_export(run, ions, tiers, vendor_rec, args, out_dir, len(peaks))

    stage_check(run, out_dir, args)
    stage_panels(run, args.method, out_dir)
    return finish(7)


if __name__ == "__main__":
    sys.exit(main())
