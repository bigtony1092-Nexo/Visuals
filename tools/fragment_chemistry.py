#!/usr/bin/env python3
"""Ask whether a set of fragment patterns carries fragmentation chemistry at all.

A pattern comparison is only worth running if both sides have signal in them. The
PAX8 page's section 3 compared Colorado's extracted patterns against a Thermo
archive and the mismatched control beat the correct pairing — because the archive
side sat at the uniform expectation on every chemical marker. It held no
fragmentation information, and a flat reference cannot be matched by anything.

This checks that before a comparison, not after. For each source it reports the
share of b/y intensity landing where chemistry says it should, against the share
expected if every backbone bond broke equally given the same observed fragments:

    y-ion share            beam-type CID/HCD favours y over b
    N-terminal to proline  the strongest rule in peptide fragmentation
    C-terminal to aspartate weaker, but real

Enrichment near 1.0x on all three means the source is indistinguishable from
uniform. It is not a reference; nothing can be matched against it.

    python3 tools/fragment_chemistry.py patterns.jsonl
    python3 tools/fragment_chemistry.py a.jsonl b.jsonl        # compare two sources
    python3 tools/fragment_chemistry.py patterns.jsonl --by-source

Input is JSON lines, one per fragment:

    {"source": "march", "sequence": "ATLTPSNTPLGR", "charge": 2,
     "fragment": "y7", "intensity": 1234.0}

`fragment` is a b or y ion with its ordinal, optionally with a charge suffix
("y7", "b4-1", "y11^2"). Anything else is counted and skipped. Intensities may be
absolute or relative; only shares within a source are used, so normalization per
spectrum does not matter here.

Exit status is 1 if any source looks uniform, so this can gate a comparison.
"""

import argparse
import json
import re
import sys
from collections import defaultdict

FRAG = re.compile(r"^([by])(\d+)")
UNIFORM_TOL = 1.25      # enrichment below this is not distinguishable from uniform


def bond_of(label, length):
    """Which backbone bond a fragment reports, 1-based, or None.

    b(n) covers the bond between residues n and n+1. y(n) covers the bond between
    residues L-n and L-n+1, so its bond index is L-n.
    """
    m = FRAG.match(label.strip())
    if not m:
        return None
    kind, n = m.group(1), int(m.group(2))
    k = n if kind == "b" else length - n
    return k if 1 <= k <= length - 1 else None


def markers(sequence, k):
    """Which chemical rules the bond at index k satisfies."""
    return {
        # The bond N-terminal to a proline: residue k+1 is P.
        "proline": sequence[k] == "P",
        # The bond C-terminal to an aspartate: residue k is D.
        "aspartate": sequence[k - 1] == "D",
    }


def summarise(rows):
    """Observed intensity share against the share expected if bonds broke equally."""
    tot = y_int = 0.0
    obs = defaultdict(float)
    cnt = defaultdict(int)
    n_frag = 0
    skipped = 0
    for seq, label, inten in rows:
        k = bond_of(label, len(seq))
        if k is None:
            skipped += 1
            continue
        n_frag += 1
        tot += inten
        if label.strip()[0] == "y":
            y_int += inten
        for name, hit in markers(seq, k).items():
            if hit:
                obs[name] += inten
                cnt[name] += 1
    if not n_frag or tot <= 0:
        return None
    out = {"fragments": n_frag, "skipped": skipped, "y_share": y_int / tot}
    for name in ("proline", "aspartate"):
        share = obs[name] / tot
        # The expectation is taken over the fragments this source actually
        # observed, which is what "if every bond broke equally" means here.
        expected = cnt[name] / n_frag
        out[name] = (share, expected, (share / expected) if expected else None)
    return out


def load(path):
    by_source = defaultdict(list)
    bad = 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
                seq = r["sequence"].strip().upper()
                inten = float(r["intensity"])
                frag = str(r["fragment"])
            except (ValueError, KeyError, AttributeError, TypeError):
                bad += 1
                continue
            if not seq or inten <= 0:
                continue
            by_source[str(r.get("source") or path)].append((seq, frag, inten))
    return by_source, bad


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="JSON-lines fragment tables")
    p.add_argument("--by-source", action="store_true",
                   help="split each file by its `source` field (the default anyway)")
    args = p.parse_args(argv)

    sources = {}
    for path in args.files:
        got, bad = load(path)
        if bad:
            print("%s: %d unparseable lines skipped" % (path, bad), file=sys.stderr)
        for name, rows in got.items():
            key = name if name not in sources else "%s (%s)" % (name, path)
            sources[key] = rows
    if not sources:
        print("No usable fragments read.", file=sys.stderr)
        return 1

    print("Share of b/y intensity on the bonds chemistry favours, against the share")
    print("expected if every backbone bond broke equally given the same fragments.")
    print()
    print("%-22s %9s %8s %22s %22s" % (
        "source", "fragments", "y-share", "proline (obs/exp)", "aspartate (obs/exp)"))
    flat, unassessable = [], []
    for name, rows in sources.items():
        s = summarise(rows)
        if s is None:
            print("%-22s %9s" % (name[:22], "none usable"))
            continue
        cells = []
        for marker in ("proline", "aspartate"):
            share, exp, enr = s[marker]
            cells.append("%.3f / %.3f = %4.2fx" % (share, exp, enr) if enr
                         else "%.3f / n/a" % share)
        print("%-22s %9d %8.3f %22s %22s" % (
            name[:22], s["fragments"], s["y_share"], cells[0], cells[1]))
        enr = s["proline"][2]
        if enr is None:
            # No proline-containing peptide, so the marker had no opportunity.
            # That is "cannot assess", not "passes": leaving it out of `flat`
            # let a literally uniform source print the all-clear and exit 0.
            unassessable.append(name)
        elif enr < UNIFORM_TOL:
            flat.append((name, enr))
    print()
    for name, enr in flat:
        print("  %s: proline enrichment %.2fx is below %.2fx — indistinguishable from"
              % (name, enr, UNIFORM_TOL))
        print("  uniform. This source carries no fragmentation chemistry, so it cannot")
        print("  serve as a reference: a wrong peptide will match it as well as the right")
        print("  one. Do not read a pattern comparison against it in either direction.")
    for name in unassessable:
        print("  %s: no proline-containing peptide, so the strongest marker had no" % name)
        print("  opportunity and this source CANNOT BE ASSESSED. That is not a pass.")
        print("  Supply peptides that can express the effect, or judge this source on")
        print("  the aspartate column and its y-share with that limitation stated.")
    if not flat and not unassessable:
        print("  Every source shows the proline effect above %.2fx, so each carries real"
              % UNIFORM_TOL)
        print("  fragmentation chemistry and a comparison between them can be read.")
    return 1 if (flat or unassessable) else 0


if __name__ == "__main__":
    sys.exit(main())
