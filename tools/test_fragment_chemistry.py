#!/usr/bin/env python3
"""Prove the chemistry check tells a real pattern from a flat one.

Two sources are planted: one where intensity concentrates N-terminal to proline
the way beam-type fragmentation actually does, and one drawn flat. The check must
credit the first and refuse the second, because a flat source cannot serve as a
reference — which is exactly how the PAX8 Thermo archive failed.

    python3 tools/test_fragment_chemistry.py

No pytest, no dependencies. Exit status is 1 if any case fails.
"""

import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "fragment_chemistry.py"
sys.path.insert(0, str(ROOT / "tools"))
import fragment_chemistry as fc  # noqa: E402

PEPTIDES = ["ATLTPSNTPLGR", "LSIDSQSSSSGPR", "GEQGLYPLPLLNSTLDDGK",
            "VQQPFNLPMDSCVATK", "QNPTMFAWEIR", "MDDSDQDSCR"]


def write(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")


def plant(source, proline_boost, y_boost, seed=3):
    rng = random.Random(seed)
    rows = []
    for seq in PEPTIDES:
        L = len(seq)
        for kind in "by":
            for n in range(1, L):
                k = n if kind == "b" else L - n
                if not 1 <= k <= L - 1:
                    continue
                inten = rng.uniform(0.5, 1.5)
                if seq[k] == "P":
                    inten *= proline_boost
                if kind == "y":
                    inten *= y_boost
                rows.append({"source": source, "sequence": seq, "charge": 2,
                             "fragment": "%s%d" % (kind, n), "intensity": inten})
    return rows


def run(*paths):
    p = subprocess.run([sys.executable, str(TOOL)] + list(paths),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def bond_mapping_is_right():
    """b(n) reports bond n; y(n) reports bond L-n."""
    bad = []
    seq = "PEPTIDEK"          # L = 8, so bonds 1..7
    for label, want in (("b1", 1), ("b7", 7), ("y1", 7), ("y7", 1), ("y3", 5),
                        ("b8", None), ("y8", None), ("a4", None), ("junk", None)):
        got = fc.bond_of(label, len(seq))
        if got != want:
            bad.append("bond_of(%r) = %r, expected %r" % (label, got, want))
    # The proline marker must fire on the bond N-terminal to P, not on P's own bond.
    s = "AAPAA"               # P is residue 3, so bond 2 is N-terminal to it
    if not fc.markers(s, 2)["proline"]:
        bad.append("proline marker missed the bond N-terminal to P")
    if fc.markers(s, 3)["proline"]:
        bad.append("proline marker fired on the wrong bond")
    # Aspartate: the bond C-terminal to D, i.e. bond k where residue k is D.
    s = "AADAA"               # D is residue 3, so bond 3 is C-terminal to it
    if not fc.markers(s, 3)["aspartate"]:
        bad.append("aspartate marker missed the bond C-terminal to D")
    return bad


def a_real_pattern_is_credited(tmp):
    f = tmp / "real.jsonl"
    write(f, plant("march", proline_boost=4.0, y_boost=2.0))
    rc, out = run(str(f))
    bad = []
    if rc != 0:
        bad.append("exit %d, expected 0 for a source with real chemistry" % rc)
    if "carries real" not in out:
        bad.append("did not credit a source showing the proline effect")
    line = [l for l in out.splitlines() if l.startswith("march")]
    if not line:
        return bad + ["no row for the source"]
    enr = float(line[0].split("=")[1].split("x")[0])
    if enr < 2.0:
        bad.append("proline enrichment %.2fx, expected well above 2 on a 4x boost" % enr)
    if " 0.6" not in line[0] and " 0.7" not in line[0]:
        bad.append("y-share did not reflect the planted 2x y boost: %r" % line[0].strip())
    return bad


def a_flat_pattern_is_refused(tmp):
    f = tmp / "flat.jsonl"
    write(f, plant("thermo_archive", proline_boost=1.0, y_boost=1.0, seed=9))
    rc, out = run(str(f))
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1 for a source indistinguishable from uniform" % rc)
    for want in ("no fragmentation chemistry", "cannot", "either direction"):
        if want not in out:
            bad.append("did not refuse it as a reference (%r missing)" % want)
    return bad


def a_real_and_a_flat_source_are_told_apart(tmp):
    a, b = tmp / "a.jsonl", tmp / "b.jsonl"
    write(a, plant("march", proline_boost=4.0, y_boost=2.0))
    write(b, plant("archive", proline_boost=1.0, y_boost=1.0, seed=9))
    rc, out = run(str(a), str(b))
    bad = []
    if rc != 1:
        bad.append("exit %d, expected 1 while one side is flat" % rc)
    if "march" not in out or "archive" not in out:
        bad.append("did not report both sources")
    if out.count("proline enrichment") != 1:
        bad.append("flagged %d sources, expected exactly the flat one"
                   % out.count("proline enrichment"))
    return bad


def unparseable_fragments_are_counted_not_guessed(tmp):
    f = tmp / "mixed.jsonl"
    rows = plant("march", 4.0, 2.0)
    rows += [{"source": "march", "sequence": "PEPTIDEK", "charge": 2,
              "fragment": frag, "intensity": 100.0}
             for frag in ("a4", "y99", "c3", "")]
    write(f, rows)
    rc, out = run(str(f))
    if "carries real" not in out:
        return ["a few unusable labels broke an otherwise good source"]
    return []


CASES = [bond_mapping_is_right, a_real_pattern_is_credited, a_flat_pattern_is_refused,
         a_real_and_a_flat_source_are_told_apart,
         unparseable_fragments_are_counted_not_guessed]


def main():
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for case in CASES:
            bad = case() if case is bond_mapping_is_right else case(tmp)
            failed += bool(bad)
            print("%s %s" % ("ok  " if not bad else "FAIL", case.__name__))
            for line in bad:
                print("     %s" % line)
    print("\n%d/%d passed" % (len(CASES) - failed, len(CASES)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
