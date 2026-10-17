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
import re
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


def a_source_with_no_proline_cannot_be_assessed(tmp):
    """No proline-containing peptide means the marker had no opportunity.

    That is not a pass. Before the fix, enrichment came back None, the source
    was left out of the flat list, and a literally uniform source printed
    "Every source shows the proline effect" and exited 0.
    """
    f = tmp / "nopro.jsonl"
    rows = []
    for seq in ("AAGSTVLK", "SDEGTVAR", "GGSSTTAAK"):     # no P anywhere
        for kind in "by":
            for n in range(1, len(seq)):
                rows.append({"source": "no_proline", "sequence": seq, "charge": 2,
                             "fragment": "%s%d" % (kind, n), "intensity": 1.0})
    write(f, rows)
    rc, out = run(str(f))
    bad = []
    if rc != 1:
        bad.append("exit %d: treated an unassessable source as a pass" % rc)
    if "CANNOT BE ASSESSED" not in out:
        bad.append("did not say the source cannot be assessed")
    if "Every source shows the proline effect" in out:
        bad.append("printed the all-clear for a source with no proline opportunity")
    return bad


def a_record_missing_its_fragment_does_not_abort(tmp):
    f = tmp / "broken.jsonl"
    rows = plant("march", 4.0, 2.0)
    rows.append({"source": "march", "sequence": "ATLTPSNTPLGR", "charge": 2,
                 "intensity": 5.0})                        # no `fragment`
    rows.append({"source": "march", "sequence": "ATLTPSNTPLGR", "charge": 2,
                 "fragment": None, "intensity": 5.0})
    write(f, rows)
    rc, out = run(str(f))
    bad = []
    if "Traceback" in out:
        bad.append("crashed on a record missing `fragment`")
    if "carries real" not in out:
        bad.append("two malformed lines lost an otherwise good source")
    if "unparseable lines skipped" not in out:
        bad.append("did not count the malformed lines")
    return bad




def bond_level_input_withholds_the_y_share(tmp):
    """b and y summed per bond: the bond markers still read, the y-share cannot.

    Real case, from the project's own FITVCDYTNPCTK profiles: intensity is
    published per cleavage position with both ion series summed. Reporting a
    y-share of 0.000 there would say the source has no y ions, when it has no
    column for them.
    """
    f = tmp / "bonds.jsonl"
    seq = "SAMPLERPEPTIDEKPAK"          # P at 8 and 16, so bonds 7 and 15
    rows = [{"source": "bondlevel", "sequence": seq, "bond": k,
             "intensity": 4.0 if seq[k] == "P" else 1.0}
            for k in range(1, len(seq))]
    write(f, rows)
    rc, out = run(str(f))
    bad = []
    line = [l for l in out.splitlines() if l.startswith("bondlevel")]
    if not line:
        return ["no row for the source"]
    if "merged" not in line[0]:
        bad.append("y-share is %r; it must read merged, not a number"
                   % line[0].split()[2])
    if " 0.000" in line[0]:
        bad.append("reported a y-share of zero for a source with no y column")
    if "not zero" not in out:
        bad.append("did not explain why the y-share is withheld")
    if "cannot be checked" not in out:
        bad.append("did not say the beam-type y-over-b rule is uncheckable here")
    # the bond markers must still work: a 4x proline boost should show
    enr = float(line[0].split("=")[1].split("x")[0])
    if enr < 2.0:
        bad.append("proline enrichment %.2fx from bond-level rows; the markers"
                   " should be unaffected by the merge" % enr)
    return bad


def a_bond_index_outside_the_peptide_is_skipped(tmp):
    f = tmp / "badbond.jsonl"
    seq = "PEPTIDEK"                     # bonds 1..7
    write(f, [{"source": "s", "sequence": seq, "bond": k, "intensity": 1.0}
              for k in (0, 7, 8, 99)])
    rc, out = run(str(f))
    bad = []
    if "none usable" not in out and "No usable" not in out:
        # bond 7 is valid, 0/8/99 are not; so one fragment should survive
        line = [l for l in out.splitlines() if l.startswith("s ")]
        if line:
            n = int(line[0].split()[1])
            if n != 1:
                bad.append("%d fragment(s) counted; only bond 7 is inside an"
                           " 8-residue peptide" % n)
    return bad


def one_marker_site_cannot_decide_either_way(tmp):
    """An enrichment over one bond IS that bond's intensity.

    This is the real FITVCDYTNPCTK case: one proline, one aspartate, twelve
    bonds. Calling the source flat from that is a verdict from n=1, and calling
    it chemical would be the same error with the opposite sign.
    """
    seq = "FITVCDYTNPCTK"
    raw = [3725720, 5827492, 3630285, 2797988, 7300418, 8855138, 9104065,
           7460752, 5041150, 4108453, 2956913, 826324]
    f = tmp / "onesite.jsonl"
    write(f, [{"source": "real", "sequence": seq, "bond": k + 1,
               "intensity": float(v)} for k, v in enumerate(raw)])
    rc, out = run(str(f))
    bad = []
    if rc == 0:
        bad.append("exit 0: a source that cannot be assessed read as a pass")
    if "CANNOT BE ASSESSED" not in out:
        bad.append("did not say the marker cannot be assessed")
    if "not a pass" not in out or "not a failure" not in out:
        bad.append("did not rule out both readings")
    if "carries no fragmentation chemistry" in out:
        bad.append("declared the source flat from a single marker site")
    if "1 marker-bearing bond" not in out:
        bad.append("did not say how many sites the ratio rests on")
    if "span an order of magnitude" not in out:
        bad.append("did not say why one bond cannot carry the ratio")
    line = [l for l in out.splitlines() if l.startswith("real")]
    if line and not re.search(r"\b1/1\b", line[0]):
        bad.append("the site count is not in the table: %r" % line[0].strip())
    return bad


def enough_sites_lets_the_verdict_through(tmp):
    """The floor must not swallow a source that can express the effect."""
    f = tmp / "many.jsonl"
    write(f, plant("many", proline_boost=4.0, y_boost=2.0))
    rc, out = run(str(f))
    bad = []
    if "CANNOT BE ASSESSED" in out:
        bad.append("the site floor blocked a source with several proline sites")
    if rc != 0:
        bad.append("exit %d on a source with real chemistry and enough sites" % rc)
    line = [l for l in out.splitlines() if l.startswith("many")]
    if line:
        counts = re.findall(r"\b(\d+)/(\d+)\b", line[0])
        if not counts:
            bad.append("no sites/peptides count in the row: %r" % line[0].strip())
            return bad
        n = int(counts[0][0])
        if n < 3:
            bad.append("only %d site(s) in the planted set; the fixture cannot"
                       " test the floor" % n)
    # and a flat source with enough sites must still be caught
    g = tmp / "flatmany.jsonl"
    write(g, plant("flatmany", proline_boost=1.0, y_boost=1.0))
    rc2, out2 = run(str(g))
    if rc2 == 0:
        bad.append("a flat source with enough sites was not refused")
    if "carries no fragmentation chemistry" not in out2:
        bad.append("a flat source with enough sites was not named flat")
    return bad


CASES = [bond_mapping_is_right,
         a_source_with_no_proline_cannot_be_assessed,
         a_record_missing_its_fragment_does_not_abort, a_real_pattern_is_credited, a_flat_pattern_is_refused,
         a_real_and_a_flat_source_are_told_apart,
         unparseable_fragments_are_counted_not_guessed,
         bond_level_input_withholds_the_y_share,
         a_bond_index_outside_the_peptide_is_skipped,
         one_marker_site_cannot_decide_either_way,
         enough_sites_lets_the_verdict_through]


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
