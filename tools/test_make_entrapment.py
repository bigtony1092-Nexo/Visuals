#!/usr/bin/env python3
"""Prove the entrapment list is derived, screened, and reproducible.

The control this list feeds exists because the tier it guards carried a figure
whose warrant lived somewhere else. A list of sequences somebody believed were
E. coli would reintroduce exactly that, so the cases here check that the list
records what it came from, that it drops peptides which really can be in the
sample, and that it refuses rather than emit one too small to certify anything.

    python3 tools/test_make_entrapment.py
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "make_entrapment.py"
sys.path.insert(0, str(ROOT / "tools"))
import proximity_control as pc  # noqa: E402

# Two small proteomes whose tryptic digests deliberately overlap in exactly one
# peptide. SHARED is flanked by K on the left and ends in R in BOTH, so trypsin
# releases the identical sequence from each -- an earlier version of this fixture
# put SHARED inside a longer run and trypsin cut through it, so neither proteome
# ever produced the peptide the case claimed they shared.
SHARED = "LSIDSQSSSSGPR"
PROBE_ORG = [
    "MAAVLTSHGKPEPTIDERSAMPLEDGTVKLLGATVLDALRAVFVDLEPTVIDEVRK",
    "MTIGGGDDSFNTFFSETGAGKYMVIQGEPGAVIRNLDIERPTYTNLNRK",
    "MSLDIEALFRVGINYQPPTVVPGGDLAKEIVHIQAGQCGNQIGAKDVNAAIATIKR",
    "MFDLMYAKQLFHPEQLITGKTLKFVTVDDGKALGIDVTEVRWVDPNDPTKR",
    "MLMTEFDAHGRIVSNASCTTNCLAPLAKGYSFTTTAERAAAK" + SHARED + "VVVK",
]
SAMPLE_ORG = ["MQQQQQNVPLYKHIIRAQLVDMLGHRDSYVGDEAQSKK" + SHARED + "HHHHHHHHK"]


def fasta(tmp, name, seqs):
    p = Path(tmp) / name
    p.write_text("".join(">p%d desc\n%s\n" % (i, s) for i, s in enumerate(seqs)))
    return p


def run(*args):
    p = subprocess.run([sys.executable, str(TOOL)] + [str(a) for a in args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def peptides_in(path):
    return [l.strip() for l in Path(path).read_text().splitlines()
            if l.strip() and not l.startswith("#")]


def case_the_list_records_what_it_came_from(tmp):
    """A derived list with no provenance is a list somebody has to trust."""
    import hashlib
    bad = []
    src = fasta(tmp, "probe1.fasta", PROBE_ORG)
    out = Path(tmp) / "e1.txt"
    rc, txt = run(src, "-o", out, "--n", 500)
    if rc not in (0, 1):
        return ["exit %d: %s" % (rc, txt)]
    head = Path(out).read_text()
    want = hashlib.sha256(Path(src).read_bytes()).hexdigest()
    if want not in head:
        bad.append("the output does not carry the source FASTA's sha256")
    if "probe1.fasta" not in head:
        bad.append("the output does not name the source FASTA")
    for mark in ("digested", "sampled", "seed"):
        if mark not in head:
            bad.append("the output does not record %r" % mark)
    if "NOT screened by mass" not in head:
        bad.append("the output does not say it is not mass-screened, which is the"
                   " screen a reader is most likely to assume happened")
    return bad


def case_the_same_inputs_and_seed_give_the_same_list(tmp):
    """Reproducible, or the provenance header is decoration."""
    bad = []
    src = fasta(tmp, "probe2.fasta", PROBE_ORG)
    a, b, c = (Path(tmp) / n for n in ("a.txt", "b.txt", "c.txt"))
    run(src, "-o", a, "--n", 6, "--seed", 0)
    run(src, "-o", b, "--n", 6, "--seed", 0)
    run(src, "-o", c, "--n", 6, "--seed", 1)
    if peptides_in(a) != peptides_in(b):
        bad.append("the same seed gave two different lists")
    if not peptides_in(a):
        bad.append("the fixture produced no peptides at all")
    elif peptides_in(a) == peptides_in(c) and len(peptides_in(a)) < len(
            [x for x in pc.digest_fasta(Path(src).read_text())]):
        bad.append("a different seed gave an identical proper subset, so the seed"
                   " is not actually sampling")
    return bad


def case_a_peptide_that_is_in_the_sample_is_not_entrapment(tmp):
    """A tryptic peptide shared residue-for-residue with the sample's proteome
    really can be there, so a hit on it is not a false positive."""
    bad = []
    src = fasta(tmp, "probe3.fasta", PROBE_ORG)
    smp = fasta(tmp, "sample3.fasta", SAMPLE_ORG)
    out = Path(tmp) / "e3.txt"
    rc, txt = run(src, "-o", out, "--exclude-fasta", smp, "--n", 500)
    peps = peptides_in(out)
    if SHARED in peps:
        bad.append("%s occurs in both proteomes and was still written as a decoy"
                   % SHARED)
    if "shared with the sample's proteome, dropped   1" not in Path(out).read_text():
        bad.append("the one shared peptide was not counted in the header")
    # And without the screen it must be kept, or the case proves nothing.
    out2 = Path(tmp) / "e3b.txt"
    run(src, "-o", out2, "--n", 500)
    if SHARED not in peptides_in(out2):
        bad.append("%s is absent even without --exclude-fasta, so the screen is"
                   " not what removed it" % SHARED)
    return bad


def case_unmassable_residues_are_dropped(tmp):
    """B, X, Z and U appear in real proteomes. A probe whose mass cannot be
    computed has no coordinate, so it would pad the denominator with a trial
    that can never match."""
    bad = []
    # The K is followed by A, not P: trypsin does not cut before proline, and an
    # earlier version of this fixture put K|P at the only internal site, so the
    # protein digested to nothing and the case tested nothing.
    src = fasta(tmp, "probe4.fasta",
                PROBE_ORG + ["MXXXBBBZZZSAMPLEUUUKAAAPEPTIDEXXXR"])
    out = Path(tmp) / "e4.txt"
    run(src, "-o", out, "--n", 500)
    peps = peptides_in(out)
    for pep in peps:
        if any(c in pep for c in "BXZU"):
            bad.append("%r carries a residue with no mass in this project" % pep)
    # And the fixture must actually have produced such peptides, or a tool with
    # no filter at all would pass this case.
    raw = pc.digest_fasta(Path(src).read_text())
    if not any(c in r for r in raw for c in "BXZU"):
        bad.append("the fixture digested to no unmassable peptide, so this case"
                   " cannot distinguish a tool that filters from one that does not")
    if "with a mass" not in Path(out).read_text():
        bad.append("the count of massable peptides is not recorded")
    return bad


def case_a_list_too_small_to_certify_says_so(tmp):
    """Below the control's own floor, and probes are lost again at run time."""
    bad = []
    src = fasta(tmp, "probe5.fasta", PROBE_ORG[:1])
    out = Path(tmp) / "e5.txt"
    rc, txt = run(src, "-o", out, "--n", 3)
    if rc != 1:
        bad.append("a 3-peptide list exited %d, not 1" % rc)
    if "cannot certify" not in txt:
        bad.append("a list too small to certify a window did not say so")
    return bad


def case_a_file_that_is_not_a_fasta_is_refused(tmp):
    bad = []
    p = Path(tmp) / "notfasta.txt"
    p.write_text("PEPTIDEK\nSAMPLEK\n")
    rc, txt = run(p, "-o", Path(tmp) / "e6.txt")
    if rc != 2:
        bad.append("a non-FASTA input exited %d, not 2" % rc)
    if "does not look like a FASTA" not in txt:
        bad.append("the refusal does not say what was wrong with the input")
    return bad


def case_everything_shared_is_refused_not_emptied(tmp):
    """If no candidate survives the screen there is no entrapment list, which is
    a refusal rather than an empty file read later as zero probes."""
    bad = []
    src = fasta(tmp, "probe7.fasta", PROBE_ORG)
    out = Path(tmp) / "e7.txt"
    rc, txt = run(src, "-o", out, "--exclude-fasta", src, "--n", 500)
    if rc != 2:
        bad.append("excluding a proteome against itself exited %d, not 2" % rc)
    if "none of them is entrapment" not in txt:
        bad.append("the refusal does not say why nothing survived")
    if out.exists():
        bad.append("an empty list was written anyway, which a later run would"
                   " read as a measurement")
    return bad


def case_the_control_reads_the_generated_file(tmp):
    """The round trip. The header is comments, and a comment read as a sequence
    would become a probe with no mass and a padded denominator."""
    bad = []
    src = fasta(tmp, "probe8.fasta", PROBE_ORG)
    out = Path(tmp) / "e8.txt"
    run(src, "-o", out, "--n", 500)
    want = peptides_in(out)
    got = pc.load_sequences(str(out))
    if sorted(got) != sorted(want):
        bad.append("proximity_control read %d sequence(s) from the generated"
                   " file, the file holds %d: %r"
                   % (len(got), len(want), sorted(set(got) ^ set(want))[:4]))
    if any(not s.isalpha() for s in got):
        bad.append("a header line was read as a sequence")
    return bad


CASES = [case_the_list_records_what_it_came_from,
         case_the_same_inputs_and_seed_give_the_same_list,
         case_a_peptide_that_is_in_the_sample_is_not_entrapment,
         case_unmassable_residues_are_dropped,
         case_a_list_too_small_to_certify_says_so,
         case_a_file_that_is_not_a_fasta_is_refused,
         case_everything_shared_is_refused_not_emptied,
         case_the_control_reads_the_generated_file]


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
