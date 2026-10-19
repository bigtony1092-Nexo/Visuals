#!/usr/bin/env python3
"""Build an entrapment peptide list from a proteome, with its provenance.

`tools/proximity_control.py` needs peptides that CANNOT be in the sample. The
one thing it must not have is a list of sequences somebody believed were E.
coli: that would be the same defect as the borrowed 33-70% this whole control
replaced -- a number, or here a sequence set, whose warrant lives somewhere
else. So this tool derives the list from a FASTA you supply and records what it
derived it from, hashes included, in the file's own header.

    python3 tools/make_entrapment.py ecoli.fasta -o entrapment.txt \\
        --exclude-fasta human.fasta --n 2000

Two screens, and the difference between them matters.

`--exclude-fasta` takes the SAMPLE's proteome and drops every probe peptide
whose exact sequence also appears in it. Those are not entrapment at all: a
tryptic peptide shared residue-for-residue between E. coli and human really can
be in a human sample, so counting a hit on one as a false positive would be
wrong. Short peptides are where this concentrates, which is also why the
default minimum length is 7.

What this does NOT do is screen by MASS. A probe whose mass happens to coincide
with a real sample peptide is the hardest and most informative probe there is:
for a ppm window, collision-with-a-real-peptide and collision-with-noise are the
same event and produce the same wrong annotation. Screening those out would
drive the measured rate toward a floor while looking like rigour. The control
excludes a probe only when it is degenerate with a listed TARGET, which is
circular for a different reason, and it reports that count.

Exit 2 if no list could be built, 1 if one was built but is too small to
certify anything, 0 otherwise.
"""

import argparse
import hashlib
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from proximity_control import digest_fasta  # noqa: E402
from validate_export import peptide_mass, MIN_CONTROL_PROBES  # noqa: E402

MIN_LEN, MAX_LEN = 7, 30
DEFAULT_N = 2000


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def n_proteins(text):
    return sum(1 for line in text.splitlines() if line.startswith(">"))


def massable(seqs):
    """Only peptides whose mass this project can compute.

    peptide_mass returns None on any residue outside the table -- B, X, Z and U
    appear in real proteomes -- and a probe with no computable mass has no
    coordinate, so it would pad the denominator with a trial that can never
    match.
    """
    return [s for s in seqs if peptide_mass(s) is not None]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta", help="proteome of an organism that CANNOT be in the sample")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--exclude-fasta",
                   help="the SAMPLE's proteome; probes whose exact sequence appears"
                        " in it are not entrapment and are dropped")
    p.add_argument("--n", type=int, default=DEFAULT_N,
                   help="how many peptides to keep (default %d)" % DEFAULT_N)
    p.add_argument("--seed", type=int, default=0,
                   help="sampling seed, recorded in the output so the list is"
                        " reproducible (default 0)")
    p.add_argument("--min-len", type=int, default=MIN_LEN)
    p.add_argument("--max-len", type=int, default=MAX_LEN)
    args = p.parse_args(argv)

    text = open(args.fasta).read()
    if not text.lstrip().startswith(">"):
        print("%s does not look like a FASTA: it does not start with '>'."
              % args.fasta, file=sys.stderr)
        return 2
    n_prot = n_proteins(text)
    peptides = [s for s in digest_fasta(text)
                if args.min_len <= len(s) <= args.max_len]
    n_digested = len(peptides)
    peptides = massable(peptides)
    n_massable = len(peptides)
    if not peptides:
        print("no tryptic peptide of %d-%d residues with a computable mass came out"
              " of %s (%d protein(s))."
              % (args.min_len, args.max_len, args.fasta, n_prot), file=sys.stderr)
        return 2

    n_shared = 0
    exclude_note = "none -- NOT screened against the sample's own proteome"
    if args.exclude_fasta:
        other = open(args.exclude_fasta).read()
        if not other.lstrip().startswith(">"):
            print("%s does not look like a FASTA." % args.exclude_fasta, file=sys.stderr)
            return 2
        shared = set(digest_fasta(other))
        before = len(peptides)
        peptides = [s for s in peptides if s not in shared]
        n_shared = before - len(peptides)
        exclude_note = "%s (sha256 %s, %d protein(s))" % (
            os.path.basename(args.exclude_fasta), sha256(args.exclude_fasta),
            n_proteins(other))
        if not peptides:
            print("every candidate peptide also occurs in %s, so none of them is"
                  " entrapment." % args.exclude_fasta, file=sys.stderr)
            return 2

    # Sorted before sampling, so the same inputs and seed give the same list on
    # any machine: dict ordering from the digest must not reach the output. This
    # is defensive -- CPython's dict preserves insertion order, so no test here
    # can distinguish its presence, and it is kept for the implementation that
    # does not.
    peptides.sort()
    rnd = random.Random(args.seed)
    keep = sorted(rnd.sample(peptides, min(args.n, len(peptides))))

    header = [
        "# Entrapment peptides: sequences that cannot be in the sample.",
        "# Generated by tools/make_entrapment.py. Do not hand-edit -- regenerate.",
        "#",
        "# source           %s" % os.path.basename(args.fasta),
        "# source sha256    %s" % sha256(args.fasta),
        "# source proteins  %d" % n_prot,
        "# excluded against %s" % exclude_note,
        "#",
        "# fully tryptic, no missed cleavages, %d-%d residues" % (args.min_len, args.max_len),
        "# digested         %d" % n_digested,
        "# with a mass      %d  (the rest carry B/X/Z/U or similar)" % n_massable,
        "# shared with the sample's proteome, dropped   %d" % n_shared,
        "# sampled          %d  (seed %d)" % (len(keep), args.seed),
        "#",
        "# NOT screened by mass. A probe whose mass coincides with a real sample",
        "# peptide is the most informative probe there is: for a ppm window that",
        "# collision and a collision with noise are the same event. Screening",
        "# them out would drive the measured rate toward a floor.",
    ]
    with open(args.out, "w") as fh:
        fh.write("\n".join(header) + "\n")
        for s in keep:
            fh.write(s + "\n")

    print("\n".join(l[2:] if l.startswith("# ") else l[1:] for l in header[3:]))
    print()
    print("wrote %d peptide(s) to %s" % (len(keep), args.out))
    if len(keep) < MIN_CONTROL_PROBES:
        print("\nThat is below the %d probes the control needs before its interval"
              " can support any claim, and probes are lost again at run time to"
              " the acquired m/z range and the charges present. This list cannot"
              " certify a window." % MIN_CONTROL_PROBES, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
