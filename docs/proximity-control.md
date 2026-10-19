# What is the false-positive rate of the `bare_proximity` window?

`bare_proximity` is the tier every ion falls to when nothing decoy-controlled
scored it: a ppm window, a charge, and enough consecutive scans. It admits more
ions than the other two tiers combined. Until now it carried one sentence
wherever it appeared — in the schema, the validator, the renderer, the run spine
and two documents — and that sentence said the same bare mechanism *passed
reversed-tryptic decoys at 33–70% elsewhere in this project*.

Two things were wrong with it, and they pull in opposite directions.

**The number was not from this repository.** All eight citation sites attributed
it to "elsewhere in this project". Nothing here derived it.

**The mechanism it indicted no longer runs.** From `export-contract.md`: *"The
retired cube pathway (`match_targets_spatial_in_cube`) is where depletion was
first described, but it no longer runs."* The 33–70% describes that pathway. The
current one is MARCH.

So the tier that admits the most ions carried a borrowed number about retired
code, and the false-positive rate of what actually runs had never been measured.
That is worse than a wrong number: it is a confidence label with no measurement
behind it.

## Reversed decoys could never have measured this window

This matters because it is where the borrowed figure came from. Reversing a
tryptic peptide permutes its residues, so it **preserves the monoisotopic mass
exactly**: a reversed decoy has the same m/z as the peptide it decoys, at every
charge. Queried through an m/z window it matches 100% by construction, not 33%
or 70%.

Whatever the 33–70% measured — fragment matching, or a pathway using more than
m/z — it was not a precursor m/z window. The same argument rules out shuffled
and permuted decoys, and it is the reason this control uses **entrapment**: real
peptides from an organism that cannot be in the sample, which have different
masses on the same peptide mass-defect lattice.

(The clustering case in `tools/test_proximity_control.py` uses that property on
purpose: four permutations of one sequence share an m/z to the last decimal, so
one peak answers all four probes.)

One thing the control deliberately does **not** do is enforce a plausibility
band. An earlier draft flagged a measured rate landing far below the geometric
expectation. That expectation assumes peak windows are spread evenly and real
peaks cluster, so no fixture separates "the probes are suppressed" from "the
peaks sit away from the probes" without building the clustering in by hand. The
expectation is still computed and printed, because it tells a reader whether
clearing the bar means anything at all — when chance alone is 0.03%, a measured
0% is unremarkable — but nothing hangs on it.

## What is measured

`tools/proximity_control.py` queries probe coordinates through the **same**
`match_one()` the annotation uses — factored out of `stage_annotate` for exactly
this reason, because a second implementation would measure a different window
and the rate would then describe code that never ran on the data.

    python3 tools/proximity_control.py peaks.jsonl targets.csv \
        --entrapment ecoli_tryptic.txt --precursor-ppm 10 --rt-tol 0.5

It is stage 3b of `tools/aria_run.py`, via `--entrapment`. Four arms, because one
number cannot hold them:

| arm | what it answers |
|---|---|
| `target` | the comparator. **Not** the rate `aria_run` prints for the match: that loop is peak-direction (for each peak, scan targets) and this is probe-direction (for each probe, find a peak). Different denominators, so the comparison is only meaningful computed in the same direction. |
| `entrapment` | the primary arm. Peptides that cannot be in the sample, at their own calculated m/z, on a real target's charge and coordinates. |
| `calc_target` | the real targets with m/z **recomputed from their sequences**. The diagnostic for the calculated-versus-observed bias below. |
| `swap` | target A's m/z on target B's retention time and drift. The only arm that asks whether rt and drift contribute anything beyond m/z. Emitted only when one of them is actually compared. |

Each arm reports `k`, `n`, a Wilson interval and a Katz ratio against the
comparator. The counts are stored, not just the rate: a rate without its
denominator cannot be given an interval or audited, and this project has already
had to correct one derivation that hung on an unstated denominator
(`replication-fdr.md`, 100 → 93).

The Katz ratio is not a nicety. **A probe rate alone cannot tell a good window
from a thin peak table** — on a sparse run every arm reports zero. It is
propagated as null, never as infinite discrimination, when either numerator is
zero.

## Three ways this could have been a check that cannot fail

Every way of breaking this control drives the measured rate toward **zero**. That
is the whole problem: a rule of the form *fail when the rate exceeds a cut* would
pass most easily when the control is broken, and would reward all of the
following.

**Calculated versus observed m/z.** `make_targets.py` joins a search to a
feature table precisely so a target carries an *observed* m/z, because a search
reports a calculated one. An entrapment probe has no observation, so its m/z
must be calculated. This project has a documented live offset — **1.570 ppm, 24
of 24 runs below zero** (`atlas-convergence.md`). Against a 10 ppm window that
offset alone sits the probe window off-centre while the target window stays
centred, and the probe rate collapses for a reason that has nothing to do with
the window. So the offset is **measured** from the targets that carry a sequence
(observed minus calculated, as a median), **applied** to every probe, and
**reported**; the `calc_target` arm then bounds what correction is left. If the
target rate falls when the same targets are queried on calculated m/z, the
entrapment rate is a floor, not an estimate, and the report says so.

**Axes that are carried but never matched.** `--rt-tol` has no default, and the
shipped `examples/pax8-targets.csv` carries `rt_lo` and `rt_hi` — neither of
which is an `rt`. In both cases a probe given "a real target's retention time"
is given a field nothing compares, and the measured quantity collapses from a
rate at an occupied address to a **global m/z-crowding rate**. Every arm
therefore reports `axes_matched`, and no arm claims a coordinate it did not
match on.

**A denominator padded with impossible probes.** A probe whose m/z falls outside
the acquired range, or whose charge no peak carries, contributes one to `n` and
zero to `k`. Those are excluded and counted, not warned about. Probe charge is
taken from the paired target rather than drawn, so a charge-distribution
mismatch cannot depress the rate at all.

A probe mass-degenerate with a listed target is also excluded: it is testing
that target, and counting it as a false positive is circular. Entrapment
sequences are **not** screened against the human proteome by mass — that would
remove the hardest probes and drive the rate to a floor while looking like
rigour. A probe colliding with a genuine co-eluting peptide *is* the mechanism by
which a real annotation is wrong.

## What the export carries, and what refuses it

- `manifest.proximity_control` — the cohort result: every arm with its counts,
  intervals, exclusions, the measured calibration offset, and the findings the
  control raised about itself. Absent or null means the rate was **not
  measured**, which is a different claim from a measured zero.
- `confidence.control` on every annotated `bare_proximity` ion — that ion's own
  **m/z stratum**, not the cohort rate. Crowding is not uniform across the m/z
  range, and one cohort number stamped on thousands of ions is the same constant
  repeated.

`decoy_controlled` stays **false** for the tier regardless. Knowing how often
the window is wrong is a property of the window, not of any one ion, and
`validate_export` is right to error on an ion that claims otherwise.

`validate_export` rejects an export when:

1. an annotated `bare_proximity` ion carries no `control`;
2. a `control` carries no `queried`/`matched`, or its rate cannot be reproduced
   from them;
3. a `control` names no `axes_matched`;
4. the manifest carries no `proximity_control` while annotated
   `bare_proximity` ions rest on it;
5. there is no target comparator in the probe direction;
6. fewer than 30 probes were queried — at n = 20, k = 0 the Wilson upper bound
   is still about 16%, which cannot certify anything;
7. no real target matched either;
8. **the entrapment interval overlaps the target interval.** A bar a decoy
   clears is not a bar, whichever way the point estimates happen to fall.

The requirement is scoped to ions that **claim an identity**. An unassigned
detection with no annotation — an unresolved isotope envelope with no determined
charge, which matched nothing — makes no claim a decoy could be wrong about. The
contract tracks those on purpose, and demanding a false-positive rate for a
match that was never asserted would be a category error.

## What the shipped target list reports

On `examples/pax8-targets.csv` with a synthetic peak table at the listed
coordinates, the control **refuses to certify**, and says why:

    axes actually compared: mz, charge
    calibration: observed - calculated = -0.79 ppm (n=45, spread 2.63)

    arm           probes     hit     rate             95% CI
    target            45      45   100.0%       92.1%-100.0%
    entrapment        17       0     0.0%         0.0%-18.4%

    Chance alone, from peak-window geometry on this run: 0.02886% --
    about 1 probe in 3,465.

    DISCRIMINATED: NO -- 17 probes is too few for the bound to mean anything

Two real facts about that list fall out of it. Retention time is carried and
never compared, because `rt_lo`/`rt_hi` are not an `rt` — so even with
`--rt-tol` the rate would be global m/z crowding. And 25 of 42 probes could not
be constructed at all — every one of them outside the m/z range the 45 listed
coordinates span, which is why only 17 were queried.

## Building the list

`tools/make_entrapment.py` derives it from a proteome you supply and records in
the file's own header what it derived it from — source filename, sha256, protein
count, peptide counts at each stage, the screen, the sample size and the seed.
A list of sequences somebody believed were *E. coli* would be the same defect as
the borrowed 33–70%: a warrant living somewhere else.

    python3 tools/make_entrapment.py ecoli.fasta -o entrapment.txt \
        --exclude-fasta human.fasta --n 2000

`--exclude-fasta` takes the **sample's** proteome and drops every probe whose
exact sequence also appears in it. Those are not entrapment at all — a tryptic
peptide shared residue-for-residue between two organisms really can be in the
sample, so counting a hit on one as a false positive would be wrong.

It does **not** screen by mass, and that is deliberate. A probe whose mass
coincides with a real sample peptide is the most informative probe there is: for
a ppm window, collision-with-a-real-peptide and collision-with-noise are the
same event and produce the same wrong annotation. Screening those out would
drive the measured rate toward a floor while looking like rigour.

It refuses a non-FASTA input, a screen that leaves nothing, and a proteome that
digests to nothing; and it exits 1 rather than 0 on a list below the control's
own 30-probe floor, because probes are lost again at run time to the acquired
m/z range and the charges present.

Sizing, measured: a 2,000-peptide list against the end-to-end fixture leaves
1,951 probes after construction and gives a 95% interval of 0.0–0.2%. At n=20,
k=0 the upper bound is still about 16%, which is why the floor exists. Note that
probes from one digest share the mass-defect lattice and several can match the
same peak, so the trials are not fully independent and the Wilson interval is
mildly anti-conservative; `distinct_peaks_matched` is recorded beside `k` so
clustering is visible, and the report says so when `k` exceeds it.

## Scale

Probe-direction matching is O(probes × peaks) written naively: 2,000 probes
against 50,000 peaks measured 44 seconds **per arm**, and a real export carries
far more peaks than that. Peaks are therefore bucketed by charge and sorted by
m/z, and each probe's candidates come from a bisect over that. The same figures
with the index: 0.03s of query at 50,000 peaks, 0.6s total at 500,000.

This is the one optimisation in the control that could do real damage, because a
candidate filter that is slightly too tight silently narrows the window being
measured — and would report the narrowing as a better false-positive rate. The
bound is exact rather than approximate (`match_one` divides by the target's m/z,
and in the probe direction the probe *is* the target), `match_one` still makes
every final decision including charge, retention time and drift, and a test
compares the indexed path against a brute-force scan over every peak —
randomised, with probes planted on, just inside and just outside the window
edge, with several peaks inside one window, and with chargeless probes that must
search every bucket. Same `k`, same distinct-peak count, same peak chosen for
every probe.

At 500,000 synthetic peaks, incidentally, 156 of 2,000 probes found one: a 7.8%
chance-collision rate. On a real dense run this number will not be small.

## Merging runs

`merge_exports.py` pools the arms across runs — the probes are independent
Bernoulli trials against one window definition, so `k` and `n` add, and the
pooled interval is correspondingly tighter. It pools **only** when every run
used the same `precursor_ppm`, `rt_tol`, `drift_tol`, `min_scans` and the same
`axes_matched`. Runs that differ on any of those measured different windows, and
their sum would be a rate for a window nobody ran: the merged manifest then
carries no control and names the runs that disagreed. A run that never measured
a rate has the same effect — a cohort rate assembled only from the runs that did
measure one would report their window as the cohort's, which is how the borrowed
33–70% got here in the first place.

Per-run strata are **not** pooled. The strata are per-run m/z quartiles and do
not align across runs, so a pooled stratum would be a different partition from
any ion's own; `per_run` keeps each run's counts and calibration instead.

## What is not established

**The rate has not been measured on real data.** Everything above is the
machinery and its refusals, exercised against planted fixtures: 18 cases in
`test_proximity_control.py`, 8 in `test_validate_export.py`, 8 in
`test_make_entrapment.py` and 2 in `test_merge_exports.py` — 280 in the suite
against a 244 baseline — each one shown to fail against a deliberately broken
implementation before it was kept, 40 mutations in all.

Four of those cases could not fail as first written, which is the defect this
whole exercise is about. An 8 ppm planted calibration offset still lands inside
a 10 ppm window, so the fixture could not tell a tool that applies the
correction from one that ignores it. A clustering case used a 200000 ppm window,
which put every probe inside its own paired target's window, excluded them all
as self-matches, and passed on zero trials against code with no clustering check
at all. And twice a fixture protein digested to nothing because trypsin does not
cut before proline: once a peptide the case claimed two proteomes shared was cut
through in both, once a `K|P` site was the only internal one. All four now assert
that the fixture produced the condition under test.

`tools/mutation_check.py` holds all 40 and runs them on demand (81 seconds; it
runs a suite per mutation, so it is not wired into the per-file suites). Three
properties make it worth having rather than a record in a commit message. A
mutation whose anchor text no longer appears exactly once is reported as a
FAILURE, not skipped — rewording a line detaches it, which happened twice while
this was being written, and detached it silently proves nothing. The baseline is
run unmutated first and a red one stops the run, because against a failing suite
every mutation is reported as caught. And every file is restored from memory in
a `finally` block and compared byte for byte afterwards.

The uniqueness check earned itself immediately: `per-ion control requirement
disabled` anchored on a line that appears twice in `validate_export.py`, and the
ad-hoc script it replaced had been mutating whichever came first.

The number the control produces on a real run is unknown until it runs on one,
and
`atlas-convergence.md` is the standing warning about what a confidently
reported number is worth when nobody asked what else could produce it.

**No entrapment list ships, so every existing export fails validation until
one exists.** `method.schema.json` has carried the default `"ecoli"` since
before this work. `make_entrapment.py` builds the list, but it needs a proteome
FASTA, and fetching one here was refused by this environment's network policy
(`rest.uniprot.org` denied). Nothing was invented to stand in for it. Supply an
*E. coli* FASTA — UniProt reference proteome UP000000625 — and the list and the
end-to-end example follow in one command.

**The organism has to suit the sample.** *E. coli* peptides are entrapment for a
HeLa lysate. They are not entrapment for a sample containing recombinant protein
expressed in *E. coli*, where that organism really can be present; such a sample
needs a different one, and `--exclude-fasta` will not save you, because it
screens sequences against the proteome you name, not against what is actually in
the tube.

**The `swap` arm is the one to watch.** On the end-to-end fixture it reports 0 of
35 — target A's m/z at target B's retention time finds nothing — which says rt
is doing real work there. On a real run with a dense peak table it may not, and
that result would say the window is narrower than its description.
