# Who can be told apart from whom

In all-ion fragmentation nothing is isolated — every precursor fragments at
once. A fragment carries its own m/z plus its **parent's** retention time and
drift. It does not carry the parent's m/z or charge. So assigning a fragment to
a parent needs the parents resolved in **(rt, drift) alone**.

Precursor m/z still matters, but indirectly, and the mechanism is worth stating
because it is easy to get backwards. At fixed charge, drift tracks m/z closely —
Spearman +0.94 at 2+, +0.93 at 3+, +0.92 at 4+ on the Colorado PAX8 runs. **That
correlation is what makes the separation work.** Two precursors apart in m/z get
a predictable, monotone drift difference, and the fragments inherit it. If drift
were orthogonal to m/z the relative distance would be random and collisions
would be *more* common, not fewer.

I initially read that 0.94 as redundancy — "a noisy copy of the first axis" —
and that was wrong. It is the mechanism.

## Relative distance, never absolute

A gap means nothing without the measurement scale beside it. The same 7.4 drift
units is 8.5 sigma on one instrument and 1.2 sigma on another, from identical
coordinates. `tools/resolvability.py` reports every separation in sigma and
refuses to run without a scale on both axes.

Getting the denominator right matters more than it sounds. My first pass divided
a 7.4-bin gap by the **population** SD across different peptides, 4.0 bins, and
concluded 1.84 sigma — "substantial overlap". That is the wrong quantity: the
spread of different peptides answers *can I predict drift from m/z*, not *can I
tell these two ions apart*. The right denominator is one precursor's own
repeated measurements:

| | |
|---|---|
| within-ion drift SD, same ion across runs | **0.87 bins** |
| 2+ at m/z 600 against 2+ at m/z 700 | gap **7.4 bins** |
| separation | **8.5 sigma** |

Cleanly resolved, not marginal. And 0.87 is an upper bound, because it includes
run-to-run calibration drift; within a single run the resolution is finer.

From the local slope, the engineering number: at fixed charge two precursors
resolve in drift once they are roughly **20–35 m/z apart** (22 at 2+, 30 at 3+,
34 at 4+, 23 at 6+).

## What the tool refuses

Three refusals, each from a failure seen in real data rather than imagined.

**No scale, no answer.** Without `--sigma-drift` and `--sigma-rt`, or
`--estimate-sigma` over data with repeated precursors, it exits rather than
reporting an absolute gap.

**A zero sigma is not a perfect instrument.** If an axis did not vary across the
repeats — a run-level constant, a coarsely rounded value — the pooled SD is zero
and dividing by it would make every pair infinitely separated and the run
flawless. It is refused with the reason. This was a live crash before the test
caught it.

**A sigma that is a large share of the axis range is a confound, not a
resolution.** Run on the PAX8 coordinates with raw minutes, `--estimate-sigma`
returns **σ_rt = 13.2 min, 19% of the whole range**, and is refused. Those runs
are 60- and 120-minute gradients mixed: the same peptide elutes at minute 28 in
one and minute 56 in the other, so the "scatter" is the gradient. `--rank-rt`
replaces rt with within-run elution rank scaled 0–100, which is MARCH's own rule
— minutes change between runs, the order does not — and the estimate becomes
usable.

## What it says about the PAX8 runs

51 distinct precursors, 1,262 observations, `--estimate-sigma --rank-rt`:

| | |
|---|---|
| σ_drift | 1.02 bins |
| σ_rt | 6.8 rank units (of 100) |
| uniquely resolved | **22 (43.1%)** |
| sharing with at least one other | **29 (56.9%)** |
| co-residents per precursor | median 1, p90 4, max 7 |

So for **one protein**, well over half the precursors share their (rt, drift)
cell with another at a 2-sigma threshold. Those fragments cannot be assigned by
coordinates alone — they are what a fragment-physics prior has to resolve, and
the percentage is the honest measure of how hard that job is.

A whole-proteome run is orders of magnitude denser. That is the binding
constraint, and it is a counting problem rather than a resolution problem: the
mechanism separates *pairs* at 8.5 sigma; what it cannot do is separate forty
things into the slots available.

Which points at the levers, in order: **reduce how many precursors share a
slot** — narrower quadrupole windows as SONAR does, or gradients that spread
elution — and then use the fragment-physics prior from
`tools/fragment_chemistry.py` to split what still collides.

## The caveats on these numbers

σ is estimated across runs, so it carries run-to-run calibration drift and the
occupancy above is **pessimistic** — the real within-run resolution is finer and
collisions fewer. Getting the optimistic bound needs peak widths from a single
run, which this data does not carry.

And the drift calibration has a known per-day offset: a leave-one-out fit on a
run's own ions lands at 1.69% median error while the previous day's calibration
lands at 7.41%, with all 618 errors negative. That offset inflates the cross-run
σ_drift here, which is part of why 0.87–1.02 bins is an upper bound rather than
the instrument's real resolving power.

## Meaning what? Sharing a cell is a candidate list, not an error

56.9% sharing an (rt, drift) cell sounds alarming and mostly is not, because a
fragment's **own m/z** must match a b/y ion of whichever parent produced it, and
two different sequences mostly have different ladder masses. Over the 29
co-resident PAX8 precursors, b and y at 1+ and 2+, 20 ppm:

| | |
|---|---|
| ladder ions examined | 4,708 |
| ions whose m/z also matches a co-resident's ladder | **909 — 19.3%** |
| per-precursor ambiguous share | median **3.8%**, p90 100%, max 100% |

So of the 56.9% sharing a cell, only about a fifth of their fragments are
genuinely ambiguous. For a typical co-resident precursor, **96% of its fragments
are assigned by their own mass** regardless of what else is in the cell.

### And the remainder splits perfectly in two

The median-3.8%-but-p90-100% shape is not noise, it is two populations. Of 36
co-resident pairs:

| pair type | n | median ladder clash |
|---|---|---|
| **nested** — one sequence contains the other | 10 | **50%** |
| unrelated sequences | 26 | **2%** |

Two unrelated peptides sharing a cell share 2% of their fragments. Essentially
solved. Two *nested* peptides share **half**:

```
MDDSDQDSCRLSIDSQSSSSGPR   vs   KMDDSDQDSCRLSIDSQSSSSGPR    52%
MDDSDQDSCRLSIDSQSSSSGPR   vs   MDDSDQDSCRLSIDSQSSSSGPRK    52%
SLSPGHTLIPSSAVTPPESPQSDSLGSTYSINGLLGIAQPGSDK  vs  …GSDKRK  100%
```

**No amount of resolution fixes that**, because they are not two things that
happen to overlap — they are the same sequence plus a residue, so they genuinely
produce the same fragments. A y-ladder from the C-terminus is *identical* between
a peptide and its N-terminally extended missed-cleavage form.

### So the lever is not where crowding suggests

**Crowding is not the enemy. Nesting is.** Narrower quadrupole windows, longer
gradients and better drift resolution all act on the 2% case, which was already
fine. They do nothing for the 50% case.

The 50% case needs a different handle: **assign on the fragments the family does
not share.** For a nested pair that means the b ions spanning the differing
terminus, and the precursor mass itself. The shared fragments are unassignable
intensity, and should be treated as such rather than counted as evidence for
either member.

This bears directly on the reversed-ladder trim in `replication-fdr.md`. Removing
reversed b(L−1) and exact coincidences strips ions near a terminus — which is
precisely where a nested family's *distinguishing* information lives. Trimming
one side there is worse than a uniform handicap: it removes the ions that do the
discriminating work.
