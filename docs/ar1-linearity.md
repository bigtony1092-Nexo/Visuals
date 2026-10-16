# Testing whether AR1 is actually loading-invariant

The contract says `ar1 = I_product / I_precursor` does not move with load, because both
terms scale with the amount injected. That is true **within the linear dynamic range**, and
only there. Detector saturation, ion suppression and dynamic-range compression all break
it, and they break it at the top of a loading series — exactly where a 50-fold range ends.

The claim was written as "loading-invariant by construction". It should read
"loading-invariant within the linear range", and the range should be measured rather than
assumed.

## The experiment already exists

The Colorado PAX8 series spans 50-fold on column with the same peptides throughout:

| On-column | Runs |
|---|---|
| 50 ng-equiv | JP213 100 ng, 0.5 µL |
| 100 ng-equiv | JP210 100 ng, 1 µL (×2 gradients) |
| 500 ng-equiv | JP210, JP211 500 ng, 1 µL |
| 1000 ng-equiv | JP212, JP213 500 ng, 2 µL |
| 2500 ng-equiv | JP214, JP215 500 ng, 5 µL |

Nothing new has to be acquired. `method.ar1-linearity.json` runs it.

## What to compute

For every ion seen across at least four load levels:

1. `ar1_sum` per run — the share of the precursor that converted to identified fragments.
2. Regress `log10(ar1_sum)` on `log10(on-column amount)`.
3. **The slope is the answer.** Zero means invariant. Negative means the precursor is
   saturating faster than its fragments, which is what detector compression looks like.

## Reading the result

| Slope (log-log) | Reading |
|---|---|
| within ±0.05 | invariant across the range tested — the contract's claim holds and is now measured |
| within ±0.05 up to some load, then negative | that breakpoint **is** the top of the linear range; quote it |
| negative throughout | ar1 is not loading-invariant here and the contract must be corrected, not footnoted |

Report the slope per ion and the distribution across ions, not a single pooled number. A
few saturating high-abundance ions can drag a pooled slope while most ions are flat — and
which ions bend is itself the useful result.

## Two controls

- **Precursor intensity alone** should rise with roughly unit slope against load. If it
  does not, the loading series is not doing what it says and nothing else from it means
  anything. Run this first.
- **A low-abundance ion** far from saturation should stay flat even if the high-abundance
  ones bend. Its slope is the floor the measurement can resolve.

## Why this is worth doing

It converts an assumption into a number, and the number is useful either way. If ar1 is
flat across 50-fold, the loading-invariance claim stops being an argument and becomes a
measurement. If it bends, you have found the linear range of your own instrument — which
every quantitative claim downstream depends on and which nobody has stated yet.
