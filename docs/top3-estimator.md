# Top3avg is a biased estimator

Top3avg takes, per peptide, the **maximum** `matched_Intensity` across every scan it hit,
then averages the top three. The maximum is the problem.

## The bias

The maximum of n draws rises with n. A peptide observed in more scans scores higher for no
abundance reason at all. Simulated, 40,000 trials per row, true intensity 100, scan-to-scan
noise σ = 15:

| Scans observed | Mean of max | Inflation |
|---|---|---|
| 2 | 108.6 | +8.5% |
| 3 | 112.7 | +12.7% |
| 5 | 117.4 | +17.4% |
| 8 | 121.3 | +21.4% |
| 12 | 124.4 | +24.5% |
| 20 | 128.0 | +28.1% |
| 30 | 130.7 | +30.6% |

The inflation itself does not matter — it would cancel if every peptide were seen the same
number of times. It is not, so it does not cancel.

**Two peptides of identical true abundance, one seen 3 times and one 20 times:**

| Estimator | Ratio measured | Ratio that is true |
|---|---|---|
| max across scans | **1.136** | 1.000 |
| mean across scans | 1.000 | 1.000 |

A 13.6% difference invented by scan count alone, before any real difference is measured.

## Why it matters here

Scan count tracks elution width, acquisition method and gradient length — none of which is
abundance. Across a cohort acquired on different gradients (60 min against 120 min, say)
the inflation differs **systematically between samples**, not randomly.

That is a candidate contributor to an instability already documented: GAPDH, the textbook
stable reference, running from the 8th to the 96th percentile across one 27-sample cohort
on this measurement, with a cross-sample rank std of 991. The bias does not explain a swing
that large on its own, but it pushes in exactly that direction and it is removable.

## The fix

Use an integrated or averaged intensity per peptide instead of the maximum.
`summed_precursor_intensity` is already in the pipeline and already in the contract, so
nothing new has to be measured.

`method.quant.peptide_aggregation` now accepts:

| Value | What it does |
|---|---|
| `max` | the current behaviour, kept so old methods still describe themselves |
| `mean` | mean across scans — unbiased with respect to scan count |
| `integrated_area` | integrated chromatographic area — the field's usual choice, and the least noise-sensitive |

The selection rule is unchanged: sort peptides descending, average the top three or fewer,
never pad. Only the per-peptide number changes.

## What to check after switching

Recompute the GAPDH percentile series under `integrated_area` and compare with the
published max-based series. If the cross-sample std falls, part of that instability was the
estimator rather than the measurement, and the published figure for rank stability should
be restated. If it does not fall, the estimator is exonerated and the instability is real —
which is worth knowing with equal certainty.

Either outcome is publishable. The current state, where the estimator's contribution is
unknown, is the only one that is not.
