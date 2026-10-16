# Verification run, 16 October 2026

Three claims recomputed from the PAX8 study's own embedded data — 1,262 Colorado PSMs
with sequence, charge, observed m/z, retention time, drift, both drift predictions and
intensity. Nothing here was taken from a report; each number was derived independently and
then compared with what was reported.

## 1. Mass identity reproduces exactly

Theoretical neutral mass computed from residue masses with carbamidomethyl fixed on every
Cys, then the residual explained as Met oxidation, Asn/Gln deamidation or an isotope pick
within 10 ppm.

| | Reported | Recomputed |
|---|---|---|
| PSMs | 1262 | 1262 |
| explained by mass | 1262 | **1262** |
| oxidised | 146 | **146** |
| deamidated | 20 | **20** |
| isotope picks | 33 | 1 (see below) |

Every PSM falls within 10 ppm of an explainable form, computed from first principles. This
is direct support for the claim that mass is computable and identical on every machine:
1,262 observed precursors, every one accounted for by sequence plus a named modification.

**The isotope-pick count does not reproduce, and the probe is at fault, not the number.**
The stored `mz` already carries the A1/A2 correction, and the ion keys carry it too for all
but one record. Testing it needs the raw reported precursor from before correction.

## 2. Drift calibration transfers per run, not per instrument

618 ions carrying drift:

| Calibration | Median absolute relative error | Reported |
|---|---|---|
| leave-one-out, fitted on the run's own ions | **1.69%** | 1.70% |
| the same instrument, previous day, applied unchanged | **7.41%** | 8.01% |

Both reproduce. The previous-day calibration is **4.4× worse**.

**The direction is stronger than reported.** Not "mostly one way" — **618 of 618 errors are
negative**, signed median −7.41%. A uniformly signed error across every ion is a systematic
offset, not drift in the colloquial sense, which raises a question nobody has asked: a
single scale factor fitted per day might recover most of that 7.41%. That would change
"the calibration does not transfer" into "it transfers with one correction", which is a
materially different claim and a cheap thing to test.

## 3. Precursor intensity barely tracks load — unexplained

Step one of the ar1 linearity test, and it has to be understood before the rest of that
test means anything. For every ion seen across at least four load levels, log10 intensity
regressed on log10 on-column amount:

| | |
|---|---|
| ions spanning ≥4 load levels | 35 |
| median log-log slope | **+0.239** |
| quartiles | −0.008 to +0.389 |
| slope > 0.5 | 3 of 35 |

A 50-fold load increase should move precursor intensity with a slope near **1.0**. The
median is **0.24**, and a quarter of these ions are flat or falling.

**Three readings remain open and this run cannot separate them:**

1. Severe detector compression at these loads, in which case the usable dynamic range is far
   narrower than the loading series suggests and every quantitative claim has a ceiling
   nobody has stated.
2. The condition labels do not reflect actual on-column amount — a dilution, split or
   injection difference not captured in the string.
3. `inten` is not a raw precursor intensity but something already normalized.

Until one of those is established this number goes no further than this document. It is
not a finding yet; it is a reason to look.

**It also gates the ar1 work.** If intensity does not track load, ar1's loading-invariance
cannot be tested against this series — the control fails before the measurement starts,
which is exactly why the control is specified first.

## What this run did not touch

No raw files, no consensus database, no Atlas, no search output. Every number above comes
from data embedded in an already-published analysis. MARCH's own result slides — HeLa_120,
HeLa_60, the water blank, the FragPipe comparison, the Amherst MYD88 study — remain
unreproduced.
