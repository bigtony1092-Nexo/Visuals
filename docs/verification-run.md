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

## The probe that closes these three gaps

Two of the three results above stop at the same wall: the published data carries one
already-corrected precursor m/z per PSM, and no record of what the instrument reported.
A field scan confirms it — the embedded fields mentioning mass or m/z are
`archive_precursors`, `colorado_mass_identity`, `explained_by_mass` and `mz`. There is no
`PrecMz` and no `CalcMz`. Probing the ion keys directly recovered **one** A1 isotope pick
of the 33 reported (`QHYPEAYASPSHTKGEQGLYPLPLLNST…` 6+, key m/z 806.91 against 806.7440
monoisotopic, +0.166 = 1.0033/6); the other 32 were corrected before anything reached the
page.

`tools/probe_consensus.py` runs the test where the uncorrected value survives:

```
python3 tools/probe_consensus.py --find --schema-only
python3 tools/probe_consensus.py --db /path/to/ARiA_Consensus.db --load-map loads.json
```

`--find` searches the usual places rather than being told the path: a `.db`,
`.sqlite` or `.sqlite3` file whose name contains `consensus` or `aria`, under the
home directory and a handful of data roots, skipping `Library`, `node_modules`,
caches and dotted directories. It opens nothing but the 16-byte SQLite header,
lists every candidate so a wrong pick is visible, and uses the largest. A
database named something else needs `--db`.

Column names are discovered rather than assumed, and `--schema-only` prints what each
table offered so a wrong guess is visible. Three tests, each skipped with a stated reason
rather than faked when its column is missing:

| | question | needs |
|---|---|---|
| **Test 2** | the ppm offset of the monoisotopic cluster, per run | a reported and a calculated m/z |
| **Test 1** | how many ions were named on an isotope peak, on calibrated deltas | the same two |
| **Test 3** | log-log slope of precursor intensity against ng loaded, per ion | an intensity column, a run column, `--load-map` |

Test 2 runs first because test 1 depends on it. **Deamidation is +0.984016 Da and an A1
isotope is +1.003355 — 0.019339 Da apart.** An uncorrected calibration offset drags every
delta upward, and on the planted test data a +3.2 ppm offset pushed 25 rows into the zone
where the two cannot be told apart. Removing the offset first dropped that to zero. A
delta is tested against known modifications before isotopes for the same reason: a
calculated m/z that omits a peptide's modification produces an offset that says nothing
about isotope selection.

### What the probe was checked against

`tools/test_probe_consensus.py` builds a database whose answers are planted and asserts
the probe recovers them. 11 of 11 cases pass. Planted against recovered:

| planted | recovered |
|---|---|
| 33 ions on A1 | 33 |
| 5 ions on A2 | 5 |
| 40 oxidised, calculated column unaware | 40, named `oxidation` at +15.9944 |
| 20 deamidated | 20, named `deamidation` at +0.9842 — **not** A1 |
| 402 monoisotopic | 402 |
| calibration +3.2 ppm | +3.20 pooled, +3.18 to +3.21 per run |
| 250 ions at slope 0.25, 250 at 1.00 | quartiles +0.249 and +0.998 |

Writing those tests found two defects in the probe. The tolerance was 0.01 Da, which is
more than half the 0.019339 Da deamidation/A1 gap, so the two windows overlapped and a
delta between them was assigned by table order rather than by distance; it is now 0.008.
And the first ambiguity warning used a 0.03 Da band — wider than the gap it was policing,
so every A1 row tripped it and the warning carried no information.

Three negative controls also pass: a table shaped like the published page (one corrected
`mz`, nothing to compare it against) is **refused** rather than treated as raw; a database
with no usable table is refused; a file that is not a database exits 2 with the reason.

### Reading the result

Compare **ions**, not rows — one ion appears once per run. The probe reports A1 alone and
any-isotope separately, because 33 could be either. If the A1 ion count lands near 33 the
published figure is confirmed and the correction is doing what it claims. If it lands far
above, more ions were named on isotope peaks than the study reported.

Test 3 is the one with consequences past this document. If a raw precursor intensity in
the consensus database reproduces the **+0.239** slope, the compression is real and the
dynamic-range ceiling belongs in the export contract. If it lands near 1.0, then `inten`
in the published JSON was never a raw precursor intensity, and reading 3 above was right.
