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
`PrecMz` and no `CalcMz`.

**Retracted:** this section previously reported that probing the ion keys recovered one A1
isotope pick, `QHYPEAYASPSHTKGEQGLYPLPLLNST…` 6+ at key m/z 806.91 against 806.7440
monoisotopic. That identification does not hold. The keys carry m/z to two decimals, which
at z=6 is a window of ±0.030 Da on the neutral mass, and inside that window A1 and
deamidation are not separable:

| | neutral offset | m/z at z=6 | key |
|---|---|---|---|
| A1 isotope | +1.0034 | 806.9112 | **806.91** |
| deamidation | +0.9840 | 806.9080 | **806.91** |

They differ by 0.0193 Da, a third of the rounding window. The authoritative list of the 33
(below) does not contain that ion, so it is the deamidation. The probe was built to keep
exactly these two apart and the key-level shortcut walked into the confound anyway — the
guard works on the data it was written for, not on a rounded summary.

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

## The probe ran on the real database, and the answer was zero — which means nothing

A session on Scott-PC ran `probe_consensus.py` (pinned at `04288c3`) against the
production database, `D:\ProcessPipelines_Latest\ARiA_Consensus.db`, about 53 GB.
Column discovery worked with no flags: `observations`, 122,031,794 rows, with
`precursor_mz`, `calc_mz`, `peptide_sequence`, `charge_state`, `sample_name` and
`precursor_intensity`. On a 2,000,000-row read:

| | result |
|---|---|
| calibration | pooled median **0.00 ppm**, per-run −2.5 to +4.4 ppm |
| monoisotopic | 1,926,966 rows / 119,837 ions |
| A1 | **0** |
| any isotope | **0** |
| near the deamidation/A1 boundary | none |
| unexplained | 73,034 rows, all ±0.01 to 0.05 Da |

**That zero is not a result.** An isotope pick sits 1.003/z Da from the
calculated m/z, which is **358 to 1003 ppm** across the charge and m/z range of
real tryptic peptides. Nothing in the sample was even 50 ppm out. So a pick
could not have been found here, and the probe searching for one was a test that
could not fail — the same defect as a control required to be reported rather
than to fail, caught earlier in `validate_export.py`.

I first attributed that to the ingest contract, and **that attribution was
wrong.** `ingest_contract.py` sets `PRECURSOR_PPM = 10.0` and is called by
`vendor_consensus_builder` for the vendor stores under `C:\VendorDBs`.
`ARiA_Consensus.db`'s `observations` is the older consensus lineage, cleaned
2026-09-26, and largely predates the contract — **it was never gated at 10 ppm**,
which is also why its 4–50 ppm rows are no contradiction. The conclusion holds
and the reason does not: something upstream, most likely a search engine that
only emits matches within tolerance, kept the errors narrow. The probe's own
guard never depended on knowing which, because it is keyed on the observed error
width rather than on a named filter. Do not cite the ingest gate as a guarantee
for this table; it does not apply to it.

This is worth being precise about, because reporting "we checked and found none"
would have been a false negative dressed as a confirmation. The published PAX8
figure of 33 picks cannot have come from this table. It must have come from
upstream of the gate: the search output (DIA-NN report parquet, FragPipe
`psm.tsv`, the Waters spectrum files) or whatever the gate quarantined.

### The re-run, and the subset that was being buried

Re-run on the production database with spread sampling, 2,000,000 rows in 740 s:
pooled calibration 0.00 ppm, **0 isotope picks**, widest error 22.8 ppm, nothing
beyond 50 ppm, and both guards fired — `THIS RESULT IS NOT EVIDENCE` and the
computed-precursor warning at **1,943,552 of 2,000,000 rows (97%)**.

A per-sample breakdown run alongside it showed what the pooled figure hid:

| samples | |
|---|---|
| 780 | no row carrying **both** `precursor_mz` and `calc_mz` — untestable for mass error at all |
| 48 | measured throughout |
| 5 | computed throughout |
| 30 | mixed |

The 97% is by *row*, and it comes from a handful of very large computed or mixed
samples — Waters UDMSE runs of 200–330k rows each. So the database is not
uniformly untestable: **48 samples are measured**, and the question could be
asked of them.

The probe was burying them. A computed row has a delta of exactly zero, so it
landed in the monoisotopic bucket, held the widest-error figure down, and turned
a potentially usable subset into a blanket "not evidence". Computed rows are now
**excluded from classification rather than counted**, in Test 1 and in the
calibration of Test 2 — where they would otherwise peg every per-run offset at
zero, which is exactly what the 0.00 ppm was. The output states the computed
share, classifies each run as measured, computed or mixed, and reports on the
measured rows alone. When nothing measured remains it says the source cannot be
tested rather than reporting a result over arithmetic.

A test plants the real shape: 3,000 computed rows in one run and 120 measured in
another, 20 of them on A1 — 96% computed by row. It fails if the probe loses
those 20, counts the computed rows as monoisotopic, or calls the result not
evidence.

### The measured subset, run on its own

The 48 measured samples were copied read-only into a scratch database — 350,419
rows — and the probe run over all of them, no sampling:

| | |
|---|---|
| still computed, excluded | 9,485 rows (3%) |
| measured | 340,934 rows |
| by run | 24 measured throughout, 24 mixed |
| calibration, measured rows only | **+0.14 ppm** over 275,181 measurements |
| isotope picks | **0** |
| monoisotopic | 275,629 rows / 12,265 ions |
| unexplained | 65,305 rows, ±0.01 to 0.05 Da |
| widest error | 28.1 ppm, none beyond 50 |

The **+0.14 ppm** replaces the earlier 0.00 ppm, which was arithmetic fitted to
itself and should not be cited anywhere. And the answer is still zero picks from
a sample whose widest error is 28.1 ppm — an order of magnitude short of the 358
ppm a pick needs. No consensus table on that machine can answer the question.

### Telling a cut distribution from one that simply ended

That left the verdict resting on a fixed 50 ppm threshold, which cannot
distinguish *filtered* from *clean*. A first attempt compared the maximum error
with the 99th percentile, on the reasoning that a cutoff leaves them close
together. **That heuristic is wrong**: for ordinary scatter the ratio depends on
how many rows were read — about 1.4 at a few hundred rows against 1.9 at a
hundred thousand — so it calls a small clean sample cut. It was replaced before
it ever ran on real data.

What does work is density at the edge. A distribution that was cut still carries
appreciable weight in its last few ppm and then stops; one that ended naturally
has thinned to a handful of rows by the time it reaches its own maximum. The
probe now reports the |ppm| quantiles, the share of rows in the top 5% of the
observed range, and which of the two it looks like. Three tests: gaussian error
rejected past 10 ppm must be called cut, unfiltered scatter must not, and a
300-row clean sample must not — the case the ratio heuristic got wrong.

### The edge test corrected the claim it was built to check

Run on the same 48-sample scratch database, all 340,934 measured rows:

| |ppm| | |
|---|---|
| median | 1.7 |
| p95 | 14.5 |
| p99 | 21.2 |
| max | 28.1 |
| beyond 50 ppm | 0 |
| in the last 5% of the range (26.7–28.1) | **59 rows, 0.02%** |

0.02% is not density at an edge. The verdict printed was *"it looks like it
ended rather than being cut"*, which **contradicts the earlier reading that a
tolerance had been applied at about 28 ppm**. No tolerance was applied there.
That correction is the whole point of having built the test: the previous
sentence was a guess, and the data disagreed with it.

### But a taper does not rescue the zero

Those two findings then sat in tension — the distribution "ended rather than
being cut", yet the probe still refused the result. The refusal is right and the
wording was what needed fixing. A smooth taper to 28 ppm rules out a cutoff **at
28 ppm**. It says nothing about a cutoff at 100 or 200, because no row in the
sample goes anywhere near one. The absence of a cliff where the data stops is not
evidence that data would have been kept far beyond it.

So the probe now states the taper narrowly, and names what is actually missing:

> **a positive control.** One row, anywhere in this source, carrying a mass error
> of a pick's magnitude — a few hundred ppm — that the pipeline nevertheless
> stored. An entrapment sequence, a deliberate misassignment, a quarantine
> record. Without one, "no picks" and "picks cannot get in" produce identical
> output, and a smooth tail at 28 ppm does not separate them.

That is the same shape of requirement as everywhere else in this repository: a
control that is reported but cannot fail is not a control. A test for something
absent needs a demonstration that the something would have been visible had it
been there.

### The circularity is confirmed, and it reaches the search outputs too

`calc_mz` in the vendor stores is built from the sequence **plus the
modifications the search engine assigned** — that is, after the search chose its
explanation of the observed mass. The MSFragger route then writes the same
computed value into both columns. So the search outputs carry the same defect as
the consensus tables, and moving upstream by one step does not escape it.

An isotope pick is a *decision* this pipeline records in an assignment, not a
residual it leaves lying in a mass column. Recovering the 33 requires either the
assignment itself — whatever field or log says "A1 was chosen here" — or the raw
spectra, where the reported precursor exists before any sequence has been
proposed. Everything between those two ends has already folded the answer into
the question.

### The quarantine may hold both the control and the answer

Each vendor store (`C:\VendorDBs\{Bruker,Thermo,Waters}_consensus.db`) has two
tables beside `observations`: **`ingest_rejected`** (file-level quarantine) and
**`ingest_row_rejects`** (rows kept out of an otherwise-accepted file by
`ingest_contract`, by reason, including the 10 ppm precursor check).

A row rejected at 10 ppm with a pick-sized error is exactly the positive control
this document asked for. And it is more than a control: **if the gate removed
pick-sized rows, the picks are in the reject pile.** That would be where to
count them from, not the accepted rows.

It turns on one thing nobody has looked at yet — whether those tables keep the
**rejected row** (its reported and calculated m/z, sequence, charge) or only
**tallies by reason**. A count records that something was dropped, not what it
was. It cannot be classified, so it cannot show that a pick-sized error would
have been caught, and it must not be read as a negative.

`tools/probe_consensus.py` Test 5 handles either shape. It finds tables whose
name suggests rejected or quarantined rows, runs the same classifier over them,
and reports:

- **picks among the rejected rows** → the positive control, and the place to
  count from. It says plainly that a zero in Test 1 then measures the gate
  rather than the pipeline.
- **rejected rows, none of them picks** → a real negative for that source,
  unlike Test 1's zero, because errors of a pick's magnitude were neither stored
  nor rejected. With a caution to confirm the gate actually ran on those files.
- **tallies only** → reported as unusable, naming what the table would have to
  keep instead.
- **no such table** → says so, and that it leaves no positive control.

Four tests cover those four shapes. Writing them also caught a defect in the
test harness itself: slicing a report on `TEST 4` ran past the end of that
section into `TEST 5`, so an ar1 case failed on the word `SKIPPED` belonging to
a different test. Sections are now bounded.

**Which sources could even produce such a row:** not Thermo. For 905 of 993
Thermo samples the reported m/z equals the calculated one, so those files cannot
generate a mass rejection at all. Only **Waters** (measured MGF `PEPMASS`) and
**Bruker** (where DIA-NN supplies an MS1 apex delta) can.

Reading those tables needs approval and has not been done.

### The thing to check upstream before anything else

Worth stating plainly, because it would make every mass test above vacuous no
matter which table it ran on. **This test assumes the calculated m/z was derived
from the sequence alone.** The PAX8 analysis's own method note says identity was
"taken from mass": carbamidomethyl assumed, then the residual explained as
oxidation, deamidation or an isotope pick, within 10 ppm — and
`explained_by_mass` is 1262 of 1262 PSMs. Every PSM was forced into an
explanation. So the ppm figures that analysis reports are residuals *after*
choosing the explanation, not raw mass errors, and its 33 isotope picks are
recorded in the assignment rather than left as 1.003 Da outliers in the data.

If a `calc_mz` column was written downstream of that choice, then reported minus
calculated is small by construction — for a pick exactly as much as for a
monoisotope — and no comparison of those two columns can ever find a pick. The
probe now says this in its own output when it refuses a result. The search
output is the right place to look not merely because it sits upstream of a
filter, but because that is where the reported precursor and a sequence-derived
theoretical mass are still independent of each other.

### Where the gate does apply, it still cannot fail

For the vendor stores, **905 of 993 Thermo samples carry a computed precursor
m/z** — `prec_mz` equals `calc_mz` exactly, via the MSFragger route. A row whose
reported m/z *is* the calculation cannot show a mass error of any size. So on
that source a mass test is not lenient, it is absent. The probe now detects this
directly: when a tenth or more of the rows have a reported m/z identical to the
calculated one, it prints `THE REPORTED m/z IS NOT MEASURED`, gives the share,
and says to find a column holding what the instrument reported. Two tests pin it
— a fully computed table must trip the warning, and the planted table with real
scatter must not.

### The ppm arithmetic, for reference

The 73,034 unexplained rows at ±0.01 to 0.05 Da work out, for ordinary tryptic
precursors, to roughly **4 to 50 ppm**:

| m/z | z | 0.01 Da | 0.05 Da |
|---|---|---|---|
| 500 | 2 | 10.0 ppm | 50.0 ppm |
| 750 | 2 | 6.7 ppm | 33.3 ppm |
| 600 | 3 | 5.6 ppm | 27.8 ppm |
| 900 | 3 | 3.7 ppm | 18.5 ppm |

Those rows predate the contract, so they are consistent. Either way 50 ppm is
still an order of magnitude short of the 358 ppm an isotope pick needs.

### What the probe now does about it

Two changes, both prompted by that run.

**It refuses to read its own zero.** Test 1 tracks the widest error from a
calculated m/z in the sample. If no isotope was found *and* no row sits beyond
`GATE_PPM` (50 ppm), the probe prints `THIS RESULT IS NOT EVIDENCE`, explains
that a pick could not have been stored, and names where to ask instead. Where
large errors do occur and the count is still zero, it says the zero carries
weight. A probe that reports a number without saying whether the number could
have come out differently is worse than no probe.

**`--limit` no longer means "the oldest rows".** A bare `LIMIT` returns the rows
SQLite reaches first, which on an append-only table of 122M rows is the earliest
ingests — so a 2M-row read sampled the oldest 1.6% of history. `--order` now
defaults to `spread`, reading equal slices from 64 points across the rowid range,
with `oldest` and `newest` available explicitly. The sample line in the output
states which was used.

Both are covered by tests: a gated table must be refused as evidence, an
unfiltered zero must be credited, and a limited read must say which rows it
covered. 17 of 17 cases pass.

### The PAX8 runs are not in the database at all

Checked read-only across every store: `observations` holds 863 distinct sample
names and **none** match CH215, JP21x or PAX8; the vendor stores match none
either. The Colorado PAX8 runs — 24 of them, Waters HDMSE — exist only upstream,
in the ApexRT/iDTs outputs and the FragPipe PSMs. **They were never ingested.**
So the premise that the raw reported precursor for PAX8 is in the consensus
database is false: the columns are there, those samples are not.

## Test 3 resolved: the load series was never a load series

The published analysis labelled its conditions `PAX8 JP210 100ng 120 min 1uL`
and so on, and the earlier run regressed intensity against the ng in that label.
The ng does not encode what went on column. Total PAX8 ion current per run, 60-
minute runs only, against both readings of the label:

| reading of "500ng" | Spearman vs run total current |
|---|---|
| total amount on column | **+0.000** |
| per µL, so on column = ng × µL | **−0.598** |

Neither works, and the second is *negative* — more nominal material, less signal.
The decisive pair: two runs on the same gradient, both labelled 100ng.

| run | volume | total PAX8 current |
|---|---|---|
| `PAX8_JP210_100ng_60_1uL` | 1.0 µL | 5.808 × 10⁷ |
| `PAX8_JP213_100ng_60_0-5uL` | 0.5 µL | 5.797 × 10⁸ |

**Exactly 10× more signal on half the volume at the same nominal amount.** JP210
through JP215 are different samples, not dilutions of one, and JP213 appears at
both 100ng/0.5µL and 500ng/2µL. There is no titration series in this data.

### The one clean comparison says the response is proportional

Hold the sample, the gradient and the injection volume fixed and only the
nominal amount changes: JP210, 60 minutes, 1 µL, 100ng against 500ng — a 5-fold
step. Per ion, best intensity in each condition, 27 ions in both:

| | |
|---|---|
| median ratio | **5.89×** against a nominal 5.0× |
| quartiles | 4.79× to 8.61× |
| implied log-log slope | **+1.102** (1.0 is proportional) |
| ions that went up | **27 of 27** |
| ions above 5× | 17 of 27 |

So **there is no detector compression**, and the earlier `+0.239` was an artifact
of regressing against a label that does not track load. Of the three readings
left open above, reading 2 — the condition labels do not reflect on-column
amount — is the one the data supports. Reading 1, severe compression, has no
support: with a real 5-fold step the response is proportional. Reading 3 is also
out: `inten` is a raw per-PSM precursor intensity, not something normalized.

**This unblocks the ar1 work.** The control that was failing was not the
detector; it was the load axis. Within a matched series intensity tracks amount,
so ar1's loading-invariance can be tested — and the dynamic-range ceiling this
document warned about does not exist at these loads.

One caveat was kept in view: 27 ions from one 100ng run against three 500ng
runs, with two ions moving 49× and 65× — more than a 5-fold step can explain.
Chasing those two corrected the headline number. See below.

### Correction: the 5.89× was inflated by taking the best of three replicates

Those 500ng figures were the **maximum** per ion across three replicate runs.
Taking the maximum of three draws and comparing it against a single draw
inflates the ratio, and this repository had already documented exactly that bias
in another place — `docs/top3-estimator.md`, where a peptide seen 20 times scores
13.6% above one seen 3 times at equal abundance. Same defect, new location.

Restricting to the 24 ions present in all four runs and taking the **median** of
the three replicates:

| | max of 3 | **median of 3** | min of 3 |
|---|---|---|---|
| median ratio | 7.40× | **5.37×** | 4.64× |
| p90 | 12.95× | 8.79× | 7.68× |
| largest | 65.45× | 23.26× | 17.75× |

So the honest figure is **5.37× against 5.0× nominal, an implied slope of
+1.044**, not the +1.102 reported earlier. Quartiles 3.44× to 8.20×, and 23 of
24 ions rose. The conclusion is unchanged and slightly stronger — the response is
proportional — but the number in the earlier section was biased upward and this
is the one to use.

### Both outliers were one irreproducible ion

`TDAFSQHHLEPLECPFER|554.01|4`, the 65× ion, spans **47,452 to 44,841,857 across
the three 500ng replicates — a factor of 945**. The 65× was max-of-3 selecting
that top replicate. On the median it moves 1.13×. It is not a load effect; the
ion is not reproducible.

That sets a floor on what this comparison can resolve. Across all 24 ions the
replicate spread is **1.69× at the median** and 4.07× at p90, so a single run at
the low end cannot distinguish a 5× step from a 3× or an 8× one per ion. It is
only the median over ions that is informative.

### The low end is unresolvable, and says so

At peptide level, collapsing charge channels, `TDAFSQHHLEPLECPFER` moves 1.67×
against a nominal 5×, and `QNPTMFAWEIR` moves 8.72×. The low one is not
interpretable: that peptide took **25.6%** of the 100ng run's entire PAX8 ion
current against 10.2% across all runs. Its denominator is high, so its ratio is
low, and with one run at 100ng there is no way to separate "under-responds" from
"measured high once". One ion fell outright —
`KMDDSDQDSCRLSIDSQSSSSGPR|665.29|4` at 0.84×.

### One thing that did check out

The data carries oxidised and unoxidised forms of the same peptide at the same
charge, 16.00 Da apart — `QNPTMFAWEIR` at 696.84 and 704.84, `MDDSDQDSCR` at
614.71 and 622.71, `VQQPFNLPMDSCVATK` at 612.30 and 617.63. Each pair responds
together: 8.8× and 8.5×, 4.3× and 5.1×, 3.6× and 3.7×. Two independent channels
for the same molecule agreeing to within 20% is a genuine internal consistency
check, and it passed without being set up as one.

## The 33 are recovered, and they were never in a mass column

`pax8_analysis.py` records an `isotope_pick` per PSM (0/1/2) during its
mass-identity step, and the published 33 comes from that field. Re-run on the
Colorado PAX8 data: **1,262 PSMs, all explained, 33 picks, every one A1**
(+1 × 1.00335), spread over 22 runs and **7 distinct ions**:

| ion | picks | of PSMs | runs |
|---|---|---|---|
| `LSIDSQSSSSGPR` 2+ | **17** | 49 | 14 |
| `LSIDSQSSSSGPRK` 2+ | 4 | 31 | |
| `DRLLAEGVCDNDTVPSVSSINR` 3+ | 3 | 31 | |
| `LLAEGVCDNDTVPSVSSINR` 3+ | 3 | 50 | |
| `VQQPFNLPMDSCVATK` 3+ | 3 | 59 | |
| `SLSPGHTLIPSSAVTPPESPQSDSLGSTYSINGLLGIAQPGSDK` 4+ | 2 | 33 | |
| `SLSPGHTLIPSSAVTPPESPQSDSLGSTYSINGLLGIAQPGSDKRK` 5+ | 1 | 15 | |

The picks are **FragPipe PSMs whose observed precursor sits +1.00335/z above the
monoisotope** — MSFragger's isotope-error allowance. So the A1 choice is the
search engine's assignment, exactly as the circularity argument predicted. The
whole mass-residual approach was looking for a decision in a place decisions are
not stored.

**What the 33 therefore measures.** It is a count of search-engine isotope-error
assignments on this protein, not a count of ARiA mis-picks. 33 of 1,262 PSMs is
2.6%, and it is concentrated: **one ion, `LSIDSQSSSSGPR` 2+, accounts for 17 of
the 33** — 35% of that ion's own PSMs. Six other ions supply the rest. This is
not a diffuse rate across the protein.

**One correction for the published page.** Its method note reads "an isotope pick
(A1/A2 chosen instead of the monoisotope)". All 33 are A1; no A2 occurred. The
"A1/A2" implies a category that is empty.

## Why no table could have held them

### Waters

Neither `ingest_rejected` nor `ingest_row_rejects` exists in the Waters store —
its ingests all come from the older CPPIS route and predate the gate. No positive
control there. The Waters CPPIS file also had **every precursor equal to
`CalcMz`**, so it is computed, not measured.

### Bruker

`ingest_row_rejects(source_file, reason, n_rows, recorded_at)` — **tallies
only**, 48 rows, mostly *"precursor outside 10 ppm of computed CalcMz"*. No
rejected row, m/z, sequence or charge is kept. That is the probe's `tallies only`
outcome: unusable as a control, and **not** a negative. `ingest_rejected` has 0
rows, so no Bruker file was quarantined wholesale — the count is not short for
that reason.

Regenerating the rejected rows from a source spec file still on disk (PXD064532
KO1): **2,613 rows / 1,154 precursors rejected at 10 ppm, all 10–30 ppm, with a
maximum neutral offset of 0.0199 Da** and **zero** within 0.02 Da of ±1× or 2×
1.00335.

That 0.0199 Da is a hard cap, and the arithmetic says it is a cap in **daltons,
not ppm** — which is what identifies it as DIA-NN's own MS1 tolerance rather than
the ingest gate, whose threshold is stated in ppm:

| reported ppm | implied neutral mass at 0.0199 Da |
|---|---|
| 10 ppm | 1,990 Da |
| 20 ppm | 995 Da |
| 30 ppm | 663 Da |

A fixed 0.02 Da window across an ordinary tryptic mass range of roughly
660–1,990 Da produces exactly the observed 10–30 ppm spread. So the Bruker route
cannot emit a pick at all: bounded upstream, the same class of vacuity as Test
1's zero.

### And the gate's own reason is circular

The rejection reason reads *"outside 10 ppm of computed CalcMz"*. The gate
compares the reported precursor against a **computed** value — so where that
computed value was derived from the search's own assignment, the gate inherits
the circularity it was being asked to break.

## Status of this line of work

**Closed.** The question — how many ions were named on an isotope peak — is
answered at 33, all A1, 7 ions, from the field that records the decision. The
mass-residual approach is retired, not merely unavailable: an isotope pick in
this pipeline is a choice recorded in an assignment, and no comparison of a
reported against a calculated m/z was ever going to find it.

Tests 1 and 2 of `probe_consensus.py` keep their value as **detectors of this
class of defect** rather than as instruments for this question. They refuse a
vacuous result four ways now — computed precursor, mass-gated sample, tally-only
quarantine, and a taper that only rules out a cutoff where the data stops.

Two routes stay open and are being fixed at the source: the Thermo MSFragger
writer is being changed to write the observed `precursor_neutral_mass` and
peak-read fragments. After that, MSFragger isotope-error PSMs will be rejected at
10 ppm — **and the store will still only count them** unless
`ingest_row_rejects` is extended to keep rows. Extending it is the one change
that would make the quarantine answer this class of question in future.

## The published page, corrected

PAX8 Across Platforms is at version 6. Three changes, and the one that matters
is not the A1/A2 wording.

**1. A `pass` that could not fail.** The page's first finding read:

> **pass** · Mass: every Colorado PSM is explained by mass within 10 ppm
> (median −1.5 ppm), once fixed carbamidomethyl and oxidation are restored.

Every PSM was forced into an explanation — `explained_by_mass` is 1262 of 1262
— so "explained" had no other possible outcome, and the −1.5 ppm is the residual
left *after* choosing carbamidomethyl, oxidation, deamidation or an A1 isotope.
It was presented as a result the method passed. It is now tagged **limit** and
says so: not an independent mass error, and not a test the method passed.

This is the same defect this repository has now caught four times in four places
— a validator requiring controls to be reported rather than to fail, a zero
isotope-pick count from a table that could not hold one, a calibration fitted to
rows whose reported m/z *was* the calculation, and now a published `pass` badge
on a tautology. Worth naming as a pattern rather than three coincidences.

**2. The A1/A2 category.** The method note said "an isotope pick (A1/A2 chosen
instead of the monoisotope)". All 33 are A1; no A2 occurred. The note now reads
A1, and adds where the picks come from: FragPipe PSMs sitting +1.00335/z above
the monoisotope under MSFragger's isotope-error allowance, recorded per PSM by
the analysis script rather than inferred from anything measured on the page.

**3. The pick count now carries its shape.** The mass note added: all 33 are A1,
on 7 ions across 22 runs, with 17 of them on `LSIDSQSSSSGPR` 2+ alone — 35% of
that ion's own PSMs. A reader seeing "33 isotope picks" would otherwise take it
for a rate across the protein. It is one ion and six stragglers.

Not changed, and worth a look by someone who can render it: the page states
"the archive stores no precursor m/z". For most archive samples that holds —
780 of 863 have no row carrying both a reported and a calculated m/z — but 48
samples are measured throughout, so the claim is too broad as written. Whether
any of those 48 carry PAX8 has not been checked.

## The fragmentation section was a broken comparison, not a negative result

Asked whether the fragment patterns really are that different, the answer turned
out to be that the question cannot be answered from this section. Recomputed from
the page's own 28 ions:

| measure | correct pairing | mismatched control | real beats control | sign test |
|---|---|---|---|---|
| Spearman of ranks | +0.119 | +0.000 | 10 of 21 | p = 1.000 |
| Cosine of vectors | **0.777** | **0.847** | **6 of 23** | **p = 0.035** |

**The control beats the test.** Not ties it — beats it, significantly. Pairing a
Colorado pattern with a *wrong* PAX8 peptide's Thermo pattern scores higher than
pairing it with the right one.

A control that reliably outperforms the thing it controls is not a weak negative.
It is an invalid null, and nothing can be read off it in either direction.

### Three explanations, two of them killed

**Thin Thermo data — refuted, and backwards.** 9 of 28 ions rest on a single
Thermo spectrum and 15 on three or fewer, so thinness was the obvious candidate.
Restricting to better-evidenced ions makes it *worse*:

| subset | ions | cosine real / control | real wins | p |
|---|---|---|---|---|
| all | 28 | 0.777 / 0.847 | 6 of 23 | 0.035 |
| ≥5 Thermo spectra | 12 | 0.777 / 0.850 | **1 of 11** | **0.012** |
| ≥10 Thermo spectra | 4 | 0.702 / 0.836 | 0 of 4 | 0.125 |
| ≥5 spectra and ≥6 shared | 10 | 0.787 / 0.840 | 1 of 10 | 0.021 |

**Cosine inflated by a smaller shared set — refuted.** The hypothesis was that a
mismatched peptide shares fewer fragment labels and cosine over fewer dimensions
is trivially higher. Within the real pairings the Spearman between shared-fragment
count and cosine is **−0.016**, and ions with ≤6 shared fragments score 0.769
against 0.796 for those with ≥10. No such effect. The rank measure, which carries
no dimension bias at all, is a plain coin flip.

**How the control's partner is chosen — untested.** If the control pairs against
a Thermo pattern that is systematically better sampled than the real partner's,
it is smoother and wins for that reason alone. The page records the real pairing's
spectrum counts but not which ion each control used, so this cannot be checked
from the published data. It is the remaining candidate and the one to check in
`pax8_analysis.py`.

### The page now says so

Version 7. The finding is **limit**, not `fail`, and the section's caveat is
retitled "Not a negative result — an invalid one", carrying the numbers above,
the two refuted explanations, and the untested one.

That is the **second** verdict on this page to change for the same reason, in
opposite directions: a `pass` on a check that could not fail, and a `fail` read
off a control that could not lose. Both were badges asserting a conclusion the
comparison could not support.

**This does not touch MARCH's `pattern` slide**, which is HeLa data from a
different study and reports 0.61 for real against 0.51 for mismatched — the real
pairing ahead, as it should be. The defect is in this PAX8 comparison, not in the
measure itself.
