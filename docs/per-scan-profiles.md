# The per-scan profile, and what it has to carry

Two questions in this project are blocked on the same missing file, and both
tools that would answer them are already written and tested.

**Does the noise fall at saturation?** A droplet of given volume carries a fixed
maximum charge, so a clamped signal should fluctuate less than a free one. That
is a real prediction and nothing measured so far can test it: a scan's total ion
current averages roughly 2 x 10^7 droplets at nanoflow, which dilutes per-droplet
clamping by the square root of that, and it sums 200 mobility bins besides. The
4-6% that survives is about 250x larger than droplet-counting noise, so it is a
whole-spray fluctuation that a per-droplet ceiling does not touch. The base peak
is no better: it is the maximum over those same bins and the ion holding it
changes scan to scan, so its scatter includes identity switching rather than one
channel's behaviour over time.

**Is suppression scan-level or channel-level?** Droplet competition depresses
every ion in a scan, dim ones included. Dead time undercounts only bright
channels. The two call for different corrections and a vendor emits one flag for
both.

Both need the same thing: one ion, followed across the scans it elutes over.

## The file

One row per ion per scan. `tools/saturation_check.py` reads this shape already,
and `tools/counting_noise.py` reads any subset with a scan index and a channel.

| column | what it is | why it matters |
|---|---|---|
| `ScanIndex` | the scan this row belongs to | the axis everything is measured along |
| `PeakGroupIndex` | which ion the row belongs to | groups the rows into one profile |
| `Counts` | **raw detector counts in this scan** | the only parameter-free expectation available |
| `Intensity` | the extracted intensity in this scan | what the correction would act on |
| `Mz`, `Drift` | the coordinate and the extraction window | so the window is known rather than assumed |
| `IsSaturated` | the vendor flag, per scan | one of the two conditioners |
| `run` | the acquisition | keeps scans from different runs apart |

## The two choices that decide whether it works

**Raw counts, not a scaled intensity.** Variance against mean is the whole test.
A normalized, smoothed or scaled intensity has no Poisson basis and the Fano
factor computed from it measures the scaling. The export in hand fails here:
`Counts` is the single value 119 in all 864 rows while `Intensity` spans 733x,
so it is not a count.

**Per scan, not per feature repeated.** The same export writes one intensity per
peak group across all three of its scans - 288 of 288 groups hold a single value.
A feature's apex repeated across its scan range is not a profile; it is one number
written three times, and every test here would be comparing a number against
itself.

## What will bounce

Both tools refuse rather than report. Worth knowing before the export is run:

- a `Counts` column that does not vary with signal
- a channel repeating its previous value in more than 20% of scans, which is a
  smoothed or forward-filled series
- peak groups holding one intensity across every scan they span
- `IsSaturated` reproduced exactly by an intensity threshold, since conditioning
  on it would then only condition on brightness
- an unpopulated flag column, which must not read as clean
- fewer than 6 scans in an ion's profile, or fewer than 60 scans in a run

## How much

One run is enough; both tests are within-run and need no cohort.

- **50-100 ions**, spanning dim to bright. All bright leaves no unsaturated
  control; all dim leaves nothing saturated.
- **10-15 scans per ion**, covering the elution rather than the apex alone. The
  estimator needs four scans outside the ones being tested to predict from.
- the ion set chosen **before** looking at intensity, not as the top N. Selecting
  on brightness selects on the quantity under test, and on planted data that
  alone contributed +1.7% of a false result.

## Where it can come from

Either the ApexRT export written per scan with real counts, or a direct
extraction from the `.raw` through `MassLynxRaw.dll`, which is how the HeLa
profiles in `examples/hela-scans.example.jsonl` were read. The direct route gives
raw counts by construction and does not depend on what the feature finder does to
intensities.
