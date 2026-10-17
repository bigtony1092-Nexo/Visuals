# The chain

Five tools, run in an order that matters, by one command:

```
python3 tools/aria_check.py \
  --export export/ --panels method.json manifest.json \
  --fragments march.jsonl search.jsonl \
  --ions ions.json --rank-rt --sigma-rt 6.8 \
  --counts counts.json
```

| | step | why it comes before the next |
|---|---|---|
| 1 | **export contract** | a violation upstream makes every later number suspect |
| 2 | **panel renderability** | no point measuring what no panel will draw |
| 3 | **chemistry** — *gate* | a flat source is not a reference; a comparison against it cannot be read either way |
| 4 | **resolvability** | the collision count bounds what any assignment can achieve |
| 5 | **ladder intervals** — *gate* | asymmetric ladders make every rate a floor |

## The one thing it enforces

**A step that did not run is never shown as one that passed.** A failed gate
stops the chain and every later step is marked `NOT REACHED`, not omitted:

```
  3 CHEMISTRY gate           GATE FAILED    a source carries no fragmentation
                                            chemistry, so a pattern comparison
                                            against it cannot be read
  4 RESOLVE   collisions     NOT REACHED    stopped by 3 CHEMISTRY gate
  5 STATS     intervals      NOT REACHED    stopped by 3 CHEMISTRY gate

  Not reached is not passed. Fix the gate and run again.
```

A step with no inputs is `skipped` and the summary says *"a skipped step has not
passed."* Exit status is 0 only when every step that had inputs ran and passed.

That rule exists because this repository has produced the opposite five times:
a validator whose controls were required to be reported rather than to fail, a
zero isotope-pick count from a table that could not hold one, a calibration
fitted to rows whose reported m/z *was* the calculation, a published `pass` badge
on an explanation every PSM was forced into, and an ingest gate that admitted
66.4M rows on a 10 ppm check against computed values. Each was a number that
looked like a result.

## The ride, on real data

Run over the PAX8 page's own two fragment sources, its 1,262 precursor
coordinates and the real replication counts, it **stopped at step 3**:

```
source                 fragments  y-share   proline (obs/exp)
colorado_fragpipe            525    0.696   0.261 / 0.128 = 2.04x
thermo_archive               404    0.483   0.093 / 0.092 = 1.01x

  thermo_archive: proline enrichment 1.01x is below 1.25x — indistinguishable
  from uniform … Do not read a pattern comparison against it in either direction.
```

That conclusion took most of a working session to reach by hand, through two
wrong explanations of my own — "thin data", then "low resolution" — before the
cause turned out to be unfiltered PSMs. **The chain reaches it in one command**,
and refuses to print the downstream numbers that would have been read off it.

With the archive excluded, all five steps run:

| | |
|---|---|
| 1 export contract | passed |
| 2 panel renderability | passed |
| 3 chemistry | passed — 2.04× proline |
| 4 resolvability | findings — **43.1% uniquely resolved, 56.9% sharing a cell** |
| 5 ladder intervals | passed — **13.49× (95% CI 5.2–35.3)** |

Step 4 reporting findings is correct: collisions are a fact to carry forward,
not a failure. It is not a gate for that reason — only a *missing measurement
scale* is fatal there, and the tool exits 2 for it, which the chain reports as
`COULD NOT RUN`.

## What the review found before any of this ran

A full review of the five tools returned **15 findings, 14 reproduced**. They are
fixed, each with a test that would have caught it. Three are worth naming because
they are the same failure the chain exists to prevent — a check that could not
fail:

**The ar1 normalization guard keyed on `(sequence, charge)` instead of on a
spectrum.** Fragments normalized per spectrum and read across three runs sum to
3.0 rather than 1.0, so the refusal never fired. Its test passed only because
the fixture used a single run. The guard built to catch this defect contained it.

**`--order spread` had no upper bound per block**, so slices ran into each other:
a 500-row table read with `--limit 5000` returned 4,654 rows, about nine copies
of each, inflating every count downstream.

**The percentiles came from a histogram capped at 400 ppm while the edge test
used the uncapped maximum.** A real isotope pick at ~1003 ppm therefore produced
`p95 400.0, p99 400.0` and the verdict *"looks like it ended rather than being
cut"* while 16.7% of rows sat at the maximum — wrong in exactly the case the
test was written for.

The rest: a flat-source gate skipped when no peptide contained proline (so a
uniform source printed the all-clear and exited 0); `--limit` emitting false
dangling-`ion_uid` errors in the documented spot-check mode; `"N-term"` parsed
character by character into residue N; four crashes on ordinary malformed input;
`isinstance(x, int)` rejecting JSON floats so symmetric ladders read as
asymmetric; an off-by-one median that reported 20 ions at 0.95 against 20 at 1.05
as a 5% systematic trim rather than as scatter; a verdict that depended on input
line order; and an exit status that consulted one test out of five.
