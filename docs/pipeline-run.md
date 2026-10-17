# The spine: raw folder to annotated ion lists

`tools/aria_run.py` is the path from a vendor raw folder to saved ion lists,
checked export and a panel plan. It exists because every piece of this
repository answered one question in isolation, and running them by hand left it
to a person to remember which answer invalidated which — the failure this
project has hit repeatedly.

    python3 tools/aria_run.py --raw CH215.raw --peaks peaks.jsonl \
        --targets targets.jsonl --out export/ --sample-id CH215 \
        --rt-tol 1.0 --drift-tol 0.02 --drop-peg --drop-z1 \
        --method examples/method.example.json \
        --sigma-drift 0.006 --sigma-rt 0.05

## The seven stages

| | stage | what it may assume |
|-|-------|--------------------|
|1| VENDOR    | nothing. It reads the folder's own files. |
|2| EXTRACT   | that peak picking happened somewhere else. |
|3| ANNOTATE  | that a coordinate match is a coordinate match. |
|4| IONLIST   | that the evidence on a record is all the evidence there is. |
|5| EXPORT    | only what was actually written. |
|6| CHECK     | nothing: it re-reads the export from disk. |
|7| PANELS    | only what the manifest declares. |

Stages 1, 2 and 4 are gates. A gate that fails stops the run, and every stage
after it is printed as `NOT REACHED` — because "we did not run it" and "it
passed" must never look alike. `run.json` in the export carries the same stage
log, plus the vendor record.

## Stage 1 is a gate because drift is not always drift

`drift` is ion mobility on HDMSE and diaPASEF, a scanning-quadrupole position on
Waters SONAR, and absent on DDA. `tools/detect_vendor.py` names the acquisition
method from `_extern.inf` or the tdf's own `GlobalMetadata` and refuses to guess,
because a tolerance keyed on the vendor alone treats a quadrupole position as a
mobility measurement.

So: an unknown `drift_kind` stops the run. An unknown **vendor** stops it too,
and `--drift-kind` cannot launder either one — it declares the drift axis, not
the vendor. A mixed-vendor folder and a folder-of-runs are separate refusals:
each run needs its own tolerances, and the method is per run.

When `--drift-kind` is used, `run.json` records `drift_kind_declared: true` and
the value the files actually gave, so a declaration never reads as a reading.

Downstream consequences, applied automatically:

- `quad_window` — drift is carried into `coord.drift` and **never** used as a
  matching dimension. `has_ims` is false. No `one_over_k0`, no `ccs`.
- `absent` — `coord.drift` is forced to null even if the peak list carried a
  number.
- `ion_mobility` — drift matches when `--drift-tol` is given. Without it, the
  run says the separating dimension sat unused rather than quietly skipping it.

## Stage 2 is the adapter boundary

Vendor peak picking does not happen in this repository. ApexRT,
`match_targets_spatial_in_cube` and the vendor SDKs live on the machine with the
data. This stage either runs `--extract 'cmd {raw} {peaks}'` or reads the
`--peaks` file that command already wrote. With neither, it stops. It does not
invent peaks.

Three failures are distinguished, because they mean different things:

- extractor exits non-zero → `COULD NOT RUN` (exit 2), with its own stderr shown.
- extractor exits **zero and writes nothing** → also `COULD NOT RUN`. Exit 0 with
  no output is a broken adapter, not an empty result.
- peak list present but empty → gate. An empty extraction is a result, but not
  one anything downstream can be read from.

A peak needs `mz`, `rt` and `charge`. JSON-lines, a JSON array or CSV with a
header all read; column aliases cover the pipeline's own CSV spellings
(`m/z`, `apex_rt`, `matched_intensity`, `1/k0`, …) so nothing is renamed in
transit.

## Stage 3: a match is a match, not an identification

Targets come from `--targets` or `--atlas`. **Not both** — a manifest that
records a consensus DB it never opened is a false audit trail, so giving both is
refused rather than resolved by precedence.

`--atlas` discovers the table and columns rather than hardcoding them (the DB has
been renamed and re-cut more than once here) and prints what it used.
`--atlas-limit` orders by m/z then charge before truncating, and says how many
rows were left behind — an unloaded target produces an ion that reads
*unassigned*, which is not the same as missing.

Matching is m/z within `--precursor-ppm` (default 10.0, the production ingest
contract's value), plus rt and drift where the target and the acquisition both
carry them. A window wider than 50 ppm is refused: an isotope pick sits 358–1003
ppm from the calculated m/z, so past that the window stops being a mass
measurement.

When several targets fall inside the window, the nearest in m/z is taken and
**said to be a tie-break**: once in the manifest warnings, and on each affected
ion's own `annotated` event and `score_kind`. Per-ion visibility is the point —
an aggregate warning lets a two-way coin flip read as a definite annotation.

### The target list selects. It does not identify.

If a peak carries its own search result *and* matches a different sequence in the
target list, the row's own sequence wins, because the score that promotes the ion
was computed for it. Exporting the target's sequence under the row's q-value
would write a decoy-controlled identification for a peptide nothing scored. The
disagreement lands on the ion's events and in the manifest warnings.

## Stage 4: tiers follow the evidence, and cannot be loosened by flag

| tier | what it requires |
|------|------------------|
| `bare_proximity` | an m/z window plus ≥2 consecutive scans. No MS2, no decoy test. |
| `pass2_confirmed` | a decoy-controlled score (`q_value`, `fdr`, `percolator_q`, …) under the cut. |
| `atlas_consensus` | agreement **across samples**. |

A single-run spine can never assign `atlas_consensus`, because one run cannot
observe agreement across samples. It does not.

A `hyperscore` or an `expect` is recorded as a score but cannot lift an ion out
of `bare_proximity`: there is no decoy model behind it.

`--min-scans` below 2 and `--q-cut` above 0.05 are **refused at the command
line**, not warned about. Those two numbers are what the tier names mean; a flag
that loosens either does not export more ions, it relabels a weaker thing as a
stronger one. A gate a flag can switch off is not a gate.

An ion that fails the scan-count requirement, the q-value cut, `--drop-z1` or
`--drop-peg` goes to `rejected.jsonl` with its reason on a `rejected`,
`z1_filtered` or `peg_filtered` event. It is not dropped silently, and it does
not sit in `ions.jsonl` under a tier whose evidence it lacks. If nothing
survives, that is a gate: there is nothing to export, check or draw.

When **every** exported ion is `bare_proximity`, the run says so — that bare
mechanism passed reversed-tryptic decoys at 33–70% elsewhere in this project.

### ion_uid

Content-addressed exactly as `schema/ion_track.schema.json` states:

    sha256(sample_id | source_file | mz@4dp | rt@3dp | drift@4dp-or-null | charge)[:16]

Fixed here and not to be changed: a change to the canonical string silently
orphans every ion list written before it. Two peak records that land on the same
address are not combined — the first wins, the later one is **dropped**, and the
run says "dropped, not combined" so nobody reads it as summed intensity.

`params_hash` on every event covers the tolerances, the cuts, the filters, the
drift semantics and the target list identity. Without the drift semantics in it,
two runs that read the same peaks under opposite interpretations of the drift
column would hash identically.

## Stage 5: the manifest declares absences

`tables` is written from what the run actually produced, so `spectra`, `traces`,
`proteins` and `fragment_hits` are **false** — this spine writes ions. Stage 7
then blanks the panels that need them, rather than failing on a missing file.

There is no `quant` block. The contract's only quant metric is `top3avg`, a
per-protein roll-up this spine does not compute. Declaring it would name a metric
nothing produced; naming `matched_intensity` would invent an enum value. Absent
is what "no quant claim" looks like in this contract.

## Stage 6: the chain, fed from what the run just wrote

`aria_check.py` re-reads the export from disk and runs the contract validator,
the chemistry gate, the collision count and the interval step.

The collision step is fed `resolve_ions.json`, written in stage 4 from the ions
themselves. Its `uid` is the **identity** (sequence and charge), not `ion_uid`:
`ion_uid` is content-addressed per detection, so each one would have exactly one
observation and the within-ion SD that resolvability needs could never be
estimated.

It is fed only when the scale can actually be estimated — at least one identity
measured three or more times, **and** at least two distinct identities. With one
identity the axis range *is* that ion's own scatter, so an estimated scale can
only come back as a confound, and there is nothing for it to collide with. One
run usually fails this, which is why `--sigma-drift` and `--sigma-rt` exist. When
it is not fed, the run says so and the step reports `skipped`, not `passed`.

## Exit status

0 when every stage that had inputs ran and passed. 1 when a gate failed or a
stage reported findings. 2 when a stage could not be run at all.

## What is still outside the spine

- **Spectra, traces and fragment hits.** The MS/MS mirror, the isotope envelope
  and the XIC panels need them and currently blank.
- **The cross-sample merge** that `atlas_consensus` requires.
- **`one_over_k0`.** On a mobility run with no 1/K0 in the peak list, every
  mobility panel plots an empty axis. The extractor has the number; it is not
  reaching the peak list yet, and the run warns about exactly that.
