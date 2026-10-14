# ARiA Visuals — export contract

The pipeline runs where the data lives. The viewer runs in a browser and cannot open
`F:/...`, cannot walk a raw folder, cannot read a SQLite file on disk and cannot launch
DIA-NN. Everything between them travels as three files. That boundary is the reason this
contract exists.

```
pipeline (python, on the box with the data)        viewer (published page)
  walk raw_dir, infer vendor/method  ─┐
  DIA-NN / pass2 / aria-spatial       ├─> manifest.json   ─┐
  consensus db, merge + cluster       │   ions.jsonl       ├─> panels
  Top3avg, normalization              ─┘   proteins.jsonl  ─┘
```

| File | Shape | Schema |
|---|---|---|
| `manifest.json` | one object | `schema/manifest.schema.json` |
| `ions.jsonl` | one ion per line | `schema/ion_track.schema.json` |
| `proteins.jsonl` | one (protein, sample) per line | `schema/protein_quant.schema.json` |
| `spectra.jsonl` | one spectrum per line, with per-peak `ar1` | `schema/spectrum.schema.json` |
| `traces.jsonl` | one XIC per line, optional | `schema/trace.schema.json` |

`manifest.tables` declares which of these were emitted, so a panel whose table is absent
blanks with a reason instead of failing on a missing file.

The viewer reads `manifest.json` **first** and decides from it which panels can honestly
render — before loading a single ion.

## 1. Panel order: the familiarity ladder

A new paradigm and an unfamiliar chart are two separate asks. Making them at once means
the reader has to grant both before they can check either, and the usual outcome is that
they grant neither.

So every panel declares a `familiarity` in `schema/panels.json`:

| Tier | What it is | Panels |
|---|---|---|
| `conventional` | reads without explanation | `msms_mirror`, `isotope_envelope`, `xic`, `intensity_vs_time`, `seqstart_vs_intensity`, `top3_percentile` |
| `bridge` | familiar axes, unfamiliar emphasis | `drift_vs_time`, `mz_vs_drift`, `stoichiometry`, `flux` |
| `novel` | asks for a new way of seeing | `ion_flight`, `ion_map_classifier`, `forensic_track` |

**`msms_mirror` is the load-bearing panel.** Head-to-tail experimental against Atlas
library, b ions blue, y ions red, unassigned grey, shared m/z axis, dot product on
square-root intensities, ppm per matched peak. It is the view in which an identification
lives or dies in every search engine's output. If geometry-based selection picked the
right ion, that is where it shows — in a frame nobody has to be taught, and one that was
not built to flatter the method.

The argument runs in that order: the conventional panels agree first, and only then does
the 3-D map get asked for, at which point it is explaining a result the reader has
already accepted rather than asserting one. `tools/plan_panels.py` prints an advisory
when a method puts a novel panel ahead of every conventional one.

This is also why the contract carries `spectra.jsonl` at all. The first draft had no
fragment data, which made the single most persuasive panel impossible to draw.

## 2. The forensic track

`ion_uid` is assigned at first detection and never reassigned:

```
ion_uid = sha256(sample_id | source_file | mz@4dp | rt@3dp | drift@4dp-or-null | charge)[:16]
```

Content-addressed, so it survives re-merge and re-clustering. `lineage.cluster_uid` and
`lineage.atlas_uid` point outward to what the ion became; `lineage.parent_ion_uids` points
back to what it came from.

`events[]` is append-only and ordered oldest first:

```
detected → z1_filtered → peg_filtered → selected → searched
         → identified → annotated → merged → clustered → atlas_member
                                                       → rejected
```

Each event carries `at`, `tool`, `tool_version`, `params_hash` and an optional
`coord_snapshot`. **A correction is a new event, never an edit to an old one.**

`coord_snapshot` is what makes recalibration visible: when RT recalibration moves an ion
0.014 min, the move is in the ledger instead of silently replacing the old number.

**The animation and the audit trail are the same object.** The flight animation is a
replay of `events[]` in order. Building the ledger is not extra work done for
bookkeeping's sake — it is what the animation reads.

### Isotopes belong to their precursor

An isotope is never an ion record of its own. The M+1 and M+2 peaks share the
precursor's RT, drift, charge, identity and `ion_uid`; only m/z and intensity differ.
Separate records would triple `ions.jsonl`, mint `ion_uid`s for things that were never
detected independently, and let one precursor's isotopes land in different clusters.

So the envelope rides on the precursor as `isotopes[]`: one entry per peak with `k`
(0 = monoisotopic, at `coord.mz`), observed `mz`, and `rel` (most intense isotope = 100).
`observed: false` marks an averagine-modelled envelope, which a panel must label as
modelled. The raw MS1 scan the envelope was read from, interferences included, stays in
`spectra.jsonl` as the `ms1_isotope` spectrum: `isotopes[]` is the assignment, the
spectrum is the evidence.

## 3. Drift is not always drift

`coord.drift_kind` is mandatory on every ion:

| `drift_kind` | Source | `drift` | Drift panels |
|---|---|---|---|
| `ion_mobility` | HDMSE, diaPASEF | real mobility | render |
| `quad_window` | Waters SONAR | quadrupole window position | **excluded** |
| `absent` | DDA | `null` | blank |

SONAR populates a drift column with real numbers that are not drift times. Plotting them
on a mobility axis produces a clean, convincing, meaningless chart. Nothing in the numbers
themselves announces the problem, which is why the gate is a required field and not a
runtime guess. Keyed the same way as `instrument_tolerance.py`'s
`drift_tolerance(vendor, acquisition_method)`.

## 4. Confidence tier is first class

| Tier | Gate | Decoy-controlled |
|---|---|---|
| `bare_proximity` | ppm window + ≥2 consecutive scans, no MS2 | **no** |
| `pass2_confirmed` | Pass2 targeting, decoy-controlled | yes |
| `atlas_consensus` | holds across samples in the Atlas | yes |

`bare_proximity` is the `direct_match_results.csv` tier from
`match_targets_spatial_in_cube()`. The same bare precursor-proximity mechanism passed
reversed-tryptic decoys at **33–70%** elsewhere in this project. It is useful for ranking
and for exploration; it is not evidence, and `tier_policy.mix_tiers` defaults to `false`
so a panel refuses rather than quietly blending it with a validated tier.

## 5. Quant: Top3avg

Exactly as the method note defines it:

1. per peptide, the **maximum** `matched_intensity` across every scan it hit — one number
   per peptide, not per scan
2. sort those maxima descending
3. average the top 3; with 2 peptides average 2, with 1 that value stands alone

`n_peptides_used` (1–3) travels with every record. **Padding to 3 is forbidden** — a
1-peptide Top3avg is a different measurement from a 3-peptide one, and a panel that
compares them without saying so is misreporting.

### Relative, not absolute

Top3avg is **relative abundance**. Absolute concentration requires a spike-in of known
amount, declared in `quant.absolute_standard`. With that field null, no panel may be
labelled as concentration — `tools/plan_panels.py` says so on every run.

### Normalization: each choice hides something

| Option | Keeps | Hides |
|---|---|---|
| `none` | everything | nothing — and nothing is comparable |
| `within_protein_max` | cross-sample comparability | real scale anomalies |
| `absolute_log10` | anomalies, five orders of magnitude | legibility on a linear axis |
| `median_scaled` | corrects loading | real global shifts |
| `reference_ratio` | anchors to a chosen protein | the reference's own instability |

Pick deliberately. In one 27-sample cohort GAPDH — the textbook stable reference — ran
from the **8th to the 96th percentile** on this measurement, with a cross-sample rank
std (991 positions) comparable to MYD88's (911). A housekeeping protein is not
automatically a stable denominator, and `reference_ratio` inherits whatever its reference
does.

`percentile` is preferred over `rank` for cross-sample work: sample protein counts range
2,331–4,235, and in a densely packed stretch of the distribution (nine proteins inside a
1.4% Top3avg span) a protein shifts several rank positions on measurement noise alone.

## 6. AR1: fragment conversion, and why it is loading-invariant

`ar1` on each assigned peak is that product ion's absolute intensity divided by its
precursor's:

```
ar1 = I_product / I_precursor
```

Scale the injection and both move together, so **ar1 does not change with load**. That is
the point of it: a difference in ar1 across runs is a difference in fragmentation, not in
how much went on column. `ar1_sum`, the total over the assigned ladder, is the share of
the precursor that converted to identified fragments — one number per ion, comparable
across runs, collision energies and loads.

Three rules keep it honest:

- **Do not normalise ar1 within a spectrum.** The precursor intensity is a constant there,
  so dividing by the largest ar1 cancels it exactly: `(I_p/I_prec) / (I_max/I_prec) =
  I_p/I_max`. That is base-peak relative intensity, which `peaks[].intensity` already
  holds. The raw ratio is the only part that carries new information.
- **`ar1_basis` is required whenever any peak carries `ar1`.** `summed_precursor_intensity`
  and `matched_intensity` are different denominators and give different numbers.
- **It must be computed upstream, from absolute intensities.** `peaks[].intensity` is
  base-peak-normalized, so ar1 cannot be reproduced from the stored spectrum.

**Where it is null, and why that matters.** ar1 needs a precursor intensity. Sources that
store fragments but no precursor intensity — the UID Atlas and `ARiA_Consensus.db` among
them — carry `ar1: null` throughout. So ar1 does not close a cross-platform gap; it
compares within sources that record precursor intensity. A panel that mixes a source
carrying ar1 with one that cannot is comparing a number against an absence.

### Ruling out a loading effect

Loading-invariance is what makes ar1 a test rather than a description. Where a difference
between runs or methods is suspected, the order is:

1. **Count what the enzyme can produce.** In a region with few cleavage sites there may be
   nothing to find; a gap there is arithmetic, not sensitivity. (PAX8: residues 288–450
   hold 3 cleavage sites against 29 in 1–287, yielding 2 observable peptides in that whole
   third.)
2. **Check whether amount moves it.** If coverage is flat across a real loading range, the
   gap is not a loading effect. (PAX8 on Colorado: 50-fold on-column range, coverage
   44.7–49.8%, zero ions past residue 287 in all ten runs.)
3. **Then compare ar1**, which holds amount constant by construction, to ask whether the
   ions that are seen fragment the same way.

Only after those does a modification hypothesis earn its cost — and the residues pick it:
N-terminal start → Met oxidation, N-terminal acetyl or a GG remnant; an S/T-rich stretch →
phosphorylation; a K-rich stretch → ubiquitin GG or acetyl.

## 7. Stoichiometry and flux

- **Stoichiometry** ratios are computed **within one run only** (`require_same_run` is
  `const: true`). `min_peptides_both_sides` guards the 1-vs-3 peptide mismatch.
- **Flux** requires a stated normalization and defaults to overlaying `total_signal`, so
  a "flux" that is really a loading change is visible as one. In the Amherst cohort, raw
  Top3avg tracked total sample signal only weakly (median r ≈ 0.23), so loading is not a
  universal explanation — but it is the first thing to rule out.

## 8. Can an ion map tell HeLa from K562?

Worth asking, and easy to answer wrongly. The defaults enforce the honest version:

- **Holdout by batch, never random.** A random split lets the model read instrument,
  column and loading, score near-perfectly, and tell you nothing about cell line.
- **`within_stratum_only`** also runs the test inside one vendor+method stratum. If
  accuracy collapses there, the signal was acquisition.
- **`min_tier: pass2_confirmed`**: a classifier fed a tier that admits decoys at 33–70%
  learns the noise.
- **≥2 batches per class.** `plan_panels.py` blanks the panel otherwise, because the
  question is unanswerable from that design — not because the code failed.

The shipped example demonstrates the trap: HeLa appears only in `amherst-A` (Bruker
diaPASEF), K562 only in `colorado-B` (Waters). Cell line is perfectly confounded with
both batch and vendor, so the panel refuses. Answering this question needs at least one
cell line run on at least two instruments, or two batches of each class on one.

## 9. The method lifecycle

`schema/method.schema.json` is the fillable form. Fill it, save it, point it at new
samples, get the same panels back.

```
fill method.json  →  save as a method  →  pipeline run (new samples)
                                       →  manifest + ions + proteins
                                       →  viewer renders the panels you picked
```

Check it before running:

```
python3 tools/plan_panels.py my_method.json                      # validate
python3 tools/plan_panels.py my_method.json export/manifest.json # feasibility
```

Per-panel requirements live in `schema/panels.json`, so a panel's needs are data, not
code. `on_requirement_unmet` defaults to `blank_with_reason`: the panel stays visible and
states why it is empty, so **a missing measurement is never mistaken for a zero**.

## 10. Decimation is disclosed

"Complete library" can be millions of ions against a few thousand renderable marks.
`manifest.decimation` records strategy, `kept`, `of` and `seed`, and
`render.show_decimation_ratio` is `const: true`. A thinned plot states what it dropped.

## 11. Depletion and interference

**Status:** the fields are defined in the schemas below. The pipeline does not produce them
yet. The Atlas reads them when present; when they are absent it computes the MS1 part
itself from `isotopes[]` and per-isotope traces, and labels that as its own estimate.

### Depletion

When the spatial matcher (`match_targets_spatial_in_cube`) annotates an ion, its spectra
go into the database and the ion **leaves the working ion list**. Each later pass searches
only what earlier passes did not explain. The remainder is the result: what is left after
depletion is what the library cannot yet account for, and how fast it shrinks pass by pass
is a measure of coverage.

In the ledger this is one more event:

```
{"stage": "depleted", "pass": 1, "explained_by": "atl-Q99836-172-z2", ...}
```

The flight animation replays it: the cloud fills, then each pass's ions fly back out
through the inlet into the database, until only the unexplained remain. An export without
`depleted` events can still be replayed by confidence tier (atlas consensus, then Pass2,
then bare proximity), and the Atlas says that is what it is doing.

### Streaming the isotope patterns with the matches

Every match already passes through the isotope envelope. Stream it alongside the match,
**scan by scan across the elution**, not as one summed envelope:

- **Expected envelope from composition, not averagine.** Once an ion is annotated its
  sequence is known, so its elemental composition and exact isotope distribution are
  known. Averagine is only for ions without an identity.
- **Residual per isotope** (`isotope_fit.residual`): observed minus expected, in
  percentage points with both envelopes scaled to a most-intense peak of 100.
- **Variation across the peak** (`isotope_fit.ratio_cv`): the coefficient of variation of
  each ratio `I_k / I_0` over the scans of the peak. A clean ion holds its ratios constant
  from leading edge to tail. A ratio that drifts means something else enters and leaves
  the window during the peak, and the scans where it departs (`interference.rt_start`,
  `rt_end`) show where.

The scan-by-scan intensities travel as **one trace per isotope** in `traces.jsonl`
(`isotope_k` = 0, 1, 2 for the same `ion_uid`), so a consumer can redraw the ratio across
the peak rather than trust a single number.

### When the interference is invisible in MS1

If the interferer begins and ends with the target, and there is not enough resolution in
m/z or drift to separate them, MS1 cannot tell them apart: the ratios hold constant across
the peak and the envelope can sit close to theory. **The only hope then is the product
ions:**

- **Unexplained fragments** (`interference.unexplained_ms2_fraction`): the share of
  experimental MS2 intensity the target's ladder does not explain.
- **`ar1_sum` falls** (`interference.ar1_sum_ratio`). `ar1 = I_product / I_precursor`.
  If the precursor intensity includes an interferer, the denominator is inflated and
  every one of the target's ar1 values drops together. Compare against the **same peptide
  and charge in runs where it is clean**, not against the library: library spectra carry
  no precursor intensity, so their ar1 is null (section 6). Loading-invariance is what
  makes the drop meaningful instead of noise.

### The fields

| Where | Field | Meaning |
|---|---|---|
| `events[]` | `stage: depleted`, `pass`, `explained_by` | which pass explained the ion, and by what |
| ion record | `isotope_fit` | `expected_basis`, `residual[]`, `rms`, `ratio_cv[]`, `n_scans` |
| ion record | `interference` | `flag`: none / ms1_envelope / ms1_profile / product_ion, plus `rt_start`, `rt_end`, `unexplained_ms2_fraction`, `ar1_sum_ratio`, `note` |
| `traces.jsonl` | `isotope_k` | one XIC per isotope, so ratios can be followed scan by scan |

`interference: null` means not assessed; `flag: none` means assessed and clean. The two
must not be conflated: an unassessed ion is not a clean one.

## 12. Open questions

- `ion_uid` hashes `source_file`. Re-converting a raw file to a new path mints new uids.
  Either freeze the path at first ingest or hash a stable file id instead.
- `params_hash` needs a home: the parameter sets it points at must be stored somewhere
  resolvable, or the hash is an audit dead end.
- `atlas_consensus` needs a definition sharp enough to compute — how many samples, at
  what agreement.
