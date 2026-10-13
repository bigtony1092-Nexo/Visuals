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

The viewer reads `manifest.json` **first** and decides from it which panels can honestly
render — before loading a single ion.

## 1. The forensic track

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

## 2. Drift is not always drift

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

## 3. Confidence tier is first class

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

## 4. Quant: Top3avg

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

## 5. Stoichiometry and flux

- **Stoichiometry** ratios are computed **within one run only** (`require_same_run` is
  `const: true`). `min_peptides_both_sides` guards the 1-vs-3 peptide mismatch.
- **Flux** requires a stated normalization and defaults to overlaying `total_signal`, so
  a "flux" that is really a loading change is visible as one. In the Amherst cohort, raw
  Top3avg tracked total sample signal only weakly (median r ≈ 0.23), so loading is not a
  universal explanation — but it is the first thing to rule out.

## 6. Can an ion map tell HeLa from K562?

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

## 7. The method lifecycle

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

## 8. Decimation is disclosed

"Complete library" can be millions of ions against a few thousand renderable marks.
`manifest.decimation` records strategy, `kept`, `of` and `seed`, and
`render.show_decimation_ratio` is `const: true`. A thinned plot states what it dropped.

## 9. Open questions

- `ion_uid` hashes `source_file`. Re-converting a raw file to a new path mints new uids.
  Either freeze the path at first ingest or hash a stable file id instead.
- `params_hash` needs a home: the parameter sets it points at must be stored somewhere
  resolvable, or the hash is an audit dead end.
- Isotopes are not yet in the contract. The flight animation draws three per precursor;
  decide whether an isotope is its own ion record or a child of one.
- `atlas_consensus` needs a definition sharp enough to compute — how many samples, at
  what agreement.
