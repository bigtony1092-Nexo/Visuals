# Visuals

Aria plotting.

Two things live here.

- **Ion Flight Atlas** — `ion-flight-atlas.html`, a standalone page that draws an
  LC-ESI-TIMS-TOF run as a flight tube, every ion placed by its own physics.
- **The export contract and its viewer** — panels built from annotated ion lists, driven
  by a saveable **method** file, fed by a contract that keeps a forensic track of every ion
  from detection to Atlas.

They share a subject and nothing else: the Atlas is self-contained and synthetic, the
contract is what real exports travel through.

## Ion Flight Atlas

`ion-flight-atlas.html` — a single self-contained page. Open it in a browser;
there is no build step and no network dependency except the Google Fonts link.

An LC–ESI–TIMS–TOF run drawn as what it physically is: a flight tube with a
hole in the front, and every ion placed in three-dimensional space by its own
physics rather than by the instrument.

### The three axes

| Axis | Quantity | What separates it |
|------|----------|-------------------|
| Down the tube | **m/z** | time of flight |
| Vertical | **1/K₀** | ion mobility drift — cross-section over charge |
| Across the tube | **RT** | retention time — hydrophobicity on the gradient |

The inlet is the origin. Each ion leaves it on the one angle those three
numbers allow, which is why the cloud forms the diagonal mobility–m/z ribbon
with its charge-state bands. Change the gradient and the ribbon changes shape;
the identifications do not.

### What it does

- **Liquid-phase panel** — the LC chromatogram, the MS1 isotope envelope of the
  precursor being traced, and the ESI capillary into the inlet. The three
  isotope peaks (red / blue / green, tallest first) close into circles sized by
  intensity, run down the capillary, and fly to their coordinate on a coloured
  dashed line.
- **Replay acquisition** — the whole sample injected as a swarm, sweeping along
  the retention-time axis as the gradient runs.
- **Flash protein** — strikes every ion of the selected protein at once, from
  the inlet outward.
- **Pick a protein, then a peptide** — its dots pulse; everything else dims to
  open rings. Unassigned features stay open rings throughout.
- **Hover or click any dot** — guide lines drop to each measurement axis, and
  the right panel draws the experimental MS2 head-to-tail against the atlas
  library entry, with the b/y coverage ladder and the dot-product score.
- **Filters** — charge state (1+ … 4+), colour by protein / isotope / charge,
  all isotopes or monoisotopic only, true or expanded isotope spacing,
  unassigned features on or off.

### The data

Synthetic, and labelled as such in the page. It is not random: masses come from
monoisotopic residue tables, retention from hydrophobicity coefficients,
cross-section from a calibrated power law in mass and charge, and drift from
Mason–Schamp. Isotope envelopes are averagine Poisson. MS2 fragments are the
real b/y series of each sequence; the experimental copy is the library spectrum
with per-peak variance, dropouts and chemical noise, and the dot product is
computed live on square-root intensities. Peptides are a mix of recognisable
tryptic sequences and sequences generated at proteome residue frequencies and
cleaved at K/R.

Everything is deterministic — the same dataset on every load.

### Loading a real export

**Load export** (or drop the files anywhere on the page) reads an export written to
`docs/export-contract.md`. `ions.jsonl` is required; `manifest.json` and `spectra.jsonl`
are optional, and other contract files are accepted and ignored. Files are recognised
by name. **Demo data** returns to the synthetic set.

What the Atlas does with it, and never hides:

- **Only real mobility goes on the vertical axis.** Ions whose `drift_kind` is not
  `ion_mobility` (SONAR, DDA) are not drawn, and the count is shown. An export with no
  mobility at all collapses to 2-D, m/z against RT.
- **Isotopes come from the precursor.** `isotopes[]` on the ion record first, then the
  raw `ms1_isotope` scan, and only failing both an averagine model, labelled *modelled*.
- **Missing spectra stay missing.** No MS2 is synthesised for a loaded export; a panel
  without data says which record is absent.
- **Thinning is disclosed.** Above 9,000 ions it draws every annotated ion plus a seeded
  sample of the rest, and says so; the manifest's own decimation ratio is shown too.
- **Sample and confidence-tier filters** appear alongside charge, so tiers need not mix.

### Depletion and interference

The method is in `docs/export-contract.md`, section 11.

- **Deplete** replays depletion. Each pass's explained ions flash, then fly back out
  through the inlet into the database, until only the unexplained remain. The HUD counts
  each pass and the remainder. Passes come from the export's `depleted` events; without
  them they are inferred from confidence tier, and the demo uses match score. The HUD
  says which. **Fill** restores the cloud.
- **Colour by Purity**: clean, MS1 envelope off, ratio drifting across the peak,
  product-ion interference, or not assessed. **Interference flags** marks flagged ions
  with a diamond in any colour mode.
- **The match panel** adds isotope fit (RMS residual against the envelope expected from
  the peptide's composition), the share of MS2 intensity the ladder leaves unexplained,
  and `ar1_sum` against the same peptide in other runs. When the export carries one trace
  per isotope, a strip plots each isotope ratio scan by scan across the peak and shades
  where it departs. The liquid-phase panel marks the expected envelope over the observed
  bars.

An export's own `isotope_fit` and `interference` always win. Without them the Atlas
estimates from what it has and labels the result an estimate. Its thresholds are first
guesses, to be tuned on real data: isotope fit above 8 points, a ratio CV above 15%,
more than 35% of MS2 intensity unexplained, or `ar1_sum` below 0.75 of the clean runs.
`ar1_sum` needs the peptide in at least three runs; with fewer it shows a dash.

## The export contract

```
docs/export-contract.md          the spec — read this first
schema/method.schema.json        the fillable method (what to plot, how to compute it)
schema/panels.json               panel catalog; each panel declares what it requires
schema/manifest.schema.json      per-export manifest (samples, counts, decimation)
schema/ion_track.schema.json     one ion, with its append-only event ledger
schema/protein_quant.schema.json Top3avg per (protein, sample)
schema/spectrum.schema.json      MS/MS and MS1 isotope peaks (the conventional panels)
schema/trace.schema.json         extracted ion chromatograms (optional)
examples/                        a filled method and a matching export
tools/plan_panels.py             says which panels will render, and why not
viewer/aria-ion-panels.html      the viewer — the six conventional panels, fed by an export
```

## Viewer

`viewer/aria-ion-panels.html` is a single self-contained page: MS/MS head-to-tail against
the Atlas library, MS1 isotope envelope, XIC, intensity vs RT, SeqStart vs intensity and
Top3avg percentile. It opens on a clearly labelled demo export and takes a real one by
drag-and-drop of `manifest.json`, `ions.jsonl`, `spectra.jsonl`, `traces.jsonl` and
`proteins.jsonl`. Each panel reports `RENDER` / `PARTIAL` / `BLANK` and states its reason
rather than vanishing. Published at https://claude.ai/artifact/EaRBSjdQ666bsfS2sVTMvx
(private until shared).

## Try it

```
python3 tools/plan_panels.py examples/method.example.json examples/manifest.example.json
```

Reports per panel: `RENDER`, `PARTIAL` (with the samples it excludes and why) or `BLANK`
(with the reason). Needs `jsonschema` for method validation; the feasibility report runs
without it.

## Panel order

Panels declare a `familiarity`: `conventional`, `bridge`, `novel`. Methods are ordered up
that ladder on purpose — the MS/MS mirror and the other views practitioners already read
come first, so geometry-based selection gets checked in a familiar frame before the 3-D
map is asked for. `plan_panels.py` flags a method that leads with a novel panel.

## Two rules worth knowing before reading the schemas

- **Drift is not always drift.** SONAR fills a drift column with quadrupole window
  positions, DDA has none. `coord.drift_kind` gates the mobility panels.
- **Top3avg is relative abundance.** Absolute concentration needs a declared spike-in
  standard. Without one, nothing is labelled as concentration.

## Working from another computer

Every computer that pushes to this project must point at the same GitHub repo:

```
https://github.com/bigtony1092-Nexo/Visuals
```

**Check where a computer is pushing.** In a terminal, inside the project folder:

```
git remote -v
```

It should print:

```
origin  https://github.com/bigtony1092-Nexo/Visuals (fetch)
origin  https://github.com/bigtony1092-Nexo/Visuals (push)
```

**If it shows a different URL,** point it here and push again:

```
git remote set-url origin https://github.com/bigtony1092-Nexo/Visuals
git push -u origin HEAD
```

**If it shows nothing,** the folder was never connected to GitHub:

```
git remote add origin https://github.com/bigtony1092-Nexo/Visuals
git push -u origin HEAD
```

The account doing the push needs write access to the repo. If GitHub refuses
with "permission denied" or "403", the repo owner adds that account under
**Settings → Collaborators** on GitHub.

**To get everyone else's work:**

```
git pull
```
