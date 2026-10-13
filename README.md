# Visuals

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
