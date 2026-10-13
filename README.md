# Visuals

Aria plotting.

Panels built from annotated ion lists, driven by a saveable **method** file, fed by an
**export contract** that keeps a forensic track of every ion from detection to Atlas.

## Layout

```
docs/export-contract.md          the spec — read this first
schema/method.schema.json        the fillable method (what to plot, how to compute it)
schema/panels.json               panel catalog; each panel declares what it requires
schema/manifest.schema.json      per-export manifest (samples, counts, decimation)
schema/ion_track.schema.json     one ion, with its append-only event ledger
schema/protein_quant.schema.json Top3avg per (protein, sample)
examples/                        a filled method and a matching export
tools/plan_panels.py             says which panels will render, and why not
```

## Try it

```
python3 tools/plan_panels.py examples/method.example.json examples/manifest.example.json
```

Reports per panel: `RENDER`, `PARTIAL` (with the samples it excludes and why) or `BLANK`
(with the reason). Needs `jsonschema` for method validation; the feasibility report runs
without it.

## Two rules worth knowing before reading the schemas

- **Drift is not always drift.** SONAR fills a drift column with quadrupole window
  positions, DDA has none. `coord.drift_kind` gates the mobility panels.
- **Top3avg is relative abundance.** Absolute concentration needs a declared spike-in
  standard. Without one, nothing is labelled as concentration.
