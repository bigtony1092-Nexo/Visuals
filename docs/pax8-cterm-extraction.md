# PAX8 C-terminal extraction

**The question.** Colorado reports no PAX8 ion past residue 287 in any of ten runs. The
Atlas reports both peptides that exist there. Did Colorado's instrument not see them, or
did its search never report them? Search output cannot tell those apart — an unreported
peptide and an absent peptide look identical in a results table. Only going to the
coordinate in the raw data separates them.

Method: `examples/method.pax8-cterm.json`.

## Why there are only two targets

PAX8 is 450 aa. Residues 288–450 hold **3 cleavage sites** against 29 in 1–287 — 1.8 per
100 aa against 10.1. Trypsin yields 4 peptides there, and two are 54 and 56 residues, past
the 6–50 window. So the entire C-terminal third offers exactly two observable peptides,
and the Atlas detects both.

| | Residues | Length | Atlas | Colorado |
|---|---|---|---|---|
| `NLSTHQTYPVVADPHSPFAIK` | 288–308 | 21 | detected | none in 10 runs |
| `FPNSSLLSSPYYYSSTSRPSAPPTTATAFDHL` | 419–450 | 32 | detected | none in 10 runs |

The 288 peptide carries the argument: fully tryptic, ordinary length, no Cys, and missed
at 2500 ng-equivalent on column. The 419 peptide is the protein's C-terminus with no
tryptic C-terminal residue, so a null there is weaker evidence.

## Coordinates

Computed from monoisotopic residue masses; carbamidomethyl applied to Cys where present.
Neither target contains Cys, so neither carries fixed-modification ambiguity.

| Peptide | Neutral | 2+ | 3+ | 4+ |
|---|---|---|---|---|
| `NLSTHQTYPVVADPHSPFAIK` | 2321.1804 | 1161.5975 | 774.7341 | — |
| `FPNSSLLSSPYYYSSTSRPSAPPTTATAFDHL` | 3461.6415 | — | 1154.8878 | 866.4177 |
| `LLAEGVCDNDTVPSVSSINR` + CAM | 2145.0372 | 1073.5259 | 716.0197 | — |
| `TDAFSQHHLEPLECPFER` + CAM | 2212.0007 | 1107.0076 | 738.3408 | — |
| `LSIDSQSSSSGPR` | 1319.6317 | 660.8231 | 440.8845 | — |

These reproduce Colorado's own recorded values for the controls (1073.53, 738.34, 660.82),
which is the check that the target coordinates are right too.

RT and drift come from each peptide's Atlas entry, mapped through the run's own
calibration. Per the drift finding, that calibration is **fitted per run** — a previous
run's calibration, even on the same instrument, lands off and off in one direction.

## The controls are not optional

Three peptides Colorado detects at n=24 are extracted alongside the targets. Without them
a null result is uninterpretable: a failed extraction and an absent peptide produce the
same empty panel. `LSIDSQSSSSGPR` is included specifically because it is serine-rich like
the C-terminal region, so it also controls for that composition.

**If the controls come back empty, stop.** The extraction is broken and the targets say
nothing.

## Reading the result

With the controls clean:

| What the targets show | Conclusion | Next |
|---|---|---|
| XIC peak, isotope envelope at 1.00335/z, fragments co-eluting | **Present in the raw data; the search did not report it** | A search-space problem — enzyme spec, length bounds, modification set |
| Precursor peak but no fragments | Present but never selected for MS2 | An acquisition problem, not a chemistry one |
| Nothing at the coordinate in any of 10 runs | **Genuinely absent** | Chemistry. Then the modification hypothesis earns its cost |
| Signal only in the 2500 ng-equivalent runs | Sensitivity after all | Revisit the loading conclusion |

Only the third row justifies a modification search, and then the residues choose it: the
288–450 stretch is **39% S+T+P and 20% serine**, so phosphorylation, not the N-terminal
modifications that do not apply to a C-terminal gap.

## Tier

`min_tier` is `bare_proximity`, deliberately. The question is presence in the raw data, not
confirmed identification, and the proximity tier is the one that answers it. Panels needing
a decoy-controlled tier will blank with their reason, which is the correct behaviour — this
method is not evidence of identity and should not be read as any.
