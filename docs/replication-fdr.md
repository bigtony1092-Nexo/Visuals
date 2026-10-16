# What replication is worth

MARCH's `replication` slide reports two rates and a count. Those numbers carry a stronger
claim than the slide makes, and the claim is about **what requiring a second run buys**.

## The slide's own numbers

| Evidence required | Real found again | Reversed found again |
|---|---|---|
| 3+ fragments | 42% | 13% |
| 4+ fragments | 50% | 7% |
| **5+ fragments** | **58%** | **4%** |

Plus: **373 real : 4 reversed** peptides found in both samples at 5+ fragments.

## Tier 1 — true with no extra assumption

The ratio of the two published rates:

| Evidence | real : decoy replication ratio |
|---|---|
| 3+ | 3.2 |
| 4+ | 7.1 |
| 5+ | **14.5** |

This is division on two numbers already on the slide. It is also the number that actually
shows the gate working, and it is not currently displayed.

**It matters because the slide's framing is slightly off.** Real "found again" *rises* as
the gate tightens — 42 → 50 → 58%. That reads as "stricter is better", but it is a
selection effect: peptides with more fragments are more abundant and replicate better
regardless of method. The decoy rate falling 13 → 7 → 4% is the part that is about the
method, and the ratio is where both appear at once.

## Tier 2 — true under the standard target-decoy assumption

Assuming the reversed set is the same size as the real set, the false discovery rate among
peptides that replicated:

```
FDR = decoy_rate / (real_rate + decoy_rate)
```

| Evidence | FDR among replicated |
|---|---|
| 3+ | 23.6% |
| 4+ | 12.3% |
| 5+ | **6.5%** |

This is the ordinary way the field computes FDR, so it needs no apology — but the
assumption should be stated on the slide rather than left implicit.

## Tier 3 — needs one thing confirmed before it is said aloud

**If** the percentages and the 373 : 4 counts share a denominator — that is, if 373 is 58%
of the real peptides that passed run 1, and 4 is 4% of the reversed peptides that passed
run 1 — then:

```
real passing run 1      373 / 0.58  =  643
reversed passing run 1    4 / 0.04  =  100

FDR after run 1 alone        100 / 643  =  15.5%
FDR after requiring both      4 / 373   =   1.07%
```

**Requiring a second run cuts estimated FDR about 14.5-fold, from 15.5% to 1.07%.**

That is the strongest statement available from this data, and it reframes replication from
a reassurance into a mechanism — the step that takes a false rate nobody would accept into
one that everybody does.

### What to check first

The derivation hangs entirely on those two denominators. Confirm, from the script that
produced the slide:

- Is 373 the count of real peptides passing 5+ in **both** runs, and 58% that count over
  the real peptides passing 5+ in **run 1**?
- Same question for 4 and 4%.
- Is the reversed set the same size as the real set?

If the percentages are over a different base — all candidates, say, rather than run-1
passers — the implied 643 and 100 are wrong and Tier 3 collapses. Tiers 1 and 2 survive
either way.

## Proposed slide wording

Tier 1 goes on the slide now. Tier 3 goes on it the moment the denominators are confirmed:

> Real peptides replicate 14.5× more often than reversed decoys at 5+ fragments.
> Requiring a second run takes the estimated false discovery rate from 15.5% to 1.07%.

One honest limit the slide already states and should keep: these are two runs on one
instrument, not biological replicates.

## Closed: the denominators, read out of the script

`studies/2026-10-12_name_grid/replication.py` re-run on the untouched
`HeLa_120_1` and `HeLa_60_1` annotated ion lists. Counts at the frozen 5+
fragment rule:

| | run 1 (HeLa_120) | run 2 (HeLa_60) | both |
|---|---|---|---|
| real | **643** | 1,261 | 373 |
| reversed | **93** | 381 | 4 |

So the two denominators are 643 and **93 — not the 100 the earlier estimate
assumed**. It is a real count, not a cap or a default.

### Three numbers on the slide were wrong

Every ratio in the table had been computed by dividing the **rounded**
percentages rather than the counts:

| bar | real again | reversed again | true ratio | slide said |
|---|---|---|---|---|
| 3+ | 922/2,190 = 42.1% | 152/1,181 = 12.9% | 3.3× | 3.2× |
| 4+ | 552/1,110 = 49.7% | 23/317 = 7.3% | 6.9× | 7.1× |
| **5+** | 373/643 = 58.0% | 4/93 = **4.3%** | **13.5×** | **14.5×** |

58 ÷ 4 is 14.5. 58.0 ÷ 4.30 is 13.5. The headline on the face of the slide —
the one its own speaker note called "the number to point at" — was inflated 7%
by rounding the denominator to one significant figure.

The derived figures: **93/643 = 14.5%** in run 1 alone, **4/373 = 1.07%**
requiring both runs, a **13.5-fold** cut. Not 15.5% and not 14.5-fold.

### Two caveats that have to travel with it

**It is not an FDR.** The decoy-to-target ratio is 1:1 by construction — every
named precursor gets exactly one reversed ladder at its own address — but this is
not a target-decoy database search. The reversed ladder is a *per-precursor
null*, so the honest term is an **estimated false rate**.

**The control is handicapped.** The reversed ladder is about **6.5% shorter**
than the real one: reversed b(L-1) and exact coincidences are removed, while the
real ladder is never trimmed. At a fixed 5-fragment bar that biases the reversed
count low by roughly 15–25% relative, which makes every figure here a **floor**:

| reversed run 1 | implied false rate |
|---|---|
| 93 as counted | 14.5% |
| ×1.15 → 107 | 16.6% |
| ×1.25 → 116 | 18.1% |

So the true run-1 rate is plausibly nearer **17–18%**, and the 1.07% is likewise
a lower bound. The asymmetry is a real defect in the control and should be fixed
at source — trim both ladders or neither — rather than carried as a caveat
forever.

### What the slide now says

MARCH is at version 13. The table carries counts alongside percentages with
ratios computed from the counts; the 5+ row reads **13.5×**; a highlighted line
on the slide face states that the reversed ladder is shorter and every figure is
a lower bound; the footnote records that this is a per-precursor null at 1:1 and
not a database search. The speaker note gives the defensible sentence:

> about 15% in one run, a lower bound; about 1% requiring both runs; a roughly
> 13-fold cut.

The "confirm those two denominators before saying it aloud" instruction is gone,
because they have been confirmed.
