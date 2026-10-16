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
