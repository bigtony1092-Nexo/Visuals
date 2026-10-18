# Does pushing more samples through make the Atlas better?

The project's premise is that Atlas quality is a function of how many samples
pass through, provided the controls at each step stay tight. That is a testable
claim, and the mass dimension can be tested right now on the 24 Colorado PAX8
runs — 1,262 PSMs, one instrument, one day.

## The run means are a tight distribution

| | |
|---|---|
| between-run SD of the per-run means | **0.643 ppm** |
| pooled within-run SD | 2.277 ppm |
| between / within | 0.283 |
| share of variance that is between-run | **7.4%** |
| range of run means | −2.873 to −0.413 ppm |

**Runs agree with each other about four times better than individual PSMs agree
within a run.** That is the shape the premise needs: the run-to-run term is the
small one, so averaging across runs is averaging a quantity that barely moves.

## And sample count buys precision at the expected rate

With run means as independent draws, the Atlas coordinate's standard error on
this dimension is `0.643 / √n`:

| runs | SE | | target | runs needed |
|---|---|---|---|---|
| 1 | 0.643 ppm | | 0.50 ppm | 2 |
| 4 | 0.322 ppm | | 0.20 ppm | 11 |
| 10 | 0.203 ppm | | 0.10 ppm | **42** |
| 24 | 0.131 ppm | | 0.05 ppm | **166** |
| 1,000 | 0.020 ppm | | 0.02 ppm | 1,035 |

So the premise holds on its own terms. Forty-odd runs put the mass coordinate
inside a tenth of a ppm. That is a real and favourable answer: this dimension
converges, and it converges cheaply.

## But precision is not accuracy, and the controls have to know the difference

| | |
|---|---|
| mean of the run means | **−1.570 ppm** |
| 95% CI on that offset | −1.841 to −1.299 |
| run means below zero | **24 of 24** |
| t | −11.96 on 23 df |

Every run is low. That is not scatter around a true value; it is a shift. And
`√n` shrinks the interval **around** the shift, never toward zero. At 1,000 runs
the coordinate would be known to 0.02 ppm and still sit 1.57 ppm off — a
confidence interval excluding the truth by about **77 standard errors**.

**More samples make a shared offset more confidently wrong.** The tightening
interval is not evidence of correctness; under a systematic offset it is the
opposite, because it removes the one visible sign that something is wrong.

This is the same signature as the drift result in `verification-run.md`: a
leave-one-out fit on a run's own ions lands at 1.69% median error while the
previous day's calibration lands at 7.41%, **with all 618 errors negative**. Two
different dimensions, one pathology.

## What this means for the controls

The premise — quality as a function of sample count, under tight per-step
controls — is sound, with one requirement that follows directly from the above:

**A control that only measures spread cannot protect the Atlas.** Scatter is the
term that averages away. The term that does not average away is a shared offset,
and it is invisible to any check on variance, repeatability or interval width.
Every such check gets *better* as the offset persists.

So the per-step controls need at least one that asks a different question:

1. **Is the mean of the run means zero?** A one-sample t against zero on the
   per-run means, not on the pooled PSMs. Here it fails at t = −11.96, which is
   the correct answer and nobody was asking.
2. **Do the per-run means all share a sign?** 24 of 24 is a sign test at
   p ≈ 1.2 × 10⁻⁷. Cheap, assumption-free, and catches exactly the failure that
   averaging hides.
3. **Does an independent reference agree?** Spread and sign tests detect an
   offset's *existence*; only an external anchor gives its size. This is the same
   requirement as export contract section 14 — a reference must be shown to carry
   signal before a comparison against it is read.

None of these is expensive. The point is that they are a different *kind* of
check from the ones already in place, and the reason the project needs them is
precisely that its premise is correct: as sample count rises, every
spread-based control reports improving health whether or not the Atlas is
converging on the right answer.

## The three checks, in code

`tools/offset_check.py` runs all three on a file of per-PSM ppm errors:

    python3 tools/offset_check.py ppm.jsonl --anchor reference.jsonl --target 0.0

It is step 6 of `tools/aria_check.py`, and it is a **gate** — not because its
finding is worse than the others, but because of the direction of its error.
A failing spread control makes the steps after it look bad. A shared offset
makes them look *good*, so a chain that carried on past it would be reporting
health computed from the pathology. "Not reached is not passed" is the only
honest state for anything downstream of an unresolved offset.

Three properties of the tool are worth stating, because each exists to stop a
specific way this check could have been written and still said nothing:

**It demonstrates the blindness rather than asserting it.** Every run injects
known offsets (0, 1, 5, 20 ppm) into the same data and reports which controls
move — not behind a flag, because a demonstration you have to ask for is one
that gets skipped. On the planted Colorado shape, **five of eight are bit-identical
at 0 ppm and at 20 ppm** — between-run SD, pooled within-run SD, range of run
means, IQR of run means, SE of the Atlas mean. Adding a constant to every run
changes no difference between runs, so this is arithmetic, not an empirical
finding; printing it as a measured result is what makes the claim checkable
instead of a sentence in a document.

One control is worse than blind. `spread / |centre|` goes 0.503 → 1.197 →
0.265 → 0.047 across those four offsets, because the centre crosses zero on
the way. It is non-monotone — wrong in both directions — and a control that
argues in the wrong direction is more dangerous than one that says nothing,
because blindness at least does not argue.

**It refuses circular input.** A ppm column computed from the same calibration
it is being used to test returns zero by construction. That already happened in
this project: 0.00 ppm over 1,994,456 rows, which read as a perfect result. The
tool rejects input that is majority exact-zero, on the grounds that a ppm of
1 × 10⁻¹² is arithmetic rather than a measurement.

**It separates existence from size.** Checks 1 and 2 can only say an offset is
there. Without `--anchor` the tool says so and does not report a magnitude —
the −1.570 ppm here is the size of the shift *relative to the mass computation
being tested*, which is not the same as the distance from truth.

## What is not established here

All 24 runs are one instrument on one day, so this separates run-to-run from
PSM-to-PSM and nothing else. **Instrument-to-instrument and day-to-day terms are
not measured** — and those are the ones likeliest to carry their own offsets. The
−1.57 ppm found here could be a property of that instrument, that day, that
calibration, or of the mass computation itself; this data cannot say which.

Running the same decomposition across instruments and across days is the natural
next test, and it is the one that would show whether a per-run calibration is
enough or whether the offsets themselves need a reference to anchor against.
