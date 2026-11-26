# Route D, first measurement: results

2026-11-25. Design frozen at `4bdead6` ([design.md](design.md)); control
stage run at `4bdead6`, main stage at `265d73a`. **A first look. No
verdict on U2 was defined and none is recorded.** AI-run; the instrument
code has not been reviewed by a person. Generator coordinates: a diagnostic
of the background order law, not an observer.

## 1. What was seen

For central events at N = 128, 192 and 256, with the order held exact:

- The displacement between two chain states keeps growing with the lag,
  from about 4 to about 11 times the one-sweep displacement.
- After the best element of the seven-parameter order-preserving conformal
  group is removed, what is left does not grow: it stays near 4.3 times
  the one-sweep displacement from a lag of `16 N` sweeps to `2 N^2` sweeps,
  at all three sizes.
- In chart units the residual at lag `2 N^2` is `12.4 / N`, `12.4 / N` and
  `12.5 / N`.

By the descriptions fixed in the design this is a **plateau** at each of
the three sizes, and **no size dependence** was found. Section 4 says what
that does not establish.

## 2. Controls (passed)

| control | result | requirement |
|---|---|---|
| fit recovery, N = 128 | 0.998, 95% interval [0.989, 1.008] | inside 1 +/- 0.05 |
| fit recovery, N = 192 | 0.994 [0.986, 1.003] | |
| fit recovery, N = 256 | 0.997 [0.991, 1.004] | |
| 1+1 chain, lag `4 N^2`, 96 orders | 1.078 [0.965, 1.192], SE 0.057 | inside 1 +/- 0.2, SE at most 0.062 |

A group element alone left every order matrix identical. The 1+1 curve was
0.57, 0.76, 0.87, 1.00 of the exact value at lags `N^2/4` to `2 N^2`: a
chain under-reports collective motion at short lags and reaches it by
about `2 N^2`.

## 3. Curves

R is the root-mean-square residual of central events divided by their
one-sweep displacement `s1`; twelve orders per size; 95% t intervals.
`s1 N` was 2.91, 2.88, 2.93.

| N | lag | total / s1 | R | R, wave coordinate | R, outer events | family relaxed |
|---|---|---|---|---|---|---|
| 128 | `N` | 4.02 | 3.78 [3.51, 4.05] | 1.70 | 5.92 | 0.21 |
| 128 | `16 N` | 5.59 | 4.37 [4.02, 4.72] | 1.91 | 7.00 | 0.84 |
| 128 | `N^2/2` | 7.38 | 4.36 [4.04, 4.69] | 2.13 | 6.95 | 1.62 |
| 128 | `2 N^2` | 10.66 | 4.27 [3.76, 4.77] | 2.15 | 6.55 | 1.18 |
| 192 | `N` | 4.46 | 4.25 [3.92, 4.57] | 1.83 | 6.82 | 0.14 |
| 192 | `16 N` | 6.04 | 4.38 [4.08, 4.68] | 2.02 | 7.28 | 0.86 |
| 192 | `N^2/2` | 8.35 | 4.43 [4.14, 4.72] | 2.06 | 7.39 | 1.02 |
| 192 | `2 N^2` | 11.38 | 4.32 [3.98, 4.65] | 1.99 | 7.50 | 1.05 |
| 256 | `N` | 4.12 | 3.92 [3.65, 4.18] | 1.80 | 7.29 | 0.09 |
| 256 | `16 N` | 5.72 | 4.19 [3.82, 4.56] | 1.93 | 7.64 | 0.64 |
| 256 | `N^2/2` | 8.61 | 4.23 [3.83, 4.64] | 1.83 | 7.45 | 1.10 |
| 256 | `2 N^2` | 11.89 | 4.26 [3.95, 4.57] | 1.99 | 7.68 | 1.10 |

All lags are in [evidence/main/summary.json](evidence/main/summary.json).

- **Outer events** (not central, not extremal) were not used in the fit.
  Their total displacement grows from about 8 to 15, 18 and 19 in the same
  unit, and their residual under the map fitted on the central events stays
  flat. So the fit has not hidden a shift of the central cluster relative
  to the outer events. It could still hide a shift common to all interior
  events: the independent check found such a shift is absorbed at 96 to
  97% and leaves only 9 to 16% in the outer residual (Section 6).
- **Family relaxed** is the squared fitted amplitude of the four stiff
  group parameters (the dilation and three boosts) over its Laplace
  prediction. It reaches about 1 by `N^2/4` and stays there within wide
  intervals (about +/- 0.5). The independent check calibrated it: fully
  relaxed draws read between 0.9 and 2.0, so this shows the stiff
  parameters near their predicted spread only to within a factor of about
  1.5 to 2. The three soft parameters (the rotation about the wave axis and
  the two tilts) have second derivatives of `-0.02` to `0.04` per event
  against 1.1 to 2.5 for the stiff ones; they keep the total displacement
  growing.
- **Pair differences.** The root-mean-square residual difference between
  central events more than 0.5 apart, over that between events less than
  0.25 apart, is 1.07 to 1.23 at every lag and size, with no trend in lag
  or N. By the design's interval rule at lag `2 N^2` this is "neither" at
  N = 128 ([1.09, 1.28]) and "rising with separation" at N = 192
  ([1.13, 1.26]) and 256 ([1.13, 1.23]). The independent check found that
  most of it is not correlation: the residual variance of an event rises
  2.4-fold from the centre to the edge of the central set, which alone
  gives a ratio close to the observed one. A weak excess covariance between
  near pairs remains.
- The residual `12.4 / N` is an average over the central set; per event it
  runs from about `9 / N` at the centre to `15 / N` at the edge.
- No fit failed and no move was skipped.

## 4. The descriptions fixed in the design, and what they do not establish

**Plateau** at every size: the ratio of R at lag `2 N^2` to R at lag
`16 N` is 0.97 [0.91, 1.04], 0.99 [0.92, 1.06], 1.02 [0.95, 1.09]; the
rule was an upper end of at most 1.10. In the design's words this is "no
residual growth seen up to `2 N^2` sweeps by this instrument". It does
**not** establish that no slower component exists, or that the chain
reaches every configuration with the given order.

**Size dependence** `R^2 = A + B N`: `B = -0.012` [-0.041, 0.018] at lag
`16 N` and `-0.008` [-0.037, 0.022] at lag `N^2/2`. The implied count-scale
coefficient is zero with an upper end of 0.39 and 0.43 (residual
`c / sqrt(N)` in chart units); the smallest value the data could have told
from zero is 0.50. Three sizes within a factor 2 make this a weak fit, as
the design said. A count-scale residual as large as `0.5 / sqrt(N)` that is
present at these lags would have been visible; a smaller one, or one that
needs more than `2 N^2` sweeps to appear, would not.

**The wave coordinate alone is not shown flat.** The same ratio for the
residual along `(T + Z) / sqrt(2)`, computed by the independent check, has
intervals [1.00, 1.23], [0.89, 1.10], [0.94, 1.20] at the three sizes and
1.06 [0.998, 1.122] pooled. By the design's rule that is "neither" at two
sizes. It is the coordinate that carries the wave phase.

**How small a late count-scale term would have shown.** From an injection
test in the independent check: a non-family field of `0.5 / sqrt(N)`
appearing after lag `16 N` would have broken the plateau at every size;
the smallest that would show is about `0.33 / sqrt(N)` per size and
`0.23 / sqrt(N)` pooled. A term already present at lag `16 N` is bounded
only by the three-size fit above.

Saturation of the residual takes about `8 N` to `16 N` sweeps, not a few:
R at lag `2 N^2` over R at lag N is 1.13 and 1.09 at N = 128 and 256.

Not established by anything here:

- that positions are fixed at the cell scale up to the family. The reading
  is one-sided: a chain can under-report, as the 1+1 control shows it does
  at short lags;
- anything at sizes above 256, or for non-central events beyond the
  out-of-sample check;
- anything about a wave, a statistic or the exponent. The heuristic step
  from "cell-scale residual" to "phase known up to a few global parameters"
  in next_question.md stays heuristic. Three of those parameters are nearly
  free at these sizes.

## 5. Run record

- Control stage: stamp `4bdead6`. Main stage: stamp `265d73a`, clean tree,
  2026-11-24T15:24Z to 16:49Z (85 minutes).
- Twelve workers, one compute thread each, runner and workers at Windows
  priority class below normal (read back as 16384 and written to the
  manifests). Checked once from outside at the start of the main stage:
  twelve workers below normal, 11.7 of 16 logical processors (73%).
- Seconds per order: 364, 1274, 3454.
- `records.jsonl` SHA-256
  `eb36653fe5d7fa9efe0e80c3191a0add1ca5bd9a9fcff2f6ad4503103b0fafd7`.
- Nothing in the design was changed after the first 3+1 number was seen.
- Two departures from the design text, both found by the independent
  check and neither affecting a reported number: at N = 192 the states
  after `256 N` and `512 N` sweeps were not stored (no lag uses them); the
  pair-difference ratio is reported as a root-mean-square ratio, which the
  design did not specify.

## 6. Independent check

[checker_report.md](checker_report.md). An AI agent in a fresh context,
given the design, the instrument and the raw records, and told not to open
this record or any summary before it had its own numbers. The first
version of Sections 3 and 4 of this record was corrected after it: the
statements on the outer events, the thermometer, the pair differences, the
wave coordinate and the saturation time.

- **Every number reproduces** from the raw records, to `1e-15`. Sixty
  recorded rows were regenerated from their seeds with the checker's own
  mask, group map and optimiser, to `1.9e-7`.
- **No defect in the instrument.** The group map agrees with an independent
  coding and preserved the order in 120 trials; the fit has no local
  minima in 96 multi-start trials and never returned the identity; the
  central and outer sets, the schedule, the unit and the random streams
  are as designed. Fresh seeds at N = 64 and 128 gave the same short-lag
  values.
- **Gaps it found,** now in Sections 3 and 4: a shift common to all
  interior events is invisible to the residual, to the out-of-sample check
  and to the thermometer; the pair-difference bins confound a variance
  gradient with correlation; the thermometer is calibrated only to a factor
  of about 1.5 to 2; the wave coordinate alone is not shown flat.
- **Its opinion:** "about `12 / N` from `16 N` to `2 N^2`" is a fair
  description (flat within 4% over three decades of lag, pooled log-log
  slope -0.001 [-0.013, 0.012]) while the total displacement doubles and is
  tracked out of sample. The strongest alternative is non-ergodicity:
  count-scale freedom reachable only through joint moves that single-event
  line moves do not make in `2 N^2` sweeps. Against it, the stiff family
  parameters, which are also collective and count-scale, do relax in the
  chain.
- **Its proposal for the next measurement:** a chain-free joint linear
  programme, with all events free and the seven family parameters pinned,
  maximising the displacement of a central event at N = 128 to 512, to see
  whether the joint slack scales like `1/N` or `N^(-1/2)`. It is a
  proposal; nothing has been designed or run.

This is AI review. It did not re-verify the chain law, run a chain above
N = 128, or test connectivity of the set of configurations with a given
order.
