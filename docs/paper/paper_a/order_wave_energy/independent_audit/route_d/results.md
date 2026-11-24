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
  flat. So the fitted map is a map of the whole configuration, and the fit
  has not hidden a common shift of the central cluster.
- **Family relaxed** is the squared fitted amplitude of the four stiff
  group parameters over its Laplace prediction. It reaches 1 by about
  `N^2/4` and stays there within wide intervals (about +/- 0.5). The three
  soft parameters (the rotation about the wave axis and the two tilts) have
  second derivatives of `-0.02` to `0.04` per event against 1.1 to 2.5 for
  the stiff ones; they keep the total displacement growing.
- **Pair differences.** Residual differences between central events more
  than 0.5 apart are 1.07 to 1.23 times those between events less than
  0.25 apart, at every lag and size, with no trend in lag or N. Nearby
  residuals are mildly correlated; nothing accumulates with separation.
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
