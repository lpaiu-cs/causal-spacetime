# Route D, first measurement: a first look, with no verdict

2026-11-25. **Frozen with the commit that adds this file under this name.**
The user approved freezing and running the second draft in this session
("동결하고 실행"). Implemented in [drift_probe.py](drift_probe.py) and
[analyse_drift.py](analyse_drift.py). Generator coordinates throughout: a
diagnostic of the background order law, not an observer and not a wave
test. AI-written. The first draft (in the history at `764ef85`) was
reviewed by an AI referee ([referee_report.md](referee_report.md)); this
design answers that review and has not itself been reviewed. The second
draft is in the history at `4220e47`; Section 9 lists the one change made
to it before this freeze.

## 1. What changed, and why there is no gate

The first draft gated a slope of the residual against `1/4` with a mixing
requirement. The referee showed that this gate cannot be trusted:

- Run on the draft's own 1+1 positive control, where a count-scale effect
  exists by construction, the mixing requirement returned "unresolved" in
  300 of 300 replays. Collective motion relaxes in about `0.6 N^2` sweeps
  there, and in 3+1 the conformal family is at most 18% relaxed at the
  gated lag. A gate that passes mixing needs chains of `4 N^2` sweeps or
  more: about 8.4 hours per order at N = 512.
- The residual is not a power law in N. Its square is a cell-scale term
  plus a count-scale term, `R^2 = A + B N`. A count-scale residual of the
  size U2 itself suggests (about `0.7 / sqrt(N)`) puts the slope on the
  decision line at these sizes, so the stated power did not hold.
- A slowly creeping residual could pass the mixing tolerance and still
  move the slope by several times the drift allowance.

So this draft drops the verdict. It measures curves and reports them with
intervals. Section 6 says in advance how the curves will be described and
what each description cannot establish. A gated design, if one is worth
its cost, comes after this look and needs its own approval.

## 2. The question (unchanged)

Given the order, how far can a central event move? Two contributions:

- **The order-preserving family.** With tips p, q and inversion
  `I(y) = y / (y . y)`, put `phi(x) = I(x - q) - I(p - q)`. The maps
  `g = phi^-1 o (lambda Lambda) o phi`, `lambda > 0`, `Lambda` orthochronous
  Lorentz, are the conformal maps of the diamond that fix both tips and
  preserve every causal relation: 7 parameters. (Referee: order matrix
  identical for 153 random elements; `phi` alone reverses the order.) Only
  the rotation about the wave axis preserves the sampling density. Of the
  other six, two (tilts of the wave axis) are nearly free at these sizes;
  the remaining four are held at about `N^(-1/2)`. Estimated family
  displacement of central events between two independent draws:
  `3.2 / sqrt(N)`.
- **The residual:** what is left after the best group element is removed.

Cell-scale residual: a fixed multiple of the single-move scale (of order
`1/N`) at every lag beyond a few sweeps. Count-scale residual: an extra
term that grows relative to that scale like `sqrt(N)`, possibly appearing
only at long lags.

## 3. Estimator, fixed here

**Chain.** The Route A chain, unchanged in law, started at the sampled
configuration (an exact draw from the conditional law given the order).
Length `4 N^2` sweeps. Seeds `SeedSequence([20261127, dimension, N, r])`.
Exact order comparison at every stored state.

**Stored states.** After `t_j = N 2^j` sweeps for `j = 0, 1, ...` while
`t_j <= 4 N^2`; after `N^2/2`, `N^2`, `2 N^2` and `4 N^2` sweeps; and one
sweep after each of these.

**Lags.** The short ladder `N 2^j` up to `64 N`, and the long lags
`N^2/2`, `N^2`, `2 N^2`.

**Event sets,** fixed from the sampled configuration and its order.
*Central*: inside the sub-diamond scaled by 0.6 about the midpoint of the
tips (as in Route A). *Outer*: not central, with at least one predecessor
and one successor among the events.

**State pairs.** For each lag t: the pairs `(0, t)` and `(t, 2 t)`; both
states are stored for every lag above. Because the chain is stationary
from its start, both are unbiased for that lag.

**Per pair (a, b).**

1. Total: mean over central events of `|x_b - x_a|^2` (Euclidean in the
   chart).
2. Fit: minimise `sum_k |g(x_a,k) - x_b,k|^2` over central events, with
   `lambda = exp(s)` and `Lambda` the exponential of the six Lorentz
   generators; `scipy.optimize.least_squares`, started at the identity,
   tolerances `1e-12`. A fit that does not report success is an error for
   that pair and is counted.
3. Residual `r_k = g(x_a,k) - x_b,k`. Mean square residual of central
   events, multiplied by `4 n / (4 n - 7)` for the seven fitted parameters
   (n central events). The same for the wave coordinate `(T + Z) / sqrt(2)`
   alone.
4. Out of sample: mean square of `g(x_a,k) - x_b,k` over outer events with
   the same fitted g, and their total mean square displacement.
5. Pair differences: mean of `|r_k - r_l|^2` over pairs of central events
   in three bins of their chart distance at state a (below 0.25, 0.25 to
   0.5, above 0.5).
6. The fitted parameters.

**Unit.** `s1^2`, per order: the mean over central events and over stored
states of the squared displacement during the one sweep after that state.
`R(t) = sqrt(mean square residual at lag t / s1^2)`, per order.

**Intervals.** Per-order values, mean over orders, two-sided 95% interval
from the t distribution with (orders - 1) degrees of freedom.

**Family thermometer.** The second-derivative matrix of the log conditional
density along the seven generators, computed without a chain at each
sampled configuration, predicts the spread of the four stiff parameters.
The ratio of the observed variance of the fitted parameters at lag `2 N^2`
to that prediction is reported as the fraction of the family that has
relaxed.

## 4. Controls with known answers

1. **Fit recovery, no chain.** At N = 128, 192, 256, 24 configurations
   each: apply a random group element that moves central events by
   `3.2 / sqrt(N)` rms, add independent noise of rms `10 / N` per event to
   both states, fit. Recovered over injected residual: 95% interval inside
   `1 +/- 0.05`. A group element alone must leave the order matrix
   identical. (The referee's own test of this gave 0.98 to 1.00.)
2. **The chain finds count-scale collective motion.** 1+1 dimensions, flat
   density, N = 64, 96 orders, chains of `8 N^2` sweeps. The exact mean
   square difference between two independent draws, averaged over events,
   is `1 / (3 (N + 1))` in the chart `(T, Z)`. The ratio of the all-event
   mean square displacement at lag `4 N^2` (pairs `(0, 4N^2)` and
   `(4N^2, 8N^2)`) to that value must have a 95% interval inside
   `1 +/- 0.2`. The whole curve against lag is reported; the referee
   measured 0.57, 0.72, 0.87 and 0.97 at lags `N^2/4` to `2 N^2`, so the
   control also shows how much a short lag under-reports.
   Sizing: relative standard deviation of the per-order value 0.47 (seen
   in development at N = 32, Section 9), 96 orders, standard error 0.048
   against the requirement `0.2 / 3.2415 = 0.062`. A larger realised
   standard error leaves the control unresolved.

If a control fails, the main stage is not run.

## 5. Campaign

CPU only; twelve single-thread workers at below-normal priority, recorded;
clean-tree preflight, commit stamp, no overwriting. Twelve orders per size.

| stage | what | sizes | estimated wall time |
|---|---|---|---|
| 0 | controls | fit: 128, 192, 256; chain: 1+1 at 64 | 15 min |
| 1 | the curves | 3+1 at 128, 192, 256, chains of `4 N^2` | 80 min |

N = 512 is left out: `4 N^2` sweeps there is about 8 hours per order.
Central events per order are about 18, 27 and 39.

## 6. How the curves will be described. None of this is a verdict on U2

Reported for each size, with intervals: `R(t)` at every lag; the same for
the wave coordinate; the out-of-sample residual in the same unit; the
pair-difference curves; total displacement of central events; fraction of
the family relaxed.

Three descriptions are fixed in advance.

- **Plateau.** At a size, the ratio of R at lag `2 N^2` to R at lag `16 N`
  has a 95% interval upper end of at most 1.10. *Cannot establish:* that no
  slower component exists, or that the chain reaches every configuration
  with the order. It is "no residual growth seen up to `2 N^2` sweeps by
  this instrument".
- **Creep.** That ratio has a lower end above 1.10. *Establishes* a lower
  bound: the residual at lag `2 N^2` is at least what was measured, in
  units of the single-move scale. *Cannot establish* its equilibrium size.
- Anything between is reported as neither.
- **Size dependence.** The least-squares fit `R^2 = A + B N` over the three
  sizes, at lag `16 N` and at lag `N^2 / 2`, with intervals for A and B,
  and the count-scale coefficient `c = s1 N sqrt(B)` it implies (residual
  about `c / sqrt(N)`), next to the smallest c the data could have told
  from zero. Three sizes within a factor 2 make this a weak fit; it is
  reported as such.

The pair-difference curves are described as flat in separation, or rising
with it, with the same interval rule on the ratio of the farthest to the
nearest bin.

## 7. Limits, stated before any data

- No verdict. A plateau is not evidence that positions are fixed at the
  cell scale; a creep is a lower bound, not a measured count-scale size.
- At lag `2 N^2` the family itself may be far from relaxed in 3+1. The
  thermometer reports how far.
- The fit is blind to a common shift of the central cluster (it keeps 3%
  of one). The out-of-sample residual is the check on that.
- Central events only, one background, one diamond, sizes 128 to 256,
  where the family displacement of central events is about 4 to 6 times
  the single-move cell radius.
- With 18 central events at N = 128 the fit estimates seven parameters from
  72 numbers.
- Whatever is seen is about positions given the order on the background.
  It says nothing about a wave, a statistic or the exponent.

## 8. What has been seen that bears on this

The Route A records and their independent check; the referee's chain-free
numbers quoted above; the referee's 1+1 relaxation table. No displacement
of central events along a 3+1 chain has been recorded or seen.

## 9. Changes to the second draft, and what was seen before this freeze

**One change.** The 1+1 chain control has 96 orders instead of 48. The
second draft sized it on an assumed spread; the development run below gave
a per-order relative standard deviation of 0.47, with which 48 orders would
have left the control underpowered by its own rule. Nothing else in the
measurement, its estimator or its descriptions was changed.

**Seen in development, before the freeze.**

- Fit recovery, four configurations per size: 1.003, 0.997, 1.006 at
  N = 128, 192, 256; order matrices identical; every fit succeeded. A group
  element is undone by the fit to `1.5e-14`.
- 1+1 chain at N = 32, six orders: ratio to the exact value 0.45, 0.57,
  0.83, 1.05, 1.19 at lags `N^2/4` to `4 N^2`; standard error 0.19 at the
  last.
- 3+1, to see that the code path runs: two orders each at N = 40 and 48
  with the full schedule. Printed: 80 state pairs, every fit succeeded, 7
  to 12 seconds per order. No displacement, residual or summary from these
  runs was printed or read.

No displacement of central events along a 3+1 chain has been seen.
