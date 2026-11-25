# Route D, U2 transport measurement. DRAFT, NOT TO BE FROZEN AS WRITTEN

2026-11-26. **Draft. Not frozen. Not implemented. Not run.** The user chose
this measurement as the next draft ("U2 수송 측정", relayed by the session
that started this one). AI-written. **An independent AI review found it
unfit to freeze as written** ([transport_referee_report.md](transport_referee_report.md);
summary in Section 11). The text below is the reviewed draft, committed at
`7131252` and unchanged apart from this paragraph and Section 11.
Generator coordinates throughout: a diagnostic of the background order law,
not an observer and not a wave test.

The cost estimate of the reviewed draft was about 6 to 7 hours of wall time
at about 75% CPU (Section 5). The referee found it consistent, at about
6.3 to 6.5 hours.

## 1. The question

U2 in [next_question.md](../next_question.md) Section 3: how finely does
the order fix the **difference of the position of two events at
separation `d`**, beyond the seven-parameter order-preserving family?

- **Count-like transport:** count increments give a posterior spread of
  the difference that grows like `sqrt(d / N)`.
- **Cell-like transport:** the spread stays at the cell scale (`~1/N`) at
  every `d`.

Route D's stop rule uses this question. Route C now depends on it too
([../route_c/RECORD.md](../route_c/RECORD.md)): an event's walls sit at
macroscopic distance along its cones, so pinning them against each other is
transport across the diamond.

## 2. What exists, and why it does not answer this

- **The first measurement and its extension**
  ([results.md](results.md), [extension_results.md](extension_results.md))
  measured the residual of **central** events after the family is fitted
  on them. It stays near `12.5 / N` up to N = 512. A count-scale term in it
  above about `0.26 / sqrt(N)` would have been seen at those lags.
  - That bounds correlated wander only **within** the central sub-diamond,
    which is about 0.6 of the diamond's scale.
  - At N = 512 it reached only lag `N^2/2`.
- **Pair differences of central events** grew from the nearest to the
  farthest bin. Most of that is a variance gradient, not correlation: an
  event's residual variance rises 2.4-fold from the centre of the central
  set to its edge.
- **Outer events,** under the map fitted on central events, had residuals
  of 6.95, 7.39, 7.45 and 8.46 one-sweep displacements at lag `N^2/2`
  (N = 128, 192, 256, 512). That is about `20 / N` rising to `25 / N` in
  chart units.
  - It mixes transport with the outer events' own looser cells.
  - No description was fixed for it.
  - It was read after the records were seen.
- **The chain under-reports collective motion at short lags.** In 1+1 it
  showed 0.57, 0.76, 0.87 and 1.00 of the exact value at lags `N^2/4`,
  `N^2/2`, `N^2` and `2 N^2`. Transport is the slowest kind of collective
  motion, so a flat reading is at risk of being an artefact of the chain.

This draft therefore adds four things: separations up to the diamond's
scale; a readout that separates correlation from a variance gradient; the
wave coordinate on its own; and a positive control in which count-scale
transport exists exactly.

## 3. Estimator, to be fixed at freeze

**Chain.** Route D's chain, unchanged in law: the Route A order-preserving
line-move chain, started at the sampled configuration, which is an exact
draw from the conditional law given the order. Seeds are
`SeedSequence([TAG, 4, N, r])` with a new tag, so no order of the earlier
runs is reused. The order is compared exactly at every stored state.

**Length and lags.**
- Chain length: `2 N^2` sweeps.
- States are stored at `0, 16 N, 32 N, N^2/2, N^2, 2 N^2`, and one sweep
  after each of these (for the unit).
- Lags `16 N`, `N^2/2` and `N^2` use the two state pairs `(0, t)` and
  `(t, 2 t)`. Lag `2 N^2` uses the single pair `(0, 2 N^2)`.

**Event sets,** fixed from the sampled configuration and its order:
- *Central*: as before.
- *Outer*: not central, with at least one predecessor and one successor.
- *Interior*: central together with outer.

**Group fit.** The seven-parameter fit of [design.md](design.md) Section 3,
now over **interior** events: unweighted least squares, started at the
identity, tolerances `1e-12`. A failed fit is counted.

**Per state pair `(a, b)`:**
1. Residuals `r_k = g(x_a,k) - x_b,k` for interior events k, plus the
   wave-coordinate component `r_k^m = (r_k,0 + r_k,1) / sqrt(2)`.
2. For every pair (i, j) of interior events:
   - the chart distance `d_ij` at state a;
   - `S_ij = |r_i - r_j|^2` and its wave-coordinate version `S_ij^m`;
   - the independence baseline `I_ij = |r_i|^2 + |r_j|^2`, and likewise
     `I_ij^m`.
3. Bins of `d`: `[0, 0.25)`, `[0.25, 0.5)`, `[0.5, 0.75)` and
   `[0.75, inf)`. Call the first "near" and the last "far". Per bin, the
   means `S(d)` and `I(d)`, and the **correlation ratio**
   `X(d) = S(d) / I(d)`:
   - `X = 1` when the residuals of the two events are uncorrelated,
     whatever their variances;
   - `X < 1` when they move together.

   So `X` removes the variance gradient that confounded the pair
   differences of the first measurement.
4. The total displacement and the fitted parameters, as before.

Per order, average the per-pair quantities over the state pairs at each
lag. `S` and `I` are reported in chart units multiplied by `N^2`, which is
the cell unit; no per-order `s1` enters the size fit (Section 6). The
one-sweep unit `s1` is still recorded so that the curves can be compared
with the earlier measurements.

**Intervals.** Per-order values, mean over orders, and a two-sided 95%
interval from the t distribution.

## 4. Controls

The main stage is run only if C1, C2 and C3 all pass.

**C1. Fit recovery, no chain.** As in design.md Section 4.1, but with the
fit over interior events. At N = 128, 256, 512, with 24 configurations
each: apply a random group element that moves central events by
`3.2 / sqrt(N)` rms, add independent noise of rms `10 / N` per event to
both states, and fit. Pass if:
- recovered over injected residual has a 95% interval inside
  `1 +/- 0.05`;
- every order matrix is unchanged by the group element alone;
- every fit succeeds.

**C2. Readout recovery with an injected transport field, no chain.**
- **Inputs.** At N = 128, 256, 512 take 24 sampled configurations each.
  To state b add a group element as in C1, cell-scale noise as in C1, and
  a Gaussian field `tau` whose structure is exactly count-like:
  `E |tau(x) - tau(y)|^2 = kappa^2 |x - y| / N`, summed over the four
  coordinates. This is Lévy's Brownian field: per coordinate, covariance
  `(kappa^2 / 8 N)(|x - o| + |y - o| - |x - y|)` about the midpoint o of
  the tips, sampled exactly by a Cholesky factor at the event positions.
  Use `kappa` in `{0.3, 1}`, with 24 field draws per configuration.
- **Expected value.** For small fields the fit acts linearly. The expected
  post-fit `S` and `I` in each bin follow without sampling: project the
  field's covariance off the seven tangent directions of the group at
  `x_a`.
- **Pass** if, for the far bin and for the near bin, recovered over
  expected has a 95% interval inside `1 +/- 0.10`.
- **Reported, not gated:** the fraction of the injected far-bin structure
  that the fit absorbs, and the description of Section 6 applied to the
  synthetic data.
- **Sizing.** Under the house rule, with a bias allowance of 0.02 the
  requirement is `SE <= 0.08 / 3.24 = 0.025`. With an assumed per-draw
  relative spread of 0.5, that needs 400 draws; `24 x 24 = 576` are
  planned per size and `kappa`. If the realised SE is larger, C2 is
  unresolved and the main stage is not run.

**C3. The chain finds count-scale transport where it exists exactly
(positive control).**
- **Setting.** 1+1 dimensions, flat density (`FLAT2`), chains of `8 N^2`
  sweeps, at N = 64 with 128 orders and N = 128 with 64 orders.
- **The exact answer.** A two-dimensional order fixes the ranks of the two
  null coordinates. Given the ranks, the coordinates are uniform order
  statistics, which is the model the frozen 1+1 control already relied on
  (its exact value `1 / (3 (N + 1))` follows from it). For two
  independent draws from the conditional law, the exact
  `E S_ij = [p_u (1 - p_u) + p_v (1 - p_v)] / (N + 2)` in the chart of
  `FLAT2`. Here `p_u` is the rank distance of i and j in `u = (T+Z)/sqrt 2`,
  divided by `N + 1`, and likewise `p_v`. For small `p` this grows
  linearly with separation: exactly the `sqrt(d / N)` transport of U2.
- **Readout.** As in Section 3, with no group fit, since the 1+1
  conditional law has no such family. The bins are scaled to the 1+1
  diamond, whose height is 1.
- **Reported at every lag:** the recovered fraction
  `F(t) = sum over far pairs of S_ij / sum over far pairs of E S_ij`, the
  same for the near bin, and the description of Section 6.
- **Pass, at N = 64 and lag `N^2`** (the main stage's primary lag):
  - the Section 6 description of the 1+1 data is "growing with d";
  - the 95% interval of `F(N^2)` has its lower end at or above 0.5.
- **Power.** The all-event value reached 0.87 at lag `N^2`, with a
  per-order relative spread of 0.42. Assume the far bin reaches 0.7 with
  a relative spread of 0.6. Then 128 orders give SE 0.037, an expected
  lower end of 0.63, and a power of about 0.9 against the line at 0.5.
- **N = 128.** Reported only. It shows whether the fraction at a fixed
  `lag / N^2` drifts with N.

What C3 does **not** show: that the 3+1 chain mixes collective transport
as fast as the 1+1 chain does. It is a necessary check, not a sufficient
one.

## 5. Campaign and cost

CPU only. Twelve single-thread workers at Windows priority class below
normal, recorded in the manifests and checked once from outside at the
start. Clean-tree preflight, commit stamp, no overwriting. The time
estimates use the measured cost of a sweep in Route D's records, about
`5.2e-5 N` seconds (N = 128 to 512).

| stage | what | per order | orders | wall time (12 workers) |
|---|---|---|---|---|
| 0 | C1, C2: chain-free, 3+1 at 128, 256, 512 | seconds | 24 configurations each | about 10 min |
| 0 | C3: 1+1, `8 N^2` sweeps, N = 64 | 76 s | 128 | about 15 min |
| 0 | C3: 1+1, `8 N^2` sweeps, N = 128 | about 10 min | 64 | about 55 min |
| 1 | main, 3+1, `2 N^2` sweeps, N = 128 | about 3 min | 24 | about 6 min |
| 1 | main, N = 256 | about 29 min | 24 | about 1 h |
| 1 | main, N = 512 | about 3.9 h | 12 | about 3.9 h |

- **Total:** about 6 to 7 hours of wall time at about 75% CPU, of which
  stage 0 is about 1.3 hours.
- Taking 24 orders at N = 512 instead of 12 adds about 3.9 hours.
- N = 1024 at `2 N^2` sweeps would take about 31 hours per order and is
  left out.
- About 71 central and 236 outer events per order are expected at
  N = 512, so about 45 000 interior pairs.

## 6. How the results will be described. None of this is a verdict on U2

Reported for each size and lag, with intervals:
- `S(d) N^2`, `I(d) N^2` and `X(d)`, in all four bins, for the full
  residual and for the wave coordinate;
- the total displacement and the fitted parameters, as before.

Three descriptions are fixed in advance.

- **Growing with d.** At a size and lag, the ratio `Y = X(far) / X(near)`
  has a 95% interval with lower end above 1.10.
  - *Establishes:* residuals of near pairs move together more than those
    of far pairs. That is the signature of transport, and in size it is a
    lower bound.
  - *Cannot establish* the equilibrium size of the transport.
- **No transport growth seen by this instrument.** The upper end of `Y`
  is at most 1.10.
  - This is the wording **whatever C3 returns**. The chain can
    under-report collective motion, and C3 shows only how much in 1+1.
  - *Cannot establish* that positions are pinned against each other at
    the cell scale.
- **Neither:** anything between.
- **Size dependence of the far bin.**
  - Fit `S_far N^2 = A + B N` by least squares over all orders at
    N = 128, 256, 512, at lag `N^2`. Do the same for the wave coordinate
    and for `I_far`.
  - The implied transport coefficient is
    `kappa = sqrt(max(B, 0) / dbar_far)`, where `dbar_far` is the mean
    separation in the far bin, so that `S_far ~ kappa^2 dbar_far / N`.
    Report it with its upper end, next to the smallest `kappa` the data
    could tell from zero.
  - Describe it as "count-scale transport found" if B's interval lies
    above zero, and "none found" otherwise.
  - The same fit on `I_far` separates loose outer cells from transport.
  - The heuristic `sqrt(d / N)` has an unstated constant of order one, so
    a bound on `kappa` is reported as a number, not as a verdict.

Following the extension's review, the bound to quote is the smallest
distinguishable `kappa`, not an upper end made small by a negative B. If
the far-bin spreads differ by size as they did before, a weighted fit is
reported beside the unweighted one. The unweighted fit is the primary.

**Sensitivity expected from the existing records.**
- **Starting point.** The far-pair difference of central events was flat
  at about 42 one-sweep units squared, or about 360 in units of `N^-2`
  (chart). Its per-order spread was about 90 to 150 at N <= 256 and about
  37 at N = 512.
- **Projection.** With 24, 24 and 12 orders the unweighted fit has
  `SE(B)` near 0.09. That gives a smallest distinguishable `kappa` of
  about 0.45 for the full residual. For the wave coordinate, whose square
  residual was about 0.22 of the full one, it is about 0.2.
- **What the heuristic predicts.** A count-like transport with the
  heuristic's unit constant has the variance of a null-coordinate
  difference equal to `d / N`. That is `S^m = 2 d / N`, or
  `kappa_m ≈ 1.4`, for the wave coordinate. It would add about `1.7 N` to
  `S_far^m N^2`: about 220 at N = 128 and 870 at N = 512, against a
  cell-scale level near 80.
- **Caveats.** These projections use central-event spreads; the interior
  spreads will be larger. The weighted fit would be sharper, since the
  N = 512 spread is smallest.

## 7. What a result would and would not mean for the routes

- **"Growing with d" with B above zero.** This is count-scale transport,
  seen as a lower bound. It would support the count-limited rows and
  Route D's `sqrt(d / N)` branch (heuristic 10/7). It would not prove
  them, since the size is a lower bound and the heuristic constant is
  unknown.
- **"No transport growth seen by this instrument" with a small `kappa`
  bound.** This removes count-scale transport above the bound at these
  lags, with the C3 fraction as the only check on the chain's reach. It
  does **not** trigger Route D's "smaller, so toward 4/3" branch: that
  branch needs a two-sided statement about the conditional law, and this
  instrument gives a one-sided one.
- Whatever is seen concerns positions given the order on the background.
  It says nothing about a wave, a statistic or the exponent, and decides
  nothing in the bracket `4/3 <= alpha_low <= alpha_high <= 2`.

## 8. Limits, stated before any data

- One-sided. The chain may not reach every configuration with the given
  order, and collective long-range modes are the slowest. C3 measures
  that only in 1+1.
- The group fit absorbs part of any smooth field, the largest scales
  first. C2 reports how much.
- Outer events include events near the tips and the boundary, whose cells
  are loose. `X` and the separate fit on `I_far` are there to keep this
  apart from transport.
- The chart's Euclidean distance is not a Lorentz-invariant separation.
  The bins are a diagnostic choice.
- The longest lag is `2 N^2`. In 1+1 the chain reported all of the
  all-event collective motion only by then, and the far-bin structure may
  take longer.
- One background, one diamond, three sizes within a factor 4, and 12
  orders at the largest.

## 9. What has been seen that bears on this

- All Route A and Route D records and their checks, including the
  outer-event residuals and the central pair-difference bins quoted in
  Section 2.
- The sizing numbers in Sections 4 to 6, computed from those records for
  this draft.
- The Route D 1+1 control records.

No transport readout (S, I, X over interior pairs), no 1+1 transport
structure and no 3+1 chain beyond the existing records has been computed
or seen. Nothing in Section 3 has been implemented.

## 10. To be confirmed by the user at approval

These follow the recommendations of the session that started this one and
have not been confirmed by the user:

1. **The positive control C3** and its pass rule (lower end of `F(N^2)` at
   least 0.5, and "growing with d" in 1+1).
2. **The flat wording** "no transport growth seen by this instrument",
   used whatever C3 returns.
3. **If any control fails, the main stage is not run**; the failure is
   reported and the user decides.
4. **Orders at N = 512:** 12 (about 3.9 hours) or 24 (about 7.8 hours).
5. **The cost:** about 6 to 7 hours of wall time at about 75% CPU, with
   12 orders at N = 512. No GPU.

## 11. Review: not to be frozen as written

[transport_referee_report.md](transport_referee_report.md). An AI referee
in a fresh context. It ran no 3+1 chain. It ran checks without a chain,
and the 1+1 chain at N = 32 with 16 orders. Its findings:

- **Both positive descriptions can appear with no transport.**
  - A cell-scale correlation over 0.2 chart units gives
    `Y = 1.17 to 1.22`: "growing with d" at every size.
  - The interior set drifts outward with N: the share of events beyond
    scale 0.9 goes from 3.5% to 8.7%. That alone gives a positive B on
    `S_far`, a spurious `kappa` of 0.56 to 1.79 depending on how loose the
    outer cells are.
  - What separates transport from both is how the near and mid-range
    covariance changes with N, read together with S.
- **C2's gate cannot fail when it should.** It is on total S, of which the
  field is 2 to 8% at `kappa = 0.3`.
- **C3 passes almost surely and says little.** "Growing with d" is vacuous
  in 1+1. The bins as written leave the far bin nearly empty. The
  premise that the order fixes the ranks is wrong but harmless: the chain
  stays in its starting realizer.
- **The sensitivity was quoted at break-even, not at 90% power.** With
  interior spreads it is no better than the existing central bound for the
  full residual. The heuristic's scale in the far bin is about `0.62 N`,
  not `1.7 N`, because U2's separation is along `u`.
- **The existing record already excludes the heuristic inside the central
  set at those lags.** A transport of the heuristic's size would leave a
  central wave-coordinate coefficient of about 0.29. The extension bounds
  it at 0.11 to 0.13, and could tell 0.16 from zero. The session that
  wrote this draft rechecked this without a chain and agrees.
- **What is left for a transport measurement:**
  - separations beyond the central set, where loose outer cells and
    composition confound;
  - lag `2 N^2` at N = 512.
- The fit itself is benign, and the cost table is consistent.

The referee lists six must-fix items (M1 to M6) and three should-fix items.
No revision has been written. Whether a revision is worth writing is the
user's decision. On the referee's projection, a revision would add little
sensitivity beyond the central bound for the full residual. It would add
mainly the larger separations, a near-null stratum for Route C, and the
longer lag at N = 512.
