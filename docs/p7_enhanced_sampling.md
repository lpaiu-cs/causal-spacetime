# P7 enhanced sampling: Wang-Landau tunneling measurement and verdict

Status: **METHOD DEVELOPMENT — NOTHING FROZEN, NO PHASE CLAIM** (2026-07-15).

## Question

The frozen P7 N=600 characterization ended in start-dependent hysteresis:
for beta >= 16 the random and bipartite starts sit in different basins with
mean-action gaps of 77-96, local Metropolis never crosses them, and no
equilibrium beta_c is extractable at any budget. The prereg therefore
requires a *validated* enhanced-sampling method before N=900/1200. This
document records the validation attempt for Wang-Landau / multicanonical
sampling, and its outcome.

## Method

`causal_spacetime_lab.positive_control.multicanonical` implements
Wang-Landau estimation of the density of states ln g(S) with a
Belardinelli-Pereyra 1/t tail, followed by fixed-weight multicanonical
sampling and canonical reweighting. Correctness is pinned to exact
enumeration at N=6 (all 720 permutations): the estimated ln g(S) and the
reweighted canonical means of S, n0 and height match the analytic Gibbs
answer at beta = 0, 0.4, 1.0, 2.0 (`tests/test_multicanonical.py`).

The walker demonstrably crosses the action barrier that local Metropolis
cannot: at N=60, eps*N=12, a single run completed 30 window traversals.
The method *works*. The question is what it costs as N grows.

The cost structure is pinned by the schedule itself: ln_f halves once per
completed round trip, so convergence to ln_f = 1e-5 costs ~17 round trips,
and the only free parameter is tau(N) = moves per round trip.

## Measurement

tau was measured at eps*N = 12 (the P7 design axis) across four run
batches (6, 25, 15, and 15 round-trip targets). The complete raw record —
every batch log and every row — is committed under
`docs/p7_tunneling_data/` (`runs_combined.csv` plus the four batch logs),
and `experiments/positive_control/p7_tau_fit.py` reproduces every number
in this section from that CSV deterministically.

All batches derived their RNG seed from the same formula
(`base + 1000*rep + n`), so a (N, rep) pair names ONE random-walk
trajectory: batches that re-ran a pair observed the *same walk* to a
different stopping point. Each walk therefore contributes exactly one
clean measurement (the one with the most completed round trips). That
leaves **16 independent clean walks** and 2 censored ones:

| N  | walks | min        | max        | mean       |
|----|-------|------------|------------|------------|
| 30 | 3     | 28,469     | 40,516     | 35,063     |
| 40 | 3     | 122,933    | 170,656    | 147,569    |
| 50 | 3     | 224,231    | 580,136    | 412,030    |
| 60 | 3     | 1,352,109  | 4,465,260  | 2,588,175  |
| 70 | 3     | 2,305,235  | 6,011,521  | 4,639,258  |
| 80 | 1     | 10,519,001 | 10,519,001 | 10,519,001 |

Censored (excluded from fits): both remaining N=80 walks stopped at the
200M-move cap with 13 and 14 of 15 traversals; the valid run-level bound
is moves/target, i.e. **tau > 13.17M for each** (moves/completed would
overstate the bound, since the budget includes the unfinished traversal).
The N=50 rep=2 walk also carries a censored 23-round-trip observation
(78.4M moves at the 25-trip target) alongside its clean 15-trip
measurement — evidence of how heavy the within-walk tail is. At N=80 the
one clean value is the *fastest* of three walks; the true tau(80) is
plausibly >= 13M.

Within-size spread is 1.4x-3.3x across walks: round-trip times are
heavy-tailed, which is why an early fit over single seeds per size
produced a spurious non-monotonic tau and a meaningless exponent of 2.95.

## Fit

Power-law fit tau ~ N^alpha over the 16 independent clean walks (every
walk one point, so within-size scatter propagates), 10,000-sample
bootstrap, intercept refit per alpha (slope and intercept are strongly
anticorrelated over this lever arm, so carrying the point-fit intercept
to the interval endpoints misstates endpoint costs by factors):

```
alpha = 5.88,   90% interval [5.44, 6.38]
```

Extrapolated cost of one converged ln_g at N=600 (17 round trips):

| alpha        | moves    | Numba kernel (135 us/move) |
|--------------|----------|----------------------------|
| 5.44 (low)   | 7.4e12   | ~278,000 core-hours        |
| 5.88 (point) | 2.2e13   | ~829,000 core-hours        |
| 6.38 (high)  | 7.8e13   | ~2,914,000 core-hours      |

Even granting a hypothetical 100x GPU ensemble speedup, the low end is
~2,800 GPU-hours for a *single* converged ln_g — before any production
sampling, seeds, or the N=900/1200 grid the prereg actually asks for.

Local slopes between adjacent sizes wobble between ~4.9 and ~9 with no
clean trend, so these data do not distinguish a steep power law from
exponential growth. The verdict does not depend on the distinction: both
put N=600 out of reach, and if the truth is exponential every number above
is an underestimate.

## Corrections after review (2026-07-15)

PR review of the first version of this document forced four corrections,
recorded here rather than silently absorbed:

1. **Duplicate-walk double counting (found while addressing the review's
   reproducibility demand).** The original headline fit pooled "25 clean
   runs", but 9 of them were re-measurements of the same trajectories
   (identical seeds across batches, different stopping points), so the
   bootstrap treated correlated points as independent. Corrected fit over
   16 independent walks: alpha 5.79 -> 5.88, interval [5.36, 6.19] ->
   [5.44, 6.38]. The verdict is unchanged and slightly strengthened.
2. **Censored bounds.** Previously reported as moves/completed (14.1M and
   15.2M); the valid run-level bound is moves/target = 13.17M for both.
3. **Interval-endpoint costs** now refit the intercept per alpha instead
   of reusing the point-fit intercept.
4. **Reproducibility.** The fit input now ships in the repository
   (`docs/p7_tunneling_data/`) with a deterministic script
   (`p7_tau_fit.py`); the first version committed only the final batch's
   12 rows while quoting a fit over 28.

Known limitation: window edge-touch flags were not recorded for these
batches (the sweep script records them from now on). A too-narrow window
truncates traversals and *underestimates* tau at that size, so the risk
runs toward understating costs — in the elimination's favour, not
against it. The tau data were also collected under the pre-normalization
1/t tail (the Belardinelli-Pereyra time is trial moves per bin; fixed in
the sampler since); the schedule's tail value does not enter the
round-trip count, and the 17-halvings cost model is unchanged, with the
normalized tail adding only a negligible ~n_bins/ln_f_final moves floor.

## Verdict

**Wang-Landau / multicanonical sampling is eliminated as the P7
enhanced-sampling method at N=600.** The 90% interval on the tunneling
exponent excludes every feasible-cost scenario, and the censoring bias
runs in the method's favour.

Standing eliminations for the N=600 barrier, all now quantified:

1. Local Metropolis — frozen P7 result (hysteresis, no crossing in 3M
   steps).
2. Naive replica exchange — beta spacing ~0.02 needed at an action gap of
   77 implies hundreds of replicas with no crossing guarantee.
3. Wang-Landau multicanonical — tunneling exponent 5.9 [5.4, 6.4] measured
   here.

## What survives

WL converges and tunnels comfortably at N <= 80 on a desktop CPU, and the
measured exponent projects single-seed costs at N=100 of ~40M moves/RT
(hours with a Numba kernel) and N=120 of ~110M moves/RT (a day). A
finite-size-scaling study over N in {60, 80, 100, 120} at eps*N = 12 is
therefore *achievable* with the code on this branch, if the program
chooses to re-scope P7's FSS axis to sizes where equilibrium sampling is
demonstrably possible. That is a scientific scope decision, not a
methods one, and is left to a future prereg freeze.

## Provenance

- Measurement scripts: `experiments/positive_control/p7_tunneling_scaling.py`,
  `experiments/positive_control/p7_wl_pilot.py`.
- Fit input and reproduction: `docs/p7_tunneling_data/runs_combined.csv`
  (all four batches, every row) with the raw batch logs alongside;
  `experiments/positive_control/p7_tau_fit.py` reproduces the exponent,
  interval, and cost table deterministically.
- Final batch: 12-run parallel sweep, 2026-07-15, 8.0 h wall on 12
  workers (Ryzen 9800X3D); earlier sequential runs 2026-07-15.
- Nothing in this document is preregistered. No beta_c, G(beta), or phase
  claim is made or implied.
