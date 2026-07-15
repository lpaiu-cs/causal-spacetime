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

tau was measured at eps*N = 12 (the P7 design axis) with 15-25 round trips
per run and 3-6 independent seeds per size. Runs that exhausted their move
budget before the target report a lower bound and are excluded from fits
(excluding them removes only *high* values at large N, so every fit below
is biased toward optimism).

Clean measurements (moves per round trip):

| N  | runs | min        | max        | mean       |
|----|------|------------|------------|------------|
| 30 | 6    | 26,306     | 40,516     | 33,861     |
| 40 | 6    | 97,238     | 170,656    | 137,567    |
| 50 | 5    | 209,585    | 666,388    | 422,412    |
| 60 | 4    | 1,171,068  | 4,465,260  | 2,233,899  |
| 70 | 3    | 2,305,235  | 6,011,521  | 4,639,258  |
| 80 | 1    | 10,519,001 | 10,519,001 | 10,519,001 |

Censored (lower bounds, excluded from fits): N=50 one run >= 3.4M;
N=80 two runs >= 14.1M and >= 15.2M (13-14 of 15 round trips completed at
the 200M-move cap). Note that at N=80 the one clean value is the *fastest*
of three seeds; the true tau(80) is plausibly ~13-15M.

Within-size spread is 1.4x-3.8x across seeds: round-trip times are
heavy-tailed, which is why early fits over single seeds per size produced
a spurious non-monotonic tau and a meaningless exponent of 2.95.

## Fit

Power-law fit tau ~ N^alpha over all 25 clean runs (every run a point, so
within-size scatter propagates), 10,000-sample bootstrap:

```
alpha = 5.79,   90% interval [5.36, 6.19]
```

Extrapolated cost of one converged ln_g at N=600 (17 round trips,
intercept refit per alpha):

| alpha        | moves    | Numba kernel (135 us/move) |
|--------------|----------|----------------------------|
| 5.36 (low)   | 5.5e12   | ~205,000 core-hours        |
| 5.79 (point) | 1.6e13   | ~613,000 core-hours        |
| 6.19 (high)  | 4.6e13   | ~1,726,000 core-hours      |

Even granting a hypothetical 100x GPU ensemble speedup, the low end is
~2,000 GPU-hours and the point estimate ~6,000 GPU-hours for a *single*
converged ln_g — before any production sampling, seeds, or the N=900/1200
grid the prereg actually asks for.

Local slopes between adjacent sizes wobble between 4.7 and 9.1 with no
clean trend, so these data do not distinguish a steep power law from
exponential growth. The verdict does not depend on the distinction: both
put N=600 out of reach, and if the truth is exponential every number above
is an underestimate.

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
3. Wang-Landau multicanonical — tunneling exponent 5.8 [5.4, 6.2] measured
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
- Raw sweep log: 12-run parallel sweep, 2026-07-15, 8.0 h wall on 12
  workers (Ryzen 9800X3D); earlier sequential runs 2026-07-15.
- Nothing in this document is preregistered. No beta_c, G(beta), or phase
  claim is made or implied.
