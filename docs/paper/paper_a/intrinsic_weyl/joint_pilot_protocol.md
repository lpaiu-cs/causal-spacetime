# Joint reconstruction: prospective development pilot v1

2026-09-22. This is a development pilot, **not** an accuracy gate or a
confirmatory result. The previous goal turn made progress: it established
local identifiability and measured an unfavorable noise cost. The present
step must test actual finite orders and remove the supplied landmark geometry.

## Estimator and input boundary

Fit inferred normal coordinates and all ten Weyl components jointly to
`S_ij = (x_i-x_j)^T eta (x_i-x_j) - C(x_i,x_j,x_i,x_j)/3`.
Use interval volume clocks `S_ij = -sqrt(n_ij/n_pq)`. The outer tips fix the
unit time axis at `(-1/2,0,0,0)` and `(1/2,0,0,0)`. The inferred tensor is
located at the midpoint of that timelike geodesic, relative to its observer.
Two inferred spatial vectors fix the remaining rotational gauge. No generator
coordinates or plane-wave tensor constraints enter the fit. Small, regular,
four-dimensional vacuum geometry remains an assumption.

Midpoint centering addresses the forward-origin expansion and the four-reference
noise bottleneck. It does not eliminate higher-order world-function terms or
the order-six volume-clock bias. Generalized least squares uses the overlapping
interval covariance measured from the SAME curved order, without a flat sample.

The blind `fit_order_archive` entrypoint accepts exactly five fields:
`selection_order`, `selected_ids`, `after`, `before`, `measurement_n`.
The last two relation arrays are packed order queries between selected
landmarks and independent measurement events. No coordinates, profiles, true
curvature, density calibration, or flat-control observations are allowed.
The entrypoint recomputes the selection and checks the saved IDs. Archives
are retained so this boundary can be independently inspected and replayed.

## Fixed development design, before stochastic data

- Use existing uniform physical-volume sampling and the guarded causal
  predicate in `window_response.py` and `p14_plane_wave.py`.
- Populations: constant vacuum plane waves with generator amplitudes
  `0`, `0.25`, `0.75`, each with two independent repetitions. The generator's
  proper tip separation is `sqrt(2)`, so the evaluation target's electric
  eigenvalues are `(-A,0,A)` in squared-tip-time units. Truth is available
  only to the evaluator. Flat populations are separate controls, not paired
  measurements and not calibration inputs to curved fits.
- Six seeds `2026092220..2026092225`, in amplitude-major, repetition-minor
  order. Independent RNG streams `[seed,0]` and `[seed,1]` supply the 1024
  selection events and 262144 measurement events. No seed removal or increase
  after results. Measurement chunk size 4096; total about 1.57 million events.
- Include the two outer tips in the selection order. Randomize interior
  labels before selection so the generator's `u` sorting is not an input.
- Select at most 32 landmarks in each of two time-proxy strata. The proxies
  are computed solely from endpoint interval counts:
  `t=(sqrt(n_pz/n_pq)-sqrt(n_zq/n_pq))/2` and
  `r^2=1/4+t^2-(sqrt(n_pz/n_pq)+sqrt(n_zq/n_pq))/2`.
  Require `0.15 <= |t| <= 0.32`, `r^2 >= 0.01`, positive endpoint clocks
  with at least two selection events, and score candidates by their number
  of clocks containing at least two selection events times `min(r^2,0.04)`.
  Peel remaining candidates with degree below eight. If fewer than sixteen
  remain, record failure. Exact score ties use input order.
- Choose all pair clocks with at least two SELECTION events, then measure
  them with the independent stream. Zero measurement counts cause failure;
  pair selection is not revised using measurement counts.
- Compute the full empirical covariance of interval membership, propagating
  it through the square root. Standardize by its diagonal and retain
  covariance eigenvalues greater than `1e-8` times the maximum. Freeze this
  cutoff. Use four deterministic optimizer initializations (seeds 0..3),
  max 600 function evaluations, and tolerances `1e-10`. Choose the flat
  initialization by its objective, then fit the ten curvature components.
- Report every eigenvalue estimate, coordinate and curvature rank, optimizer
  status, boundary violations, nominal component covariance, run time, and
  failure. A successful optimizer is NOT a scientific success criterion.
  The covariance is a local model approximation, not validated 95% coverage.
- Commit code/protocol before data; run from a clean worktree; bind each
  archive and result to hashes and a reachable commit. Preserve earlier
  diagnostics. Never replace this pilot by a post-result tuned version.

## Required checks and interpretation

Before data, check analytic derivatives by finite differences, recovery of
zero curvature from exact flat clocks with no coordinate input, and agreement
of rectangular order queries with the existing full causal predicate.

This pilot can establish that an actual order-only pipeline executes and can
reveal failures or encouraging estimates. It cannot by itself prove accuracy,
coverage, arbitrary-profile generalization, nonvacuum validity, or a global
curvature theorem. A later accuracy claim must freeze the actual estimand,
95% two-sided construction, bias treatment, and sample size powered to 90%
before independent confirmation data. The full active goal remains open.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_joint_pilot.py --check
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_joint_pilot.py
```
