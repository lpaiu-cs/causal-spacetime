# P1: Geometry-Dilution Epsilon Sweep

Status: DRAFT v1 (not frozen). Builds on the frozen PC-V1 positive control
(`docs/prereg/pc_v1_positive_control.md`), which validated the response-profile
pipeline as a latent-geometry discriminator (H-SENS and H-SPEC both supported).
This document lives in version control. Nothing under `outputs/` is ever an
input to a decision defined here.

## 1. Motivation

PC-V1 established the two endpoints: measured Minkowski-geometric order passes
the fixed gates and recovers true spatial order (epsilon = 0), while
density-matched geometry-free order blocks (epsilon = 1). PC-V1 cannot say
whether the transition is graded or a cliff, or where it sits. P1 makes
geometry the manipulated variable and measures the dose-response curve between
the endpoints. This is the experiment that tests the working hypothesis in its
proper conditional form: metric representation emerges to the extent that
causal order is geometrically consistent.

## 2. Hypotheses

- H-MONO (monotone degradation): as the geometry-dilution fraction epsilon
  increases from 0 to 1 at fixed relation density, the recovery of true
  spatial order degrades monotonically.
- H-THRESH (identifiable transition): there is an identifiable epsilon* at
  which recovery crosses from geometric to non-geometric, reported with its
  across-seed spread. P1 does not prereg a specific functional form (cliff vs
  smooth); it estimates and reports the crossing.
- H-LAG (embeddability outlasts geometry, secondary): pure embeddability
  (held-out violation) degrades at a larger epsilon than geometry recovery
  (truth-order error), i.e., there is a "false-pass" window where profiles
  still embed in 1D but no longer encode true space. Calibration showed this;
  P1 quantifies it as a secondary result.

Only degradation OF RECOVERY is interpretable. A monotone curve with an
identifiable transition supports the conditional hypothesis; it is not a
spacetime-emergence claim. A flat or non-monotone curve would falsify the
graded-emergence reading at this scale.

## 3. Design overview

Reuse the frozen PC-V1 instrument unchanged: scene construction, bracket-echo
measurement, parallax dissimilarity, leak-free split, the fit policy
(RepresentabilityFitPolicy, >= 1500 steps / 5 restarts / dims {1,2,3}), and
the no-silent-clamp assertion. The only new element is the epsilon-diluted
order generator (Section 4). The frozen PC-V1 gate thresholds are imported as
a fixed reference, not re-derived.

| Stage | Purpose | Seeds | Thresholds/tests |
| --- | --- | --- | --- |
| P1-A | calibration, exploratory | 0-9 | proposes the H-MONO/H-THRESH test constants |
| P1-B | confirmatory | 300-319 | frozen P1 constants + frozen PC-V1 gates |

Seed lists are exact and disjoint from each other and from all PC-V1 stages
(calibration reuses 0-9, which is acceptable for the exploratory stage;
confirmatory uses fresh seeds 300-319). Invalid (epsilon, seed) cells are
recorded, never silently replaced.

## 4. Epsilon-diluted order generator

`causal_spacetime_lab.positive_control.epsilon_sweep.build_epsilon_scene`.
For dilution fraction epsilon and a scene's geometric causal order:

- epsilon = 0: the pristine geometric order (identical to PC-V1 structured).
- epsilon > 0: take the geometric covering relation (transitive reduction),
  keep each non-chain covering edge with probability 1 - epsilon, add random
  time-respecting edges at a bisection-tuned probability, always keep
  chain-internal covering edges, then transitively close.
- Density hold: the random-edge probability is bisected per (epsilon, seed) so
  the achieved post-closure relation density among time-ordered pairs matches
  the geometric density (target tolerance 0.02). Achieved density is recorded
  for every cell; the sweep is only interpretable where density is held.
- Targets are re-selected under each diluted order by the PC-V1 policy
  (two-sided bracketed by all reference chains, band-limited, <= 40).

Rewiring covering edges (not closed edges) makes dilution effective: removing
a covering edge genuinely removes reachability, whereas closed-edge removal is
absorbed by redundant paths.

Epsilon grid (primary): {0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0} (8 points).
Endpoints epsilon = 0 and epsilon = 1 are anchors whose behavior is already
known from PC-V1 and must reproduce (epsilon = 0 passes gates and recovers x;
epsilon = 1 blocks).

## 5. Metrics (fixed here)

Per (epsilon, seed), at the gate dim d = 1 unless noted:

- achieved_density (must stay within 0.02 of geometric; cells outside are
  report-only, excluded from the monotonicity test).
- truth_order_error: sign-discordance between fitted pair distances and true
  |x_i - x_j| (PC-V1 metric m5). PRIMARY recovery metric.
- coord_x_rank_corr: |Spearman| between the fitted 1D coordinate and true x
  (reported; sign is free by reflection).
- heldout_violation, restart_order_disagreement: as in PC-V1 (for H-LAG and
  the frozen-gate readout).
- gate_pass: the frozen PC-V1 decision (G1 heldout <= 0.05 AND G4 truth <=
  0.15 AND G3 stability <= 0.15) at each cell. G2 (null gap) is optional and
  report-only in P1 (the column-shuffle null is epsilon-independent).

## 6. Preregistered tests and threshold-setting rule (P1-A -> frozen)

Hard floors (fixed now):

- HF1: at epsilon = 0 the median truth_order_error must be <= 0.15 and median
  heldout <= 0.05 (reproduce the PC-V1 endpoint). Else the sweep harness is
  inconsistent with PC-V1; stop.
- HF2: at epsilon = 1 the median truth_order_error must be >= 0.40 (reproduce
  the geometry-free endpoint). Else the dilution does not reach a geometry-free
  order; stop.

Mechanical test constants from P1-A distributions:

- H-MONO test: per seed, Spearman(epsilon, truth_order_error) over the grid
  (density-held cells only). H-MONO supported for a seed if rho_seed >=
  rho_min, where rho_min = min(0.85, round_down(p10_A(rho_seed), 0.05)).
- H-THRESH estimate: epsilon* = the linearly interpolated epsilon at which
  truth_order_error first crosses the midpoint level
  (truth@eps0 + truth@eps1) / 2. Report the across-seed median and IQR. No
  functional form is assumed.
- H-LAG estimate: epsilon_heldout* (crossing of heldout past the frozen 0.05
  gate) minus epsilon_truth* (crossing of truth past the frozen 0.15 gate).
  Report the median gap; H-LAG supported if the median gap > 0.

## 7. Decision rules (P1-B confirmatory)

- H-MONO supported: >= 16 of 20 valid seeds have rho_seed >= rho_min.
- H-THRESH supported: epsilon* is estimable (both endpoints reproduce HF1/HF2)
  for >= 16 of 20 seeds; report its distribution. (This is a reporting
  criterion, not a pass/fail of a physical claim.)
- H-LAG supported: median (epsilon_heldout* - epsilon_truth*) > 0 across valid
  seeds.
- Endpoint reproduction: epsilon = 0 must pass the frozen PC-V1 gate and
  epsilon = 1 must block, on >= 18 of 20 seeds, or the harness is
  inconsistent with PC-V1 and results are void.

## 8. Stop rules and repair policy

- Instrument repair (generator, grid, metrics) is allowed ONLY before the P1
  freeze, via a versioned amendment that re-runs P1-A on fresh seeds. The
  frozen PC-V1 instrument (scene/echo/dissimilarity/fit) must NOT be modified
  by P1; if P1 needs a change there, it is a new PC-V(n) with its own
  positive-control re-validation.
- No threshold retuning after the P1 freeze. At most 2 repair cycles.
- Never: dropping (epsilon, seed) cells after seeing outcomes, changing the
  epsilon grid post hoc, applying tests to report-only (non-density-held)
  cells.

## 9. Freeze and provenance

Freezing P1 mirrors PC-V1 Section 12: the maintainer commits P1-A summary and
the mechanical test constants under `docs/prereg/frozen/p1_test_constants.json`
with the calibration commit hash; the runner refuses P1-B without it. Every
output row carries scene digests, epsilon, achieved density, seed, stage, code
version, and the requested-vs-executed fit budget.

## 10. Claim boundary

- A monotone, density-held degradation curve with an identifiable epsilon*
  shows that the validated discriminator responds to graded geometric
  consistency, not to density or to a single cliff artifact. It licenses
  reading metric representation as a graded effective property of order
  consistency in this controlled family.
- It is NOT a proof of spacetime emergence, a continuum limit, or dynamics. It
  is a controlled dose-response inside a fixed generator family at fixed
  finite scale.
- The H-LAG "false-pass" window is a cautionary methodological result:
  embeddability alone (without the truth gate) would over-report geometry.

## 11. Deviations log

(empty — append dated entries; deviations do not retroactively change recorded
decisions)
