# P6a: Near-critical N=600 characterization

Status: **PROTOCOL FROZEN, no outcome prediction** (2026-07-14). Constants are
in `docs/prereg/frozen/p6_near_critical_constants.json`.

## Question

What phases and reconstruction verdicts occur at beta in {12, 16, 24}, inside
the P5 reconnaissance bracket `(8, 32)`? This arm refines the transition at
N=600 and supplies the first P7 finite-size-scaling grid. It is explicitly a
characterization arm: a continuum state may legitimately pass, and a layered
or crystal state may block structurally or at the numerical gates.

## Frozen chains

The ensemble is the P5 2D-order ensemble at N=600 and eps=0.02. Each beta has
three 3M-step chains: two seeded random starts and one exact bipartite start.
Burn-in is 60%; nominally 48 post-burn samples are retained every 25,000
steps. Matching the validated sampler, a scheduled point with a no-op proposal
is omitted; at least 45 samples are required. The first, middle, and last
retained states from each chain are judged by the unchanged P3/P5 instrument.
All nine chains use the trajectory-equivalent accelerated replay, with exact
RNG proposal streams and periodic exact action resynchronization.

## Recorded quantities

Every retained sample records action, n0/n1/n2, Myrheim-Meyer dimension, and
height. Per chain, report acceptance, integrated autocorrelation time and ESS
for action/n0/height, and a phase label using the already used P5 reference
rules. The three instrument snapshots record structural status, heldout,
null-gap, truth error, and frozen gate pass.

The IAT estimator uses Geyer's initial-positive-pair sequence: consecutive
autocorrelations are summed as `rho[2k] + rho[2k+1]`, and truncation occurs at
the first non-positive pair. With only 48 retained samples it is a screening
diagnostic, not a precision ESS claim. Dual-start disagreement is reported
directly and is not averaged away.

## Decision semantics

There is no confirmatory PASS/FAIL for this arm and no post hoc gate. The
scientific output is the observed phase/reconstruction map plus chain-quality
caveats. These data may motivate a later P7 freeze but do not themselves fix
P7's geometry-order-parameter threshold or transition hypothesis.

## Deviations log

- On 2026-07-14, PR review identified that the archived aggregation had
  truncated IAT at the first non-positive single lag rather than the first
  non-positive paired sum. The estimator was corrected to Geyer's
  initial-positive-pair sequence and summaries were regenerated from the
  unchanged retained chains. No MCMC chain or instrument snapshot was rerun.
  Phase labels and reconstruction results were unchanged; the downstream P7
  beta=14 minimum ESS changed enough to alter its screening verdict, as
  recorded in `p7_geometry_fss.md`.

## Outcome

Executed 2026-07-14. This characterization does not produce a confirmatory
PASS/FAIL. All nine chains met the frozen sample-count requirement and all
scheduled outputs are archived under `docs/prereg/frozen/`.

- At beta = 12, all three starts were classified continuum, all 9/9
  instrument snapshots passed, and the smallest reported ESS was 31.8. This
  is the only tested point with clean dual-start agreement and no material
  mixing warning.
- At beta = 16, both random starts remained continuum and passed all 6/6
  snapshots. The bipartite start ended with a continuum mean label, but only
  its last snapshot became structurally measurable (1/3 pass); its action,
  n0, and height ESS were only 2.69, 3.17, and 2.94. The apparent phase
  agreement is therefore not evidence of equilibration.
- At beta = 24, both random starts remained continuum and passed 6/6
  snapshots, while the bipartite start remained intermediate and was
  structurally blocked in 3/3 snapshots. Acceptance was 0.0116. This is
  direct finite-run hysteresis / basin dependence, not an equilibrium
  estimate of beta_c.

The N = 600 refinement therefore establishes a usable lower anchor at
beta = 12 and detects severe metastability by beta = 16, with explicit
dual-start disagreement at beta = 24. It does **not** locate an equilibrium
critical beta inside `(12, 24)`. P7 must retain dual starts, treat convergence
as a prerequisite for transition inference, and use replica exchange or an
equivalent mixing intervention before making a thermodynamic coincidence
claim.
