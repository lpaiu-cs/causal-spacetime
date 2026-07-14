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

The IAT estimator sums positive sample autocorrelations until the first
non-positive lag. With only 48 retained samples it is a screening diagnostic,
not a precision ESS claim. Dual-start disagreement is reported directly and
is not averaged away.

## Decision semantics

There is no confirmatory PASS/FAIL for this arm and no post hoc gate. The
scientific output is the observed phase/reconstruction map plus chain-quality
caveats. These data may motivate a later P7 freeze but do not themselves fix
P7's geometry-order-parameter threshold or transition hypothesis.

## Deviations log

(empty)
