# P2: 2+1D robustness and dimension selection

Status: DRAFT v1 (not frozen). Builds on the frozen PC-V1 discriminator and P1
dose-response, extending them to 2+1D. This document lives in version control;
nothing under `outputs/` is ever a decision input.

## 1. Motivation

PC-V1 and P1 are 1+1D: there is one spatial dimension, so the effective
embedding dimension is trivially 1 and dimension selection is untested. In 2+1D
there are two spatial dimensions, so a pipeline that genuinely recovers spatial
structure must (a) still pass and recover true position on geometric order, (b)
require dimension 2 — one embedding dimension must be insufficient — and (c)
still block on geometry-free order. P2 tests these. The novel test is H-DIM:
the effective dimension is 2, not 1.

## 2. Hypotheses

- H-SENS-2D: on 2+1D Minkowski-sprinkled causal sets, measured bracket-width
  parallax profiles pass the gate at embedding dimension d = 2 and recover the
  true 2D spatial order of targets.
- H-DIM: the effective spatial dimension is 2. Dimension 1 is insufficient —
  its truth-order error exceeds the truth gate — while dimension 2 is within
  it; dimension 3 does not improve truth recovery beyond tolerance (2 already
  saturates recovery).
- H-SPEC-2D: density-matched geometry-free order blocks at d = 2 (held-out
  above the gate, or structural failure to form a constraint pool).

Only the pass-at-2 / fail-at-1 / block-on-geometry-free pattern is
interpretable. It is not a spacetime-emergence claim.

## 3. Design

Reuse the frozen PC-V1 measurement/dissimilarity/fit pipeline unchanged; the
only new element is 2+1D scene construction (Section 4). Fit at d = 1, 2, 3;
the gate dimension is 2. P2-A calibration (seeds 0-9) proposes thresholds
mechanically; P2-B confirmatory (fresh seeds 400-419) applies the frozen
thresholds. Invalid seeds are recorded, never silently replaced.

## 4. 2+1D scene

`causal_spacetime_lab.positive_control.scene_2d.build_scene_2plus1d`:

- Bulk: `sprinkle_minkowski_causal_diamond(n = 2600, spacetime_dim = 3,
  T = 2.0)`.
- References: K = 8 stationary observer chains on a ring of radius 0.25 in the
  (x, y) plane (non-collinear, so cross-observer parallax fixes a 2D position),
  96 ticks each over t in [-0.7, 0.7], all inside the diamond.
- Causal order: `causal_matrix_minkowski` (dt > 0 and dt^2 >= ||dx||^2).
- Targets: bulk events with |t| <= 0.10 and radius sqrt(x^2 + y^2) <= 0.22 that
  are two-sided bracketed by all 8 chains; >= 30 required, subsampled to <= 44.
- Provenance: sha256 digests of events and causal matrix on every row.

Bracket-width measurement, parallax dissimilarity, margin, leak-free
pair-level split, and the fit policy (>= 1500 steps, 5 restarts) are the frozen
PC-V1 pieces, all dimension-agnostic. Truth recovery is scored against the true
2D coordinates (x, y).

## 5. Metrics

Per seed, at each dim d in {1, 2, 3}: held-out violation and truth-order error
(sign discordance of fitted pair distances vs true 2D pair distances). Plus the
geometry-free control's d = 2 held-out (or a structural-block flag).

## 6. Threshold-setting rule (P2-A -> frozen)

Hard floors (fixed now):

- HF1: Stage A median d = 2 truth-order error <= 0.20 (2D recovery works) and
  median d = 2 held-out <= 0.10.
- HF2: Stage A minimum d = 1 truth-order error > the proposed truth gate
  (1D genuinely underfits, so H-DIM is testable).

Mechanical thresholds from Stage A distributions (amended, deviation D1 — gates
are placed at the MIDPOINT between the pass cluster and the fail cluster so both
sides keep margin, since resolution does not move the d=2-truth floor):

- gate_truth = min(0.20, round(0.5 * (max_A(d2 truth) + min_A(d1 truth)), 0.05))
- gate_heldout = min(0.10, round(0.5 * (max_A(d2 held-out) +
  min_A(control d2 held-out)), 0.05))

## 7. Decision rules (P2-B confirmatory)

- H-SENS-2D supported: >= 16 of 20 valid seeds have d = 2 held-out <=
  gate_heldout AND d = 2 truth <= gate_truth.
- H-DIM supported: >= 16 of 20 valid seeds have d = 1 truth > gate_truth
  (1D underfits) AND d = 2 truth <= gate_truth (2D suffices); additionally
  report the fraction with d = 2 truth <= d = 3 truth + 0.05 (recovery
  saturates at 2).
- H-SPEC-2D supported: >= 16 of 20 valid seeds block (control d = 2 held-out >
  gate_heldout, or structural pool failure).

## 8. Stop rules and repair policy

Instrument repair (scene parameters, gate construction) is allowed only before
the P2 freeze, via a versioned amendment re-running P2-A on fresh seeds. The
frozen PC-V1 measurement/fit pipeline must not be modified by P2. No threshold
retuning after freeze; at most 2 repair cycles; no dropping seeds after seeing
outcomes.

## 9. Freeze and provenance

Freezing mirrors PC-V1 Section 12: commit the P2-A summary and mechanical
thresholds to `docs/prereg/frozen/p2_test_constants.json` with the calibration
commit hash; P2-B refuses to run without it. Every row carries scene digests,
seed, stage, and code version.

## 10. Claim boundary

Pass-at-2 / fail-at-1 / block-on-geometry-free demonstrates that the validated
discriminator recovers the correct spatial dimension and 2D structure of a
2+1D causal order, robust beyond the 1+1D case. It is not a spacetime-emergence
claim, a continuum-limit claim, or a dynamics claim; the geometry is sprinkled
in and recovered.

## 11. Deviations log

- D1 (2026-07-09, midpoint gate placement). The first P2-A calibration
  (seeds 0-9, code b941be0) passed all three hypotheses and both hard floors,
  but the original truth gate (p90 of d=2 truth, rounded) landed at 0.10 and
  hugged the d=2 cluster (d=2 truth 0.089-0.117), leaving one seed above it and
  no margin. A resolution sweep (chains 8->12, N 2600->3200, ticks 96->128)
  did NOT move the d=2-truth floor (mean ~0.097 throughout): the ~0.10 floor is
  a near-degenerate-pair effect in the compact 2D target band, not a resolution
  limit. Margin is therefore secured by gate PLACEMENT, not resolution: both
  gates are placed at the midpoint between the pass cluster and the fail cluster
  (Section 6). This is a pre-freeze gate-construction refinement; the frozen
  PC-V1 measurement/fit pipeline is untouched and no scene parameter changed.
  Thresholds derive only from the post-refinement re-calibration.
