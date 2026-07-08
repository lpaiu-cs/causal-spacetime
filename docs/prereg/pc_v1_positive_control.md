# PC-V1: Positive-Control Preregistration For Response-Profile Representability

Status: FROZEN v1 (2026-07-08). Calibration provenance: Stage A at commit
9162e8e, seeds 0-9. Frozen thresholds: `docs/prereg/frozen/pc_v1_thresholds.json`.
Instrument repairs D1/D2 (Section 14) were applied before this freeze; the
frozen thresholds derive only from the post-repair calibration. After this
point, thresholds are locked (Sections 10-12) and Stage B/C are confirmatory.
This document lives in version control. Nothing under `outputs/` is ever an
input to a decision defined here.

## 1. Background and motivation

The 2026-07 program review found that the earlier confirmatory pipeline
(M29-M60) never demonstrated that its carry-forward gates are passable by any
input: production profiles were synthetic rank tables rather than measured
order data, no positive control existed, the fit budget was silently reduced,
and one blocking metric was saturated by construction. This experiment
restarts the representability program from the verified M17 foundation and
asks the missing question first.

The pipeline under test is: sprinkled 1+1D Minkowski causal set -> supplied
observer tick chains -> order-level bracket-width echo profiles -> pairwise
response-profile dissimilarity -> ordinal-embedding representability
diagnostics -> fixed gates.

## 2. Hypotheses

- H-SENS (sensitivity): on causal sets sprinkled from 1+1D Minkowski geometry,
  measured response profiles admit a stable low-dimensional ordinal
  representation that passes the fixed gates, and the fitted d=1 coordinates
  agree with the true spatial order of targets.
- H-SPEC (specificity): on matched geometry-free order data (time-respecting
  random partial orders) and on column-shuffled profiles, the same pipeline
  with the same frozen gates blocks.

A pipeline that passes H-SENS and H-SPEC is a validated discriminator for
latent spatial structure in discrete causal data. Passing does not prove
spacetime emergence; blocking on geometry-free data is the expected null
outcome, not a physics discovery. Only the pass/block CONTRAST is
interpretable.

## 3. Design overview

Three stages. Stage A calibrates the instrument and proposes thresholds.
Stage B is the confirmatory sensitivity run on fresh seeds with frozen
thresholds. Stage C is the confirmatory specificity run.

| Stage | Purpose | Seeds | Data | Thresholds |
| --- | --- | --- | --- | --- |
| A | calibration, exploratory | 0-9 | geometric + all controls | none (proposes) |
| B | confirmatory sensitivity | 100-119 | geometric | frozen from A |
| C | confirmatory specificity | 200-209 | geometry-free + shuffled | frozen from A |

Seed lists are exact and disjoint. A seed that fails scene-validity
preconditions (Section 4) is recorded as invalid and reported; it is never
silently replaced.

## 4. Data generation (primary configuration)

All quantities in natural units, c = 1. Generator:
`causal_spacetime_lab.positive_control.scene.build_positive_control_scene`.

- Bulk: `sprinkle_1p1_causal_diamond(n=900, T=2.0, seed=s)` (verified module).
- Reference protocols (supplied observer structure, per axiom A4): K = 6
  stationary tick chains at `x0 in linspace(-0.25, 0.25, 6)`, each with 96
  tick events uniform in `t in [-0.7, 0.7]`, inserted into the event set.
  Clock labels used downstream are tick INDICES (0..95); clock values are
  never used. All chain events must lie inside the diamond (hard check).
- Causal order: `causal_matrix_1p1` on the combined event set (J+ convention,
  atol 1e-12, as verified).
- Targets: bulk events with `|t| <= 0.10` and `|x| <= 0.25` that are two-sided
  bracketed (predecessor and successor tick) on all 6 chains; deterministic
  subsample to at most 40 targets. Scene validity requires >= 30 targets.
- Provenance: every scene records sha256 digests of the event array and the
  causal matrix plus the full config; these appear in every output row.

Rationale for the region: with this geometry every chain brackets every
target (worst case tau_return = 0.60 < 0.70), so the primary configuration
isolates representability from missing-data effects. Tick resolution
(delta_tau/2 ~ 0.0074) is below the typical nearest-target spacing (~0.0125).

Secondary configuration (report-only in PC-V1): two outer chains replaced by
boosted chains (beta = +/-0.3), and a wider target band that induces missing
brackets. No gate is applied to the secondary configuration.

## 5. Measurement protocol (order-level echo)

For target j and chain r, with tick indices as clocks:

- `p = max{ i : tick_i precedes j }`, `s = min{ i : j precedes tick_i }`
  (strict causal order; `find_radar_ticks_from_order`, verified M5 module).
- Bracket-width echo delay `W[j, r] = s - p` (integer tick units). Missing if
  either side is absent.

`W` is the passive same-protocol echo observable (legacy W_rank semantics),
now measured from geometric causal order instead of synthesized. Continuum
expectation: `W[j, r] ~ 2 |x_j - x0_r| / delta_tau + O(1)`.

A single measurement protocol is used for all columns of one profile matrix;
columns differ only by reference position (protocol-invariance rule from the
M40 audit, kept).

## 6. Dissimilarity, split, constraints

- Parallax profile (amended, see deviation D1): each target's bracket widths
  are centered across the references that reach it,
  `P[j, r] = W[j, r] - mean_r' W[j, r']`. This removes the shared per-target
  common mode (a global/temporal scalar that any time-respecting order
  injects) and keeps only cross-observer parallax, the part that carries
  spatial information. It operationalizes the underdetermination principle
  that a single shared scalar across observers is not a distance structure.
- Profile dissimilarity: `D[i, j] = RMS over common reachable columns of
  (P[i, .] - P[j, .])`, defined only if >= 4 common columns (primary: all 6).
- Constraint margin: `delta = 25th percentile of positive |D(i,j) - D(k,l)|`
  over 20,000 deterministic probe quadruples. The margin is derived from the
  measured resolution, not hand-tuned.
- Pair-level held-out split: target pairs are assigned train/held-out with
  probability 0.8/0.2 by a seeded deterministic draw. A train constraint uses
  two train pairs; a held-out constraint uses two held-out pairs; mixed
  quadruples are discarded. This prevents the column-split leakage of the
  legacy pipeline.
- Counts: 4000 train / 1000 held-out constraints of the form
  `D(i,j) + delta < D(k,l)`.

## 7. Fit policy (instrument integrity)

- Fitter: `fit_ordinal_embedding_gradient_descent` (verified M12 module),
  hinge margin default, learning rate 0.05, `batch_size = None` (full batch).
- Budget: dims {1, 2, 3} x 1500 steps x 5 restarts (best-of). The runner
  passes these values explicitly and asserts the diagnostics echo them; any
  silent clamp is an instrument defect and stops the run.
- Budget floor justification (measured during skeleton validation,
  2026-07-08, smoke configuration): with identical near-perfect input signal
  (Spearman rho ~0.99 between profile dissimilarity and true |dx|), a
  200-step/2-restart fit produced heldout violation 0.33 (reads as null)
  while 800/3 produced 0.04 and 1500/5 produced 0.05 with truth-order error
  0.09. An underpowered optimizer is indistinguishable from absent signal;
  the production budget must never be reduced below 800 steps / 3 restarts.
- Restart stability: 4 additional independent single-restart fits; stability
  is `pairwise_order_stability` (sign-discordance over sampled pair-pairs).
- FORBIDDEN METRIC: positional argsort mismatch
  (`mean(argsort(d_a) != argsort(d_b))`) is banned as a stability or
  agreement metric. It saturates at ~1 for dense continuous distances under
  jitter and produced the unpassable-by-construction criterion in M52-M59.
- Convergence is operationalized as restart-order stability (below), not as a
  loss-trajectory heuristic.

## 8. Metrics (all formulas fixed here)

Per scene and embedding dimension d:

- m1 `train_violation`, m2 `heldout_violation`: fraction of (train/held-out)
  constraints violated by the fitted coordinates.
- m3 `generalization_gap = m2 - m1`.
- m4 `restart_order_disagreement`: mean pairwise sign-discordance across the
  4 stability fits (5000 pair-pair comparisons, seeded).
- m5 `truth_distance_order_error` (geometric scenes only): sign-discordance
  between fitted pair distances and true `|x_i - x_j|` pair distances
  (`embedding_distance_order_error` against the true x coordinates).
- m6 `null_gap = heldout_violation(column-shuffled) - heldout_violation
  (structured)` at the same dimension, same seed, same budget.
- m7 `dimension_preference`: heldout_violation(d=1) <= heldout_violation(d=2)
  + 0.02 (1+1D expectation; report-only in PC-V1).
- m8 `symmetry_shift`: |heldout_violation(relabeled) - heldout_violation
  (structured)| under a consistent target relabeling (label-equivariance
  sanity; report-only, expected ~0).

Gate metrics are evaluated at d = 1 (declared physical expectation for 1+1D);
all dims are reported.

## 9. Controls

- C-NULL (destructive, marginal-preserving): each profile column permuted
  independently across targets; destroys cross-column consistency, preserves
  per-column marginals. Used for m6.
- C-FLIP (constraint-side flip, report-only): each constraint's sides swapped
  with probability 0.5; heldout_violation must move to ~0.5.
- C-SYM (symmetry, non-destructive): consistent relabeling; metrics must be
  unchanged (m8).
- C-RANDOM-ORDER (Stage C): time-respecting random partial order (transitive
  percolation) whose pre-closure edge probability is tuned by bisection so
  the achieved post-closure relation density among time-ordered pairs matches
  the geometric scene's causal density (amended, see deviation D2);
  chain-internal edges preserved; transitively closed. A single shared
  uniform draw makes density monotone in the edge probability, so the
  bisection is exact and deterministic. Achieved density and the tuned edge
  probability are recorded. Target selection policy is re-applied under the
  random order. This is the epsilon = 1 endpoint; the graded epsilon sweep is
  deferred to the next preregistration (P1).

## 10. Gates and threshold-setting rule

Hard floors are fixed now, before any run:

- F1: if Stage A median heldout_violation(structured, d=1) > 0.35, the
  instrument is declared inadequate. Stop; no threshold may be created above
  0.35.
- F2: if Stage A effect size (Cohen's d of heldout_violation between
  structured and C-NULL across seeds) < 2.0, apply the single allowed
  escalation (Section 11) or stop.

Threshold construction (mechanical, from Stage A distributions):

- G1 heldout_violation <= min(0.35, round_up(p90_A(structured), 0.05))
- G2 null_gap >= max(0.10, round_down(p10_A(null_gap), 0.05))
- G3 restart_order_disagreement <= min(0.30, round_up(p90_A, 0.05))
- G4 truth_distance_order_error <= round_up(p90_A, 0.05)  (sensitivity only)
- G5 scene validity + no-clamp assertion + constraint counts as specified.

Decision rules (per configuration):

- Stage B PASS (H-SENS supported): >= 16 of 20 valid seeds satisfy G1-G5.
- Stage C BLOCK (H-SPEC supported): for each control family, >= 8 of 10 valid
  seeds fail G1 or G2. G4 is not applied to controls.
- Any other outcome: the corresponding hypothesis is unsupported; report and
  go to Section 11.

## 11. Stop rules, escalation, repair policy

- Escalation (at most ONE, only if F2 triggers in Stage A): n=1800 events,
  128 ticks per chain, rerun Stage A with seeds 50-59. No other parameter may
  change.
- Instrument repair (metric definitions, fit budget, split logic) is allowed
  ONLY before the freeze, or after a failed Stage B via a versioned amendment
  (PC-V2) that reruns Stage A on new seeds. Repairs never move thresholds
  directly; they re-derive them through Section 10 on fresh calibration data.
- Threshold retuning after the freeze is forbidden. At most 2 repair cycles
  (PC-V2, PC-V3); if Stage B still fails, the representability approach at
  this scale is declared inadequate and the negative result is reported with
  the calibration record and a power analysis.
- Never: dropping seeds after looking at outcomes, swapping metrics post hoc,
  applying gates to report-only quantities, reusing seed ranges across
  stages.

## 12. Freeze procedure and provenance

- Freezing PC-V1 means: (1) this file is amended from DRAFT to FROZEN in a
  dedicated commit by the maintainer; (2) Stage A summary and the mechanical
  threshold table are committed as `docs/prereg/frozen/pc_v1_thresholds.json`
  together with the Stage A summary CSV under `docs/prereg/frozen/`; (3) the
  commit hash of the freeze is recorded inside the thresholds file.
- `pc02` (Stage B) and `pc03` (Stage C) refuse to run unless the frozen
  thresholds file exists. Their decision registries are committed under
  `docs/prereg/frozen/` as well. Raw per-run artifacts may live under
  `outputs/` but are never decision inputs.
- Every output row carries: scene digests, config hash, seed, stage,
  condition, code version (git describe), and the requested-vs-executed fit
  budget.

## 13. Claim boundary

- Stage B pass = the gates are jointly satisfiable by measured geometric
  order data at this scale. Not a spacetime-emergence claim.
- Stage C block = the gates reject matched geometry-free order data. Not a
  falsification of any physical theory.
- Pass + block = the pipeline is a validated discriminator, enabling the P1
  epsilon-sweep (geometry as manipulated variable) and any future application
  to dynamically generated order data.
- A Stage B fail after two repair cycles is a negative result about THIS
  instrument at THIS scale, with the calibration record quantifying why.

## 14. Deviations log

All entries below are PRE-FREEZE instrument repairs (permitted by Section 11,
which allows metric/split/control changes before the freeze). They were
prompted by a first Stage A calibration run (seeds 0-9, code 4a2d0fa) whose
specificity controls exposed two instrument defects. Thresholds were not set
from that run; Stage A is re-run after these repairs and thresholds derive
only from the re-run. Deviations do not retroactively change recorded
decisions.

- D1 (2026-07-08, parallax dissimilarity). The first calibration showed that
  a geometry-free order, once its density was made comparable (D2), still
  passed the embeddability gate: its bracket-width profiles were dominated by
  a per-target shared scalar (cross-chain column correlation ~0.7 vs ~0.0 for
  geometric scenes) that embeds in 1D without encoding space (fitted 1D
  coordinate correlated with x at ~0.1-0.5 vs ~0.98 for geometric scenes).
  Fix: compute dissimilarity over reference-centered parallax profiles
  (Section 6) so only cross-observer disagreement counts as spatial signal.
  After the fix, geometric scenes still recover x (corr ~0.98, heldout ~0.02)
  while density-matched random orders block (heldout ~0.22-0.30).

- D2 (2026-07-08, density-matched control). The first calibration's
  C-RANDOM-ORDER matched PRE-closure edge density, but transitive closure
  percolated to a near-complete order (density ~0.997), an unfair degenerate
  foil. Fix: tune the pre-closure edge probability by bisection to match the
  geometric POST-closure density (Section 9). The transitive-closure helper
  was also switched to a float32 BLAS matmul (identical result, ~140x faster)
  to keep the per-scene bisection affordable.

## 15. Confirmatory outcome (post-freeze factual record)

This section records the result of executing the frozen Stage B/C decisions.
It changes no design or threshold. Registries and per-seed CSVs are under
`docs/prereg/frozen/`.

- Frozen thresholds (commit b77f588, calibration provenance 9162e8e): heldout
  <= 0.05, null_gap >= 0.15, restart_order_disagreement <= 0.15,
  truth_order_error <= 0.15, at gate dim d=1.
- Stage B (sensitivity, seeds 100-119): H-SENS SUPPORTED. 20/20 valid seeds
  pass all four gates (rule >= 16/20); every individual gate passed 20/20,
  including the tight stability gate.
- Stage C (specificity, seeds 200-209): H-SPEC SUPPORTED. Both control
  families block 10/10 (rule >= 8/10): density-matched random_order (9
  evaluated-and-blocked + 1 scene-invalid structural block) and
  column_shuffled_geometric (10 evaluated-and-blocked).
- Conclusion: on fresh confirmatory seeds under frozen thresholds, the
  pipeline passes on measured Minkowski-geometric order and blocks on matched
  geometry-free order, with clean non-overlapping separation. PC-V1 validates
  the response-profile representability pipeline as a latent-geometry
  discriminator at this scale. This is not a spacetime-emergence claim; it
  licenses the P1 epsilon-sweep (geometry as manipulated variable) and
  application to dynamically generated order data.
