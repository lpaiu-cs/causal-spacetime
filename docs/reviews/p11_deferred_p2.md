# Deferred P2 review items (P11 campaign)

Tracking document for P2-severity review findings deliberately
deferred during the P11 merge cycle (user policy 2026-07-27: P1 and
data-level findings handled in-round; P2 wording/artifact-polish
findings batched here). Each item carries its source, the finding,
and the intended fix. Fixes land on this branch; the PR stays draft
until they do.

## 1. Publish the calibrated per-rung variance bound

- Source: PR #28,
  [comment 3654551205](https://github.com/lpaiu-cs/causal-spacetime/pull/28#discussion_r3654551205)
  on `docs/prereg/frozen/p11/p11_pilot_summary.json`.
- Finding: the pilot artifact's `per_rung.variance_bound_95` records
  the NOMINAL Bonett bound, while `power.s2_90` is computed from the
  CALIBRATED bounds (`z_used` after the bootstrap coverage check), so
  the artifact's per-rung field does not reproduce the power input.
- Intended fix: `run_pilot` / `run_pilot_b` record the calibrated
  per-rung bound (and keep the nominal one under a distinct name);
  regenerate + refreeze both pilot artifacts at the fix commit (gate
  numbers must be bit-identical — the bound values feeding power are
  unchanged, only the published fields gain the calibrated copy).

## 2. Correct the Stage B projected-wall narrative

- Source: PR #30,
  [comment 3654994020](https://github.com/lpaiu-cs/causal-spacetime/pull/30#discussion_r3654994020)
  on `docs/prereg/p11_continuum_metric.md` (Section 9.4).
- Finding: the narrative says "Projected Stage B wall 1.5 s", but the
  frozen artifact records `projected_stage_b_hours = 0.0010651...`
  = 3.83 s — the certificate-enabled regeneration tripled per-sample
  cost and the narrative kept the pre-certificate number.
- Intended fix: quote the artifact (3.8 s) with the standard
  correction marker, same class as the g4 misread fixed in 9f97328
  (narrative numbers must be read from artifacts, never from console
  tails or memory).
