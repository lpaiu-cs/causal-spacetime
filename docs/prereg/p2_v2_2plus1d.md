# P2-v2: 2+1D robustness, remediated

Status: FROZEN v1 (2026-07-09). Calibration provenance: P2v2-A at commit
aee54e0, seeds 0-9 (0/10 scene-invalid). Frozen thresholds:
`docs/prereg/frozen/p2_v2_test_constants.json` (gate_truth 0.15, gate_heldout
0.10; coverage floor V>=18, >=80% of valid). After this point the thresholds
are locked and P2v2-B is confirmatory. A preregistered remediation of P2
(frozen, `docs/prereg/p2_2plus1d.md`), whose confirmatory H-SENS-2D missed the
16/20 bar for a scene-generation reason (20% 2+1D scene-invalidity against a
fixed /20 denominator; P2 §12). P2-v2 does not retune any P2 threshold; it
applies the two remedies P2 §12 named and re-runs the full calibrate -> freeze
-> confirm cycle on fresh seeds. This document lives in version control.

## 1. What changed from P2 (and why)

Exactly two upstream changes, both named in P2 §12; nothing else (the PC-V1
measurement/fit pipeline and the midpoint gate-placement rule are unchanged):

1. **Reduce scene-invalidity by scene size.** N is raised 2600 -> 3600. A
   20-seed scan (seeds 0-19) gives 0/20 scene-invalid at N = 3600 (40-44
   targets per scene) versus 2-4/20 at N = 2600, while the d=1-underfits /
   d=2-recovers separation is preserved (d=1 truth ~0.22-0.24, d=2 truth
   ~0.10-0.11). No other scene parameter changes (band 0.22, ring 0.25, 8
   chains, 96 ticks).
2. **Coverage-aware denominator (adopt the P1 rule).** Scene-invalid seeds are
   recorded and excluded from the denominator; the experiment additionally
   requires a coverage floor so invalidity cannot silently shrink the sample.

## 2. Hypotheses

Identical to P2: H-SENS-2D (2+1D geometric passes at d = 2 and recovers true 2D
order), H-DIM (effective dimension is 2: d = 1 underfits, d = 2 suffices),
H-SPEC-2D (density-matched geometry-free order blocks).

## 3. Design

Reuse the frozen PC-V1 pipeline and the P2 measurement/fit code (`sweep_seed`);
only the scene N and the decision denominator change. Fit d = 1, 2, 3; gate
dimension 2. P2v2-A calibration seeds 0-9 propose thresholds by the midpoint
rule; P2v2-B confirmatory FRESH seeds 500-519 apply the frozen thresholds.

## 4. 2+1D scene (v2)

`build_scene_2plus1d` with `Scene2DConfig(n_events = 3600)`; all other fields as
P2 (diamond T = 2.0, 8 stationary ring chains at radius 0.25, 96 ticks, target
band |t| <= 0.10 and radius <= 0.22, >= 30 targets, subsampled to <= 44).

## 5. Metrics

As P2: per seed and dim d in {1,2,3}, held-out violation and truth-order error
against true 2D coordinates; the geometry-free control's d = 2 held-out (or a
structural-block flag).

## 6. Threshold-setting rule (P2v2-A -> frozen)

Hard floors (fixed now): HF1 median d = 2 truth <= 0.20 and median d = 2
held-out <= 0.10; HF2 min d = 1 truth > the proposed truth gate.

Midpoint gate placement (unchanged from P2 deviation D1):

- gate_truth = min(0.20, round(0.5 * (max_A(d2 truth) + min_A(d1 truth)), 0.05))
- gate_heldout = min(0.10, round(0.5 * (max_A(d2 held-out) +
  min_A(control d2 held-out)), 0.05))

## 7. Decision rules (P2v2-B confirmatory), coverage-aware

Let V = number of scene-valid seeds of the 20. Report V and the scene-invalid
seed list.

- Coverage floor: V >= 18 (scene-invalidity <= 10%). If V < 18 the run is
  under-covered and inconclusive (report; do not claim support).
- H-SENS-2D supported: >= 0.80 * V valid seeds have d = 2 held-out <=
  gate_heldout AND d = 2 truth <= gate_truth.
- H-DIM supported: >= 0.80 * V valid seeds have d = 1 truth > gate_truth AND
  d = 2 truth <= gate_truth; also report the saturation fraction
  (d = 2 truth <= d = 3 truth + 0.05).
- H-SPEC-2D supported: >= 0.80 * V valid seeds block (control d = 2 held-out >
  gate_heldout, or structural pool failure).

The coverage floor plus the 80%-of-valid rule together prevent invalidity from
either shrinking the denominator unfairly (P2's failure mode) or hiding a real
sensitivity failure.

## 8. Stop rules and repair policy

Instrument repair only before the P2-v2 freeze, via a versioned amendment
re-running P2v2-A on fresh seeds. The frozen PC-V1 pipeline and the P2 midpoint
gate rule are unchanged. No threshold retuning after freeze; at most one further
remediation (P2-v3) if H-SENS-2D still misses, justified by root cause; no
dropping seeds after seeing outcomes.

## 9. Freeze and provenance

Mirror PC-V1 Section 12: commit the P2v2-A summary and thresholds to
`docs/prereg/frozen/p2_v2_test_constants.json` with the calibration commit
hash; P2v2-B refuses to run without it. Every row carries scene digests, seed,
stage, and code version.

## 10. Claim boundary

As P2: pass-at-2 / fail-at-1 / block-on-geometry-free shows the validated
discriminator recovers the correct spatial dimension and 2D structure of a 2+1D
causal order. Not a spacetime-emergence, continuum-limit, or dynamics claim.
P2-v2 is a scene-size and denominator remediation, not a change to what is
claimed.

## 11. Deviations log

(empty — no pre-freeze repair was needed; the P2v2-A calibration passed cleanly.)

## 12. Confirmatory outcome (post-freeze factual record)

Result of executing the frozen P2v2-B decisions (P2-v2 constants at aee54e0;
PC-V1 gates at 9162e8e), fresh seeds 500-519. Registry and per-seed CSV under
`docs/prereg/frozen/`. Changes no threshold or rule.

- **All three hypotheses SUPPORTED, 20/20 valid seeds** (0 scene-invalid;
  coverage floor met).
  - H-SENS-2D: 20/20 (d=2 held-out 0.026-0.072 all <= 0.10; d=2 truth
    0.090-0.116 all <= 0.15).
  - H-DIM: 20/20 (d=1 truth 0.203-0.276 all > 0.15 so 1D underfits; d=2 truth
    <= 0.15 so 2D suffices; recovery saturates at d=2 on 20/20).
  - H-SPEC-2D: 20/20 (every geometry-free control blocks).
- The two P2 remedies worked: raising N to 3600 eliminated scene-invalidity
  (0/20), and the seed with a borderline d=2 held-out in P2 (0.133) did not
  recur (max d=2 held-out 0.072 here). The coverage-aware denominator was not
  even needed (V = 20) but is in force.

Conclusion: P2-v2 confirms the full 2+1D robustness result — the validated
discriminator passes on measured 2+1D geometric order, recovers true 2D
structure, selects the correct effective spatial dimension (2, not 1), and
blocks on matched geometry-free order, all on fresh confirmatory seeds under
frozen thresholds. This completes the H-SENS-2D confirmation that P2 missed for
a scene-generation reason. It remains a controlled result: not a
spacetime-emergence, continuum-limit, or dynamics claim.
