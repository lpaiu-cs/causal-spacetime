# P8: 3+1D robustness and dimension selection

Status: **STAGE A RUN; NO GATE PLACED; NOT FROZEN** (2026-07-25).
Stage A produced no gate because the pass and fail clusters overlap --
the outcome Section 6 named in advance as admissible. P8-B has not run
and cannot until a versioned amendment re-calibrates on fresh seeds.
See Section 12.

This document is the preregistration. Sections 1-11 were written before
Stage A ran and are unedited since; Section 12 is the factual record of
what Stage A produced. No gate value appears anywhere, because none was
placed.
The scene was chosen from an exploratory sweep whose numbers are recorded
in Section 4 and are not permitted to become thresholds.

## 1. Motivation

Every empirical arm of this programme lives in one or two spatial
dimensions: PC-V1 and P1 in 1+1D, P2/P2-v2 in 2+1D. The theory has since
outrun them. T1/G4c states an identifiability law for every spatial
dimension and pins the physical case exactly — in 3+1, five observers and
19 targets determine the scene up to congruence — and A1 states the
operating conditions that go with it. None of that has been put to a
discriminator on measured causal sets in the dimension it is about.

P8 asks the P2 question one dimension up: does the validated instrument
still pass on 3+1D Minkowski-sprinkled causal order, still select the
correct spatial dimension, and still block on matched geometry-free
order?

It is worth saying in advance what a negative result would mean, because
the design work below already shows this is not a formality. If H-DIM
fails, the finding is that the instrument's *dimension selection* weakens
as dimension grows — a real result about the instrument, reportable as
such, and not a defect in the arm.

## 2. Hypotheses

- **H-SENS-3D.** The discriminator passes on 3+1D Minkowski-sprinkled
  causal sets: at the gate dimension `d = 3`, held-out violation and
  truth-order error both stay under their gates.
- **H-DIM-3D.** It selects the correct spatial dimension: `d = 2`
  underfits (truth error above gate) while `d = 3` does not, and recovery
  saturates at 3 rather than continuing to improve at `d = 4`.
- **H-SPEC-3D.** It blocks on a density-matched geometry-free control.
- **H-OBS (secondary, descriptive).** The observer count matters in the
  direction the conditioning study predicts: an arm at `R = 8` separates
  the dimensions less well than the primary `R = 12`. Reported as a
  measured contrast, not gated.

## 3. Design

Reuse the frozen PC-V1 measurement/dissimilarity/fit pipeline unchanged.
The only new element is 3+1D scene construction (Section 4). Fit at
`d = 2, 3, 4`; the gate dimension is 3.

Stage A calibration proposes thresholds mechanically and is then frozen;
Stage B applies them to fresh seeds. Invalid seeds are recorded, never
silently replaced.

Two deliberate differences from P2, both logged in Section 11:

- **Seed hygiene.** Scene selection used seeds 0-3. Stage A therefore
  uses seeds **100-119** rather than P2's 0-9, so no seed that informed
  the instrument also sets its gates. Stage B uses **400-419**.
- **Stage A is twenty seeds, not ten.** The pass and fail clusters are
  measurably closer in 3+1D than in 2+1D (Section 4), so the gate is
  being placed in a narrower gap and deserves a better-resolved estimate
  of where the clusters sit.

## 4. 3+1D scene, and how it was chosen

`causal_spacetime_lab.positive_control.scene_3d.build_scene_3plus1d`:

- Bulk: `sprinkle_minkowski_causal_diamond(n = 7000, spacetime_dim = 4,
  T = 2.0)`.
- References: `R = 12` stationary observer chains on a Fibonacci sphere
  of radius 0.25, 96 ticks each over `t` in `[-0.7, 0.7]`, all inside the
  diamond. The builder asserts the shell has affine rank 3 — a coplanar
  observer set is exactly the degenerate case for multilateration, so it
  is checked rather than assumed.
- Causal order: `causal_matrix_minkowski`.
- Targets: bulk events with `|t| <= 0.10` and `||x|| <= 0.30`, two-sided
  bracketed by all 12 chains; `>= 30` required, subsampled to `<= 44`.
- Provenance: sha256 digests of events and causal matrix on every row.

**Exploratory sweep behind those choices (seeds 0-3; no gate may be
derived from these numbers).** Transplanting P2's layout unchanged leaves
almost no room between the correct dimension and one short of it. The
separation between the `d = 3` and `d = 2` truth-error clusters:

| `R` | 6 | 8 | 12 | 16 | 20 |
|---|---|---|---|---|---|
| gap | 0.015 | 0.016 | 0.035 | 0.038 | 0.031 |

against **0.074** for P2 in its own dimension. There is a step of about
two between `R <= 8` and `R >= 12` and a turnover by 20. Four seeds
cannot separate 12 from 16 and no such claim is made; 12 is the cheap end
of the plateau. It is also where the conditioning study independently put
the margin optimum for three spatial dimensions, which is corroboration
and not evidence — that study measured a different quantity.

Two levers that did not work, recorded so they are not retried:

- Doubling ticks per chain to 192: separation `0.035 -> 0.034`. The
  readout resolution is not what binds here, consistent with A1 finding
  the instrument already comfortable at 96.
- Raising the target cap to 80: separation `0.035 -> 0.016`, and every
  error worse. The frozen fit policy spends a fixed constraint budget, so
  more targets divide it more thinly. **Target count is a lever in the
  wrong direction under a frozen policy.**

Even at `R = 12` the expected gate margin is roughly half P2's. That is
stated here, before the fact, so that a tight Stage A result is not
reported later as a surprise.

**A second thin margin, on the other gate.** A three-seed development run
put `d = 3` held-out at 0.058-0.087 against P2's `d = 2` range of
0.026-0.072, and the midpoint rule then wants a held-out gate of about
0.12 — which the standing ceiling of 0.10 pulls back down. So the
held-out gate is expected to sit *at the ceiling*, with only about 0.013
of pass-side margin.

The ceiling is deliberately not raised. It is a programme-wide standard
for acceptable held-out violation, and loosening it because the harder
case is harder would be exactly backwards. The consequence is that
**H-SENS-3D may fail**, and if it does the finding is a real one and
should be reported as it stands: held-out violation in 3+1D exceeds the
ceiling the programme holds every other arm to.

## 5. Metrics

Per seed, at each `d` in `{2, 3, 4}`: held-out violation rate and
truth-order error (sign discordance of fitted pair distances against true
3D pair distances). Plus the geometry-free control's `d = 3` held-out, or
a structural-block flag.

## 6. Threshold-setting rule (P8-A -> frozen)

Identical in form to P2's, so the two arms remain comparable:

- `gate_truth` is placed at the midpoint between the pass cluster
  (`d = 3` truth over valid Stage A seeds) and the fail cluster (the
  smaller of `d = 2` truth and the control), so both sides keep margin.
- `gate_heldout` is placed at the midpoint between the `d = 3` held-out
  cluster and the control held-out cluster.
- Both are computed mechanically from the Stage A table by
  `p8_3plus1d.py --stage a`. No hand adjustment. Re-running at the frozen
  commit must reproduce them.
- `pass_min = 16`, `denominator = 20`, as P2.

If the two clusters overlap, no gate is placed, and P8 is reported as a
failure to establish the instrument in 3+1D. That outcome is admissible
and is not grounds for changing the scene after the fact.

## 7. Decision rules (P8-B confirmatory)

- **H-SENS-3D supported**: `>= 16` of 20 valid seeds have `d = 3`
  held-out `<= gate_heldout` AND `d = 3` truth `<= gate_truth`.
- **H-DIM-3D supported**: `>= 16` of 20 valid seeds have `d = 2` truth
  `> gate_truth` AND `d = 3` truth `<= gate_truth`; additionally report
  the fraction with `d = 3` truth `<= d = 4` truth + 0.05 (recovery
  saturates at 3).
- **H-SPEC-3D supported**: `>= 16` of 20 valid seeds block (control
  `d = 3` held-out `> gate_heldout`, or structural pool failure).
- **H-OBS**: descriptive only. The `R = 8` arm is run on the same Stage B
  seeds and its separation reported beside the primary. No gate, no
  claim of significance.

## 8. Stop rules and repair policy

Scene parameters and gate construction may be repaired only *before* the
P8 freeze, via a versioned amendment that re-runs Stage A on fresh seeds.
After the freeze, nothing about the scene, the gates or the decision
rules may change; deviations are logged in Section 11 and reported.

A Stage B seed that raises `SceneValidityError` is recorded as invalid
with its reason and excluded from the denominator; the denominator is not
topped up from further seeds.

## 9. Freeze and provenance

To be completed at freeze: frozen commit, PC-V1 frozen commit, Stage A
table and summary under `docs/prereg/frozen/`, `p8_test_constants.json`
carrying the gates and their provenance block.

## 10. Claim boundary

- P8 concerns the *discriminator*, not the identifiability theory. A pass
  does not verify G4c, and G4c's `R >= d + 2` threshold is not what is
  being tested: the arm runs at 12 observers, far above 5, and its
  targets are selected for bracketing rather than for rigidity.
- The truth-order error is scored against coordinates the builder knows.
  That is a validation instrument, not something the pipeline has access
  to.
- Nothing here speaks to continuum emergence, to quantum dynamics, or to
  causal sets not built by sprinkling.
- The exploratory numbers in Section 4 are descriptive characterization
  of the instrument at design time. They are not results, and they are
  not gates.

## 11. Deviations log

- **D1.** Stage A seeds are 100-119 rather than P2's 0-9, and there are
  twenty rather than ten. Reason: seeds 0-3 informed the scene choice, so
  reusing them to set gates would let the instrument grade itself; and
  the narrower cluster gap in 3+1D warrants a better-resolved estimate.
- **D2.** Twelve observer chains rather than P2's eight, and a target
  band radius of 0.30 rather than 0.22. Both are Section 4 choices, made
  before freeze, from the exploratory sweep recorded there.
- **D3.** The gate rounding grid is 0.01 rather than P2's 0.05. P2's grid
  is comparable to P8's entire cluster gap, so rounding could carry the
  gate out of the interval it is meant to sit inside. A grid has to be
  small against what it discretises.
- **D4.** Seeds 100-102 — three of the twenty Stage A seeds — were run
  once during development, as a smoke test that the geometry-free control
  arm builds at all in 3+1D (it initially did not; the density helper
  returns a triple). Recorded because those seeds are no longer wholly
  unseen. The scene was already fixed from seeds 0-3 at that point and
  nothing was changed in response, but the fact belongs on the record
  rather than in a commit message.

## 12. Stage A outcome (factual record)

**Run 2026-07-25 at commit `9f6e495`, seeds 100-119, all 20 valid.**
**No gate was placed. P8 is not frozen and P8-B has not run.**

Section 6 named this outcome in advance and it is the one that occurred:
both cluster pairs overlap at their extremes, so the midpoint rule has no
interval to place a gate in.

| | pass cluster | fail cluster | overlap |
|---|---|---|---|
| truth | `d = 3`: 0.1097 - **0.1395** | `d = 2`: **0.1328** - 0.1998 | 0.0066 |
| held-out | `d = 3`: 0.0270 - **0.0990** | control: **0.0950** - 0.2490 | 0.0040 |

Artifacts: `docs/prereg/frozen/p8_stage_a_calibration.csv` and
`p8_stage_a_summary.json`. There is deliberately **no
`p8_test_constants.json`** — its absence is the record, and the P8-B
guard refuses to run without it.

### What kind of failure this is

The instrument is not indifferent to dimension. **`d = 3` truth error is
below `d = 2` truth error in 20 of 20 seeds**, without exception, and the
medians are cleanly apart (0.1237 against 0.1712). What fails is the
*decision rule*, which places one absolute gate across seeds and
therefore needs the two populations to separate — not merely the
within-seed comparison to come out right. In 2+1D the populations did
separate; in 3+1D they interleave at the edges. Three seeds carry it:
105 (`d2 = 0.153`, `d3 = 0.136`), 111 (`d3 = 0.139`), 114
(`d2 = 0.133`).

This is recorded as a description of the failure, not as a proposal. A
paired or within-seed rule would very likely pass on this data, and that
is exactly why adopting one now is not available: it would be a decision
rule chosen after seeing the numbers it is applied to. A different rule
needs its own preregistration and its own fresh seeds.

### What is and is not concluded

- **Not concluded:** that the discriminator fails in 3+1D. It orders
  every seed correctly.
- **Concluded:** that the P2 decision procedure, transplanted to 3+1D at
  this configuration, does not yield a gate. The design note in Section 4
  predicted a margin about half P2's from four exploratory seeds; twenty
  seeds show the margin is not merely thin but absent under a min/max
  cluster rule.
- **Also on the record:** the held-out clusters overlap too, at 0.0040.
  The ceiling of 0.10 was never reached — `d = 3` held-out topped out at
  0.0990 — so H-SENS-3D's anticipated failure mode did not materialise;
  a different one did.

### Permitted next steps

Section 8 allows exactly one repair route before a freeze: a **versioned
amendment**, stating what changes and why, re-running Stage A on **fresh
seeds**. Candidate amendments, none of which may be adopted without that:

1. More observers (`R = 16`), the only lever the design sweep showed
   still moving in the right direction.
2. A within-seed decision rule, preregistered as such, with its own
   calibration seeds.
3. Reporting P8 as a negative result at this configuration and stopping.

The choice is a decision for the programme, not a repair to be made
quietly.
