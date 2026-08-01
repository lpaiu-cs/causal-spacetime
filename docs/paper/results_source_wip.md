# Results source of truth (WIP) — P8 through P13

**Status: WORK IN PROGRESS. Not a manuscript, and deliberately not a
closing account.** The programme is continuing (the curvature track
closed 2026-08-01 and handed over a `d >= 3` lever — see
`docs/prereg/p14_weyl_curvature.md`), so this file exists only so that a
future write-up quotes numbers from ONE place, each traceable to a frozen
artifact. Nothing here interprets beyond what the frozen records already
state; where a record carries a scope or caveat, it is reproduced, not
smoothed.

Rule inherited from the records: every number below is read from an
artifact, never from a console or from memory (the defect class that
recurred three times in P11 — see `p11_continuum_metric.md` 9.6).
Measured timings are quoted to one significant figure for the same
reason.

## Fixed-size emergence (what the order alone gives)

| result | claim | record |
|---|---|---|
| P8 | 3+1D arm, informative no-gate | `docs/prereg/p8_3plus1d.md` |
| P9 | dimension SELECTED from a 3+1D order, 20/20/20 within-seed, SUPPORTED | `docs/prereg/p9_paired_dimension_rule.md` |
| T1 G4c | rigidity threshold `R >= d + 2`, `n = d^2 + 3d + 1`; 3+1D pinned at R=5, n=19; measured d=1..6 over 20 preregistered predictions, proved up to one finite-dimensional hypothesis | `docs/theory/t1_g4c_proof.md` |

## The continuum question (two instrument families, opposite outcomes)

**P10 — the discriminator family: PROCEDURAL CLOSURE.** Three frozen
gates failed as frozen. The mechanism reading is scoped: no
improvement is detected in any arm, while "measurably degrades" is a
frozen-sample result that replicates in sign, ordering and dense-tail
direction but NOT in CI separation under a fully stream-separated
replication (0 of 16 pointwise bins separate; the sparse extreme tail
inverts). Gate failure is a stopping rule operating, not evidence for
the null. Record: `docs/prereg/p10_continuum_limit.md` 9.13, 9.10.

**P11 — the metric family: three gates, all IMPROVES.** All artifacts
at one stamp, `535a5b7`.

| stage | quantity | Delta (600 -> 2400) | 95% CI | verdict |
|---|---|---|---|---|
| A | proper time, chain estimator | -0.187 | [-0.289, -0.091] | IMPROVES |
| B | spacelike distance, dual-box chain | -0.213 | [-0.303, -0.122] | IMPROVES |
| C | coordinates, metric-only trilateration | -0.276 | [-0.384, -0.171] | IMPROVES |

Derived targets: `Delta*_A = Delta*_B = Delta*_C = -0.2007` dex, all
from the same longest-chain law. Labelled slope checks: A `-0.311`,
B `-0.354`, C `-0.458` against a predicted `-1/3`, the last steeper by
prediction stated in advance (the coordinate fit compounds three
measurements). Median Stage C coordinate error falls
`0.156 -> 0.107 -> 0.089` across the ladder.

Scope carried from the records, not to be dropped in any write-up:
- claims bind to this estimator family over this pre-asymptotic
  ladder (`m ~ 24-240`), never to the continuum limit as such;
- inputs are relabeling-invariant order quantities PLUS the count as
  density calibration — order fixes conformal structure, the count
  fixes scale ("order + number = geometry");
- one Stage C sample in 48 is realizer-assisted; removing it gives
  `Delta = -0.294`, CI `[-0.396, -0.190]`, and its inclusion is the
  conservative direction;
- Stage B's box membership is order data (conjugate-order
  betweenness), certified per scored pair: 215 of 216 certified, the
  one exception in the middle rung the gate never uses.

## Curvature (does the same instrument read CURVED geometry, and how far)

Two experiments, both closed 2026-07-31, each with its own frozen stamp.
`dS_2` throughout, Ricci scalar `R = 2` in units `ell = 1`.

**P12 Stage A-12 — the instrument reads curved proper time.** Stamp
`6d1d1bf`. `Delta = -0.104`, CI `[-0.177, -0.041]`, **IMPROVES**, at
`tau/ell ~ 0.3`. Flat verdict was unavailable at this variance
(`flat_available: false`), which re-scopes what a non-IMPROVES outcome
would have meant; no selection caveat.

**P12 Stage B-12 — the instrument recovers the curvature itself.** Stamp
`5252ee7`, `n = 1117` per rung per arm, operating point `tau/ell = 1.5`,
so `R tau^2 = 4.5`.

| gate | quantity | value | 95% CI | outcome |
|---|---|---|---|---|
| (i) rate | `Delta_B` over `m` 18.87 -> 75.77 | `-0.2637` | `[-0.3037, -0.2233]` | IMPROVES |
| (ii) recovery | top-rung `\|R tau^2 hat - R tau^2\| / R tau^2` | `0.1708` | `[0.0784, 0.2644]` | passes `<= 0.25` |

Verdict **RECOVERS-CURVATURE**. The derived target was `Delta*_B =
-0.2508`, propagated from the Poisson and chain rates rather than
borrowed from P11's `-0.2007`; the campaign sits `0.63 sigma` from it.

Per-rung recovery, all three undershooting `1 - Q` by 7 to 17% — the
residual bias mismatch the design priced in advance:

| rung | `Q_hat` | `R tau^2 hat` | rel error | 95% CI |
|---|---|---|---|---|
| 600 | 0.92382 | 4.1540 | 0.0769 | `[0.0038, 0.2410]` |
| 1200 | 0.93228 | 3.6386 | 0.1914 | `[0.0705, 0.3138]` |
| 2400 | 0.93074 | 3.7316 | **0.1708** | `[0.0784, 0.2644]` |

**P13 Stage A-13C — how far the flat-normalized reading goes.** Stamp
`33e371f`, `n = 1731` curved / `1650` twin, ladder `tau/ell` in
`{0.3, 0.6, 1.0, 1.5}`. Endpoint contrast `Delta = +0.0072`, CI
`[-0.00039, +0.01481]`, against an equivalence margin `delta_eq = 0.02`.
Verdict **CURVATURE-ROBUST** with **CONTROL-CLEAN** and no control
caveat. The reading does not notice curvature out to `R tau^2 = 4.5`,
i.e. through a 25-fold growth in `(tau/ell)^2`.

Scope carried from the records, not to be dropped in any write-up:

- **The verdict word is not an accuracy figure** (P12 §10.9). Any
  sentence quoting `0.1708` must carry `[0.0784, 0.2644]` with it: the
  frozen rule is a POINT comparison and the interval **crosses** the
  0.25 threshold. Realized headroom `1.67 sigma` against `1.92 sigma`
  disclosed pre-run. §9.3 declined to re-score on an interval rule after
  seeing the interval; the rule for future designs is now fixed in
  `AGENTS.md`, and it is not retroactive.
- **Both records confine themselves to the instrument.** In 1+1D every
  metric is conformally flat, so curvature reaches the causal order only
  through the volume this estimator reads — the "number" half of P11's
  "order + number = geometry". The "order fixes conformal structure"
  half has never been tested at finite density, because every ensemble
  so far had flat conformal structure. Testing it needs non-vanishing
  Weyl, which needs spacetime dimension >= 4; at the continuum level the
  order-signal's existence is a theorem (causal order determines the
  conformal class), so the open question is finite-density
  detectability. See `docs/prereg/p14_weyl_curvature.md`.
- **Selection caveat: YES on both.** P12 Stage B replaced two curved
  seeds at rung 600 and one each at 1200 and 2400 from reserve slots
  (twin arm skipped none); identities are published in the artifact. The
  caveat also inherits from P13 Stage A-13C, which carries its own.
- **The cross-arm gate passed with a real chance of having tripped.**
  Spurious-trip probability `0.154 / 0.064 / 0.001` per rung, computed
  against each rung's own frozen offset (`+0.55% / -0.48% / -0.15%`).
  About one campaign in six would have halted at the bottom rung on
  correctly calibrated arms. The tolerance was not widened for it.
- **P13's `CURVATURE-ROBUST` is an ARM difference, not an identified
  geometric one** (its §13.2), and the twin matches the curved arm's
  residual to `+0.00007`.

## What a future manuscript still needs

Not gaps in the record — gaps in the writing: a figure programme, the
P10/P11 contrast stated as one narrative (why a fixed-relative-
resolution instrument cannot answer a question a counting-measure
instrument can), positioning against the causal-set literature
(`docs/theory/t1_literature_positioning.md` is the source), and a
decision on whether the curvature track is written up now or held until
the `d >= 3` question has an answer.
