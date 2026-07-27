# Results source of truth (WIP) — P8 through P11

**Status: WORK IN PROGRESS. Not a manuscript, and deliberately not a
closing account.** The programme is continuing (P12 curvature is in
design), so this file exists only so that a future write-up quotes
numbers from ONE place, each traceable to a frozen artifact. Nothing
here interprets beyond what the frozen records already state; where a
record carries a scope or caveat, it is reproduced, not smoothed.

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

## What a future manuscript still needs

Not gaps in the record — gaps in the writing: a figure programme, the
P10/P11 contrast stated as one narrative (why a fixed-relative-
resolution instrument cannot answer a question a counting-measure
instrument can), positioning against the causal-set literature
(`docs/theory/t1_literature_positioning.md` is the source), and
whatever P12 returns.
