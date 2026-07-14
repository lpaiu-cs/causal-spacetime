# P7: Geometry order parameter and finite-size scaling

Status: **G DEFINITION FROZEN; N=600 GRID NOT YET RUN** (2026-07-14).

## Question

Does loss of instrument-level reconstructable geometry coincide with the
thermodynamic crystallization transition of the smeared-action 2D-orders
ensemble, or does one precede the other? P7 treats the reconstruction verdict
as a continuous observable while retaining action, abundance, height, and
mixing diagnostics.

## Frozen geometry score

For a numerically valid P3/P5 reconstruction, define normalized margins

```text
m_h = (0.10 - heldout) / 0.10
m_n = (null_gap - 0.10) / 0.10
m_t = (0.40 - truth_error) / 0.40
G   = clip(0.5 + min(m_h, m_n, m_t), 0, 1).
```

A structural block has `G = 0`. Thus `G >= 0.5` is exactly equivalent to
passing all three frozen gates; no new fitted decision threshold is
introduced. True coordinates are used only to evaluate the order-intrinsic
fit, never as an input to reconstruction.

## G calibration outcome

The definition was evaluated once on previously frozen P1/P3/P5/P6 data
before any new P7 chain was run. All prespecified requirements were met.

- P3 sprinkled controls: 10/10 above 0.5; median G = 0.977.
- P5 uniform controls: 10/10 above 0.5; median G = 1.000.
- P5 beta=2 and beta=8 continuum configurations: 18/18 above 0.5.
- P5 beta=32 crystal: 4/4 structural zeros.
- P3 column-shuffle: the heldout-only upper bound on G is 0 for 10/10.
- P1 dilution analogue: G decreases with epsilon for all 20 seeds; median
  Spearman rho = -0.939.
- P6 layered hard negatives: median G = 0, but 11/160 are at or above 0.5.
  G is therefore a continuous form of the validated gate verdict, not a
  complete or uniquely sufficient manifoldlikeness observable.

The frozen calibration rows and summary are
`p7_stage_a_geometry_score.csv` and
`p7_stage_a_geometry_score_summary.json` under `docs/prereg/frozen/`.

## Inference constraint from P6a

At N=600, beta=12 mixed cleanly across starts. At beta=16 the bipartite-start
chain had ESS only 2--3 despite an eventual continuum mean label, and at
beta=24 random and bipartite starts remained in different basins. Therefore
phase labels alone cannot identify beta_c. P7 must report each start
separately, make convergence a prerequisite for equilibrium inference, and
not interpret a pooled two-start histogram as an equilibrium double peak.

## Confirmatory status

No H-G-MONO, H-COINCIDE, or H-HYST prediction is frozen yet. Those hypotheses
may be fixed only after the N=600 characterization grid and a viable mixing
strategy are established. N=900/1200 confirmatory work must use fresh seeds
and may not proceed by silently averaging non-converged starts.

## Deviations log

(empty)
