# P7: Geometry order parameter and finite-size scaling

Status: **N=600 CHARACTERIZATION COMPLETE; CONFIRMATORY FSS NOT FROZEN**
(2026-07-14).

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

## Frozen N=600 characterization grid

The grid is beta in `{12, 14, 16, 18, 20, 22, 24, 28}` at N=600 and
epsilon=0.02. Existing P6a chains at 12, 16, and 24 are reused without
rerunning. At each added beta, run three 3M-step chains under the exact replay
sampler: two seeded random starts and one exact bipartite start, with 60%
burn-in and nominally 48 retained samples every 25,000 steps. Evaluate G on
the first, middle, and last retained states.

For every chain report acceptance, phase label, action/n0/height IAT and ESS,
and the three G values. Also report start-separated mean action, n0
susceptibility `var(n0)/N`, centered n0 Binder cumulant, height, and G. The
screening rule marks a beta mixing-adequate only when all three starts share a
phase label and the minimum ESS over action/n0/height and chains is at least
20. This is a warning rule, not a proof of convergence. Non-adequate starts
must not be pooled into an equilibrium histogram or transition estimate.

The added betas are characterization-only and carry no expected phase or G
verdict. Constants and seeds are frozen in
`p7_n600_stage_a_constants.json` before executing any added chain.

## Confirmatory status

No H-G-MONO, H-COINCIDE, or H-HYST prediction is frozen yet. Those hypotheses
may be fixed only after the N=600 characterization grid and a viable mixing
strategy are established. N=900/1200 confirmatory work must use fresh seeds
and may not proceed by silently averaging non-converged starts.

## N=600 characterization outcome

All 15 added chains completed, and the three P6a beta points were reused as
frozen. Only beta = 12 and 14 passed the prespecified mixing screen.

- Beta 12 and 14: all starts reached the continuum label. Minimum ESS was
  31.4 and 29.4, respectively. Mean G was 0.98--0.99 for random starts and
  0.99/0.81 for bipartite starts.
- Beta 16: all starts had a continuum mean label, but the bipartite chain had
  action ESS 2.34, mean action -25.7 versus +1.0 for random starts, and mean
  G 0.26 versus 0.95. This point is not mixed.
- Beta 18--24: random starts stayed continuum with mean G 0.95--1.00, while
  bipartite starts stayed intermediate with G = 0. Their start-separated
  action gaps grew from 76.8 to 95.9. These are explicit finite-run
  hysteresis points, not equilibrium mixtures.
- Beta 28: all starts were labelled intermediate and had G = 0, but the
  random and bipartite means still differed by 32.4 in action and 9,773 in
  n0. Random-start action ESS was only 3.50--3.97. Matching coarse labels do
  not establish convergence.

Operationally, the bipartite basin loses reconstructable geometry between
beta 14 and 16, phase-label disagreement begins by beta 18, and the random
basin loses geometry between beta 24 and 28. This wide, start-dependent
window is the principal Stage A observation. It prevents an equilibrium
estimate of either `beta_c^geo` or `beta_c^thermo`, so H-G-MONO and
H-COINCIDE remain unfrozen.

Naive coarse replica exchange is not an adequate automatic fallback here.
At beta 18 the observed start-separated action gap is about 77; a direct
cross-basin exchange with target acceptance 0.2 would require a local beta
spacing of order `-log(0.2)/77 ~= 0.02`, far denser than this grid and still
not a guarantee of barrier crossing. Before N=900/1200, P7 therefore needs a
validated enhanced-sampling method (for example multicanonical/umbrella
sampling with overlap checks, or a demonstrated sufficiently dense replica
ladder). Scaling non-converged local Metropolis chains would add compute but
not resolve the scientific question.

Frozen outputs include all added raw chains and instrument snapshots,
start-separated chain/beta summaries, the JSON registry, and the
`p7_n600_stage_a_observables` figure.

## Deviations log

(empty)
