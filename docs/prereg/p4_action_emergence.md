# P4: action-weighted emergence in the restricted 2D-orders ensemble

Status: FROZEN v1 (2026-07-10). Gates and predictions in
`docs/prereg/frozen/p4_test_constants.json`, frozen before the confirmatory
run on seeds disjoint from all calibration/exploration seeds.

## 1. Question

P3 established that a geometry-free growth dynamics (transitive percolation)
does not produce geometry. E1 (exploratory) established that weighting the
UNRESTRICTED poset ensemble by the verified 2D Benincasa-Dowker action does
not either: the raw action's low-action region is the complete bipartite
crystal (S = -2(N^2/2 - N)), not the manifoldlike region, and the smeared
(nonlocal) action does not change this structurally -- the bipartite order
has S_eps = eps N (2 - eps N) exactly, below the sprinkled value whenever
eps N >> 2 (i.e. whenever the nonlocality scale still resolves the system).

P4 asks the remaining controlled question: in the RESTRICTED ensemble of 2D
orders (permutations; the intersection of two total orders), with the
verified smeared action at fixed eps, does a continuum (manifoldlike) phase
exist at finite beta, and does it give way to the crystal phase at large
beta? The restriction kills the Kleitman-Rothschild entropy problem but does
NOT exclude the crystal geometrically (the bipartite order is the
permutation (k..1, 2k..k+1)) -- so the continuum/crystal competition is
non-trivial and the outcome is not built in.

## 2. Instruments (validated before freezing)

- Action: `positive_control/action.py`, conventions verified against
  Dowker-Glaser 2013 (arXiv:1305.2588). Tests: eps=1 reduces to the raw
  action exactly; bipartite closed form; flat-sprinkling value ~ 0 with
  Sorkin fluctuation damping.
- Sampler: Metropolis over permutations (`positive_control/two_orders.py`).
  Tests + stage A: exact Gibbs at N=6 (TV = 0.018, all 720 states visited);
  beta=0 reproduces the uniform ensemble, which is verified to coincide with
  sprinkled 2D diamonds (identity deviation 0.17%).
- Observables: interval abundances (n0, n1, n2), height, mean action.
  MM dimension is recorded but EXCLUDED from gates: crystals also sit at
  MM ~ 2 (E1/E2 lesson); it cannot separate the phases.

## 3. Design

N = 100, eps = 0.12 (eps N = 12: nonlocality scale well inside the system),
beta in {0, 1, 2, 3, 4, 6}, dual start per beta (R = seeded random
permutation, X = bipartite crystal), seeds 100-104 per (beta, start): 60
chains of 400k Metropolis steps, second half sampled every 1000.

Reference profile (frozen, stage A, sprinkled N=100 seeds 0-29):
(R, n0, n1, n2, height) = (2486.7, 324.4, 228.6, 181.3, 17.0).

## 4. Frozen classification gates

Per chain (means over samples):
- continuum: |n0/324.43 - 1| <= 0.5 AND height >= 12
- crystal: (n1 + n2) <= 0.5 x 409.97 AND height <= 9
- otherwise "other".

Per beta (10 chains):
- "continuum" iff all chains continuum AND |<S>_R - <S>_X| <= 5
- "crystal" iff all chains crystal
- else "transition/hysteresis".

## 5. Frozen predictions

beta 0, 1, 2 -> continuum; beta 3 -> transition/hysteresis; beta 4, 6 ->
crystal. Experiment PASS iff (i) continuum at {0,1,2}, (ii) crystal at 6,
(iii) non-continuum at 3 and 4. Predictions are sharper than the pass
criteria and are scored separately (beta=4 rests on a ~7% margin in the
exploratory data and may honestly land in transition/hysteresis).

## 6. Interpretation rules

- PASS = a continuum phase of the action-weighted restricted ensemble exists
  at finite beta and crystallizes beyond a transition in (2, 4): emergence
  requires BOTH the ensemble restriction (entropy control; P3/E1 show its
  necessity) AND bounded coupling (crystallization shows beta cannot be
  large). Not a claim that geometry is dynamically generated from nothing:
  at beta=0 the restricted ensemble is already sprinkling-equivalent; the
  claim is about survival and destruction of geometry under action
  weighting.
- Any continuum verdict at beta >= 4 or crystal verdict at beta <= 1 would
  contradict the exploratory picture and would trigger a deviations-log
  investigation (mixing failure vs genuine physics) before any claim.
- Slow mixing at high beta (acceptance ~ 0.003-0.05) is expected; the
  high-beta claim is crystallinity of the reached states (profile-based),
  not full equilibration. Hysteresis (start-gap > 5) at beta=3 is the
  first-order-like signature seen in exploration.

## 7. Deviations log

- (pre-freeze) order_height initially assumed topologically-sorted indices;
  wrong on sprinkled matrices (reference height 8.8 vs MCMC 16.8 exposed the
  bug). Fixed via ancestor-count topological order + labelling-invariance
  test BEFORE freezing; stage A re-run (reference height 17.0).

## 8. Confirmatory outcome

(to be recorded after the stage-B run; changes no gate or rule)
