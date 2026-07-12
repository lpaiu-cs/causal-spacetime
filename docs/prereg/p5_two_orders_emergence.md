# P5: does the action-weighted continuum phase pass the frozen geometry discriminator?

Status: FROZEN v1 (2026-07-10). Constants and gates in
`docs/prereg/frozen/p5_test_constants.json`, frozen before any confirmatory
chain was run.

## 1. Question

P4 establishes (profile-level) that the smeared-action-weighted 2D-orders
ensemble has a continuum phase at beta < beta_c and crystallizes beyond.
P5 asks the sharper, instrument-level question this program was built for:
do equilibrium samples of the continuum phase PASS the validated
order-intrinsic geometry discriminator -- i.e. can spacetime geometry
actually be reconstructed from them, order-intrinsically, and does the
crystal phase correctly FAIL? To our knowledge phases of causal-set
statistical ensembles have previously been characterized by macroscopic
observables (action, abundances); judging them with a validated
reconstruction instrument is the novel step.

## 2. Instruments

The FROZEN P3 discriminator protocol, verbatim (`p3_dynamics.analyze_order`:
6 disjoint chains >= 25 ticks, up to 44 two-sided bracketed targets, min 20;
order-intrinsic bracket profiles; ordinal embedding; held-out quadruplet
violation; column-shuffle null; truth = distance-order error against true
coordinates). 2D orders supply true lightcone coordinates (t, x) =
(i + pi_i, i - pi_i), so the truth gate is exact. Scale requirement measured:
N = 600 is the smallest scale where uniform 2D orders pass structurally
(6/6 seeds; N = 400 gives 0/6). Calibration at N = 600 (P5-A, seeds 0-9):
10/10 valid, heldout 0.016-0.084, null_gap 0.106-0.322, truth 0.10-0.21.
Gates inherited unchanged from P3: heldout <= 0.10, null_gap >= 0.10,
truth <= 0.40.

## 3. Ensemble and sampler

2D orders at N = 600, smeared action at eps = 0.02 (eps N = 12 as in P4, so
the crystal's exact action advantage S_bip = eps N (2 - eps N) = -120 is
held fixed across scales). Sampler: mcmc_2d_order_fast, validated to
reproduce the reference sampler's trajectory exactly; the reference sampler
is validated against exact Gibbs enumeration. Recon (400k steps, dual
start): beta = 2 and beta = 8 are continuum (a bipartite start MELTS back to
the reference profile), beta = 32 is not (hysteretic, layered drift);
beta_c in (8, 32).

## 4. Frozen design

- Continuum betas {2, 8}: seeds 100-102, random start, 3M steps, burn 0.6,
  one configuration collected every 400k steps -> 3 configs/chain,
  9 configs/beta.
- Crystal control beta = 32: seeds 100-101, bipartite start, 1M steps,
  2 configs/chain.
- Every collected configuration is judged by the frozen discriminator with
  its own true coordinates.

## 5. Frozen expectations

- beta = 2 and beta = 8: "pass" -- >= 75% of the 9 configs clear all three
  gates.
- beta = 32: "block" -- 0 configs clear the gates (structural block counts
  as non-pass).

PASS of the experiment = all three expectations met. Interpretation: the
continuum phase of an action-weighted causal-set ensemble supports genuine
order-intrinsic geometry reconstruction, and the discriminator correctly
rejects the crystal phase. A near-miss pattern (e.g. beta = 8 configs
failing null_gap marginally) would be reported as a quantitative bound on
how far action weighting preserves reconstructable geometry -- itself a
result; no gate will be moved post hoc.

## 6. Deviations log

(empty)

## 7. Confirmatory outcome

Executed 2026-07-10/11 (chains detached from the authoring session; wrap-up
2026-07-12). Registry and per-config CSV under `docs/prereg/frozen/`.
**Experiment PASS -- all three frozen expectations met, unanimously.**

- beta = 2: **9/9 configs pass** all gates. Median heldout 0.031, null_gap
  0.226, truth 0.140. Chain <S> ~ +1.7 to +2.1, acceptance 0.98.
- beta = 8: **9/9 configs pass**. Median heldout 0.039, null_gap 0.239,
  truth 0.133. Chain <S> ~ +1.6, acceptance 0.92.
- beta = 32 (crystal control): **4/4 structural block** ("only 0 chains" --
  the crystal states cannot furnish even one 25-tick chain), 0 gate passes.
  Chain S ~ -106, acceptance 0.006.

The continuum-phase configurations are statistically indistinguishable from
the uniform-ensemble calibration (P5-A median heldout 0.046, truth 0.144):
action weighting up to beta = 8 -- an order of magnitude in beta, reaching
within the decade of the transition at beta_c in (8, 32) -- fully preserves
order-intrinsically reconstructable geometry, with true lightcone position
recovered to truth ~ 0.13-0.14 (chance 0.5).

Conclusion: equilibrium samples of the smeared-action-weighted 2D-orders
ensemble in its continuum phase PASS the validated geometry discriminator --
geometry in this ensemble is not merely profile-level but reconstructable,
and the discriminator correctly rejects the crystal phase. Together with
P3 (dynamics alone: no geometry), E1/de-risk (action alone, unrestricted:
no geometry, exact bipartite obstruction), and P4 (restricted + action:
continuum phase exists, crystallizes at beta_c), this completes the
program's emergence chain with an instrument-level positive result.
