# Paper B strengthening levers

Status: **revision decision memo, 2026-07-17**.

## Decision

The best next scientific investment is an **independent, family-held-out,
diagnostic-matched benchmark whose generator cells are selected without
consulting the instrument margin**. It directly addresses the largest remaining
validity gap in P6b. More local P7 sampling and a headline KPZ-universality
claim are not efficient Paper B critical-path work.

## Ranked levers

| Rank | Lever | Type | Payoff | Cost | Recommendation |
| ---: | --- | --- | --- | --- | --- |
| 0 | Correct P6b, integrate P7, narrow T1 claims | manuscript + deterministic audit | removes validity errors | low | completed in v0.7 |
| 1 | Family-stratified, cluster-level P6b sensitivity | existing-data post-hoc analysis | exposes heterogeneity and sampling uncertainty | low | do before submission |
| 2 | Independent diagnostic-matched hard-negative holdout | new preregistered experiment | resolves P6b selection bias and external-validity gap | medium | best new experiment |
| 3 | P5 low-beta dual-start replication | new preregistered experiment | permits a stronger equilibrium interpretation at beta 2/8 | medium-high | only if “equilibrium” stays central |
| 4 | 2+1D dilution dose-response | new preregistered experiment | tests whether the core H-LAG result generalizes beyond 1+1D | medium-high | strong follow-up |
| 5 | Observer-congruence holdout | new robustness experiment | tests dependence on observer placement/layout | medium | useful secondary lever |
| 6 | Heavier manifoldlikeness comparators | analysis + possible new computation | answers reviewer baseline questions | medium | pair with lever 2 |
| 7 | T1 wandering-class study | separate theory/numerics project | could support a real universality claim | high | remove from Paper B critical path |
| 8 | P7 equilibrium FSS | currently blocked program | very high if solved | prohibitive | do not resume without a new method |

## 1. Existing-data sensitivity analysis

The pooled corrected order-only AUC is 0.96796 versus height 0.96676, but
family-pair orderings reverse. P1 positives against P6 negatives favor height
(0.938 versus 0.913), while P5 positives against P6 negatives favor the
order-only margin (0.997 versus 0.968). P6a cells were also selected at Stage A
partly for instrument blocking. A single pooled row-level AUC therefore hides
both heterogeneity and selection.

Freeze an analysis plan before adding more summaries:

- report each source and normalization group separately;
- report leave-one-family-out and equal-family-weight macro AUC;
- use paired `delta AUC = instrument - height`;
- bootstrap independent units, not rows: seed for P1/P3/P6 and MCMC chain for
  P5;
- treat the three P5 configurations from a chain as correlated;
- include a simple multivariate cheap baseline combining height, MM dimension,
  and abundance;
- label all results post-hoc sensitivity analyses, not confirmatory evidence.

Success here is not necessarily a positive delta. A precise near-tie would
support the narrower and useful conclusion already adopted in v0.7: height is
a strong cheap screen, while the instrument supplies a validated
reconstructability verdict.

## 2. Independent diagnostic-matched holdout

This is the highest-value new experiment because it removes the
instrument-conditioned family selection in P6a/P6b.

### Design constraints

- Choose generator families and parameter cells using chain operability and
  cheap diagnostics only; never inspect instrument margins during selection.
- Match negative cells to geometric references on N, relation density, height,
  MM dimension, and interval-abundance profile.
- Freeze generator code, cell list, seeds, evaluator hash, row exclusions, and
  stopping rule before confirmatory execution.
- Separate generator-development seeds, calibration seeds, and confirmatory
  seeds; hold out at least one entire generator family.
- Include fresh geometric positives in every reference group, so a detector
  cannot win merely by rejecting everything.
- Make the gate-complete order-only margin primary. Report truth-assisted
  margin as secondary.
- Compare against each scalar and a frozen multivariate cheap baseline.
- Make family-balanced, cluster-level paired delta AUC the primary statistic.

### Decision criterion

A defensible superiority claim would require a positive lower confidence bound
for order-only delta AUC against both height and the multivariate cheap
baseline on held-out families. If that bar is missed, retain the
protocol-meaning claim and drop ranking superiority entirely.

## 3. Low-beta equilibrium replication

P5 currently supports configuration-level verdicts: three random-start chains
per beta, three stored post-burn-in configurations per chain, all passing at
beta 2 and 8. It does not certify equilibrium.

If equilibrium language is important, preregister a separate qualification:

- random and bipartite starts at beta 2 and 8;
- chains, not saved configurations, as primary units;
- pre-frozen start-agreement, IAT/ESS, and trace-stationarity criteria for
  action, n0, height, and geometry score;
- sampler qualification decided before geometry outcomes are opened;
- enough independent chains to estimate between-chain uncertainty;
- trace-level or sufficiently granular archived artifacts;
- unchanged reconstruction gates.

Without this replication, “random-start post-burn-in configurations” is the
correct and sufficient description.

## 4. Generality levers

### 2+1D dilution

Port P1's held-density geometry dilution to the remediated P2-v2 scene and
freeze dimension selection, truth recovery, and H-LAG crossing behavior. This
tests the paper's central methodological result rather than merely another
endpoint. It is a stronger generality lever than adding more 1+1D seeds.

### Observer congruence

Hold out observer layouts: positions, unequal spacings, boosted/mixed chains,
and target hull coverage. Keep generator and gates fixed. The existing
selector probe varies chain extraction but does not establish robustness to
the observer congruence itself.

### Heavier comparators

Spacelike-distance and homology diagnostics are useful reviewer-response
baselines, but should be evaluated on the independent holdout. Adding them only
to the instrument-conditioned P6b pool cannot cure its selection bias.

## 5. Theory lever

The order-only transverse exponent -0.168 is a seven-density point estimate
numerically consistent with -1/6. A real wandering-class result requires a
separate design with multiple system sizes, anchor separations, independent
replicates, uncertainty on slopes, finite-size corrections, competing exponent
models, and preferably distributional or scaling-collapse checks. That is a
substantial standalone theory/numerics project. Paper B is stronger when it
uses only the current protocol-dependence and wandering-dominance results.

## 6. Blocked P7 lever

P7 cannot currently supply an N=600 equilibrium transition or FSS:

- local Metropolis does not cross the beta>=16 basin barrier;
- the replica-ladder arithmetic is prohibitive;
- measured Wang-Landau tunneling cost grows too rapidly;
- the sampler-feasible small-N window and instrument-operable window do not
  overlap.

Reopen P7 only after either a qualitatively different sampler demonstrably
tunnels at N>=500, or a new small-N geometry observable passes its own positive
control. The latter would be a new instrument and a related, not identical,
research question.

## Submission sequence

1. Keep the v0.7 correction and claim-boundary changes.
2. Run the existing-data family/cluster sensitivity analysis.
3. Complete the Paper A evidence package and rebuild both manuscripts.
4. If one new experiment is affordable, run the independent matched holdout.
5. Add low-beta dual-start replication only if equilibrium language is worth
   the cost.
6. Move KPZ-universality and P7 FSS to follow-up projects.
