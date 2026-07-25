# P6b: Same-data diagnostic comparison

Status: **FROZEN v1** (2026-07-14), before aggregate comparison. Constants are
in `docs/prereg/frozen/p6b_test_constants.json`.

## Question and inputs

What does the frozen order-intrinsic reconstruction discriminator add over
Myrheim-Meyer dimension, interval abundances, and order height on identical
orders? Inputs are the frozen P1 epsilon sweep, P3 transitive-percolation
orders, P5 continuum/crystal configurations, and P6 layered hard negatives.
The deterministic generators and recorded seeds reproduce every order. P5
MCMC regeneration is analysis-only and does not create a new sample.

## Labels

Primary ROC labels are fixed as follows. Geometric: P1 epsilon <= 0.15 and P5
beta in {2, 8}. Non-geometric: P1 epsilon >= 0.90, every P3 cell, P5 beta=32,
and every confirmed P6 layered cell. Intermediate P1 cells are report-only for
ROC and enter the rank and H-LAG analyses.

## Diagnostics

- Frozen instrument margin: the minimum normalized distance inside its
  applicable P1 or P3/P5/P6 gates. Structural blocks receive -5.
- MM distance: absolute standardized distance from the matching geometric
  reference's MM dimension.
- Interval-abundance distance: RMS standardized distance in
  `log1p(n0), log1p(n1), log1p(n2)`.
- Height distance: absolute standardized distance in `log1p(height)`.

Because P1 contains supplied chains, P3 has N=1500, and P5/P6 have N=600,
normalization is performed separately against deterministic geometric
references: P1 pristine scenes seeds 0-19, sprinkled 1+1D orders seeds 0-19,
and uniform N=600 2D orders seeds 1000-1019. This normalization is fixed
before confirmatory metrics and is not fitted to class labels.

## Metrics

1. ROC/AUC on the frozen endpoint/class labels. Higher instrument margin and
   lower cheap-diagnostic distance mean more geometric.
2. Per-seed P1 Spearman correlation of each diagnostic with epsilon, and with
   frozen truth-order error; report the median and range over covered seeds.
3. H-LAG overlap: among P1 cells where heldout <= 0.05 but truth error > 0.15,
   report the fraction each cheap diagnostic calls geometric. A cheap
   diagnostic passes when its distance is no larger than the 95th percentile
   of its own reference-distance distribution.

ROC/AUC uses average ranks for ties. Spearman uses average ranks. No optional
spacelike-distance or homology diagnostic enters the primary comparison.

## Execution policy

Raw regenerated observables are cached by source. P5 is sharded by beta and
seed because exact deterministic regeneration of the frozen 3M-step chains is
long-running. Aggregation refuses to run unless all frozen source cells and
all reference groups are present. No missing cell is silently dropped.

The optional Numba replay is an implementation accelerator, not a new chain:
it pre-generates the exact conditional RNG stream used by the validated NumPy
sampler and applies the same accepted/rejected permutation swaps. An automated
equivalence test requires acceptance rates and every sampled permutation to
match the validated sampler exactly; each replayed P5 shard additionally must
reproduce its frozen mean action before it is accepted.

## Deviations log

(empty)

## Outcome

Executed on the complete frozen input set: 442 orders, of which 58 geometric
and 304 non-geometric orders enter the primary ROC comparison; 80 intermediate
P1 cells are rank/overlap only. All eight replayed P5 chains reproduced their
frozen mean actions before their 22 configurations were admitted.

| Diagnostic | ROC AUC | Median P1 rho vs epsilon | H-LAG false-pass overlap |
|---|---:|---:|---:|
| Frozen instrument margin | **0.993** | 0.976 | n/a |
| MM dimension distance | 0.939 | 0.012 | 25/27 (0.926) |
| Interval-abundance distance | 0.933 | 0.786 | 0/27 (0.000) |
| Height distance | 0.967 | **0.994** | 0/27 (0.000) |

The discriminator has the strongest overall class separation, but the result
does **not** support claiming that cheap diagnostics are generally inadequate.
Normalized height is highly competitive and tracks the P1 dose response even
slightly more monotonically; abundance and height both reject all frozen P1
H-LAG cells under their preregistered reference cutoffs. MM dimension is the
clear weak diagnostic for graded degradation and false-passes 25/27 H-LAG
cells, consistent with the known crystalline-MM degeneracy.

The added value of the instrument is therefore narrower and operational: it
directly establishes whether observer-derived profiles support held-out
ordinal reconstruction and, where truth labels exist, whether the recovered
order is correct. It improves mixed-class AUC by 0.027 over height, while
height remains an effective low-cost screening observable. Paper B should
present both facts and avoid a blanket superiority claim.

## Post-analysis correction (2026-07-17; frozen files unchanged)

The frozen P6b proxy was reproduced exactly, but its P1 formula in
`p6b_test_constants.json` omitted the preregistered restart-stability gate.
P1's declared decision is the conjunction of held-out violation, truth-order
error, and restart-order disagreement. This is a specification inconsistency,
not a discrepancy in the executed rows, so the historical constants, scored
rows, and summary remain unchanged.

A deterministic gate-complete correction joins
`p6b_scored_rows.csv` to `p1_stage_b_epsilon_sweep.csv` on `(seed, epsilon)`.
The dated correction, source hashes, and exact values are in
`docs/paper/paper_b/figures/p6b_margin_correction.json`. The corrected full
margin has AUC 0.9898. Removing only the truth-coordinate term while retaining
all order-only gates gives AUC 0.9680, effectively tied in point estimate with
height at 0.9668. Adding stability changes no labelled pass/block decision but
does change ranks. Accordingly, the original 0.9934 value is reported only as
the historical frozen proxy, and the earlier claim that the aggregate ranking
does not depend on truth assistance is withdrawn.

In the 27 P1 H-LAG cells, the corrected order-only margin is nonnegative for
26, while the full margin's truth gate defines the window. This is direct
evidence that truth recovery is load-bearing rather than a dispensable
evaluation-only term.

The aggregate also mixes generator families and reference normalizations.
Within-family and family-pair comparisons can reverse the apparent
order-only/height ordering, while repeated configurations from a chain or seed
are correlated. These results are descriptive point estimates, not an
uncertainty-qualified proof of general diagnostic superiority.
