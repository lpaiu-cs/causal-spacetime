# Representation Stability And Null Models

Milestone 13 strengthens the ordinal embedding layer by asking whether a fitted
low-dimensional metric representation is a stable compression of structured
order data rather than an artifact of a finite optimizer.

This remains an order-first diagnostic. The primitive data are order relations,
not meters, ratios, or a metric tensor. A metric embedding is treated as a
candidate effective representation when it preserves held-out order constraints,
is stable across independent constraint subsets, and outperforms simple null
models.

## Why Milestone 12 Was Not Enough

Milestone 12 showed that ordinal embeddings can be fitted to controlled
distance-order constraints. A low training loss alone is not enough. A flexible
optimizer may partially fit inconsistent or random constraints, especially when
there are many adjustable coordinates.

Milestone 13 therefore adds:

- held-out order validation,
- shuffled-side and fully random null-model baselines,
- controlled noisy constraints,
- stability across independent subsets of the same order data,
- observer-derived order comparisons against null baselines,
- complexity curves across candidate embedding dimensions.

## Held-Out Order Validation

Held-out validation separates fitting from testing. The embedding is fitted on
training order constraints and evaluated on constraints it did not see during
optimization. Structured geometric order data should generalize better than
shuffled or random constraints if a low-dimensional representation captures real
regularity in the order data.

This is still a finite diagnostic. Good held-out performance is evidence of
stable representation quality in the controlled model, not a representation
theorem.

## Null-Model Baselines

Null models are baselines, not physical alternatives. They answer a practical
question:

```text
Would the same optimizer appear successful on order data whose direction has
been shuffled or randomized?
```

Random constraints can sometimes be partially fit, but they should not show the
same held-out performance, stability, or low-dimensional compression behavior
as structured geometric constraints.

## Embedding Stability

If the same structured order data is sampled in independent subsets, the fitted
embeddings should become more similar as the number of constraints grows. The
project checks this using Procrustes-aligned coordinate stability and pairwise
distance-order disagreement between embeddings.

Stability matters because a useful effective metric representation should not
depend sensitively on one arbitrary finite sample of comparisons.

## Observer-Derived Distance Order

Observer-derived distance order is generated from reconstructed oriented radar
coordinates, using hidden coordinates only for validation. Milestone 13 compares
these constraints against shuffled and random baselines. The question is whether
the observer protocol produces order data that behaves more like structured
geometric data than like null-model data.

## Complexity Curves

A low-dimensional metric representation is useful only if it compresses order
data. The complexity curve compares held-out violation against candidate
embedding dimension and a simple dimension penalty. This heuristic is not an
intrinsic-dimension theorem; it is a finite check for low-complexity
representability.

## Interpretation

These experiments support the idea of metric geometry as an effective
calibrated representation when structured order data has stable low-dimensional
regularity. They do not show that distance order alone gives full metric
geometry, derive space, prove spacetime emergence, or make the optimizer a
physical principle.

## Simultaneity-Sliced Baselines

Milestone 14 refines the observer-derived baseline by generating distance-order
constraints only within radar-time bins. This tests whether weaker performance
in all-event observer-derived constraints was partly due to mixing events from
different observer-time slices. Null-model comparisons remain important: a
sliced observer-derived order should still be compared against shuffled and
random constraints before treating it as a useful representation source.
