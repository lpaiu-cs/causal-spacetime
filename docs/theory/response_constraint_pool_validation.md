# Response Constraint Pool Validation

Milestone 30 validates response-comparison constraint pools before any future
representability experiment. A pool is a finite set of quadruplets:

```text
(i, j, k, l)
```

meaning that target pair `(i, j)` has smaller response-profile dissimilarity
than target pair `(k, l)` under a declared admissible comparison protocol.
This is a response-comparison relation, not a spatial-distance relation.

## High-Margin Constraints

A high-margin constraint requires the protocol dissimilarity gap to exceed a
chosen threshold. Larger margins usually reduce the number of available
constraints. They can also improve held-out agreement by avoiding near-ties.

## Held-Out Protocol Agreement

Held-out validation builds a pool from one subset of protocol columns and
evaluates it on another subset. Agreement, inversion, and unresolved fractions
are recorded. This checks protocol-column stability without fitting any
embedding.

## Bootstrap Stability

Bootstrap validation resamples protocol columns and evaluates the fixed pool.
The stable constraint fraction records how many constraints agree in enough
bootstrap samples. This is a finite-sample robustness diagnostic.

## Null-Baseline Separation

Null validation compares the structured profile with shuffled or
marginal-preserving profiles. A useful pool should outperform simple nulls,
but null-baseline separation remains a pre-metric diagnostic. It does not show
that the constraints are spatial distances.

## Coverage Diagnostics

Target coverage reports how many target events appear in the pool. Pair-node
coverage treats unordered target pairs as nodes and reports how many of those
nodes participate. Sparse coverage is a warning that later representation
tests would be underdetermined by the selected pool.

## Validation Gates

A validation gate combines minimum constraint count, evaluable fraction,
agreement fraction, inversion threshold, unresolved threshold, null-separation
threshold, and bootstrap-confidence threshold. Passing validation gates makes
a pool eligible for future representability experiments only under explicitly
stated assumptions. Passing gates does not make constraints spatial distances.

This milestone does not perform ordinal embedding and does not validate
distance-order constraints.

Milestone 31 exports handoff manifests for pools that satisfy predeclared
criteria. These manifests are frozen input specifications for future
experiments. They do not contain fitted embeddings and do not imply distance
or geometry.

Milestone 32 is the first fitting consumer of those manifests. It fits latent
ordinal representation models from the frozen train split and evaluates the
frozen held-out split. This does not alter the validation gates, and it does
not infer spatial distance or metric geometry.
