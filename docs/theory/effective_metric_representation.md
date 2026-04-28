# Effective Metric Representation

Milestone 12 studies when order data admits a useful low-dimensional metric
representation. This is an order-first question:

Milestone 18 places this question inside the state-change causal order
hierarchy. The metric representation layer attempts to compress causal trigger
order, observer-chain order, and observer-slice distance order after
calibration. It is not part of the primitive state-changing event definition.

```text
Given order relations, when is a metric representation a stable low-complexity
compression of those relations?
```

The metric representation is not treated as fundamental. It is an effective
description that becomes useful when causal order, observer-relative distance
order, calibration, measure information, and atlas consistency satisfy suitable
conditions.

## Order Preservation Is Weaker Than Metric Accuracy

Distance order can say:

```text
d(a,b) < d(c,d)
```

without assigning a meter value or a ratio. A representation may preserve many
such comparisons while still having substantial coordinate RMSE. Conversely, a
small metric error is a stronger statement than preserving ordinal comparisons.

Milestone 12 therefore reports order-preservation loss, hinge loss, violation
rates, and Procrustes-aligned RMSE as separate diagnostics.

## Monotone Transformations And Ratios

Positive monotone transformations preserve order. They need not preserve
ratios. This is why ratios are not primitive in the order-first framing.

Calibration-stabilized ratios require additional restrictions, such as:

- equal-step clock constraints,
- repeated calibrated processes,
- concatenation rules,
- supplied density or measure,
- dynamics selecting an affine or metric representation class.

## Ordinal Embedding

Ordinal embedding tries to find coordinates whose pairwise distances satisfy
quadruplet constraints:

```text
d(i,j) < d(k,l)
```

The implementation in this repository is a simple finite diagnostic solver
using a quadruplet hinge loss. It is not a production optimizer and not a
representation theorem.

## Low-Dimensional Representation Is Nontrivial

Not every distance order admits a useful low-dimensional metric representation.
Some constraint sets contain cycles or noise. Some candidate numeric distance
matrices violate triangle inequalities or fail Euclidean embedding diagnostics.

Low-dimensional metric geometry is therefore a representability condition, not
an automatic consequence of distance order alone.

## Relation To Observer Protocols

Observer-relative distance order can be produced by radar protocols and
oriented charts. Milestone 12 tests whether those order relations support a 1D
effective spatial embedding in controlled 1+1D validation data.

The hidden coordinates are used only for controlled validation. The embedding
uses order constraints generated from reconstructed observer-relative distance
data.

## Lorentzian Atlas Is Stronger

A single ordinal embedding is weaker than a Lorentzian atlas. An atlas requires
multiple observer protocols, overlap, transition maps, invariant consistency,
and calibration across charts. Metric representation of one distance order is
only one layer of that broader structure.

## Interpretation

These experiments support the order-first reconstruction program by clarifying
when metric geometry can act as a low-complexity compression of order data.
They do not show that distance order alone gives full metric geometry, derive
meters, prove spacetime emergence, or imply that all order structures are
embeddable.

## Held-Out And Null-Model Stability

Milestone 13 adds the next diagnostic layer. A finite ordinal embedding is more
informative when it:

- preserves held-out order constraints not used during fitting,
- outperforms shuffled-side and random null constraints,
- remains stable across independent subsets of the same structured constraints,
- degrades in controlled ways under flipped or noisy comparisons,
- shows clearer low-dimensional compression for structured order than for null
  order.

These tests guard against reading too much into optimizer training loss. They
still remain finite diagnostics, not proofs that metric geometry is fundamental
or that all observer-relative distance orders are embeddable.

## Simultaneity-Sliced Order

Milestone 14 adds a domain restriction for observer-derived spatial order:
distance comparisons should be made within an observer-selected radar-time
slice. This avoids treating time-separated events as if their pair separations
were purely spatial. Same-slice distance order is still weaker than metric
distance and remains dependent on the observer and slice protocol.

Milestone 15 clarifies why a single global coordinate across time is an
additional representation choice. Slice-local embeddings can preserve
same-slice distance order while independent per-slice gauges change global
alignment, same-position judgments, and velocity estimates. Transport and
anchors are therefore part of the effective representation layer.

Milestone 16 adds a metric-free relational diagnostic. Ordinal shape histories
can remain invariant under slice gauges even when coordinate velocities change.
A metric representation is still required for quantitative dynamics, but some
weaker relational evolution can be expressed at the order level.

Milestone 17 studies whether persistent labels are available, partial, or
inferred. Identity matching is another representation-layer hypothesis. A
low-dimensional metric representation of a relational history is meaningful
only relative to the persistence assumptions used to assemble that history.
