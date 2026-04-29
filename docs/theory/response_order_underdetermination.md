# Response-Order Underdetermination

Milestone 28 clarifies what a single-reference echo-response signature can
support. A response signature orders targets by their recovered `D_echo` values
relative to one chosen reference protocol. That is a scalar response order.

## Scalar Response Order

For targets `e_i`, a fixed reference chain, and a fixed emission protocol:

```text
D_i = D_echo(e_i; k)
```

The response-order signature records whether `D_i < D_j`, `D_i = D_j`, or the
comparison is unresolved. This is a scalar preorder over target events.

A single reference response signature is not a distance-order structure. It
does not provide target-target pair comparisons, ratios, coordinates, metric
scale, or spatial adjacency.

## Underdetermination Example

The same scalar target order:

```text
e_1 < e_2 < e_3 < e_4
```

can be represented by many synthetic 1D layouts:

```text
[1, 2, 3, 4]
[1, 2, 100, 101]
[1, 10, 11, 100]
```

All preserve the same scalar response order. They can disagree on synthetic
target-target pair ordering. These layouts are representation witnesses only;
they are not physical coordinates.

## Multi-Reference Response Profiles

A multi-reference response profile records several protocol columns for the
same target set:

```text
target x protocol -> D_echo or unresolved
```

This can distinguish targets that tie under one protocol. It is richer than a
single scalar order, but it remains pre-metric. It does not provide pairwise
target comparison protocols by itself.

## Additional Structure

Before ordinal embedding or metric-representation tests are meaningful, the
project must specify what structure supplies target-target pair comparisons or
equivalent constraints. Possible future inputs include explicit pairwise
comparison protocols, multiple calibrated reference protocols, measure or
density assumptions, atlas consistency, and held-out validation against null
baselines.

Milestone 28 therefore treats scalar response-order representability as a weak
precondition. It is not metric reconstruction and not spatial geometry.
