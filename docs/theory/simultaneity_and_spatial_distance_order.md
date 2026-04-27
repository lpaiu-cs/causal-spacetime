# Simultaneity And Spatial Distance Order

Milestone 14 refines the order-first spatial layer. Spatial distance order is
not treated as absolute. It is observer-relative and depends on a
slice-selection protocol.

The core distinction is:

```text
causal order -> primitive temporal order
observer protocol + causal order -> radar tick brackets
radar tick brackets -> radar-time order / radar-time bins
same radar-time slice -> spatial distance-order comparisons
```

## Why Slices Matter

Pairwise distance comparisons over all accessible events can mix spatial
separation with time separation. In a relativistic setting, "spatial distance"
is not an observer-independent relation across arbitrary events. It is defined
relative to an observer protocol and a simultaneity or slice-selection rule.

For this project, that rule is built from observer radar data rather than from
hidden coordinates.

## Radar-Time Bins

For each target event, the causal-order protocol finds:

```text
predecessor tick = latest observer tick preceding the target
successor tick   = earliest observer tick succeeding the target
```

The order-level radar-time rank is:

```text
predecessor_tick_position + successor_tick_position
```

The order-level radar-distance rank is:

```text
successor_tick_position - predecessor_tick_position
```

These ranks use only the causal matrix and observer tick order. They do not use
numeric clock labels or hidden target coordinates. Radar-time bins group nearby
radar-time ranks into observer-relative slices.

## Same-Slice Spatial Comparison

Within a slice, pairwise comparisons can be interpreted as observer-relative
spatial distance-order comparisons. This is still not an absolute distance and
not a metric by itself. It is an order relation that may later admit an
effective metric representation if it is stable, calibrated, and consistent
with other observer protocols.

## Interpretation

Milestone 14 tests whether observer-derived distance order becomes cleaner when
restricted to same-slice comparisons. This supports the order-first thesis by
separating temporal order, slice selection, spatial distance order, and metric
representation.

The milestone does not prove spacetime emergence, derive space, make spatial
distance observer-independent, or show that causal order alone gives metric
geometry.

## Cross-Slice Limitation

Milestone 15 adds the next separation. Same-slice distance order is not a
global spatial coordinate through time. Without a transport, anchor,
persistence, or calibration rule, cross-slice predicates such as same position,
same direction, velocity, constant velocity, and spatial evolution are
undefined. With a chosen rule, they become transport-relative statements.

Milestone 16 shows that persistence plus same-slice distance order can still
define weaker relational histories. The comparison is not "same position across
slices"; it is whether persistent object-pair distance-order relations change
from one slice to another.

Milestone 17 asks what happens when persistence is not supplied. Same-slice
distance order can constrain identity matchings, but it does not determine
them uniquely in general. Relational histories built from inferred identities
are persistence-hypothesis dependent.
