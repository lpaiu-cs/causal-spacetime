# Core Definitions

This project studies whether spacetime-like quantities can be operationally
reconstructed from primitive causal or information-accessibility structure plus
additional scale and observer information.

## Event Set

Let `E` be a set of events. Events are treated as primitive elements in the
theory-facing formulation. Coordinates, if present in simulations, are hidden
validation data rather than primitive structure.

## Primitive Relation

The primitive accessibility relation is written:

```text
p ≺ q
```

and is read as "`q` is causally or informationally accessible from `p`". In
causal-set simulations this is implemented as a causal order.

## Strict Partial Order

A strict partial order is a relation that is:

- irreflexive: not `p ≺ p`,
- transitive: if `p ≺ q` and `q ≺ r`, then `p ≺ r`.

The simulations use strict causal precedence matrices.

## Local Finiteness

Local finiteness requires every Alexandrov interval to contain finitely many
events:

```text
|I(p, q)| < infinity
```

for all `p ≺ q`.

## Alexandrov Interval

The Alexandrov interval between `p` and `q` is:

```text
I(p, q) = { r in E : p ≺ r ≺ q }
```

Interval cardinality is the main causal-order statistic used for timelike
volume reconstruction.

## Counting Measure And Event Density

The counting measure assigns volume-like information by event cardinality. A
density `rho` converts counts into estimated continuum volumes:

```text
V_hat = count / rho
```

Density is additional scale information. It is not determined by causal order
alone in the current simulations.

## Observer Chain

An observer chain is an ordered subset of events interpreted as an observer's
worldline. Observer protocols, such as radar decomposition, use this structure
to define operational time and distance assignments.

In discrete radar reconstruction, the observer chain is supplied with strictly
increasing clock labels. These labels are protocol structure, not information
derived from causal order alone.

## Null Accessibility

Null accessibility describes boundary accessibility relations analogous to
light-cone structure in a Lorentzian continuum. The code currently validates
null and causal relations in known Minkowski intervals.

## Radar Decomposition

Radar decomposition assigns an event a time and distance relative to an
observer protocol:

```text
radar_time = (tau_plus + tau_minus) / 2
radar_distance = (tau_plus - tau_minus) / 2
```

where `tau_minus` and `tau_plus` are observer emission and reception times.
In the causal-order-based protocol, these are reconstructed as the latest
observer tick preceding the target and the earliest observer tick following the
target.

## Reconstruction Map

A reconstruction map is a procedure:

```text
(E, ≺, counting measure, observer protocol) -> reconstructed quantities
```

Examples include proper-time estimates, dimension estimates, radar coordinates,
and exploratory spatial-distance proxies.
