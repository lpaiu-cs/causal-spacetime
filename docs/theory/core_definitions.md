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

In the order-first reformulation, this relation is the primitive temporal
structure. Duration values are representation-layer quantities introduced only
after clocks, density, calibration, or dynamics are supplied.

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

## Measure Encoding

Measure information can be supplied as weights on events or encoded
statistically by the sampling process. In Milestone 10, events are sprinkled
uniformly with respect to conformal physical volume, so the coordinate event
density is proportional to the supplied `Omega^2` profile.

This counting measure is still additional structure. A global constant scale is
not identifiable from normalized event positions without an absolute density or
volume convention.

## Coarse Graining By Thinning

Random thinning keeps each event with probability `p`. If the original density
is `rho`, the expected thinned density is:

```text
rho_thinned = p * rho
```

Density rescaling is part of the reconstruction protocol. Using the original
density after thinning is a deliberate failure case.

## Conformal Ambiguity

Positive conformal rescaling preserves causal order:

```text
ds_Omega^2 = Omega(t, x)^2 ds^2
```

with `Omega > 0`. The light cones are unchanged, but physical volume and
proper-time scale generally change. This is why causal order by itself
determines conformal/light-cone structure rather than the full metric scale.

In 1+1D conformal validation models:

```text
dV_Omega = Omega^2 dt dx
d tau_Omega = Omega dt
```

The conformal factor or equivalent measure is additional structure.

## Observer Chain

An observer chain is an ordered subset of events interpreted as an observer's
worldline. Observer protocols, such as radar decomposition, use this structure
to define operational time and distance assignments.

In discrete radar reconstruction, the observer chain is supplied with strictly
increasing clock labels. These labels are protocol structure, not information
derived from causal order alone.

The order along an observer chain is denoted `≺_O`.

## Observer-Relative Distance Order

An observer-relative distance order compares targets or pairs relative to a
supplied observer protocol:

```text
A ≺^d_O B
(a,b) ≺^d_O (c,d)
```

The first expression means that target `A` is closer than target `B` according
to observer protocol `O`. The second means that pair `a-b` is closer than pair
`c-d` relative to that protocol.

This is weaker than metric distance. It does not by itself provide meters,
ratios, signed coordinates, curvature, or a unique metric tensor.

## Metric Representation Layer

A metric representation layer assigns numerical quantities to order data:

```text
durations, distances, coordinates, metric tensors, curvature values
```

These objects can be extremely useful effective descriptions, but they are not
primitive in the order-first framing. Their stability depends on calibration,
measure/density input, dynamics, and observer-atlas consistency.

## Effective Metric Geometry

Effective metric geometry is a compact representation of regular order
structures. It is appropriate when causal order, observer-relative distance
order, measure information, and transition maps satisfy suitable consistency
and representability conditions.

## Orientation Reference

An orientation reference is additional observer-protocol structure used to
distinguish mirrored spatial positions. In the current 1+1D implementation, it
is represented by a synchronized comoving beacon chain at a known rest-frame
separation from the primary observer chain.

The beacon separation is supplied protocol data. It is not recovered from
causal order alone.

## Observer Atlas

An observer atlas is a collection of observer protocols whose reconstructed
charts overlap on some events. Atlas validation asks whether transition maps on
overlaps are mutually consistent.

In Milestone 7, the controlled transition maps are affine Lorentz/Poincare maps
between oriented inertial charts:

```text
chart_B ~= L(beta_AB) chart_A + translation_AB
```

Observer origins and translations are supplied or fitted protocol-level
structure, not quantities derived from causal order alone.

## Reconstruction-Inaccessibility Boundary

A reconstruction-inaccessibility boundary is a boundary beyond which a supplied
observer protocol cannot assign the relevant operational coordinates. In
Milestone 8, two-way radar accessibility requires both `tau_minus` and
`tau_plus` from the observer chain.

The Rindler horizon in flat 1+1D Minkowski spacetime is used as a controlled
horizon analogue. It is not a black hole model.

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
(E, ≺, counting measure, observer protocol, orientation reference)
  -> reconstructed quantities
```

Examples include proper-time estimates, dimension estimates, radar coordinates,
oriented coordinates, coordinate maps between observer protocols, and
exploratory spatial-distance proxies.
