# Horizon Inaccessibility

## Controlled Rindler Analogue

Milestone 8 uses a Rindler observer in flat 1+1D Minkowski spacetime. This is a
controlled horizon analogue for accelerated observers. It is not a black hole
simulation, and it does not use curved spacetime.

For a uniformly accelerated observer, only events in the appropriate Rindler
wedge have two-way radar accessibility in the infinite-duration idealization.
The horizon is therefore treated here as a boundary of observer-dependent
causal reconstruction.

## Two-Way Radar Accessibility

The reconstruction protocol uses:

- the causal matrix,
- observer-chain indices,
- observer proper-time clock labels.

For a target event, two-way radar accessibility requires both:

```text
observer tick before target -> tau_minus
observer tick after target  -> tau_plus
```

If either tick is missing, that event is inaccessible to this radar
reconstruction protocol. This does not mean the event does not exist; it means
the supplied observer protocol cannot reconstruct radar coordinates for it.

## Infinite Wedge Versus Finite Chain

The analytic Rindler wedge describes ideal infinite-duration accessibility. A
finite simulated observer chain has a shorter clock range, so some events
inside the wedge still lack both radar ticks.

Milestone 8 reports both notions:

- infinite-wedge accessibility,
- finite-chain coverage accessibility.

This separation prevents finite protocol duration from being mistaken for the
true Rindler horizon boundary.

## Interpretation

The experiment tests whether a reconstruction-inaccessibility boundary can be
detected operationally from causal order and a supplied observer clock protocol
in a known flat spacetime.

Order-first interpretation: the boundary is first a limitation of
observer-relative accessibility order. Rindler radar coordinates are a
representation of the accessible region, not primitive evidence that metric
geometry has been derived.

It does not prove spacetime emergence. It does not prove black hole physics. It
does not show that causal order alone derives horizons.

Milestone 9 separately addresses conformal ambiguity. The Rindler causal wedge
is a light-cone accessibility statement; physical volume and clock scale still
depend on measure or conformal-factor information.

Milestone 10 studies how such measure information can be represented by
event-density statistics and how reconstruction behaves under thinning. This
is separate from the Rindler accessibility boundary: thinning changes sampling
density, not the ideal causal wedge.

Milestone 11 keeps the same distinction: accessibility order can be studied
before assigning metric distances, ratios, or curvature-like representation
values.

Milestone 12 focuses on embeddability of distance-order data rather than
horizon accessibility. Both remain representation-layer diagnostics in
controlled models, not derivations of metric geometry as fundamental.

Milestone 13 adds held-out and null-model validation for effective metric
representations. It is orthogonal to the Rindler horizon experiment: horizon
accessibility tests observer-dependent causal access, while representation
stability tests whether distance-order data supports a robust low-dimensional
metric compression.

Milestone 14 adds radar-time slice selection for observer-relative spatial
order. In horizon-limited settings, the same distinction would matter: first
identify which events are accessible, then define observer-relative slices
inside the accessible region before making spatial comparisons.

Milestone 15 adds a separate cross-slice identification issue. Even inside an
accessible Rindler region, same-position, same-direction, velocity, and spatial
evolution statements require transport, anchors, persistence, calibration, or
dynamics. Horizon accessibility does not by itself supply those structures.

Milestone 16's relational histories could be applied inside an accessible
region if persistent object labels and slice-local distance order are supplied.
This would remain a relational diagnostic, not a horizon dynamics model.

Milestone 17 clarifies that those labels are another layer. Horizon
accessibility, slice-local distance order, and persistence matching are
separate questions. A Rindler accessibility boundary does not by itself define
object identity across slices.
