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

It does not prove spacetime emergence. It does not prove black hole physics. It
does not show that causal order alone derives horizons.

Milestone 9 separately addresses conformal ambiguity. The Rindler causal wedge
is a light-cone accessibility statement; physical volume and clock scale still
depend on measure or conformal-factor information.
