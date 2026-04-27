# Observer Atlas Consistency

## Why Atlases Matter

An isolated observer protocol can assign operational coordinates to accessible
events, but a spacetime-like reconstruction needs more structure than one
chart. Overlapping reconstructed charts should have transition maps that are
mutually consistent on shared events.

Milestone 7 tests a small observer atlas in controlled 1+1D Minkowski
intervals. The atlas consists of multiple oriented observer protocols applied
to the same causal order.

## Supplied Protocol Structure

Each protocol includes:

- a primary observer chain,
- a synchronized beacon chain,
- clock labels along both chains,
- a supplied beacon separation,
- a supplied affine origin event in the lab validation frame,
- an inertial boost parameter used to construct the controlled protocol chains.

These are protocol inputs. They are not derived from causal order alone.

## Transition Maps

In the flat validation setting, the expected transition from chart `A` to chart
`B` is an affine Lorentz/Poincare map:

```text
chart_B ~= L(beta_AB) chart_A + translation_AB
```

The reconstruction experiment fits this map from overlapping reconstructed
chart coordinates. Hidden Minkowski coordinates are used only to construct the
controlled protocol chains and to validate expected values.

## Invariant Interval Agreement

For pairs of events visible in two charts, the experiment compares:

```text
s^2 = dt^2 - dx^2
```

computed from each reconstructed chart. Agreement of this invariant is a
chart-consistency diagnostic. It is not a proof that a manifold has been
derived.

## Loop Consistency

Atlas coherence also requires transition maps to compose consistently. The
experiment compares the fitted direct map `A -> C` with the composed map:

```text
A -> B -> C
```

and records beta and translation composition errors.

## Interpretation

Observer-atlas consistency is a stronger validation target than a single
observer reconstruction. It checks overlapping reconstructed charts, affine
transition maps, invariant interval agreement, and loop consistency in a known
flat model.

The result remains conservative. It does not prove spacetime emergence, it does
not show that causal order alone gives a manifold, and it does not establish
uniqueness of the atlas.

Milestone 8 studies a different consistency layer: horizon-limited
reconstruction for one accelerated observer protocol. It uses Rindler
accessibility as a controlled flat-spacetime analogue rather than adding curved
spacetime dynamics.

Milestone 9 adds a limitation that applies to any atlas layer: transition maps
and causal order do not by themselves fix conformal volume scale. A measure or
density field remains additional structure.

Milestone 10 extends that limitation into a stability test. If measure is
encoded by event density, local relative measure shape can be estimated
statistically, but absolute scale and thinning rescaling remain explicit
additional inputs to the reconstruction protocol.
