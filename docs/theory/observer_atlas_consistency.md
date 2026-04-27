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

In the order-first reformulation, an atlas is an effective metric
representation of compatible order data. Transition maps and invariant
intervals are not primitive; they are consistency conditions indicating that
multiple observer-relative order structures admit a shared low-dimensional
Lorentzian representation.

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

Milestone 11 adds explicit representability diagnostics. These emphasize that
not every distance order or chart overlap can be assumed to admit a useful
metric representation.

Milestone 12 implements ordinal embedding diagnostics for that
representability question. A successful single embedding is still weaker than a
Lorentzian atlas, because an atlas also requires calibrated overlap maps and
composition consistency.

Milestone 13 adds stability checks for the embedding layer: held-out
constraints, null-model baselines, and independent subset comparisons. These
diagnostics help identify whether a chart-like metric representation is a
stable compression of structured order data rather than an optimizer artifact.

Milestone 14 emphasizes that spatial distance order inside each chart requires
a slice-selection rule. Atlas consistency is stronger than same-slice order,
but the same conceptual separation applies: observer protocol, simultaneity
protocol, and metric transition maps are distinct representation layers.

Milestone 15 adds transport as another atlas-level assumption. Comparing
spatial structures across chart time slices requires cross-slice
identification. Different transport choices can preserve slice-local order
while changing global alignment and velocity judgments.

Milestone 16 identifies atlas-independent relational content inside a supplied
persistent-object history. Pair-distance order histories are weaker than chart
transition maps, but they can test transport-gauge invariant relational change.

Milestone 17 makes the supplied-history assumption explicit. If object labels
are absent, atlas-independent relational histories require a persistence
hypothesis. Symmetry or crossing can leave several matchings compatible with
the same slice-local order data.
