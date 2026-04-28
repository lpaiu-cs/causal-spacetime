# Cross-Slice Identification And Transport

Milestone 15 clarifies a missing layer in the order-first program. Same-slice
distance order gives slice-local spatial structure. It does not by itself
identify points, directions, or objects across different observer-time slices.

Milestone 18 keeps transport in the effective/protocol layer. The primitive
state-change trigger order does not define same position, same direction,
velocity, or metric spatial evolution across observer slices.

The corrected interpretation is:

```text
observer-relative, slice-relative distance order
  + chosen transport / anchors / persistence / calibration
  -> transport-relative cross-slice statements
```

Without such additional structure, questions such as "is this the same
position?", "is this the same direction?", "is this object's velocity
constant?", or "has space changed over time?" are not false. They are
undefined.

## Slice-Local Structure

Milestone 14 restricts spatial comparisons to observer-defined radar-time
bins. This avoids treating all accessible events as if they shared one
simultaneity slice. Within each bin, distance-order comparisons can be tested
as observer-relative spatial order.

That local order still has gauge freedom. Each slice can be independently
translated, reflected, and rescaled without changing the ordering of distances
within that slice. Therefore a global spatial coordinate across time is
underdetermined by same-slice distance order alone.

## Transport Is Additional Structure

Transport is a protocol-dependent identification rule. It is not a restoration
of absolute space and is not derived from causal order alone in these
experiments.

Examples of additional structures that can restrict transport choices include:

- supplied anchor worldlines,
- persistent object identities,
- beacon chains or orientation references,
- gyroscope-like direction references,
- matter configurations,
- repeated calibration processes,
- dynamical laws.

Once such a rule is chosen, same-position, same-direction, velocity, constant
velocity, and spatial evolution become transport-relative statements.

## Anchors And Persistence

Milestone 15 includes controlled anchor and persistence diagnostics. Synthetic
anchors are fixed-position events supplied for each slice. Persistent objects
are supplied identity tracks across slices. They are protocol or dynamics
inputs, not consequences of bare causal order.

Anchor-constrained transport can improve global coordinate alignment relative
to no-transport or random per-slice gauge baselines in controlled validations.
Noisy anchors degrade derived cross-slice quantities. This is expected: those
quantities depend on the supplied transport information.

## What This Does Not Show

These diagnostics do not prove spacetime emergence, derive object identity, or
make velocity primitive. They show where cross-slice statements enter the
representation layer and how their conclusions depend on additional protocol,
anchor, persistence, or calibration structure.

## Relation To Relational Invariants

Milestone 16 adds a weaker class of cross-slice statements. If object
persistence is supplied, pair-distance order histories can be compared across
slices without same-position transport. These histories are transport-gauge
invariant under per-slice affine/reflection transformations, but they still do
not define velocity or absolute position.

Milestone 17 makes the persistence assumption explicit. Object persistence is
not derived from causal order alone here. When labels are absent, identity
matching is a supplied or inferred hypothesis; symmetric, noisy, or crossing
histories can leave multiple compatible matchings.
