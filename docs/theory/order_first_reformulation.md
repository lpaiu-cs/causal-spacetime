# Order-First Reformulation

Milestone 11 reframes the project as an order-first reconstruction program.
The point is not that metric geometry is unnecessary. The point is that metric
geometry is treated as a stable, calibrated representation layer for more
primitive order structures.

## Core Thesis

```text
Time is not primitive duration but causal order.
Space is not primitive length but observer-relative distance order.
Metric geometry is a stable calibrated representation of order structures.
```

The fundamental ingredients in this framing are:

- `E`: an event set,
- `≺`: causal or accessibility order on events,
- `O`: an observer protocol or observer chain,
- `≺_O`: order along the observer chain,
- `A ≺^d_O B`: target `A` is closer than target `B` relative to `O`,
- `(a,b) ≺^d_O (c,d)`: pair `a-b` is closer than pair `c-d` relative to `O`.

Seconds, meters, ratios, metric tensors, curvature values, and quantitative
dynamical parameters are not taken as primitive in this theory-facing
formulation. They are representation-layer quantities that may become stable
and useful when order structures satisfy appropriate calibration, consistency,
and representability conditions.

## Why Speed Language Was Rejected

The phrase "time = information-transfer speed" is rejected because speed is
already a ratio of distance and time in ordinary physics. Using it as a
primitive explanation would be circular and would confuse signal speed with
causal order.

The order-first replacement is:

```text
primitive temporal structure = causal/accessibility order
```

Duration values are later representation choices. They require clocks,
calibration, density, repeated processes, or dynamics.

## Primitive Space As Distance Order

The spatial primitive is not a length value. It is an observer-relative ordinal
relation such as:

```text
A ≺^d_O B
```

meaning that `A` is closer than `B` according to an observer protocol. For
pairwise spatial structure, the primitive comparison can be:

```text
(a,b) ≺^d_O (c,d)
```

meaning that the separation of `a` and `b` is smaller than the separation of
`c` and `d` in the observer's reconstructed chart or operational protocol.

This is weaker than metric distance. It does not supply lengths, ratios,
coordinates, curvature, or a unique manifold by itself.

## Metric Geometry As Representation Layer

Metric geometry remains extremely useful. A Lorentzian atlas compresses large
amounts of order data into coordinates, intervals, volumes, and transformation
rules. The order-first thesis is that this geometry is an effective
representation when the order structures are sufficiently regular.

The project therefore asks representation questions:

- Does the causal order behave like a low-dimensional Lorentzian order?
- Does observer-relative distance order admit a useful coordinate
  representation?
- Do overlapping observer protocols have consistent transition maps?
- Are ratios stable under calibration or repeated processes?
- Does thinning preserve reconstructed structure after density rescaling?

## Ratio Stability

Ratios are not primitive in bare order data. A monotone transformation can
preserve all ordinal relations while changing ratios. Ratio stability requires
additional restrictions, such as:

- equal-step clock calibration,
- concatenation rules,
- repeated physical processes,
- supplied density or measure,
- dynamics that selects an affine or metric representation class.

Milestone 11 includes deterministic checks showing that arbitrary monotone
maps preserve order but not equal-step or ratio structure.

## Metric Representability Is Nontrivial

Not every distance order admits a useful low-dimensional metric
representation. A finite distance-order relation may contain cycles, and a
candidate numeric distance matrix may fail triangle inequalities or Euclidean
embedding diagnostics.

The current code only implements lightweight finite diagnostics. These are not
a complete solution to ordinal metric representability. They are sanity checks
for the gap between:

```text
distance order
```

and:

```text
metric geometry
```

## Reinterpreting Milestones 1-10

Milestones 1-3 tested timelike metric reconstruction in known 1+1D Minkowski
intervals. Under the order-first thesis, these are representation-layer checks:
causal interval counts plus supplied density can support a duration
representation.

Milestone 4 tested dimension as an order-statistical observable. This is closer
to the order-first program because it uses causal-order statistics before
assigning full metric quantities.

Milestones 5-6 supplied observer protocols and orientation references. These
can now be read as procedures for deriving observer-relative distance orders
and, with extra calibration, signed coordinate representations.

Milestone 7 tested observer-atlas consistency. In order-first terms, it checks
whether multiple calibrated observer protocols admit compatible affine
Lorentz/Poincare representation maps on overlap.

Milestone 8 tested reconstruction-inaccessibility for a Rindler observer. This
is an order/protocol accessibility boundary in a controlled flat-spacetime
model.

Milestones 9-10 clarified that measure and density are additional structure.
They determine which metric-scale representation is available, and how stable
that representation is under physical-volume sprinkling and thinning.

Milestone 11 adds explicit ordinal diagnostics so that future work can separate
order preservation from metric-value accuracy.

Milestone 12 adds ordinal embedding as a finite diagnostic for effective metric
representation. It asks whether a rich set of distance-order comparisons can be
compressed into a low-dimensional coordinate system, and how that compression
fails under insufficient constraints, wrong dimension, or noisy comparisons.
This strengthens the order-first framing: metric geometry is useful when it is
a stable low-complexity representation of order data, not because it is assumed
as the primitive substrate.

Milestone 13 adds held-out validation and null-model baselines. It asks whether
structured order constraints generalize better than shuffled or random
constraints, whether fitted embeddings are stable across independent subsets,
and whether observer-derived order behaves more like structured geometric data
than like null order. This further separates useful effective representation
from optimizer artifacts.

Milestone 14 adds the missing slice-selection layer for spatial order.
Observer-relative distance order is evaluated within radar-time bins derived
from causal order and observer tick order. This keeps spatial comparison
relative to both an observer protocol and a simultaneity protocol.

Milestone 15 adds the cross-slice identification layer. Same-slice order alone
does not define same position, same direction, velocity, constant velocity, or
spatial evolution. These become transport-relative statements only after a
transport, anchor, persistence, calibration, or dynamics protocol is supplied.
