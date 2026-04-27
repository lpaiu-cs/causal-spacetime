# Observer Protocols

## Why Causal Order Alone Is Not Enough

Causal order can identify which events are accessible from which other events.
It does not, by itself, choose an observer, assign clock readings, or define an
observer-dependent spatial decomposition.

For radar reconstruction, the observer chain and its clock labels are additional
protocol structure.

A single observer chain gives an unsigned distance. It does not define a
left-right orientation or a signed spatial coordinate.

## Observer Chain

An observer chain is a totally ordered sequence of events:

```text
o_0 ≺ o_1 ≺ ... ≺ o_n
```

with clock labels:

```text
tau_0 < tau_1 < ... < tau_n
```

The labels represent an operational clock along the observer chain. In the
current controlled validation experiment, the chain is supplied externally.

## Recovering Radar Ticks From Order

For a target event `e`, the causal-order-based radar protocol finds:

```text
tau_minus = latest tau_k such that o_k ≺ e
tau_plus  = earliest tau_k such that e ≺ o_k
```

This uses only:

- the causal relation matrix,
- observer-chain indices,
- observer clock labels.

The target event's hidden coordinates are not used in the reconstruction.

## Radar Decomposition

When both ticks are available:

```text
radar_time = (tau_plus + tau_minus) / 2
radar_distance = (tau_plus - tau_minus) / 2
```

with natural units `c = 1`.

This is an operational spatial decomposition relative to an observer protocol.
It is not a claim that causal order alone gives spatial distance.

In the order-first reformulation, the first object reconstructed by a radar
protocol can be an observer-relative distance order, such as which target
returns earlier to the observer. Numeric radar distance is a calibrated
representation of that order, not a primitive spatial fact.

## Orientation Reference

An oriented observer protocol can add a synchronized comoving beacon chain at a
known rest-frame separation `a`. Radar distances to the primary chain and beacon
chain can then define a signed coordinate:

```text
x_hat = (R0^2 - Ra^2 + a^2) / (2a)
```

The beacon separation is additional reference structure. Without it, the
single-observer reflection degeneracy remains.

## Validation Scope

Milestone 5 validates observer-chain radar reconstruction in known 1+1D
Minkowski causal diamonds. Hidden coordinates are used only to compare the
reconstructed radar time and radar distance against the analytic stationary
observer values.

The result supports one layer of the mathematical reconstruction program:

```text
observer chain + causal accessibility -> operational radar decomposition
```

It does not prove spacetime emergence.

Milestone 6 extends this validation to two-chain radar reconstruction and
Lorentz-map recovery between oriented inertial protocols.

Milestone 7 uses multiple oriented protocols as an observer atlas. Atlas
validation asks whether reconstructed chart overlaps admit affine
Lorentz/Poincare transition maps, agree on invariant intervals, and compose
consistently around simple loops.

Milestone 8 applies the observer-chain radar protocol to a uniformly
accelerated Rindler observer in flat spacetime. Missing `tau_minus` or
`tau_plus` is treated as observer-protocol reconstruction inaccessibility, with
finite chain coverage reported separately from ideal wedge accessibility.

Milestone 9 clarifies that observer protocols do not remove conformal
ambiguity by themselves. Clock rates and volume scale depend on supplied
measure or conformal-factor information.

Milestone 10 keeps that distinction explicit when measure is encoded in event
statistics. Observer protocols still require clock labels and reference
structure, while count-based volume reconstruction requires a supplied density
scale and appropriate rescaling under thinning.

Milestone 11 adds radar-return order diagnostics that use only causal order and
observer tick order. Numeric clock labels are not needed for that ordinal
comparison, but metric distances and ratios still require calibration.

Milestone 12 then asks whether observer-relative distance-order data can be
embedded into a low-dimensional effective metric representation. The observer
protocol supplies the order data; the embedding is a diagnostic compression of
that data, not a primitive spatial substrate.

Milestone 13 compares observer-derived distance-order constraints against
shuffled and random null baselines. The goal is to check whether the supplied
observer protocol produces structured order data that generalizes and supports
a stable representation better than null-model order. Hidden coordinates remain
validation data, not reconstruction input.

Milestone 14 adds a simultaneity protocol to this observer layer. Radar-time
ranks and bins are reconstructed from tick brackets, then spatial
distance-order comparisons are restricted to same-slice pairs. This keeps
spatial order observer-relative and slice-protocol dependent.
