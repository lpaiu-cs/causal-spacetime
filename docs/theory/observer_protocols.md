# Observer Protocols

## Why Causal Order Alone Is Not Enough

Causal order can identify which events are accessible from which other events.
It does not, by itself, choose an observer, assign clock readings, or define an
observer-dependent spatial decomposition.

For radar reconstruction, the observer chain and its clock labels are additional
protocol structure.

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

