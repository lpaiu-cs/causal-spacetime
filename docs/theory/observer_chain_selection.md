# Observer-Chain Selection

Milestone 20 studies observer-like chain selection inside finite
state-change trigger networks.

This is an order-theoretic diagnostic for possible observer protocols. It is
not metric reconstruction, radar distance reconstruction, or a derivation of a
unique observer.

## Observer-Like Chain

An observer-like chain is a finite sequence of events

```text
O = (o_0, o_1, ..., o_m)
```

such that `o_i ≺_T o_j` for `i < j` in the transitive causal trigger order.
The chain supplies an ordered protocol candidate. It does not by itself supply
seconds, meters, velocities, or a calibrated physical clock.

## Candidate Sources

Milestone 20 compares several candidate sources:

- local-system chains: state-change sequences belonging to one local system,
- greedy order-only chains: heuristic walks through the causal order,
- longest order chains: dynamic-programming chain candidates in the finite
  transitive order,
- random order chains: baselines for selection ambiguity.

Local-system chains use event metadata. Greedy, longest, and random order
chains use only the order matrix once the finite network has been generated.

## Chain Quality Diagnostics

The diagnostics are finite, heuristic, and representation-facing:

- comparability coverage: the fraction of events comparable to at least one
  chain event or lying on the chain,
- two-sided bracketing: the fraction of non-chain events with at least one
  chain event before and one chain event after,
- interval profile: the causal interval cardinalities between adjacent chain
  elements,
- interval regularity: the coefficient of variation of those interval
  cardinalities,
- local-system purity: the fraction of chain events coming from the most common
  local system.

Interval profiles are order-statistical regularity diagnostics. They are not
seconds and are not a calibrated clock.

## Selection Ambiguity

Different chains can have similar coverage, bracketing, and regularity scores.
The top-score gap and pairwise chain overlap are therefore recorded as
ambiguity diagnostics. A high score means that a chain is useful for the chosen
finite diagnostics, not that it is the true physical observer.

## Relation To Milestone 19

Milestone 19 checked that finite state-change trigger networks behave as
locally finite strict partial orders. Milestone 20 adds a first observer-like
protocol selection layer on top of that finite order. It deliberately stops
before radar distance reconstruction, finite-speed spatial geometry, metric
representation, quantum amplitudes, or curved-spacetime dynamics.
