# Reference-Chain Selection

Milestone 20 studies reference-chain utility inside finite state-change
trigger networks. A reference chain is any causal chain selected as a possible
reference backbone for later order-level protocols.

The utility score does not measure whether a chain is an observer. It measures
how useful the chain is for finite diagnostics such as coverage, two-sided
bracketing, and interval-profile regularity.

## Reference Chain

A reference chain is a finite sequence

```text
R = (r_0, r_1, ..., r_m)
```

with `r_i ≺_T r_j` for `i < j`. The chain gives ordered reference positions.
It does not provide clock labels, seconds, meters, velocities, or metric
coordinates.

## Candidate Sources

Milestone 20 compares:

- local-system chains from one local system's state-change history,
- greedy order-only chain candidates,
- longest order-only chain candidates,
- random chain baselines.

Local-system chains use event metadata. Order-only chain candidates use only
the transitive trigger order.

## Utility Diagnostics

The finite diagnostics are:

- comparability coverage: the fraction of events comparable to the chain,
- two-sided bracketing: the fraction of non-chain events with a chain event
  before and after,
- interval profile: causal interval cardinalities between adjacent chain
  elements,
- interval-profile regularity: variation in those interval cardinalities,
- local-system purity: how much of the chain comes from one local system.

None of these gives calibrated time, metric distance, or spatial geometry.

## Protocol-Reference Choice Dependence

Different chains can have similar utility scores while covering or bracketing
different events. The top-score gap and pairwise chain overlap diagnose this
protocol-reference choice dependence.

The highest-utility candidate is not a physical-observer claim. It is only a
chosen reference protocol for a specified finite diagnostic.

## Relation To Milestones 19 And 20

Milestone 19 introduced finite state-change trigger networks and checked
strict partial-order behavior. Milestone 20 added reference-chain utility
diagnostics for choosing candidate reference backbones inside those finite
orders. Milestone 21 then uses selected reference chains to compute
order-level bracket ranks. Milestone 22 fixes an emission position on a
selected reference chain and computes same-emission echo-return positions and
echo-delay ranks. Those echo diagnostics remain reference-protocol dependent
and do not add metric calibration.

Milestone 23 plants controlled echo-response motifs relative to selected
reference chains. This checks whether the same echo-order protocol recovers
known planted echo-delay ranks in finite DAGs and how visibility changes when
a different reference protocol is chosen.

Milestone 24 shows that shortcut classification is also reference-protocol
dependent. A shortcut return relative to one reference chain can be missing,
later, or differently ranked relative to another. This is protocol-reference
dependence in the order-level diagnostic, not a conflict between physical
observers.
