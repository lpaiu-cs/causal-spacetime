# State-Change Causal Order

Milestone 18 refactors the theory-facing layer around a narrower primitive:
locally finite state-changing events ordered by causal trigger relations. The
goal is to keep the reconstruction program conservative while reducing the
appearance that every important structure is inserted independently.

Milestone 19 implements a minimal finite toy model for this primitive layer.
The toy model stores state-changing events in a finite trigger graph and uses
transitive closure for order queries. It does not add metric coordinates or
global physical time.

Milestone 20 adds reference-chain selection diagnostics in those finite trigger
networks. Local-system chains and order-only chain candidates can be ranked by
coverage, two-sided bracketing, interval-profile regularity, and
protocol-reference choice diagnostics. This is a reference-protocol diagnostic,
not an observer-reality claim or a metric reconstruction.

Milestone 21 then computes order-level predecessor and successor brackets from
selected reference chains. The resulting radar-time ranks and bracket-width
ranks are rank-level diagnostics, not calibrated clock or metric distance
quantities.

Milestone 22 adds fixed-emission echo-order diagnostics. After choosing a
reference chain and an emission position, the code asks which target events
are after that emission and return to a later reference-chain position. The
echo-return position and echo-delay rank remain order-level quantities.

Milestone 23 adds controlled echo-response motifs. These motifs deliberately
insert trigger paths with known planted echo-delay ranks so the echo-order
protocol can be checked against finite DAG validation labels.

## State-Changing Events

A state-changing event is written:

```text
e_{alpha,n} = (alpha, n, s_alpha^n -> s_alpha^{n+1})
```

Here `alpha` denotes a local system, node, or physical subsystem. The index `n`
orders state changes along that subsystem. The event is the transition itself:
the change from state `s_alpha^n` to state `s_alpha^{n+1}`. Events are
primitive in this theory-facing layer. Coordinates, metric durations, and
lengths are not part of the primitive event definition.

## Causal Trigger Relation

The causal trigger relation is written:

```text
e_i ≺_T e_j
```

It means that the state change at `e_i`, or information/trigger emitted by
`e_i`, enables or contributes to the state change at `e_j`. This relation is an
order-theoretic abstraction. It does not assign a metric duration, speed, or
spatial separation by itself.

## Strict Partial Order

The trigger relation is assumed to be a strict partial order:

```text
not (e ≺_T e)
```

and

```text
e_i ≺_T e_j and e_j ≺_T e_k => e_i ≺_T e_k
```

Irreflexivity excludes an event triggering itself. Transitivity records
indirect trigger dependence. These axioms define causal ordering, not a global
clock.

## Local Finiteness

For comparable events `a ≺_T b`, define the interval:

```text
I(a,b) = { e in E : a ≺_T e ≺_T b }
```

Local finiteness requires:

```text
|I(a,b)| < infinity
```

This is a finite-interval condition. It does not imply global discrete time
slices, a preferred update tick, or a regular lattice.

## Observer Chains

An observer chain is a selected ordered sequence of events:

```text
O = (o_0, o_1, ..., o_m)
o_i ≺_T o_j for i < j
```

Observer time is the causal order along such a chain. Seconds are a calibration
layer placed on top of the chain, not primitive structure. The observer chain
itself is also protocol structure; it is not automatically derived from the
bare trigger order in the current simulations.

## Reference-Chain Candidates

A finite state-change network may contain many valid chains. Milestone 20
therefore treats reference-chain selection as a diagnostic problem:

- local-system chains use the stored local-system event metadata,
- greedy and longest order chains use the transitive trigger order,
- random order chains provide baselines.

Candidate quality is summarized by comparability coverage, two-sided
bracketing, adjacent interval cardinality profiles, local-system purity, and
top-score gaps. These quantities help identify useful reference protocol
candidates inside a finite order. They do not supply calibrated time or select
a physically privileged observer by themselves.

## Reference-Chain Brackets

For an event `e` and a reference chain `R`, define order-level brackets:

```text
p_R(e) = max { i : r_i ≺_T e }
q_R(e) = min { j : e ≺_T r_j }
```

When both exist, `e` is two-sided accessible relative to `R`. Milestone 21
uses these positions to define `T_rank = p_R + q_R` and `W_rank = q_R - p_R`.
`W_rank` is a bracket-width rank, not a metric spatial distance.

## Same-Emission Echo Order

For a selected reference chain `R` and fixed emission position `k`, define the
echo-return position:

```text
return_R(e; k) = min { j : r_k ≺_T e ≺_T r_j and j > k }
```

The echo-delay rank is:

```text
delay_R(e; k) = return_R(e; k) - k
```

This gives a same-emission echo-order relation among reachable targets by
comparing delay ranks. It is a reference-protocol-dependent rank diagnostic,
not a spatial metric, calibrated time, or physical-distance assignment.

## Controlled Echo-Response Motifs

A controlled motif plants trigger structure:

```text
r_k -> ... -> e -> ... -> r_{k+d}
```

The label `d` is the planted echo-delay rank. Recovery compares this label with
the order-level echo-delay rank computed from the transitive trigger order. If
additional trigger structure gives an earlier return to the reference chain,
the result is recorded as a shortcut return or background interference.

The planted label is a controlled validation tag. It does not represent metric
distance, calibrated duration, or a finite-speed spatial model.

## Radar Brackets

For an event `e` and observer chain `O`, define the order-level radar brackets:

```text
p_O(e) = max { i : o_i ≺_T e }
q_O(e) = min { j : e ≺_T o_j }
```

If both brackets exist, define rank-level quantities:

```text
T_rank(e) = p_O(e) + q_O(e)
R_rank(e) = q_O(e) - p_O(e)
```

These use only observer tick order and causal trigger order. They do not use
metric coordinates, seconds, or meters.

## Observer Slices

Given bin width `w`, an observer slice is:

```text
S_{O,k}^{(w)} = { e : floor(T_rank(e) / w) = k }
```

This is an observer-relative slice-selection protocol. It is not a global time
slice and does not select an absolute simultaneity surface.

## Observer-Slice-Relative Distance Order

Within an observer slice `S`, spatial comparison is order-level and
observer-relative:

```text
(a,b) ≺^d_{O,S} (c,d)
```

This means that, within observer `O`'s slice `S`, pair `a-b` is closer than
pair `c-d`. Such order can be induced by finite-speed information exchange,
radar-return order, or another specified observer protocol. It remains weaker
than metric distance.

## Metric Representation

A metric representation of a slice assigns coordinates:

```text
phi_{O,S}: S -> X
```

with a metric `d_X` satisfying the order constraint:

```text
(a,b) ≺^d_{O,S} (c,d)
  => d_X(phi(a), phi(b)) < d_X(phi(c), phi(d))
```

This representation is effective. It can be useful when stable, calibrated,
and compatible across observer protocols, but it is not primitive. Seconds,
meters, ratios, velocity, curvature values, and quantitative dynamics belong to
the representation layer.

## Echo Shortcut Classification

Milestone 24 remains inside the finite state-change trigger-network layer. It
adds return-spectrum diagnostics for controlled echo-response motifs and
classifies exact recovery, shortcut returns, and missing or late returns.
Shortcut returns are earlier causal paths in the transitive trigger order.
They are not metric quantities and do not add spatial geometry.
