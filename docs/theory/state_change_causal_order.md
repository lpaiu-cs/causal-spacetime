# State-Change Causal Order

Milestone 18 refactors the theory-facing layer around a narrower primitive:
locally finite state-changing events ordered by causal trigger relations. The
goal is to keep the reconstruction program conservative while reducing the
appearance that every important structure is inserted independently.

Milestone 19 implements a minimal finite toy model for this primitive layer.
The toy model stores state-changing events in a finite trigger graph and uses
transitive closure for order queries. It does not add metric coordinates or
global physical time.

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
