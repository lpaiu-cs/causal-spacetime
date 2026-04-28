# State-Change Echo Order

Milestone 22 extends reference-chain diagnostics from passive brackets to a
fixed-emission protocol. The setting is still a finite state-change trigger
network with a chosen reference chain. No clock labels, metric coordinates, or
calibrated units are introduced.

## Same-Emission Echo Protocol

Let

```text
R = (r_0, r_1, ..., r_m)
```

be a selected reference chain. Fix an emission reference position `k`. A target
event `e` participates in the same-emission echo diagnostic when:

```text
r_k ≺_T e
e ≺_T r_j for some j > k
```

The first condition says that the target is reachable after the chosen
reference event in the trigger order. The second says that the target has a
later return to the reference chain.

## Echo Return Position

The echo-return position is:

```text
return_R(e; k) = min { j : e ≺_T r_j and j > k }
```

This is a position along the chosen reference chain. It is not a clock reading.

## Echo-Delay Rank

The echo-delay rank is:

```text
delay_R(e; k) = return_R(e; k) - k
```

This is an order-level rank. It is not physical distance, metric radar
distance, or a calibrated duration.

For two reachable targets, the same-emission echo-order relation is:

```text
e_a ≺^{echo}_{R,k} e_b
  iff delay_R(e_a; k) < delay_R(e_b; k)
```

Ties indicate finite rank resolution rather than a quantitative equality of
physical quantities.

## Reference-Protocol Dependence

Echo order depends on both:

- the selected reference chain,
- the selected emission reference position.

Different high-utility reference chains can produce different reachable target
sets and different return-rank orderings. Different emission positions along
the same chain can also change reachability and tie structure. This is
reference-protocol dependence, not an inconsistency.

## Relation To Milestone 21

Milestone 21 computed passive predecessor and successor brackets:

```text
p_R(e) = max { i : r_i ≺_T e }
q_R(e) = min { j : e ≺_T r_j }
```

Milestone 22 instead fixes a reference emission position first and asks which
targets are reachable after that position and return later to the same
reference chain. The resulting echo-return position and echo-delay rank are
rank/order-level diagnostics for causal response structure relative to a
chosen reference chain.

## Limits

This milestone does not define spatial geometry, reconstruct metric radar
distance, implement a finite-speed spatial model, assign seconds or meters,
or introduce speeds, velocities, or metric coordinates. It provides finite
diagnostics for how causal response structure is ordered by a selected
reference protocol.

## Controlled Motif Validation

Milestone 23 adds controlled echo-response motifs as validation cases for this
protocol. A motif inserts trigger paths so a target event is after a selected
reference emission and before a planted return position. The planted
echo-delay rank is a known finite-DAG label. The recovered echo-delay rank is
then compared against that label.

Shortcut returns are recorded when additional trigger structure makes the
target return earlier than the planted label. These are background causal
interference in the order, not metric perturbations.
