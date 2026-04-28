# Controlled Echo-Response Motifs

Milestone 23 adds controlled echo-response motifs to finite state-change
trigger networks. The purpose is to validate the reference-chain echo-order
protocol in networks where the intended response-depth labels are known by
construction.

## Motif Definition

Given a selected reference chain

```text
R = (r_0, r_1, ..., r_m)
```

a controlled echo-response motif chooses:

- a planted emission position `k`,
- a planted delay rank `d`,
- a target event `e`,
- an intended planted return position `k + d`.

The inserted trigger structure guarantees:

```text
r_k ≺_T e ≺_T r_{k+d}
```

The planted delay rank is a validation label for the controlled finite DAG. It
is not distance, not time, and not a speed.

## Motif Insertion

The implementation inserts ordinary state-changing events and trigger edges:

```text
r_k -> outward path -> target -> return path -> r_{k+d}
```

The outward and return paths can have zero or more intermediate events. The
existing reference-chain events are not duplicated. Event ids are storage
identifiers for the finite graph and are not physical time labels.

## Recovery

The reference-chain echo-order protocol recovers:

```text
return_R(e; k) = min { j : e ≺_T r_j and j > k }
delay_R(e; k) = return_R(e; k) - k
```

In a clean motif network, the recovered echo-delay rank should match the
planted echo-delay rank. If additional trigger structure creates an earlier
return, the recovered rank can be smaller than the planted label. Milestone 23
records this as a shortcut return or background interference.

## What This Validates

The motif experiments test whether the echo-order code can recover controlled
order-level response-depth motifs, distinguish exact recovery from shortcut
returns, and report rank-resolution limits when many targets share a small set
of delay ranks.

This is not metric reconstruction, calibrated radar distance, finite-speed
spatial geometry, quantum mechanics, curved spacetime, or a claim that a
reference chain is a physical observer. It is a controlled order-level response
diagnostic.
