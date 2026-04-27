# Relational Spatial Evolution

Milestone 16 separates two kinds of cross-slice statements.

Transport-dependent predicates include same position, same direction,
coordinate velocity, coordinate acceleration, and metric spatial evolution.
They require a chosen transport, anchor, persistence, calibration, or dynamics
protocol.

Some weaker relational-evolution statements do not require same-position
transport. They do require supplied object identity or persistence. Given
persistent labels and slice-local distance order, one can compare pair-distance
order histories:

```text
slice S_i: d(A, B) < d(C, D)
slice S_j: d(A, B) > d(C, D)
```

This says that an ordinal relation among persistent object pairs changed. It
does not identify absolute positions across slices and does not define
velocity.

## Pair-Distance Order History

A pair-distance order history records, for each slice, the order of distances
among persistent object pairs. The resulting ordinal shape signature is
invariant under per-slice translation, reflection, and positive scale changes.
Positive monotone transforms of slice-local distance values also preserve the
signature.

This makes the signature a candidate transport-gauge invariant: it is
persistence-dependent but transport-independent.

## Object Persistence Is Supplied

The repository does not derive persistent object identity from causal order.
Object labels are additional structure in these diagnostics. The point is to
test what becomes definable after persistence is supplied, while still avoiding
absolute same-position transport.

Milestone 17 studies the case where those labels are absent, partial, or
inferred by a matching criterion. Without persistence, pair-distance order
history is undefined. With an inferred matching, relational evolution becomes
hypothesis-dependent rather than absolute.

## Weaker Than Metric Dynamics

Relational spatial evolution is weaker than metric dynamics. It can say that
an ordinal shape changed, or that a pair-order relation remained stable. It
cannot assign velocities, accelerations, forces, or metric rates of change
without calibration and transport.

These experiments therefore support the order-first program by clarifying which
cross-slice content can be expressed relationally and which still belongs to
the effective metric representation layer. They do not prove spacetime
emergence or derive dynamics.
