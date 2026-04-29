# Response-Order Representability Preconditions

Milestone 27 adds scalar ordinal representability diagnostics for stable
response-order cores. This is a precondition diagnostic, not metric
representability.

## Scalar Ordinal Representability

Given a response-order sign matrix `Omega`, a nonzero pair imposes an ordinal
constraint:

```text
Omega_ij = -1  means target i has smaller D_echo than target j
Omega_ij =  1  means target i has larger D_echo than target j
```

Tied or unresolved pairs impose no scalar-rank constraint. The nonzero signs
define a directed relation. If that directed relation is acyclic, it admits a
topological rank representation.

## Topological Rank Representation

A topological rank representation assigns scalar ordinal ranks so every
nonzero response-order edge points from lower rank to higher rank. This checks
whether the stable response-order core is internally consistent as a scalar
ordinal relation.

The diagnostic reports:

- nonzero pair count,
- tied or unresolved pair count,
- directed 3-cycle count,
- cycle presence,
- scalar representability,
- rank span.

## Limits

Scalar ordinal representability is not metric representability. It does not
assign seconds, meters, speeds, velocities, coordinates, or calibrated
distances. It also does not prove that a response-order signature is geometry.

A scalar-rank-representable stable core can be used as input for future
representability tests. Those future tests would still need explicit
calibration or stronger physical protocol structure before making metric
claims.

## Milestone 28 Clarification

Scalar ordinal representability is only the first rung of a representability
ladder. It shows that the response-order relation can be represented by
topological scalar ranks. It does not provide pairwise target comparison data.

Milestone 28 therefore treats single-reference scalar response order as
underdetermined with respect to target-target pair ordering. Multi-reference
response profiles are richer inputs, but they remain pre-metric unless a
separate pairwise comparison or calibration protocol is supplied.
