# Persistence And Identity Matching

Milestone 17 adds the next missing layer after relational spatial evolution:
object persistence. Milestone 16 assumed persistent object labels. Milestone 17
tests what remains undefined when those labels are absent, and what becomes
defined only after a persistence hypothesis is supplied or inferred.

Milestone 18 keeps persistence outside the primitive layer. The state-change
causal trigger order supplies event ordering, not object identity across
observer slices. Persistence remains supplied, inferred, or dynamically
constrained structure.

Object identity is additional structure in this project. Causal order,
observer protocols, and slice-local distance order do not by themselves say
that an object in slice `S_i` is the same object as one in slice `S_j`.
Without supplied persistence, matching, anchors, partial labels, or dynamics,
cross-slice object identity is undefined rather than false.

## Why Relational Histories Require Persistence

A pair-distance order history compares persistent object pairs across slices:

```text
slice S_i: d(A, B) < d(C, D)
slice S_j: d(A, B) > d(C, D)
```

This statement assumes that `A`, `B`, `C`, and `D` have been identified across
slices. Without that identity relation, the pair `AB` in one slice cannot be
compared with `AB` in another slice. Relational histories are therefore
persistence-relative.

## Persistence Hypotheses

When labels are not supplied, a matching rule can be used as a persistence
hypothesis. Milestone 17 implements a finite relational-continuity diagnostic:
enumerate candidate permutations between adjacent unlabeled slices and rank
them by how well they preserve slice-local pair-distance order signatures.

This is a hypothesis-selection diagnostic, not a derivation of object identity.
The best matching according to one criterion should not be described as the
true identity relation. Hidden labels are used only for validation in
controlled experiments.

## Matching Ambiguity

Symmetric or noisy configurations can admit multiple equally good identity
matchings. In such cases, the ambiguity gap between the best and second-best
matching can be small or zero. This is expected: slice-local relational order
does not contain enough information to distinguish all persistence hypotheses.

Crossing histories provide another failure mode. When objects pass through
configurations with similar relational signatures, a continuity criterion may
swap identities or become ambiguous. This shows why persistence may require
additional anchors, partial labels, matter dynamics, or smoothness assumptions.

## Partial Labels And Anchors

Partial identity labels or anchors restrict the admissible matching set. In the
Milestone 17 experiments, fixed-point constraints are supplied from controlled
validation labels and then used to filter candidate permutations. More supplied
identity information should reduce ambiguity, but this information remains
additional structure.

## Hypothesis-Dependent Relational History

Different persistence hypotheses can produce different pair-distance order
histories and different relational-evolution claims. Relational history is
therefore not absolute unless persistence is part of the protocol or dynamics.
This remains weaker than velocity or metric dynamics, which also require
transport and calibration.

Milestone 17 supports the order-first program by making the identity layer
explicit. It does not claim that causal order alone tracks objects, that object
identity is primitive, or that the simulations prove spacetime emergence.
