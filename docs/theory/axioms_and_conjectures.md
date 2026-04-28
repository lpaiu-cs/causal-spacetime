# Axioms And Conjectures

This file sketches a conservative theory-facing hypothesis. It is not a theorem
and is not established by the simulations.

## Candidate Axioms

A0. Events are primitive state-changing events:

```text
e_{alpha,n} = (alpha, n, s_alpha^n -> s_alpha^{n+1})
```

A1. Causal trigger order `≺_T` is represented by a strict partial order.

A1a. Primitive temporal structure is causal trigger order, not a duration
value.

A2. Local finiteness gives discrete interval cardinalities.

A2a. Local finiteness does not imply global discrete time slices, a preferred
global update tick, or absolute discrete spatial cells.

A3. Counting measure supplies volume-like information.

A3a. Metric-scale reconstruction depends on a supplied measure, density, or
equivalent volume structure; causal order alone leaves conformal scale
underdetermined.

A3b. A measure may be represented either by supplied weights or by a specified
sampling/counting process, but the representation is still additional structure.

A4. Observer chains allow operational time and space decomposition.

A4c. Primitive spatial structure is observer-relative distance order, not an
absolute length value.

A4a. Signed spatial coordinates require orientation or reference structure
beyond a single observer chain.

A4b. Coordinate atlases require overlap and transition-map consistency between
multiple observer protocols.

A5. Lorentzian continuum behavior is a limiting property, not assumed at finite
`N`.

A6. Finite signal speed alone is insufficient; statistical Lorentz
compatibility or equivalent structure is required.

A7. Metric tensors, units, ratios, and curvature values belong to an effective
metric representation layer, not the primitive order layer.

A8. Quantum compatibility, if pursued, requires an additional amplitude,
Hilbert-space, or equivalent probabilistic layer. It is not supplied by the
current state-change order axioms.

## Candidate Conjectures

C1. Causal order determines conformal structure in suitable continuum limits.

C2. Causal order plus counting measure reconstructs volume-scaled geometry.

C2a. Positive conformal rescalings should preserve causal order while changing
physical volume and clock scale, so conformal factors require additional
measure information.

C2b. If events are sampled uniformly with respect to a supplied physical-volume
measure, local relative measure shape may be estimated statistically from event
counts, while global constant scale remains underdetermined without an
absolute density convention.

C2c. Reconstruction procedures should be stable under random thinning when the
density is rescaled appropriately, with larger finite-sampling noise at lower
retention probability.

C3. Timelike duration is recoverable from chain and interval structure under
appropriate density assumptions.

C4. Dimension is recoverable from order statistics in suitable Alexandrov
intervals.

C5. Observer-dependent spatial distance can be reconstructed from radar or
overlap protocols under additional assumptions.

C5e. Observer-relative distance order may admit a stable low-dimensional metric
representation only under additional consistency and representability
conditions.

C5f. Ratio stability should arise only after calibration, concatenation,
repeated processes, or dynamics restrict the admissible order-preserving
representations.

C5g. Rich, mutually consistent observer-relative distance-order data may admit
low-dimensional ordinal embeddings that serve as effective metric
representations, but noisy or inconsistent order data may not.

C5h. Structured geometric or observer-derived order constraints should
outperform shuffled, random, or inconsistent null-model constraints under
held-out validation and subset-stability diagnostics if they support a useful
effective metric representation.

C5i. Observer-relative spatial distance order should be defined only after a
slice-selection protocol. Radar-time bins derived from causal order and
observer tick order provide one controlled finite protocol for same-slice
spatial comparisons.

C5j. Cross-slice predicates such as same-position, same-direction, velocity,
constant velocity, and spatial evolution are undefined without a transport,
anchor, persistence, calibration, or dynamics rule.

C5k. Given a supplied transport or anchor protocol, cross-slice predicates
become protocol-relative and should be tested for stability under finite
sampling and noisy transport information.

C5l. Given supplied persistent object labels, pair-distance order histories
can define transport-gauge relational spatial evolution without identifying
same positions across slices.

C5m. Relational shape changes should remain invariant under independent
per-slice affine/reflection gauges, while coordinate velocity should not.

C5n. Without supplied persistence labels, anchors, dynamics, or a chosen
matching criterion, cross-slice object identity and pair-distance order
histories are undefined.

C5o. Relational-continuity matching can constrain persistence hypotheses in
asymmetric low-motion cases, but symmetric, noisy, or crossing configurations
may admit multiple equally compatible matchings.

C5a. Given a suitable observer chain with clock labels, radar time and radar
distance can be approximated from causal accessibility relations, with
resolution controlled by clock tick density and accessibility coverage.

C5b. Given a suitable two-chain oriented observer protocol with supplied beacon
separation, signed 1+1D coordinates can be approximated from causal
accessibility relations and clock labels.

C5c. Coordinate maps between compatible oriented inertial protocols should
approach Lorentz maps in suitable flat continuum limits.

C5d. Multiple compatible oriented protocols should form transition maps with
approximately consistent composition on chart overlaps in suitable flat
continuum limits.

C6. Horizons correspond to reconstruction-inaccessibility boundaries for
restricted observers.

C6a. In controlled flat-spacetime analogues, two-way radar accessibility for an
accelerated observer should identify the appropriate Rindler wedge up to finite
observer-chain coverage and clock-resolution effects.

These conjectures are targets for mathematical analysis and controlled
simulation. They are not presented as proven by this repository.
