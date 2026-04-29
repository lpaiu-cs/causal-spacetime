# State-Change Toy Model

Milestone 19 implements the first concrete toy model for the state-change
causal trigger primitive. The model is a finite directed acyclic trigger graph
whose vertices are state-changing events:

```text
e_{alpha,n} = (alpha, n, s_alpha^n -> s_alpha^{n+1})
```

Edges are immediate causal trigger relations. Their transitive closure is used
for causal-order queries.

## What The Model Implements

The model stores:

- state-changing events,
- local system chains,
- local successor trigger edges,
- external trigger edges,
- initial seed markers for initial events.

It validates finite order-theoretic properties:

- acyclicity of the immediate trigger graph,
- irreflexivity and transitivity of the transitive closure,
- finite causal intervals,
- local chain ordering inside each system.

## What Event IDs Mean

`event_id` and construction order are finite-storage and topological
bookkeeping. They are not physical global time. They are used only so the
finite graph can be generated and indexed deterministically.

## Local Chains Are Not Global Slices

Each local system has a chain of state changes. A local chain can be treated as
a reference protocol for diagnostics, but it is not a universal clock, a global
slice, or a calibrated observer protocol.

Local finiteness means intervals contain finitely many events. It does not
mean global synchronous ticking.

## Out Of Scope

Milestone 19 does not reconstruct metric geometry, extract observer chains from
raw networks, implement finite-speed spatial geometry, add quantum amplitudes,
or model curved spacetime. It is a finite DAG diagnostic for the primitive
state-change order layer.

The toy model does not derive relativity or quantum mechanics. It only checks
that a minimal state-change trigger network behaves as a locally finite strict
partial order.

## Milestone 20 Extension

Milestone 20 keeps the same finite trigger-network setting and adds
reference-chain selection diagnostics. Local-system chains, greedy
order-only chains, longest order chains, and random baseline chains are ranked
by coverage, two-sided bracketing, interval profiles, and ambiguity measures.

These diagnostics identify candidate reference protocols for later
order-level tests. They do not add metric coordinates, metric radar distance,
clock calibration, finite-speed spatial geometry, or a unique observer
selection rule.

## Milestone 21 Extension

Milestone 21 adds order-level bracket diagnostics for selected reference
chains. It records predecessor and successor reference positions, two-sided
accessibility, radar-time rank, bracket-width rank, and rank slices. These are
rank-level diagnostics only.

## Milestone 22 Extension

Milestone 22 adds same-emission echo-order diagnostics. Given a selected
reference chain and a fixed emission position, the code records which target
events are after that reference event and return to a later reference-chain
position. Echo-return position and echo-delay rank are order-level diagnostics
only. They do not reconstruct metric radar distance, define physical distance,
or calibrate a clock.

## Milestone 23 Extension

Milestone 23 inserts controlled echo-response motifs into finite trigger
networks. A motif plants a target event with a known echo-delay rank relative
to a selected reference chain. Clean networks test exact recovery; background
interference tests shortcut returns. The planted rank is a validation label
for the finite order, not a metric quantity.

## Milestone 24 Extension

Milestone 24 classifies shortcut returns using return spectra. It records all
later reference-chain return positions for a motif target and separates
targeted shortcut-return stress tests from generic acyclic background edge
perturbations. This is a finite DAG interference diagnostic, not a metric or
geometry layer.

## Milestone 25 Extension

Milestone 25 tests return-spectrum stability under finite coarse-graining.
Closure-preserving event thinning hides intermediate events without deleting
their reachability effect. Immediate-edge thinning removes trigger edges and
can delete causal paths. Reference-chain subsampling coarsens reference ranks
and can reduce order resolution.

## Milestone 26 Extension

Stable echo-response order signatures build on the same finite trigger-network
toy model. They use controlled motif targets and reference-chain protocols to
compare response-rank order across protocol variants. The toy model still does
not introduce metric coordinates, calibrated clocks, speed, or spatial
geometry.
