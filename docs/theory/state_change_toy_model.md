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
an observer-like ordered protocol for diagnostics, but it is not a universal
clock, a global slice, or a metric observer with calibrated seconds.

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
observer-like chain selection diagnostics. Local-system chains, greedy
order-only chains, longest order chains, and random baseline chains are ranked
by coverage, two-sided bracketing, interval profiles, and ambiguity measures.

These diagnostics identify candidate observer protocols for later operational
tests. They do not add metric coordinates, radar distance, physical clock
calibration, finite-speed spatial geometry, or a unique observer selection
rule.
