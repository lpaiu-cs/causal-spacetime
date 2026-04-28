# Validation Plan

## Milestone 18: State-Change Causal Order Refactor

Milestone 18 is validated primarily by documentation and language checks. The
theory-facing docs should consistently use the hierarchy:

```text
locally finite state-changing events
  -> causal trigger order
  -> observer-chain time
  -> observer-slice-relative distance order
  -> calibrated effective metric representation
```

The validation target is conceptual consistency, not a new numerical
experiment. The language guard checks that theory docs avoid risky overclaim
phrases except in explicit rejected-language contexts.

## Milestone 19: Minimal State-Change Toy Model

Milestone 19 validates a finite DAG diagnostic for state-changing events and
causal trigger relations. Exact tests check a hand-coded network. The toy-model
experiment checks generated networks for irreflexive transitive closure and
finite intervals. The reference-chain diagnostic checks whether one local
system chain is totally ordered by the closure without treating it as a global
clock.

## Milestone 20: Reference-Chain Selection

Milestone 20 validates finite reference-chain selection criteria. Exact
tests check chain validity, comparability masks, two-sided bracketing, bracket
positions, quality reports, and ranking output on a hand-coded state-change
network.

The numerical experiments compare local-system chains, greedy order chains,
longest order chains, and random baselines. They track coverage, bracketing,
interval-profile regularity, local-system purity, top-score gaps, and
candidate overlap. The validation target is reference-chain utility and
protocol-reference choice dependence, not metric reconstruction or clock
calibration.

## Milestone 21: Reference-Chain Bracket Diagnostics

Milestone 21 validates order-level brackets induced by selected reference
chains. Exact tests check predecessor and successor reference positions,
accessible masks, radar-time ranks, bracket-width ranks, rank slices, and
coverage summaries.

The experiments compare bracket accessibility and rank slices across
reference-chain sources, test pairwise bracket-rank dependence among
high-utility reference chains, and measure how trigger density changes
two-sided accessibility. The validation target is what a chosen reference
protocol can order-access, not metric radar distance or calibrated time.

## Milestone 22: Same-Emission Echo Order

Milestone 22 validates fixed-emission echo-order diagnostics relative to
selected reference chains. Exact tests check echo-return positions, echo-delay
ranks, reachable masks, echo-order matrices, and summary fields on a
hand-coded finite trigger network.

The experiments compare same-emission echo reachability across reference-chain
sources, test reference-protocol dependence among high-utility chains, measure
emission-position sensitivity, and study trigger-density effects on
reachability and rank resolution. The validation target is order-level causal
response structure, not physical distance, calibrated time, or metric radar
reconstruction.

## Milestone 23: Controlled Echo-Response Motifs

Milestone 23 validates the echo-order protocol against planted causal trigger
motifs. Exact checks insert motifs on a minimal reference backbone and verify
that the recovered echo-delay rank matches the planted validation label.

The numerical experiments separate clean recovery, background shortcut
interference, motif-density tie structure, and reference-choice visibility.
The validation target is controlled order-level response depth. Planted delay
rank is not physical distance, and shortcut returns are background causal
interference rather than errors in a metric representation.

## Milestone 24: Echo Shortcut And Interference Classification

Milestone 24 validates return-spectrum and shortcut-classification utilities.
Exact sanity checks verify that a planted motif has the expected return
spectrum, that an injected earlier return changes the recovered delay rank,
and that the immediate graph remains acyclic.

The numerical experiments separate targeted shortcut-return injection from
generic acyclic background edge perturbations. Targeted injection sweeps test
how exact recovery and planted-order agreement degrade under controlled
earlier returns. Background perturbation sweeps test whether generic extra
trigger edges create new shortcuts or deepen existing ones. Additional
experiments measure motif path-length robustness and reference-protocol
dependence. The validation target is causal interference classification before
any stronger distance-order interpretation.

## Milestone 25: Return-Spectrum Stability Under Coarse-Graining

Milestone 25 validates coarse-graining utilities for echo return spectra.
Exact sanity checks verify protected event selection, closure-preserving
restriction, motif remapping, reference-chain subsampling, expected coarse
delay ranks, and acyclic immediate-edge thinning.

The numerical experiments separate three cases. Closure-preserving event
thinning should preserve return spectra among retained reference and target
events. Immediate-edge thinning can delete path information and create missing
or delayed recoveries. Reference-chain subsampling reduces rank resolution,
creates ties, and shifts expected coarse return ranks. Shortcut
classification is then compared across these coarse-graining choices.

## Why The Full-Diamond Check Is Limited

The original endpoint reconstruction experiment uses the full generated causal
diamond between two boundary events. If the same total event count is used both
to define the global density and to estimate the full interval volume, then the
interval-count proper-time estimate can become nearly tautological:

```text
rho = N / V_full
tau_hat = sqrt((N / rho) / eta_2)
```

For the same full interval, this algebra mostly returns the input scale. The
experiment remains a useful plumbing check, but it is not strong validation of
an operational reconstruction procedure.

## Internal-Pair Reconstruction

The non-tautological validation experiment generates a larger known 1+1D
Minkowski causal diamond and then samples many internal timelike event pairs.
For each internal pair, the code measures the Alexandrov interval cardinality
between those two events and compares the reconstructed proper time against the
hidden Minkowski value.

This is more informative because the density is supplied once from the known
generated region, while each tested sub-interval has an independently varying
size and location.

## What Event Density Contributes

Event density is additional scale information:

```text
rho = N_events / V_region
```

It converts a causal interval cardinality into a volume estimate. Without this
scale, causal order can identify relational structure but not an absolute metric
scale in these simulations.

## What Causal Order Contributes

Causal order determines which sampled events lie between a timelike pair:

```text
i precedes k and k precedes j
```

The interval count is then used as a finite-sampling proxy for continuum
Alexandrov volume. The validation checks whether this proxy behaves as expected
in known 1+1D Minkowski spacetime.

## Longest Chains

The longest-chain diagnostic uses the 1+1D asymptotic relation:

```text
L ~ sqrt(2 * rho) * tau
```

The implementation documents whether endpoint events are included in the chain
length. The current estimator subtracts two endpoint contributions when the
chain length includes both selected endpoints. This is a conservative convention
for the existing dynamic-programming implementation, but longest-chain estimates
are finite-size sensitive and should not be treated as precisely calibrated.

## Scientific Scope

These experiments validate reconstruction procedures inside a known spacetime
model. They do not show that causal order alone determines metric scale, and
they do not prove a new theory of spacetime. The project should continue to
report failure modes, finite-size effects, boundary effects, and assumptions
explicitly.

## Milestone 3: Statistical Calibration

The independent probe-pair experiment improves the statistical interpretation of
the interval-count estimator. Support events are sprinkled once in the global
diamond, while probe endpoint pairs are sampled independently and are not
inserted into the support set. This makes the interval count for each probe pair
cleaner to compare with finite-sampling models.

For Poisson sprinkling, the interval count has mean:

```text
lambda = rho * tau^2 / 2
```

For the current fixed-count support sprinkling, the closer comparison is:

```text
K ~ Binomial(N_support, V_interval / V_global)
```

The validation question is whether observed RMSE is broadly consistent with
these sampling-noise predictions. A mismatch would point to bias, boundary
effects, implementation mistakes, or convention issues. A match is only a
consistency check in known 1+1D Minkowski spacetime.

## Milestone 4: Dimension Reconstruction

Dimension reconstruction tests a different kind of claim from timelike
proper-time estimation. Instead of assuming the dimension and estimating a scale,
the Myrheim-Meyer estimator asks whether the causal-order relation fraction in a
flat Alexandrov interval is consistent with a spacetime dimension.

The experiment generates known D-dimensional Minkowski intervals for `D = 2, 3,
4`, computes:

```text
r_obs = (# ordered causal relations) / (N * (N - 1))
```

and inverts the theoretical Myrheim-Meyer curve:

```text
r(D) = Gamma(D + 1) * Gamma(D / 2) / (4 * Gamma(3D / 2))
```

This is closer to the mathematical reconstruction program because dimension is
treated as an order-statistical observable. The assumptions remain substantial:
flat intervals, known sprinkling model, finite sampling, and a continuum formula
used for interpretation.

## Milestone 5: Observer-Chain Radar Validation

The discrete radar experiment validates a causal-order-based radar protocol in
known 1+1D Minkowski intervals. It supplies an observer chain and clock labels,
then reconstructs target-event radar coordinates from order relations alone:

```text
latest observer tick before target -> tau_minus
earliest observer tick after target -> tau_plus
```

The validation checks whether radar time and radar distance approach the hidden
stationary-observer values as clock tick density increases. The observer chain
and clock labels are additional protocol structure, so this is not a test of
causal order alone.

## Milestone 6: Orientation And Lorentz-Map Validation

The oriented radar experiment validates a two-chain observer protocol in known
1+1D Minkowski intervals. A single observer chain supplies only unsigned radar
distance, so mirrored events `(t, x)` and `(t, -x)` remain degenerate. The
reflection-degeneracy diagnostic records that failure mode directly.

The two-chain protocol supplies a synchronized comoving beacon chain with known
rest-frame separation. This orientation reference is additional protocol
structure. Given radar distances to the primary and beacon chains, the protocol
reconstructs a signed coordinate and compares it against hidden validation
coordinates.

The Lorentz-map recovery experiment then compares two oriented protocols: one
stationary and one moving inertially. In the controlled validation setting, the
map between their reconstructed coordinates should approach the known Lorentz
map as observer tick density increases. This is a test of protocol consistency
inside known special relativity, not a proof that causal order alone determines
signed space or coordinate transformations.

## Milestone 7: Observer Atlas Consistency

The observer-atlas experiment uses three supplied oriented inertial protocols
with different boosts and affine origins. Each protocol reconstructs chart
coordinates for support events from the same causal matrix, observer-chain
indices, beacon-chain indices, and clock labels.

For overlapping chart pairs, the validation fits:

```text
target_chart ~= L(beta_relative) source_chart + translation
```

and checks fitted beta error, transition residual RMSE, accessible overlap
fraction, and invariant interval disagreement. For the chart loop
`A -> B -> C`, it compares the composed transition map with the direct `A -> C`
fit.

This moves closer to atlas structure because it tests overlapping charts and
transition-map consistency. The origins, clocks, beacon separations, and
orientation references remain supplied protocol data, and the validation is
still performed inside known 1+1D Minkowski intervals.

## Milestone 8: Rindler Horizon Validation

The Rindler experiment tests reconstruction-inaccessibility for an accelerated
observer in flat 1+1D Minkowski spacetime. It supplies a uniformly accelerated
observer chain and proper-time clock labels, then reconstructs two-way radar
accessibility from causal order alone relative to that protocol.

The validation explicitly separates:

```text
ideal Rindler wedge accessibility
finite observer-chain coverage
reconstructed two-way radar accessibility
```

The expected result is that reconstructed accessibility mostly matches
finite-chain coverage, while events outside the Rindler wedge generally lack
two-way radar reconstruction. False positives should remain low, and
radar-coordinate errors inside finite coverage should decrease as tick
resolution increases.

This is a controlled flat-spacetime horizon analogue, not a black hole
simulation and not a proof that causal order alone derives horizons.

## Milestone 9: Conformal Ambiguity

The conformal ambiguity experiment validates an important limitation. Positive
conformal rescaling preserves causal order, so causal-order statistics such as
relation fraction and Myrheim-Meyer dimension do not fix metric scale or
physical volume.

The constant-scale check uses the same causal matrix while changing:

```text
proper time by a factor alpha
volume by a factor alpha^2
```

The weighted conformal volume experiment then supplies local `Omega^2` measure
weights and compares weighted interval-volume estimates with unweighted
coordinate-volume estimates. This tests measure-dependent reconstruction in a
controlled 1+1D conformal toy model.

The validation should be read as a limitation statement: causal order determines
light-cone/conformal structure, while density or volume measure is additional
metric-scale structure.

## Milestone 10: Measure Encoding And Thinning Stability

Milestone 10 tests how measure information behaves when it is represented by
counting statistics. In the physical-volume sprinkling experiment, support
events are sampled with probability density proportional to `Omega^2` in
coordinate space. The causal order is still the flat causal order, but the event
distribution encodes the supplied conformal measure model.

The validation compares two density conventions:

```text
rho_physical = N / V_global_Omega
rho_coordinate = N / V_global_coordinate
```

When events are sampled from conformal physical volume, counts divided by
`rho_physical` should reconstruct conformal physical interval volumes. Counts
divided by coordinate density are a deliberate failure baseline for non-flat
physical volume.

The local profile experiment bins events by coordinate time and compares the
normalized density shape to the normalized `Omega^2` shape. This can recover
relative local measure statistically, but it cannot identify an overall
constant conformal scale from normalized event positions alone.

The thinning experiment keeps each event with probability `p`. Reconstruction
should remain approximately stable when density is rescaled to `p * rho`, while
using the original density after thinning should produce a scale error. This is
a coarse-graining stability check, not proof of continuum convergence.

## Milestone 11: Order-First Validation Layer

Milestone 11 changes the theory-facing interpretation from metric-first to
order-first. Earlier experiments mostly reconstructed metric representation
values in controlled settings. The new validation layer asks whether ordinal
relations themselves are preserved:

```text
causal order -> primitive temporal order
observer protocol + causal order -> observer-relative distance order
```

The radar-return experiment reconstructs a same-emission return-tick order from
causal order and observer tick order without numeric clock labels. The oriented
chart experiment compares pairwise spatial distance order, which is weaker than
coordinate RMSE but closer to the order-first thesis.

The scale and calibration experiments show that positive monotone
transformations preserve order while changing ratios. Ratios become stable only
after additional calibration or dynamics restricts the admissible
representations.

The metric-representability diagnostics show that a distance order or candidate
distance matrix is not automatically a useful metric geometry. Cycles,
triangle-inequality failures, and embedding diagnostics are separate
representability checks.

## Milestone 12: Ordinal Embedding Validation

Milestone 12 tests effective metric representation directly. It samples
quadruplet distance-order constraints and fits low-dimensional coordinates
using a simple finite ordinal-embedding optimizer.

The validation separates several questions:

```text
Does the embedding satisfy the sampled order constraints?
Does increasing constraints improve a consistent representation?
Does ordinal stress drop at or above an effective dimension?
How much do noisy or incomplete comparisons degrade the representation?
Can observer-derived distance order support a 1D effective spatial chart?
```

The hidden coordinates are used only to generate controlled constraints or to
validate fitted embeddings. Ordinal embedding is treated as a finite diagnostic
for representability, not a proof that all distance orders are metric or that
metric geometry is fundamental.

## Milestone 13: Held-Out Validation And Null Baselines

Milestone 13 adds stronger baselines for the ordinal embedding layer. The core
validation question is whether structured order data supports a stable
low-complexity representation better than shuffled, random, noisy, or
inconsistent constraints.

The validation separates:

```text
training constraint satisfaction
held-out order generalization
stability across independent constraint subsets
observer-derived order quality relative to null baselines
complexity curves across candidate embedding dimension
```

Structured constraints should generally show lower held-out violation, better
subset stability, and cleaner low-dimensional compression than null-model
constraints. These are finite diagnostics for effective metric representation,
not physical proof that distance order alone determines geometry.

## Milestone 14: Simultaneity-Sliced Spatial Order

Milestone 14 tests whether observer-derived spatial distance order improves
when comparisons are restricted to observer-selected radar-time slices. The
slice labels are reconstructed from causal order and observer tick order using
radar tick brackets; hidden coordinates are used only for validation.

The validation separates:

```text
radar-time order recovery from tick brackets
same-slice distance-order preservation
sliced observer-derived order versus null baselines
slice-width sensitivity
```

This keeps spatial distance explicitly observer-relative and slice-protocol
dependent. It does not claim that spatial distance is observer-independent or
that causal order alone determines metric geometry.

## Milestone 15: Cross-Slice Identification And Transport

Milestone 15 validates the statement that cross-slice predicates are undefined
without additional identification structure. The experiments check that
same-slice constraints decompose into slice-local components, that slice-local
embeddings can be validated independently, that random per-slice gauges
preserve same-slice order while changing cross-slice judgments, that supplied
anchors constrain transport, and that noisy transport destabilizes derived
spatial-evolution statements.

## Milestone 16: Relational Spatial Evolution

Milestone 16 validates transport-gauge relational statements. The experiments
check the predicate-definability table, synthetic pair-distance order histories,
gauge invariance under independent per-slice affine/reflection transformations,
observer-derived slice histories for supplied persistent objects, and the
contrast between relational invariants and transport-dependent velocities.

## Milestone 17: Persistence Ambiguity And Identity Matching

Milestone 17 validates the next layer: object persistence. The tests verify
that object identity and pair-distance histories are reported as undefined
without persistence. The experiments compare identity-matching hypotheses in
symmetric, small-motion, partial-label, crossing, and deliberately alternative
histories. Hidden labels are used only for controlled validation of matching
accuracy and track consistency.

Expected behavior is conservative: asymmetric low-motion histories may be
matched by relational continuity, while symmetric or crossing cases can remain
ambiguous or fail. Partial labels should reduce ambiguity because they are
additional supplied structure.
