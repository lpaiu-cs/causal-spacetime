# Validation Plan

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
