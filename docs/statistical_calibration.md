# Statistical Calibration

## Why Milestone 2 Was Useful

Milestone 2 sampled internal timelike event pairs from a single sprinkled causal
set and reconstructed each pair's proper time from Alexandrov interval
cardinality. This is operationally natural: the endpoints and interval elements
are all part of the same finite causal order.

That experiment is more informative than reconstructing only the full generated
diamond, but it still mixes endpoint sampling and interval-count sampling. The
observed error includes fluctuations in where the endpoint events happened to
land as well as fluctuations in how many events lie between them.

## Why Independent Probe Pairs Help

Milestone 3 adds independent probe endpoint pairs. The procedure is:

```text
1. sprinkle N support events in a known 1+1D Minkowski causal diamond,
2. independently sample many hidden probe timelike endpoint pairs,
3. do not insert the probe endpoints into the support set,
4. count support events in each probe pair's Alexandrov interval,
5. estimate proper time from interval count and supplied global density.
```

This separates support-event counting noise from endpoint sampling. The hidden
probe coordinates are used for validation and for selecting timelike probes, not
as input to the causal interval count.

## Poisson And Fixed-N Models

For a 1+1D timelike interval:

```text
V_interval = tau^2 / 2
```

If support events are a Poisson process with density `rho`, then the interval
count has mean:

```text
lambda = rho * tau^2 / 2
```

The project often uses fixed-count sprinkling inside a global diamond. In that
case, for a global diamond of duration `T`:

```text
K ~ Binomial(N_support, V_interval / V_global)
V_global = T^2 / 2
```

The Poisson formulas remain useful approximations, especially when the probed
interval is small compared with the global region. The binomial model is the
closer finite-N comparison for the current fixed-count support sprinkling.

## Delta-Method Error Scale

The interval-count estimator is:

```text
tau_hat = sqrt(2K / rho)
```

For Poisson fluctuations, the leading delta-method absolute standard deviation
is approximately:

```text
std(tau_hat) ~= 1 / sqrt(2 rho)
```

The relative standard deviation scales like:

```text
std(tau_hat) / tau ~= 1 / (sqrt(2 rho) * tau)
```

This is why small intervals have large relative error. Even if the absolute
error scale is modest, dividing by a small true proper time amplifies relative
error.

## Interpretation

Matching the expected finite-sampling error scale is a consistency check for the
operational reconstruction procedure in known 1+1D Minkowski spacetime. It is
not evidence that causal order alone determines metric scale. Event density is
still supplied as additional scale information, and agreement with sampling
noise does not prove a new spacetime theory.

## Relation To Dimension Reconstruction

Milestone 4 extends calibration beyond timelike scale by estimating dimension
from the causal-order relation fraction. Statistical scatter still matters:
finite `N` changes the observed relation fraction and therefore the inverted
dimension estimate.

The dimension experiment is theory-facing because it treats dimension as an
order-statistical observable. It remains a controlled validation inside known
flat causal intervals, not a proof that physical spacetime dimension has been
derived from information alone.

## Relation To Oriented Radar

Milestone 6 adds a different validation target: orientation and coordinate-map
compatibility. The statistical-calibration experiments focus on interval counts
and finite-sampling noise, while oriented radar focuses on finite observer-clock
resolution and accessibility coverage.

The signed-coordinate reconstruction uses a supplied beacon chain and known
beacon separation. Those protocol choices are additional structure, just as
event density is additional scale information in interval-count reconstruction.

## Relation To Atlas Consistency

Milestone 7 uses different diagnostics again: affine transition-map residuals,
relative-beta errors, invariant interval disagreement, accessible overlap, and
loop consistency. These are not Poisson interval-count errors. They are
finite-protocol-resolution and overlap diagnostics for reconstructed observer
charts.

## Relation To Rindler Horizon Validation

Milestone 8 is primarily an accessibility and finite-coverage validation. Its
error sources are observer clock resolution, finite chain duration, and the
observer-dependent Rindler wedge, rather than interval-count sampling noise.
The experiment therefore reports confusion rates and radar-coordinate errors
relative to finite-chain coverage.

## Relation To Conformal Measure Dependence

Milestone 9 separates counting noise from measure ambiguity. Uniform coordinate
sprinkling plus unweighted counts estimates coordinate volume. Conformal
physical volume requires additional `Omega^2` weights or an equivalent density
field. The weighted estimator can reduce bias when the supplied measure matches
the validation profile, but the measure is not inferred from causal order.

## Relation To Measure Encoding And Coarse Graining

Milestone 10 keeps the same finite-sampling perspective but changes where the
measure enters. Instead of attaching weights to coordinate-uniform events, it
sprinkles events uniformly with respect to conformal physical volume. Then the
coordinate event density is a statistical representation of the supplied
measure model.

The relevant density scale is still supplied. For constant conformal rescaling,
normalized event positions are unchanged, so an absolute scale cannot be
recovered from positions alone. For random thinning, the expected density is
multiplied by the keep probability. Correct density rescaling should preserve
volume estimates up to increased sampling noise; uncorrected density is a
deliberate failure baseline.

## Relation To Order-First Diagnostics

Milestone 11 separates ordinal error from metric-value error. An order
inversion/error rate asks whether pairwise comparisons are preserved, not
whether a reconstructed value has small RMSE in seconds or meters. This is
useful for the order-first thesis because metric values are treated as
representation-layer objects.

Finite tick resolution can preserve coarse order only up to ties. Calibration
and repeated processes are needed before order data supports stable ratios or
metric units. These diagnostics complement, rather than replace, the sampling
noise models used for interval-count reconstruction.

## Relation To Ordinal Embedding

Milestone 12 introduces another finite-error layer: ordinal embedding stress.
Constraint violation rate and hinge loss quantify how well a low-dimensional
coordinate representation satisfies sampled distance-order comparisons. These
are not Poisson interval-count errors, but they play a similar diagnostic role:
they separate finite sampling, incomplete comparisons, noisy constraints, and
representability limits from stronger metric claims.

## Relation To Held-Out And Null-Model Validation

Milestone 13 adds generalization and baseline diagnostics for ordinal
embedding. Training violation alone can be misleading because a finite optimizer
may partially fit random or inconsistent constraints. Held-out order validation
checks whether fitted coordinates satisfy independent constraints from the same
order source.

Null models such as shuffled-side constraints, fully random quadruplets, and
controlled flips estimate how much apparent success can come from optimizer
flexibility rather than structured order data. Subset-stability diagnostics ask
whether independent constraint samples lead to similar embeddings. These are
finite validation tools, not proofs of metric emergence.

## Relation To Simultaneity-Sliced Distance Order

Milestone 14 adds a slice-selection error layer. Spatial distance-order
comparisons are evaluated within radar-time bins derived from observer tick
brackets. Narrow bins can reduce time mixing but may produce too few pair
comparisons; wide bins increase pair counts but may mix observer-time slices.

The relevant uncertainties are finite tick resolution, finite event sampling,
and the chosen bin width. These are representation diagnostics, not evidence
that spatial distance is observer-independent or primitive.
