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
