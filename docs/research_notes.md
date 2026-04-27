# Research Notes

## Milestone 1: Timelike Reconstruction In 1+1D

This milestone works in 1+1D Minkowski spacetime with natural units `c = 1`.
Events are represented as `(t, x)`.

The causal relation is:

```text
p precedes q iff dt > 0 and dt^2 >= dx^2
```

where `dt = q.t - p.t` and `dx = q.x - p.x`.

The experiment samples events inside the causal diamond bounded by:

```text
p = (-T/2, 0)
q = ( T/2, 0)
```

The diamond condition is:

```text
abs(x) <= T/2 - abs(t)
```

For this 1+1D interval, the continuum volume is:

```text
V = tau^2 / 2
```

so an event-density estimate `rho` and an interval count `N` give:

```text
tau_hat = sqrt((N / rho) / (1/2))
```

The density is not derived from causal order alone. It is extra metric-scale
structure supplied to the reconstruction procedure.

## Interpretation

The current experiment is a controlled check that causal interval cardinality
and longest chains behave consistently in a known spacetime. It is not a claim
that causal order alone fixes all metric quantities, nor that finite signal
speed by itself implies Lorentzian spacetime.

## Milestone 2: Internal Timelike Pair Validation

The original full-diamond endpoint experiment is useful as a sanity check, but
it can be close to tautological if the same full-diamond event count is used to
set the density and then reconstruct the full-diamond proper time.

The milestone 2 validation experiment instead:

```text
1. generates one larger 1+1D Minkowski causal diamond,
2. treats global event density as supplied scale information,
3. samples many internal timelike pairs using the causal order,
4. counts each pair's Alexandrov interval cardinality,
5. compares reconstructed proper time to the hidden Minkowski value.
```

For an internal pair `(i, j)`, causal order contributes the set of sampled
events `k` satisfying:

```text
i precedes k and k precedes j
```

Event density contributes the scale needed to convert that count into a volume
estimate. In 1+1D, the operational reconstruction being tested is:

```text
tau_hat = sqrt((interval_count / rho) / (1/2))
```

Longest-chain estimates are also reported with the asymptotic 1+1D relation:

```text
L ~ sqrt(2 * rho) * tau
```

The implementation marks this normalization as finite-size sensitive. These
experiments are validation in known 1+1D Minkowski spacetime, not a proof of a
new spacetime theory.

## Milestone 3: Probe-Pair Statistical Calibration

Milestone 3 separates support-event sampling from endpoint selection. The
support events are sprinkled in a known 1+1D causal diamond, while probe
timelike endpoint pairs are sampled independently and are not inserted into the
support set.

For each probe pair, the code counts support events in the probe pair's
Alexandrov interval and uses:

```text
tau_hat = sqrt(2K / rho)
```

where `rho` is supplied from the known global support density. The hidden probe
coordinates are used to validate the result and to compute expected
finite-sampling error scales.

For Poisson support-event fluctuations:

```text
E[K] = rho * tau^2 / 2
std(tau_hat) ~= 1 / sqrt(2 rho)
```

For the fixed-count sprinkling used by the code, the interval count is better
approximated as binomial with:

```text
p = V_interval / V_global
K ~ Binomial(N_support, p)
```

The probe-pair experiment compares observed RMSE against both Poisson and
fixed-N binomial delta-method predictions. This helps distinguish ordinary
finite-sampling noise from possible bias, boundary effects, estimator mistakes,
or convention errors. It does not imply that causal order alone determines
metric scale.

## Milestone 4: Dimension As An Order-Statistical Observable

Milestone 4 adds D-dimensional flat Alexandrov intervals and the
Myrheim-Meyer dimension estimator. The validation question is whether dimension
can be estimated from causal-order relation statistics in controlled known
intervals.

For a finite causal matrix `C`, the ordered relation fraction is:

```text
r_obs = (# ordered relations) / (N * (N - 1))
```

The Myrheim-Meyer curve used here is:

```text
r(D) = Gamma(D + 1) * Gamma(D / 2) / (4 * Gamma(3D / 2))
```

The experiment inverts this curve by bisection to estimate `D`. This tests
dimension as an order-statistical observable. It does not prove spacetime
emergence, and it still assumes the events were generated inside known flat
Minkowski Alexandrov intervals.

The theory-facing hypothesis is a mathematical reconstruction program:

```text
primitive causal/information-accessibility structure
  + counting measure / event density
  + observer protocol
  -> operational time, distance, dimension, and spacetime-like geometry
```

The project separates established mathematics, implemented numerical models,
reconstruction hypotheses, speculation, and open proof obligations in
`docs/theory/`.

## Milestone 5: Observer-Chain Radar Reconstruction

Milestone 5 moves radar decomposition away from direct coordinate formulas. A
stationary observer chain with clock labels is supplied, and the reconstruction
uses only:

```text
causal matrix + observer chain indices + observer clock labels
```

For a target event, the protocol identifies:

```text
tau_minus = latest observer clock tick preceding the target
tau_plus  = earliest observer clock tick following the target
```

Then:

```text
radar_time = (tau_plus + tau_minus) / 2
radar_distance = (tau_plus - tau_minus) / 2
```

Hidden coordinates are used only for validation against the known stationary
observer values `t` and `abs(x)`. This supports the observer-dependent spatial
decomposition layer of the reconstruction program, while keeping the observer
chain and clock labels explicit as additional protocol structure.

## Radar Coordinates

For the stationary observer worldline `O(tau) = (tau, 0)`, the radar
decomposition of an event `(t, x)` is:

```text
tau_minus = t - abs(x)
tau_plus  = t + abs(x)
radar_time = (tau_plus + tau_minus) / 2
radar_distance = (tau_plus - tau_minus) / 2
```

For an inertial observer moving with velocity `v`, the event is first Lorentz
boosted into that observer's rest frame:

```text
t' = gamma * (t - v x)
x' = gamma * (x - v t)
```

The same stationary radar decomposition is then applied to `(t', x')`. This is
an operational coordinate protocol relative to a chosen observer.

## Lorentz Length Contraction

For a rod at rest in `S'` with endpoints `x'=0` and `x'=L0`, the transform to a
lab frame `S` where the rod moves at velocity `beta` is:

```text
t = gamma * (t' + beta x')
x = gamma * (x' + beta t')
```

A lab-frame length measurement compares endpoint events with the same lab time.
At lab time `t=0`, the selected rest-frame endpoint events are:

```text
left:  t'=0,        x'=0
right: t'=-beta L0, x'=L0
```

Transforming these two events gives a lab-frame separation:

```text
L = L0 / gamma
```

This is the usual special-relativistic length contraction as an event-selection
result, not a claim about a new spacetime model.

## Finite-Speed Lattice Counterexample

A regular 1+1D lattice graph can impose a finite propagation rule:

```text
(t, x) -> (t + 1, x - 1)
(t, x) -> (t + 1, x + 1)
```

This graph has a finite causal cone, but it also has preferred lattice
directions, a preferred time slicing, and a parity structure. Comparing the
lattice cone with a Poisson sprinkle in a continuum causal cone is a simple
counterexample to the claim that finite signal speed alone implies
Lorentz-invariant spacetime behavior.

The comparison is intentionally modest. It shows that extra structure matters;
it does not claim that a Poisson sprinkle is a complete physical model.

## Exploratory Spacelike-Distance Proxies

Spacelike distance is harder to reconstruct from causal order than timelike
proper time in this toy setting. A timelike interval has a direct Alexandrov
interval cardinality, while a spacelike pair has no causal interval between the
two events themselves.

The current exploratory implementation records three simple proxy counts for
spacelike-separated event pairs:

```text
common future overlap count
common past overlap count
minimal enclosing Alexandrov interval count
```

These quantities are compared against the true 1+1D Minkowski spacelike
separation:

```text
sqrt(dx^2 - dt^2)
```

The proxies are strongly affected by the finite sampling region, event density,
and boundary placement. They are useful diagnostics for failure modes and
possible correlations, but they are not validated estimators of spacelike
distance.
