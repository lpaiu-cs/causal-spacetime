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
