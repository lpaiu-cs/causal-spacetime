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

