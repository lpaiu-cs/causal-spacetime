# Conformal Ambiguity And Measure

## Causal Order And Conformal Structure

In a Lorentzian continuum, positive conformal rescaling preserves null cones:

```text
ds_Omega^2 = Omega(t, x)^2 * (dt^2 - dx^2)
```

for `Omega(t, x) > 0`. The causal order is unchanged because multiplying the
metric by a positive factor does not change which directions are timelike,
null, or spacelike.

This means causal order can encode light-cone or conformal structure, but it
does not by itself fix metric scale, physical volume, or clock rates.

## Measure As Additional Structure

In 1+1D, the conformally rescaled volume element is:

```text
dV_Omega = Omega(t, x)^2 dt dx
```

For a clock at fixed `x = 0`, the proper-time element is:

```text
d tau_Omega = Omega(t, 0) dt
```

Therefore two models with the same causal order can assign different physical
volumes and different clock scales. A measure, density, or conformal weight
field is additional structure.

## Simplest Non-Identifiability Example

For constant `Omega = alpha`, causal order is unchanged, while:

```text
proper time scales as alpha
volume scales as alpha^2
```

No causal-order statistic alone distinguishes `alpha = 1` from `alpha = 2`.

## Weighted Volume Reconstruction

Milestone 9 uses support events sprinkled uniformly in coordinate volume. The
unweighted count estimates coordinate volume. A supplied conformal measure
weight estimates conformal physical volume:

```text
V_hat_Omega = (1 / rho_coordinate) * sum_{events in interval} Omega(event)^2
```

The weights are supplied measure information. The experiment does not derive
the conformal factor from causal order.

## Interpretation

This milestone strengthens the reconstruction program by making a limitation
explicit:

```text
causal/accessibility order -> conformal/light-cone structure
causal/accessibility order + measure/density -> volume-scaled reconstruction
```

It is a controlled conformal toy model, not a general-relativistic dynamics
simulation and not proof of spacetime emergence.
