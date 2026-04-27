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

In the order-first reformulation, conformal ambiguity is a statement about the
gap between primitive causal order and metric representation. Causal order can
remain fixed while metric tensors, duration values, and volume values change.

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

## Measure Encoding By Event Distribution

Milestone 10 tests a complementary encoding. Support events are sprinkled
uniformly with respect to conformal physical volume, so coordinate event
density is proportional to the supplied `Omega^2` profile. In that setting,
unweighted counts estimate conformal physical volume only if the correct
physical density scale is supplied.

For a constant conformal scale, the normalized coordinate distribution is the
same as the flat profile. This is a simple reminder that absolute constant
scale remains underdetermined without an external density or volume convention.

Random thinning is treated as a coarse-graining check. If each event is kept
with probability `p`, density must be rescaled by `p`; otherwise interval
volumes are systematically underestimated.

Milestone 11 adds the related point that even ratios are representation-layer
objects. Positive monotone transformations can preserve order while changing
ratio structure, so metric scale and ratio stability require additional
calibration or dynamics.

Milestone 12 uses ordinal embedding to test another representation condition:
even after order data is supplied, a low-dimensional metric representation must
earn its usefulness by satisfying many distance-order comparisons with low
loss. This does not remove conformal or measure ambiguity.

## Interpretation

This milestone strengthens the reconstruction program by making a limitation
explicit:

```text
causal/accessibility order -> conformal/light-cone structure
causal/accessibility order + measure/density -> volume-scaled reconstruction
```

It is a controlled conformal toy model, not a general-relativistic dynamics
simulation and not proof of spacetime emergence.
