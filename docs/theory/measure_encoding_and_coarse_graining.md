# Measure Encoding And Coarse Graining

Milestone 9 used supplied `Omega^2` weights on coordinate-uniform support
events. Milestone 10 tests a complementary representation: events are sprinkled
uniformly with respect to conformal physical volume, so the counting measure
itself carries the supplied measure model.

Milestone 11 reinterprets this as a representation-layer issue. Counting
measure and density help stabilize metric-volume representations, but they are
not part of bare causal order or bare observer-relative distance order.

## Physical-Volume Sprinkling

For a positive 1+1D conformal factor,

```text
ds_Omega^2 = Omega(t, x)^2 * (dt^2 - dx^2)
dV_Omega = Omega(t, x)^2 dt dx
```

the causal order is unchanged, but the coordinate density of events sampled
uniformly with respect to `dV_Omega` is proportional to `Omega^2`. In this
controlled model, unweighted interval counts can estimate conformal physical
volume if the physical density scale is supplied:

```text
V_hat_Omega = count / rho_physical
rho_physical = N / V_global_Omega
```

The sampling process is additional measure structure. It is not derived from
causal order alone.

## Local Relative Measure Shape

A nonconstant conformal profile changes the relative density of events in
coordinate space. Binning events by coordinate time gives a statistical
estimate of the local relative density shape. After normalization by the mean,
this can be compared to the normalized `Omega^2` profile.

This only estimates relative shape. A constant conformal scale multiplies all
physical volumes by the same factor while leaving normalized event positions
unchanged, so the absolute constant scale remains underdetermined without an
absolute density or volume convention.

## Thinning And Density Rescaling

Random thinning is a simple coarse-graining test. If events are retained with
probability `p`, the expected density changes as:

```text
rho_thinned = p * rho_original
```

Interval-volume reconstruction should remain approximately stable when the
density is rescaled correctly, with larger finite-sampling noise at lower keep
probability. Reusing the original density after thinning produces a predictable
scale error.

## Interpretation

These checks support the reconstruction program by clarifying how measure may
be represented:

```text
supplied weights on events
or
sampling/counting measure plus supplied density scale
```

Coarse-graining stability is necessary for a serious continuum reconstruction
program, but it is not sufficient. These controlled conformal toy models do not
derive the conformal factor, prove spacetime emergence, or replace general
relativity.

The order-first program therefore treats measure encoding as one condition
under which an effective metric representation may be useful. It does not imply
that event distribution alone determines an absolute metric scale.

Milestone 12 studies a different representability condition: whether
distance-order constraints can be compressed into low-dimensional coordinates.
Measure stability and ordinal embeddability are complementary requirements for
an effective metric representation.

Milestone 13 adds stability and null-model checks for that embeddability layer.
Held-out validation and subset stability help distinguish structured
low-dimensional order data from optimizer fits to shuffled, random, or noisy
constraints. These tests complement, but do not replace, measure and density
assumptions.

Milestone 14 adds simultaneity-sliced distance order. It addresses which
comparisons should count as spatial for a chosen observer protocol, while
measure encoding addresses how counting statistics and density represent
volume scale. Both are additional structure beyond causal order alone.
