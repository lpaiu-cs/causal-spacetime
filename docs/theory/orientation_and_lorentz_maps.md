# Orientation And Lorentz Maps

## Single-Observer Reflection Degeneracy

A single observer-chain radar protocol reconstructs radar time and an unsigned
radar distance. Relative to a stationary observer at `x = 0`, the events
`(t, x)` and `(t, -x)` have the same emission and reception clock readings.
They therefore have the same radar time and radar distance.

This degeneracy is expected. It is not a defect in the protocol, and it is not a
claim that causal order alone supplies a signed spatial coordinate.

## Oriented Observer Protocol

Signed space requires additional reference structure. Milestone 6 uses a
two-chain radar reconstruction protocol:

```text
primary observer chain: x' = 0
beacon chain:           x' = a
```

The two chains are comoving and synchronized in the protocol rest frame. The
beacon separation `a` is supplied protocol data. It is not reconstructed from
causal order alone in the current implementation.

For a target event, the protocol reconstructs the unsigned radar distances
`R0` and `Ra` to the primary and beacon chains. The signed coordinate relative
to the oriented protocol is then:

```text
x_hat = (R0^2 - Ra^2 + a^2) / (2a)
```

This is a controlled 1+1D construction. It tests whether a supplied orientation
reference can remove the reflection degeneracy.

## Lorentz-Map Recovery

Two oriented inertial protocols can assign two coordinate decompositions to the
same target events. In known 1+1D Minkowski validation data, the compatibility
map between a lab protocol and a moving protocol should be Lorentzian:

```text
t' = gamma * (t - beta x)
x' = gamma * (x - beta t)
```

Milestone 6 fits `beta` from reconstructed coordinates and compares the fitted
map with the known hidden boost. The hidden coordinates are used only for
validation and for reporting error.

Milestone 7 generalizes this to affine protocols with different origins. In
that setting, the expected transition map includes a translation:

```text
chart_B ~= L(beta_AB) chart_A + translation_AB
```

This is an observer-atlas consistency test, not a derivation of translations or
origins from causal order alone.

## Interpretation

This experiment supports the observer-dependent spatial decomposition layer of
the mathematical reconstruction program:

```text
observer chain + beacon orientation + causal accessibility
  -> operational signed coordinates
  -> compatibility map between observer protocols
```

It does not prove spacetime emergence. It does not show that causal order alone
gives signed space, and it does not derive Lorentz transformations without
additional protocol structure.

Milestone 8 does not extend the Lorentz-map fitting layer. Instead, it uses a
single accelerated observer protocol to study reconstruction-inaccessibility at
a Rindler horizon analogue in flat spacetime.

Milestone 9 addresses a different ambiguity: even when causal order and
observer protocols are fixed, positive conformal rescaling can change physical
clock scale and volume unless measure information is supplied.

Milestone 10 tests that measure information can also be represented by a
physical-volume sprinkling distribution and by density rescaling under
thinning. This does not change the need for explicit orientation/reference
protocols when reconstructing signed coordinates.
