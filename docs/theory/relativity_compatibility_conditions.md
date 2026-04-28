# Relativity Compatibility Conditions

Milestone 18 does not make Lorentzian geometry primitive. It treats
Lorentz/Poincare-compatible observer atlases as target effective structures.

## Structures To Avoid

The state-change order program must avoid adding:

- global discrete update slices,
- preferred global time,
- regular lattice preferred frames,
- absolute rest frames,
- absolute discrete space cells.

Local finiteness is compatible with an order-theoretic program only if it is
not implemented as a regular global lattice with a preferred frame.

## Observer-Relative Decomposition

Observer time is order along an observer chain. Radar brackets and radar-time
bins define observer-relative slices. Distance order is then interpreted within
those slices. This is not a preferred global slicing and not absolute space.

Milestone 20 adds finite diagnostics for selecting reference-chain candidates.
Milestone 21 adds order-level brackets from those selected reference chains.
These diagnostics can compare local-system chains, order-only chains, and
random baselines, but they do not define calibrated time, metric distances, or
a preferred observer frame.

## Effective Atlas Condition

Given observer protocols `O` and `O'`, reconstructed charts should satisfy, on
overlaps:

```text
phi_O' o phi_O^{-1} ~= Lorentz/Poincare map
```

This is an effective representation condition. It is not a primitive axiom and
not an assumed background. The controlled atlas experiments test whether
overlap maps become approximately Lorentz/Poincare-compatible under supplied
observer protocols, clocks, orientation references, and calibration.

The desired outcome is a Lorentz-compatible observer atlas as an emergent
representation layer, not a claim that the primitive order itself already
contains metric coordinates.
