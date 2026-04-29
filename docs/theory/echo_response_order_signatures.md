# Echo Response Order Signatures

Milestone 26 moves from individual motif recovery to response-order
signatures over populations of controlled echo-response motif targets.
Individual recovered echo-delay ranks are useful, but they do not by
themselves show whether a shared reference protocol supports a stable ordinal
response structure.

## Response-Order Signature

For motif targets:

```text
M = {e_1, ..., e_n}
```

and a fixed reference-chain emission protocol, let:

```text
D_i = D_echo(e_i; k)
```

where `D_i` is the recovered echo-delay rank for target `e_i`, if one is
reachable. The response-order sign matrix is:

```text
Omega_ij = -1 if D_i < D_j
            0 if D_i = D_j or unresolved
            1 if D_i > D_j
```

This matrix is the echo-response order signature. The response-order signature
is an ordinal diagnostic, not a metric space.

## Ties And Rank-Resolution Loss

Ties occur when two targets have the same recovered echo-delay rank or when at
least one target is unresolved under the chosen protocol. Ties are not errors
by themselves. They record rank-resolution loss in the finite reference
protocol.

Reference-chain subsampling can merge distinct fine-grained ranks into the same
coarse reference rank. Shortcut returns can also change pairwise order by
moving one target to an earlier recovered return.

## Protocol Variants

For protocol variants `c in C`, each variant induces its own sign matrix:

```text
Omega^(c)
```

Variants may include closure-preserving event thinning, immediate-edge
thinning, reference-chain subsampling, or targeted shortcut injection. These
variants test stability of the ordinal response structure, not calibrated
distance or clock behavior.

## Stable Response-Order Core

A pair `(i, j)` is stable when `Omega^(c)_ij` agrees for enough variants `c`.
The stable response-order core is the set of pairwise signs that survive the
chosen protocol variants above a specified agreement threshold.

Stable response order is a candidate input for future representability tests,
not a proof of geometry. A later representability test could ask whether the
stable ordinal response structure admits a compact scalar or low-dimensional
representation. That would still not be metric reconstruction unless
calibration, scale, and stronger physical protocols are added.

## Scope

Milestone 26 does not interpret response-order signatures as spatial distance.
It does not reconstruct metric geometry, calibrated radar distance, speed,
velocity, or finite-speed spatial geometry. It only asks whether controlled
echo-response motif populations produce ordinal response structures that remain
stable under specified finite protocol variations.
