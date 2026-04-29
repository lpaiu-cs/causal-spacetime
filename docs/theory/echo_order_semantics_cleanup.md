# Echo-Order Semantics Cleanup

Milestone 27 fixes terminology before any representability diagnostics are
applied to echo-response order signatures. The correction is semantic, not a
new physical model.

## Passive Brackets Versus Fixed-Emission Echo

Passive reference-chain brackets use:

```text
p_R(e) = max { i : r_i ≺_T e }
q_R(e) = min { j : e ≺_T r_j }
```

The passive bracket-width rank is:

```text
W_rank(e) = q_R(e) - p_R(e)
```

`W_rank` belongs to the passive bracket diagnostic. It is not the same object
as a fixed-emission echo-delay rank.

## Deprecated Terminology

`R_rank` is deprecated in the state-change echo layer because it conflated
passive bracket-width diagnostics with fixed-emission echo diagnostics. Use
`W_rank` only for passive bracket width and `D_echo` only for fixed-emission
echo-delay rank.

For a chosen emission position `k`, the fixed-emission echo protocol uses:

```text
D_echo(e; k) = min S_full(e; k)
```

where `D_echo` is the earliest-return delay rank. It is not distance,
duration, speed, or a calibrated time label.

## Return Spectrum Types

The full transitive return spectrum is:

```text
S_full(e; k) = { j - k : e ≺_T r_j, j > k }
```

When the full transitive closure and full reference chain are used, this
spectrum is generally a suffix from the earliest return through later
reference positions, because the reference chain is itself ordered.

The retained reference spectrum is:

```text
S_retained(e; k)
```

It reports only returns to retained or subsampled reference ticks. It can be
sparse because intermediate reference positions may not be present.

The immediate-edge return spectrum is:

```text
S_immediate(e; k)
```

It reports direct immediate trigger edges from the target to later reference
events. It can also be sparse and should not be confused with the full
transitive spectrum.

## Shortcuts And Gated Protocols

A shortcut return occurs when:

```text
D_echo(e; k) < D_plant(e; k)
```

This is an earlier causal return selected by the earliest-return rule. It is
not a bug and not metric noise.

A gated echo rule may instead choose the first return at or after a
predeclared gate:

```text
min { d in S_full(e; k) : d >= gate }
```

That is a separate protocol selected before evaluation. It is not a post-hoc
adjustment after seeing shortcuts.

## Scope

Stable response-order signatures are ordinal diagnostics. They are not
distance orders, spatial representations, or metric reconstructions. Milestone
27 only clarifies which rank is being used and which return spectrum is being
queried.

## Milestone 28 Follow-Up

Milestone 28 uses these cleaned semantics to show that a single-reference
response signature remains a scalar preorder. Even when `D_echo` values are
scalar-rank representable, they do not supply target-target pair comparisons.
That additional structure belongs to later levels of the representability
ladder.
