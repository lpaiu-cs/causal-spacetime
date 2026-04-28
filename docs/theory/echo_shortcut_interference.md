# Echo Shortcut And Interference Classification

Milestone 24 classifies shortcut returns in controlled echo-response motifs.
The setting is still a finite state-change trigger DAG with a selected
reference chain. No metric coordinates, calibrated clock labels, or spatial
geometry are introduced.

## Return Spectrum

For a reference chain

```text
R = (r_0, r_1, ..., r_m)
```

and a fixed emission position `k`, a controlled motif supplies a planted delay
rank:

```text
D_plant(e; k) = d
```

This means the motif was inserted so that the intended return position is
`k + d`. The planted value is a validation label for the finite DAG.

The echo protocol recovers the earliest available return:

```text
D_echo(e; k) = min { j - k : e ≺_T r_j, j > k }
```

The full return spectrum is:

```text
S_R(e; k) = { j - k : e ≺_T r_j, j > k }
```

The spectrum records all later reference-chain positions reachable from the
target. It does not require that the target was reached from the emission; that
reachability is checked separately.

## Classifications

Exact recovery occurs when:

```text
D_echo = D_plant
```

A shortcut return occurs when:

```text
D_echo < D_plant
```

Missing or late recovery occurs when:

```text
D_echo is missing
```

or:

```text
D_echo > D_plant
```

The shortcut depth is the rank difference:

```text
D_plant - D_echo
```

It is a diagnostic of how much earlier the earliest return occurred relative
to the planted validation label.

## Interference Interpretation

Shortcut returns are causal-order interference. They arise because the
transitive order contains an earlier return path than the planted motif path.
The earliest-return rule is part of the same-emission echo protocol, so the
recovered delay rank is intentionally the first available return in the order.

Milestone 24 separates two stress tests:

- naturally occurring background shortcut paths from random or existing
  trigger structure,
- deliberately injected shortcut-return paths used to test robustness.

Those cases should not be conflated. The first asks how generic background
trigger structure alters motif visibility. The second deliberately inserts
earlier returns to measure exact-recovery degradation and rank-resolution loss.

Return-spectrum analysis is therefore a prerequisite before any stronger
distance-order interpretation is attempted. In this milestone, planted delay
rank and recovered delay rank remain order-level response-depth diagnostics.

## Coarse-Graining Stability

Milestone 25 applies coarse-graining to the same return-spectrum
classification. Closure-preserving event thinning should preserve spectra
among retained reference and target events. Immediate-edge thinning can remove
or alter return paths. Reference-chain subsampling can hide shortcuts or merge
fine ranks into coarse ties. These changes are classified as coarse-graining
artifacts of the finite protocol.
