# Echo Return-Spectrum Coarse-Graining

Milestone 25 studies return-spectrum stability under controlled
coarse-graining and subsampling of finite state-change trigger networks. The
goal is to separate stable order-level echo diagnostics from artifacts of event
resolution, edge resolution, or reference-chain resolution.

## Return-Spectrum Stability

For a motif target `e` and reference emission position `k`, the return
spectrum is the set of later reference-rank returns:

```text
S_R(e; k) = { j - k : e ≺_T r_j, j > k }
```

Return-spectrum stability asks whether this set, its earliest return, and its
shortcut classification survive a controlled reduction of the finite
description.

The diagnostics include:

- spectrum Jaccard agreement,
- earliest-return stability,
- classification stability,
- rank-resolution loss,
- tied target-pair fraction.

## Closure-Preserving Event Thinning

Closure-preserving event thinning retains selected events and restricts the
original transitive closure to those retained events. Hidden intermediate
events are integrated out. Causal reachability among retained events is
therefore preserved by construction.

This tests whether motif targets and reference events can keep the same
return-spectrum diagnostics after unobserved intermediate events are removed
from the representation.

## Immediate-Edge Thinning

Immediate-edge thinning removes selected trigger edges before recomputing the
transitive closure. This can delete causal path information. A motif that was
previously reachable can become missing or delayed if its trigger path is not
protected.

This is different from closure-preserving event thinning. One hides
intermediate events while preserving their reachability effect; the other
changes the immediate trigger graph.

## Reference-Chain Subsampling

Reference-chain subsampling keeps only selected positions on the reference
chain. The resulting coarse reference rank can merge distinct return ranks into
the same bin or shift the expected coarse return to the next retained
reference position.

This reduces rank resolution and can create ties. It is not a clock
calibration procedure.

## Shortcut Classification

Shortcut status can change under coarse reference protocols. A shortcut that is
visible at a fine reference resolution may be hidden, merged with the planted
return, or shifted to a different coarse rank after reference-chain
subsampling.

Milestone 25 therefore treats shortcut classifications as protocol- and
coarse-graining-dependent finite diagnostics. No metric distance, calibrated
time, speed model, or spatial geometry is inferred.
