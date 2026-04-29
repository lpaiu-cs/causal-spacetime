# Admissible Pairwise Response Protocols

Milestone 29 defines admissible pairwise response-profile comparison protocols
for finite response profiles. These protocols are chosen finite diagnostics,
not physical measurements.

## Protocol Definition

An admissible protocol specifies:

- a response-profile dissimilarity method,
- a missing-data policy,
- a minimum common-protocol requirement,
- whether rank gaps are normalized by protocol-column range,
- whether reachability mismatch is penalized.

The output is a pre-metric pairwise response comparison over target pairs.
It can produce response-comparison constraints, but those constraints are not
spatial distance-order constraints.

## Missing-Data Policies

`common_reachable` compares targets only on protocol columns where both are
reachable. It is conservative about rank gaps but can discard sparse pairs.

`require_all_reachable` keeps only pairs reachable in all protocol columns. It
is stricter and can sharply reduce valid-pair coverage.

`penalize_mismatch` keeps more pairs by adding a reachability mismatch term.
This changes the interpretation and must be declared before analysis.

## Null Baselines

Null baselines are required before any ordinal embedding attempt. They include:

- shuffled delay ranks within each protocol,
- shuffled reachability masks,
- independently permuted target profiles,
- random profiles with the same per-column marginals.

These baselines preserve selected summary structure while disrupting target
profile alignment. They help test whether a comparison protocol sees
nontrivial profile structure, not whether geometry has been recovered.

## Scope

No embedding is performed in Milestone 29. No calibrated time, spatial
distance, metric scale, speed, or geometry is inferred.
