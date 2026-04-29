# Response Representability Ladder

Milestone 28 defines a representability ladder for response-order diagnostics.
Each level states its input, output, additional assumptions, and limits.

## Levels

1. `scalar_response_rank`

Input: response-order signature.
Output: topological scalar ranks.
Additional assumption: acyclicity.
Does not imply: target-target distance order or geometry.

2. `stable_scalar_response_core`

Input: multiple protocol variants.
Output: robust scalar order core.
Additional assumption: protocol stability threshold.
Does not imply: spatial distance.

3. `multi_reference_response_profile`

Input: multiple reference chains or emission positions.
Output: response-profile embedding candidate.
Additional assumption: aligned target identities and protocol labels.
Does not imply: metric space.

4. `pairwise_distance_order`

Input: explicit target-target pair comparisons.
Output: distance-order relation.
Additional assumption: a pairwise comparison protocol.
Does not imply: ratios or a metric tensor.

5. `ordinal_embedding_candidate`

Input: pairwise distance-order constraints.
Output: low-dimensional coordinates preserving many constraints.
Additional assumptions: stable constraints, held-out validation, and null
baselines.
Does not imply: calibrated metric geometry.

6. `calibrated_metric_representation`

Input: ordinal embedding plus calibration, measure, or atlas consistency.
Output: effective metric representation.
Additional assumptions: calibration, scale, consistency, and dynamics-like
constraints.
Does not imply: a fundamental metric field or GR dynamics.

## Interpretation

The ladder prevents a category error: a scalar response order is not enough for
metric interpretation. Multi-reference profiles are richer but still
pre-metric. Metric geometry requires calibration, measure, consistency
conditions, and validation beyond response order.
