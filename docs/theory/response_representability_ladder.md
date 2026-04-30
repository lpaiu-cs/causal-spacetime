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

Input: multiple reference chains inside one fixed measurement protocol.
Output: protocol-invariant multi-reference response profile.
Additional assumption: aligned target identities and explicit measurement
protocol metadata.
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

Milestone 29 fills in the transition from `multi_reference_response_profile`
to `pairwise_distance_order` by defining pre-metric pairwise response-profile
comparison protocols. These protocols do not themselves reach the
pairwise-distance-order rung; they only define response-comparison constraints
that future work may test against stronger assumptions.

Milestone 30 adds validation gates for those response-comparison constraint
pools. Held-out protocol agreement, bootstrap stability, null-baseline
separation, and coverage are necessary pre-embedding diagnostics. They still
do not move the project to ordinal embedding or calibrated metric
representation.

Milestone 31 records successful or failed validation outcomes in handoff
manifests. A handoff manifest is a frozen input specification for future
representability experiments, not a new rung of metric interpretation.

Milestone 32 adds latent ordinal representation diagnostics that consume
handoff manifests. This checks whether frozen response-comparison constraints
can be fit by finite latent variables, but the fitted coordinates are not
spatial coordinates and the exercise does not reconstruct metric geometry.

Milestone 33 adds preregistered family comparison across frozen manifests.

Milestone 40 audits response-profile protocol invariance before v3 execution.
This is a pre-execution design correction, not threshold retuning. It does not
change M34, M38, or M39 decisions. It does not generate v3 manifests. It does
not run stress tests. It does not fit new representation models. A single
response profile must not mix measurement protocols.

A response profile used for pairwise response-profile dissimilarity may vary
reference chains inside one fixed measurement protocol. If emission, gate,
echo rule, spectrum type, subsampling, normalization, missing policy, tie
policy, or margin policy varies, those variants must form separate profile
families or be explicitly marked exploratory/report-only.
This is still inside the latent ordinal representation diagnostic level: it
compares structured families, failed or ineligible controls, destructive
nulls, symmetry controls, and marginal-preserving nulls without retuning
thresholds. It does not add a calibrated metric level to the ladder.

Milestone 34 adds a carry-forward gate for future stress-test eligibility.
This remains a pre-stress-test administrative and diagnostic step; it does not
advance the ladder to calibrated geometry.
