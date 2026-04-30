# V3 Protocol Metadata In Manifests

Every production v3 handoff manifest must include measurement protocol metadata,
response profile metadata, and handoff provenance metadata. The serialized
manifest records `measurement_protocol_id`, `measurement_protocol_hash`,
`reference_set_id`, `reference_chain_ids`, profile invariance status,
admissibility for pairwise response-profile dissimilarity, provenance type,
design digest, and any top-down template identifiers.

Production pairwise response-profile dissimilarity must use only
protocol-invariant profiles. Protocol-mixed and underspecified profiles are not
production-admissible. Report-only controls may carry explicit metadata while
remaining ineligible for production pairwise dissimilarity and later
carry-forward decisions.

Protocol metadata, profile metadata, and provenance metadata are admissibility
metadata, not physical interpretation. Metadata does not imply geometry, metric
scale, or spatial coordinates.

## Milestone 42 Boundary

Milestone 42 applies fixed carry-forward criteria to the M41 v3 protocol bundle.
Protocol metadata, profile metadata, and provenance metadata are preconditions
for admissibility. They do not imply carry-forward. M42 does not regenerate
manifests, rerun fits, retune thresholds, or run stress tests. Carry-forward is
stress-test eligibility, not geometry.
