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

