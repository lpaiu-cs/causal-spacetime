# Handoff Manifest Provenance

Milestone 41 supports bottom-up, top-down, and hybrid handoff provenance. The
allowed `handoff_provenance_type` values are:

- `bottom_up_profile_derived`: constraints are generated from a concrete
  response profile by a fixed algorithm.
- `top_down_preregistered_template`: the handoff is specified by a
  preregistered theoretical or protocol template before response-profile
  instantiation.
- `hybrid_template_instantiated_from_profile`: a preregistered top-down template
  determines family, protocol, admissible targets, comparison method, margin
  policy, and constraint schema, while concrete constraints are instantiated
  from generated response-profile data.
- `report_only_control`: the manifest exists for accounting or failed-control
  tracking and is not eligible for production pairwise dissimilarity or
  carry-forward.

Top-down handoff is allowed only when preregistered and frozen before
evaluation. A top-down handoff is not evidence by itself and is not evidence of
geometry. It must record a design source, design digest, frozen milestone,
allowed data dependencies, forbidden data dependencies, constraint-selection
policy, template identifier, and evaluation boundary.

M41 supports bottom-up, top-down, and hybrid handoff provenance. M41 does not
evaluate carry-forward decisions. M41 does not run stress tests. M41 does not
retune thresholds. M41 does not infer metric geometry.

## Milestone 42 Boundary

Milestone 42 applies fixed carry-forward criteria to the M41 v3 protocol bundle.
Top-down/hybrid provenance is allowed only as preregistered manifest provenance,
not as evidence. M42 does not regenerate manifests, rerun fits, retune
thresholds, or run stress tests. Carry-forward is stress-test eligibility, not
geometry.
