# V3 Measured Failure Taxonomy

Milestone 43 uses the following measured-failure taxonomy for the v3 protocol
blocked-decision audit.

`heldout_failure`: mean held-out violation above the fixed threshold.

`generalization_gap_failure`: held-out minus train gap above the fixed
threshold.

`destructive_null_gap_failure`: destructive nulls are not sufficiently worse
than the structured fit.

`symmetry_control_gap_failure`: symmetry-control behavior differs too much from
structured behavior.

`restart_instability`: held-out violation varies too much across optimizer
restarts.

`latent_order_instability`: fitted latent pair-order varies too much across
restarts.

`profile_resolution_failure`: the response profile has too many ties, too
little strict order, or weak pair resolution.

`reachability_failure`: common reachability or pair-node coverage is weak.

`provenance_precondition_failure`: metadata or provenance is invalid. In M42
this category should be absent for v3 structured families, but it remains a
separate category for audit completeness.

`control_family_blocking`: a report-only or failed-control family is
intentionally not eligible.

Diagnostic completeness is not success. Protocol invariance is not success.
Provenance integrity is not success. Carry-forward failure is not a proof that
the research hypothesis is false.

## Milestone 44 Boundary

Milestone 44 executes the preregistered v4 protocol manifest-generation run
after this measured-failure audit. It does not evaluate carry-forward
decisions. It does not run stress tests. It does not retune thresholds. It does
not infer metric geometry. Diagnostic-complete v4 output is not carry-forward
success.
