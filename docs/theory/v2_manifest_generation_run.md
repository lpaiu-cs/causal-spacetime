# V2 Manifest Generation Run

Milestone 37 executes the preregistered v2 manifest-generation run described
by the Milestone 36 remediation plan. The planned v2 families become
production v2 handoff manifest families, and their JSON artifacts are written
under `outputs/manifests_v2/` rather than replacing the earlier manifest
directory.

The purpose is diagnostic-complete output production. Milestone 37 executes
the v2 manifest-generation run, computes the required v2 diagnostics, and
exports a bundle that can be read by a future carry-forward evaluation. It
does not evaluate carry-forward decisions. It does not run stress tests. It
does not infer metric geometry.

The generated handoff manifests contain response-comparison constraints,
validation summaries, target identifiers, margins, and frozen train/held-out
constraint splits. They do not contain latent fitted variables. The fitted
diagnostics used to populate the v2 output bundle are finite diagnostic
measurements, not physical spatial coordinates.

Carry-forward evaluation is deferred to Milestone 38. Until a later
carry-forward or provisional status exists, future stress tests remain
forbidden by the stop-rule logic from Milestones 34 and 35.

Diagnostic-complete does not mean successful. A complete output bundle means
that the required metric columns are present for later evaluation; it does
not mean that any family passes future criteria.

Milestone 38 applies fixed criteria to v2 outputs. It does not regenerate
manifests, rerun fits, retune thresholds, or run stress tests.
Carry-forward is stress-test eligibility, not geometry.

Milestone 39 reads the v2 decision outputs as blocked results. It audits root
causes and preregisters v3 design work, but it does not change the M37 v2
generation outputs.

Milestone 40 audits response-profile protocol invariance before v3 execution.
It does not change v2 generation outputs, M38 decisions, or M39 blocking
analysis. It does not generate v3 manifests. It does not run stress tests. It
does not fit new representation models.
