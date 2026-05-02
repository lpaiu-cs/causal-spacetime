# V4 Execution Boundary

Milestone 43 does not execute v4. It creates a no-execution remediation plan and
exports a planned-only v4 preregistration.

Milestone 44 executes v4 only from the Milestone 43 preregistered spec. That
execution generates production manifests under `outputs/manifests_v4/` and
produces a diagnostic-complete v4 bundle without applying carry-forward
decisions in the same milestone.

Milestone 45 should apply the fixed carry-forward criteria to v4 outputs if and
only if a future v4 generation milestone has produced the required bundle.

Stress tests remain impossible unless a later carry-forward or provisional
family exists. Blocked v3 families remain blocked under the M42 decision layer.

Milestone 43 does not retune thresholds. It does not run stress tests. It does
not rerun v3 fits. It does not infer metric geometry.

Milestone 44 executes the preregistered v4 protocol manifest-generation run. It
does not evaluate carry-forward decisions. It does not run stress tests. It does
not retune thresholds. It does not infer metric geometry. Diagnostic-complete v4
output is not carry-forward success.

Milestone 45 applies fixed carry-forward criteria to the M44 v4 protocol bundle.
It does not regenerate manifests. It does not rerun fits. It does not retune
thresholds. It does not run stress tests. Carry-forward is stress-test
eligibility, not geometry. Diagnostic-complete v4 output is not carry-forward
success. Top-down/hybrid provenance is admissibility metadata, not evidence.
