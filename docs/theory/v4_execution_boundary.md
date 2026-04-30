# V4 Execution Boundary

Milestone 43 does not execute v4. It creates a no-execution remediation plan and
exports a planned-only v4 preregistration.

A future Milestone 44 may execute v4 only from the Milestone 43 preregistered
spec. That future execution would need to generate production manifests under
`outputs/manifests_v4/` and produce a diagnostic-complete v4 bundle without
applying carry-forward decisions in the same milestone.

Milestone 45 should apply the fixed carry-forward criteria to v4 outputs if and
only if a future v4 generation milestone has produced the required bundle.

Stress tests remain impossible unless a later carry-forward or provisional
family exists. Blocked v3 families remain blocked under the M42 decision layer.

Milestone 43 does not retune thresholds. It does not run stress tests. It does
not rerun v3 fits. It does not infer metric geometry.
