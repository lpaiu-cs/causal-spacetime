# V4 Protocol Manifest-Generation Run

Milestone 44 executes the preregistered v4 protocol manifest-generation run.
It reads `outputs/remediation/v4_protocol_preregistration_spec_m43.json` and
writes production v4 handoff manifests under `outputs/manifests_v4/`.

M44 retains the fixed M34 criteria for later evaluation. It retains
protocol-invariant profile construction, parameter-complete measurement
protocols, provenance-aware handoff manifests, failed-control accounting, and
report-only controls.

M44 is not a carry-forward decision layer. It does not evaluate carry-forward
decisions. It does not run stress tests. It does not retune thresholds. It does
not infer metric geometry. Diagnostic-complete v4 output is not carry-forward
success.

M45 is the future v4 carry-forward evaluation milestone.

## Milestone 45 Decision Layer

Milestone 45 applies fixed carry-forward criteria to the M44 v4 protocol bundle.
It does not regenerate manifests. It does not rerun fits. It does not retune
thresholds. It does not run stress tests. Carry-forward is stress-test
eligibility, not geometry. Diagnostic-complete v4 output is not carry-forward
success. Top-down/hybrid provenance is admissibility metadata, not evidence.

## Milestone 46 Follow-Up

Milestone 46 audits why v4 protocol families remain blocked. It does not change
M45 decisions. It does not retune thresholds. It does not run stress tests. It
does not execute v5. V5 design is planned-only. Repeated remediation is an
overfitting risk unless justified by root-cause analysis.
