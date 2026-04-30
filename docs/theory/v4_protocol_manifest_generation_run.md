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

