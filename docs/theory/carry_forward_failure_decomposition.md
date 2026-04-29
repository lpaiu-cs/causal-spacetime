# Carry-Forward Failure Decomposition

Milestone 35 decomposes failure modes after the carry-forward stop condition.
Milestone 34 found no current manifest family eligible for future stress tests,
so the next conservative step is a carry-forward failure decomposition and
diagnostic completeness audit.

This is not threshold retuning. It is not a stress test. It does not fit new
representation models, change the Milestone 34 carry-forward registry, or
reinterpret blocked families as eligible.

## Failure Categories

A measured failure is an available diagnostic value that fails a fixed
criterion. Examples include a held-out violation above the preregistered
threshold, a generalization gap above threshold, or a null-taxonomy gap below
threshold.

A missing diagnostic is a required value that was not present in the
cross-family output bundle. Missing metrics are not success and not the same
as measured failure. They are reported separately as missing-metric warnings
or blocking missing metrics, depending on the criterion.

Hard criteria can block carry-forward on their own. Non-hard criteria can
contribute to provisional or blocked status. Warning criteria identify
diagnostic incompleteness or weak evidence without being treated as direct
measured failures. Accounting criteria ensure failed-control and ineligible
families remain visible.

## Stop Condition

If no family is `carry_forward` and no family is explicitly `provisional`,
future stress tests are not allowed. The correct mode is report-only
continuation: decompose the decision, audit diagnostic completeness, and
preserve all blocked, report-only, and failed-control rows.

Threshold sensitivity is not retuning. It is a fixed report showing how
decisions would vary under declared variants; it is not used to choose a new
threshold after seeing failure.

## Root-Cause Categories

Milestone 35 reports root-cause categories such as heldout failure,
generalization-gap failure, null-separation failure, coverage failure, restart
instability, latent-order instability, missing output, missing metric,
accounting control, and ineligible control.

These categories explain the failure mode of the diagnostic pipeline. They do
not say whether any physical theory is true or false.

## Upstream Remediation Boundary

An upstream remediation design is a future proposal for improving manifest
generation or diagnostic reporting. It is not execution and not post-hoc
rescue. Future remediation must be preregistered before it changes data
generation, target selection, comparison protocols, thresholds, or manifest
families.

Pure reporting fixes, such as computing missing coverage metrics in a future
run, can be identified here. They still do not change the Milestone 34
decision.

Milestone 36 produces a preregistered remediation plan from this decomposition.
The plan is not executed in Milestone 36. It does not retune thresholds. It
does not run stress tests. It does not fit new representation models. Future
execution requires a new preregistered run.
