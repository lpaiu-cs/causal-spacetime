# Preregistered Upstream Remediation Plan

Milestone 36 preregisters a remediation plan. It follows from the Milestone 35
finding that no manifest family currently passes carry-forward and that the
output bundle is diagnostically incomplete.

The plan is not executed in Milestone 36. It does not retune thresholds. It
does not run stress tests. It does not fit new representation models. Future
execution requires a new preregistered run.

## Remediation Types

Measured-failure remediation targets diagnostics that were present and failed
fixed criteria, such as held-out violation or generalization gap. These actions
are upstream design proposals for future profile or manifest generation.

Missing-metric remediation targets diagnostics that were unavailable in the
current family-level output bundle. Missing metrics are not treated as
successful, and they are not the same as measured failures.

Reporting-completeness remediation targets absent or incomplete output files,
schema rows, and audit rows. It improves future observability but does not
change any current decision.

Future data-generation remediation proposes new manifest-family designs. Any
proposal that changes data generation, target inclusion, protocol columns,
constraint sampling, or manifest-family selection requires preregistration
before execution.

## Execution Boundary

No stress tests are allowed while no `carry_forward` or explicitly
`provisional` family exists. Planned v2 families are not current results. They
are future candidate specifications that must still report failed and
ineligible families if executed later.

Diagnostic-complete manifest design means that required metrics are
predeclared and produced. It does not imply metric structure or physical
interpretation.

## Milestone 37 Boundary

Milestone 37 executes the preregistered v2 manifest-generation run derived
from this plan. The plan itself remains the preregistration artifact; M37
does not rewrite it as a result. It produces diagnostic-complete v2 outputs,
but it does not evaluate carry-forward decisions, does not run stress tests,
and does not fit any model for the purpose of improving a criterion.
