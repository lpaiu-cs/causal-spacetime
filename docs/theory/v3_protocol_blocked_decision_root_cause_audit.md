# V3 Protocol Blocked-Decision Root-Cause Audit

Milestone 43 audits why v3 protocol families remain blocked. It follows from
Milestone 42 because no v3 protocol family received `carry_forward` or
`provisional` status. The M42 decision layer applied the fixed M34 criteria to
the M41 protocol-invariant, provenance-aware v3 bundle.

Milestone 43 does not change M42 decisions. It does not retune thresholds. It
does not run stress tests. It does not execute v4. V4 design is planned-only and
requires a later execution milestone.

The audit separates measured metric failures from precondition failures and
control-family outcomes. In M42 the required bundle inputs were present and the
v3 protocol/profile/provenance preconditions passed, so the structured-family
blocking is analyzed primarily as measured metric failure.

The root-cause audit keeps blocked families visible. It reports held-out
failure, generalization-gap failure, null taxonomy weakness, symmetry-control
failure, restart instability, latent-order instability, coverage weakness, and
family/protocol-specific weaknesses. It also records provenance-mode drilldown
and manifest-level variance without treating top-down or hybrid provenance as
evidence.

Stress-test design remains stopped unless a later milestone produces at least
one carry-forward or provisional family under the fixed criteria.

