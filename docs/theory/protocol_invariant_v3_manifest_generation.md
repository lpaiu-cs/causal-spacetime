# Protocol-Invariant V3 Manifest Generation

Milestone 41 executes the protocol-invariant patched v3 manifest-generation run
defined by the Milestone 40 patched preregistration. It reads the patched v3
family specifications, checks that structured families use a fixed measurement
protocol, and writes production v3 handoff manifests under
`outputs/manifests_v3/`.

Each structured v3 family is protocol-invariant. Within one response profile,
only `reference_chain_ids` vary. If echo rule, emission rule, gate rule,
reference-subsampling rule, spectrum type, normalization rule, missing-data
policy, tie policy, margin policy, or any declared rule parameter differs, the
variant belongs to a separate family or a report-only context.

Milestone 41 supports bottom-up, top-down, and hybrid handoff provenance. The
default structured v3 mode is a hybrid template-instantiated handoff: the family,
protocol, comparison method, margin policy, admissible targets, and constraint
schema are preregistered and frozen, while concrete constraints are instantiated
from a protocol-invariant response profile.

M41 executes the protocol-invariant patched v3 manifest-generation run. M41
does not evaluate carry-forward decisions. M41 does not run stress tests. M41
does not retune thresholds. M41 does not infer metric geometry. The resulting
diagnostic-complete v3 output bundle is not carry-forward evaluated; that
decision boundary is reserved for a later milestone.

## Milestone 42 Boundary

Milestone 42 applies fixed carry-forward criteria to the M41 v3 protocol bundle.
It does not regenerate manifests. It does not rerun fits. It does not retune
thresholds. It does not run stress tests. Carry-forward is stress-test
eligibility, not geometry. Top-down/hybrid provenance is allowed only as
preregistered manifest provenance, not as evidence.
