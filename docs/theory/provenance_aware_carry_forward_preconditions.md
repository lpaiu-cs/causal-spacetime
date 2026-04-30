# Provenance-Aware Carry-Forward Preconditions

Milestone 42 applies fixed carry-forward criteria to the M41 v3 protocol bundle,
but first checks provenance-aware admissibility preconditions. Preconditions are
not performance thresholds and do not replace the fixed M34 criteria.

Production v3 manifests must include `measurement_protocol` metadata,
`profile_metadata`, and `handoff_provenance`. Structured families must be
`protocol_invariant`, parameter-complete, admissible for production pairwise
response-profile dissimilarity, and backed by valid bottom-up, top-down, or
hybrid provenance. Report-only controls must remain ineligible.

Invalid provenance blocks production carry-forward eligibility because the
manifest cannot be treated as an admissible frozen handoff. This is an
admissibility decision, not a metric result.

Diagnostic completeness does not imply carry-forward. Protocol invariance does
not imply carry-forward. Provenance integrity does not imply carry-forward.
Top-down/hybrid provenance is allowed only as preregistered manifest provenance,
not as evidence.

M42 does not regenerate manifests. M42 does not rerun fits. M42 does not retune
thresholds. M42 does not run stress tests. Carry-forward is stress-test
eligibility, not geometry.

## Milestone 43 Follow-Up

Milestone 43 audits why v3 protocol families remain blocked after those
preconditions have been checked. It does not change M42 decisions. It does not
retune thresholds. It does not run stress tests. It does not execute v4. V4
design is planned-only and requires a later execution milestone.
