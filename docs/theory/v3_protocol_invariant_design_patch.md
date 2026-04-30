# V3 Protocol-Invariant Design Patch

Milestone 40 patches v3 design semantics before execution. The v3 design was
planned but not executed in Milestone 39, so Milestone 40 can correct the
planned family definitions without generating v3 manifests.

This is a pre-execution design correction, not threshold retuning. It does not
change M34, M38, or M39 decisions. It does not run stress tests. It does not fit
new representation models.

## Patch Rule

Each structured v3 family must be protocol-invariant. Different measurement
protocols must become separate v3 families. A single response profile must not
mix measurement protocols.

The patched v3 families remain planned-only. V3 production manifests will be
generated only in a later milestone.

## Planned Families

The patch separates earliest full, earliest retained-reference, gated full,
coverage-enriched full, combined earliest full, and immediate-edge report-only
contexts. Failed controls remain report-only controls with explicit metadata.

These planned v3 families are not current results. They do not imply future
carry-forward status and do not infer metric geometry.
