# V3 Execution Boundary

Milestone 39 defines a future execution boundary for v3. It may specify what a
future Milestone 40 is allowed to generate, but it must not generate those
artifacts itself.

Milestone 40 does not execute the preregistered v3 manifest-generation run.
It audits response-profile protocol invariance and patches the planned v3
family semantics before any production v3 manifests exist.

A later milestone may execute the protocol-invariant patched v3
manifest-generation run if the patched preregistration remains fixed. Such a
run would write production manifests under `outputs/manifests_v3/` and
produce a diagnostic-complete v3 output bundle.

Milestone 39 must not execute v3. It must not generate v3 production
manifests, run representation fits, apply carry-forward decisions, or run
stress tests. It must not alter the Milestone 34 fixed criteria.

Carry-forward evaluation for v3 should be separated into a later milestone.
No stress test is allowed before a later fixed-criteria evaluation produces a
carry-forward or provisional family.

Milestone 40 is a pre-execution design correction, not threshold retuning. It
does not change M34, M38, or M39 decisions. It does not generate v3 manifests.
It does not run stress tests. It does not fit new representation models. A
single response profile must not mix measurement protocols.
## Milestone 41 Boundary

M41 executes the protocol-invariant patched v3 manifest-generation run and writes
production v3 handoff manifests under `outputs/manifests_v3/`. M41 does not
evaluate carry-forward decisions. M41 does not run stress tests. M41 does not
retune thresholds. M41 does not infer metric geometry.
