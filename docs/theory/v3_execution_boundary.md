# V3 Execution Boundary

Milestone 39 defines a future execution boundary for v3. It may specify what a
future Milestone 40 is allowed to generate, but it must not generate those
artifacts itself.

Milestone 40 may execute the preregistered v3 manifest-generation run if the
v3 preregistration remains fixed. Such a run would write production manifests
under `outputs/manifests_v3/` and produce a diagnostic-complete v3 output
bundle.

Milestone 39 must not execute v3. It must not generate v3 production
manifests, run representation fits, apply carry-forward decisions, or run
stress tests. It must not alter the Milestone 34 fixed criteria.

Carry-forward evaluation for v3 should be separated into a later milestone.
No stress test is allowed before a later fixed-criteria evaluation produces a
carry-forward or provisional family.

