# V2 Stress-Test Handoff Decision

The v2 carry-forward registry determines whether future stress-test design is
allowed. The registry contains family-level decisions only. It does not
contain embeddings, fitted coordinates, or metric interpretations.

Milestone 38 evaluates the v2 diagnostic-complete bundle. It does not
regenerate manifests. It does not rerun fits. It does not retune thresholds.
It does not run stress tests. Carry-forward is stress-test eligibility, not
geometry.

If at least one v2 family is `carry_forward`, later milestones may design
future stress tests for those families only. If only `provisional` families
exist, future stress tests may be designed only with explicit caveats. If no
`carry_forward` or `provisional` family exists, future stress tests must stop
or remain report-only.

The handoff plan lists allowed future uses by family. Blocked, report-only,
and failed-control families are not future stress-test inputs.

Future stress tests still cannot claim physical distance, metric
reconstruction, a GR limit, or quantum structure. They would remain finite
diagnostics over frozen response-comparison artifacts.

Milestone 39 keeps this boundary closed. It audits blocked v2 families and
preregisters v3 design work, but it does not run stress tests.

Milestone 40 keeps this boundary closed as well. It is a pre-execution design
correction for response-profile protocol invariance. It does not generate v3
manifests, fit new representation models, retune thresholds, or run stress
tests.
