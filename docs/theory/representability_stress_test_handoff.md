# Representability Stress-Test Handoff

The future stress-test handoff boundary separates Milestone 34 carry-forward
decisions from later experiments. A family may be carried forward only after
the fixed cross-family robustness criteria are applied to the frozen
family-level diagnostics.

Future stress tests may examine:

- perturbation robustness,
- held-out manifest transfer,
- harder null baselines,
- constraint sparsity,
- optimizer stability.

These stress tests still cannot claim spatial geometry, physical metric
distance, calibrated scale, a GR limit, or quantum structure. They would remain
representability diagnostics over response-comparison constraints unless
additional calibration and interpretation assumptions are added.

The handoff stop rules are:

- if no family meets carry-forward criteria, future stress tests must report
  that and stop or run only report-only controls,
- provisional families may be stress-tested only with explicit caveats,
- blocked families are report-only,
- failed-control families stay in accounting and are not stress-test inputs.

Milestone 34 defines cross-family carry-forward decisions. Carry-forward is
stress-test eligibility, not geometry. Provisional and blocked families remain
reported. Threshold sensitivity is not threshold retuning. If no family passes,
later stress tests should stop or run only report-only controls.

Milestone 35 makes the stop condition explicit. It decomposes measured
failures and missing diagnostics, audits completeness, and records upstream
remediation designs without changing thresholds or running stress tests.
Missing metrics are not success and not the same as measured failure. Future
remediation must be preregistered.

Milestone 36 produces a preregistered remediation plan but does not execute it.
No stress tests are allowed until carry-forward or provisional families exist.
Planned v2 families are future specifications, not current stress-test inputs.

Milestone 37 executes the v2 manifest-generation run but does not cross the
stress-test boundary. The v2 outputs are prepared for a future carry-forward
evaluation, and stress tests remain unavailable until that later evaluation
allows them.

Milestone 38 performs that carry-forward evaluation for v2 families. If no v2
family is carry-forward or provisional, the future stress-test handoff remains
stopped or report-only.

Milestone 39 does not reopen the handoff. It audits why the v2 handoff stopped
and preregisters a v3 design, while leaving stress tests unavailable.

Milestone 40 also does not reopen the handoff. It audits response-profile
protocol invariance before v3 execution and patches v3 design semantics only.
No stress tests are allowed, no v3 manifests are generated, no new
representation models are fit, and no carry-forward decisions are changed.
## Milestone 41 Stress-Test Boundary

M41 executes the protocol-invariant patched v3 manifest-generation run and
produces a diagnostic-complete v3 output bundle. It does not authorize stress
tests. Stress tests remain blocked unless a later carry-forward evaluation finds
at least one carry_forward or explicitly provisional v3 family.

## Milestone 42 Boundary

Milestone 42 applies fixed carry-forward criteria to the M41 v3 protocol bundle
and writes a future stress-test handoff plan. It does not run stress tests. It
does not regenerate manifests, rerun fits, or retune thresholds. Carry-forward
is stress-test eligibility, not geometry.

## Milestone 43 Stop Rule

Milestone 43 audits why v3 protocol families remain blocked after the M42
decision layer. It does not change M42 decisions. It does not retune thresholds.
It does not run stress tests. It does not execute v4. V4 design is planned-only
and requires a later execution milestone.
