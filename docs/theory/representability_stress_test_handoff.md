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
