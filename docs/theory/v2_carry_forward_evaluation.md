# V2 Carry-Forward Evaluation

Milestone 38 evaluates the v2 diagnostic-complete bundle. It applies the
fixed cross-family robustness criteria from Milestone 34 to the v2 output
bundle produced by Milestone 37.

Diagnostic completeness is necessary for this evaluation, but it is not
sufficient for carry-forward. A family can have all required metrics present
and still fail fixed robustness criteria.

V2 family decisions use the same decision vocabulary as Milestone 34:

- `carry_forward`: eligible for future stress-test design.
- `provisional`: eligible only for limited future stress-test design with
  explicit caveats.
- `blocked`: not eligible for future stress tests and retained in reports.
- `report_only`: retained for accounting only.
- `failed_control`: retained as a failed-control family.

Carry-forward is stress-test eligibility, not geometry. It does not imply
physical coordinates, calibrated scale, or metric reconstruction.
No metric geometry is inferred.

Milestone 38 is a decision layer. It does not regenerate manifests. It does
not rerun fits. It does not retune thresholds. It does not run stress tests.
It reads the M37 v2 metrics, writes v2 decisions, exports a v2 registry, and
builds a future stress-test handoff plan.

Blocked, provisional, report-only, and failed-control families remain visible
in accounting. Missing or incomplete v2 bundle inputs are reported explicitly
rather than silently ignored.
