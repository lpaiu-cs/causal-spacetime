# Cross-Family Robustness Criteria

Milestone 34 defines carry-forward eligibility for future stress tests. It
follows the Milestone 33 manifest-family comparison and asks which families
are worth carrying into later representability stress tests under fixed
diagnostic thresholds.

The decisions are:

- `carry_forward`: the family passes the preregistered hard criteria and the
  available non-hard robustness checks.
- `provisional`: the family passes the hard criteria but has limited non-hard
  warnings or failures.
- `blocked`: the family fails a hard criterion or too many non-hard criteria.
- `report_only`: the family is retained for accounting but is not eligible for
  future stress tests.
- `failed_control`: the family is a failed-control family and remains visible
  as a control.

Carry-forward is stress-test eligibility, not representation success.
Carry-forward eligibility is not evidence of metric geometry. A carried family
only becomes an allowed input to later perturbation, transfer, null, sparsity,
or optimizer-stability diagnostics.

Target-label permutation remains a symmetry control. It should not be grouped
with destructive nulls because consistent relabeling preserves the
response-comparison structure. Destructive nulls include shuffled sides and
random constraints. Marginal-preserving nulls preserve selected profile
marginals and are reported separately.

Stricter-threshold failures must be reported. Failed, ineligible, blocked, and
provisional families must remain visible in reports. No thresholds are retuned
after seeing fit outcomes, and threshold sensitivity is not threshold retuning.

Milestone 34 performs no new fitting, embedding, or stress testing. It reads
Milestone 33 family-level outputs, aggregates fixed diagnostics, writes family
decisions, and exports a registry for future work.

Milestone 35 decomposes failure modes after the carry-forward stop condition.
It does not retune thresholds. It does not run stress tests. Missing metrics
are not success and not the same as measured failure. Future remediation must
be preregistered.

Milestone 36 preregisters a remediation plan for diagnostic-complete future
manifest generation. It does not execute remediation, fit new representation
models, or change the fixed carry-forward thresholds.
