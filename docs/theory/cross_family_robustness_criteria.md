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

Milestone 37 does not apply these criteria to v2 outputs. The v2 bundle is
prepared for a future Milestone 38 evaluation, where the fixed criteria can be
read without changing the M37 generation results.

Milestone 38 applies fixed criteria to v2 outputs. It does not regenerate
manifests, rerun fits, retune thresholds, or run stress tests. Carry-forward
is stress-test eligibility, not geometry.

Milestone 39 separates structural blocking from measured blocking in the v2
results. Structural counterfactuals are report-only, and v3 design is
preregistered but not executed.

Milestone 40 does not apply robustness criteria. It audits response-profile
protocol invariance before v3 execution and patches planned v3 family
semantics. The patch does not change M34, M38, or M39 decisions, does not
generate v3 manifests, does not run stress tests, and does not fit new
representation models.
## Milestone 41 Non-Decision Boundary

M41 produces a diagnostic-complete v3 output bundle but does not apply the fixed
cross-family robustness criteria. M41 does not evaluate carry-forward decisions.
M41 does not run stress tests. M41 does not retune thresholds. M41 does not infer
metric geometry. A future milestone must apply the unchanged fixed criteria to
the M41 outputs.

## Milestone 42 Boundary

Milestone 42 applies fixed carry-forward criteria to the M41 v3 protocol bundle.
The M34 thresholds remain fixed. M42 does not regenerate manifests, rerun fits,
retune thresholds, or run stress tests. Carry-forward is stress-test
eligibility, not geometry. Top-down/hybrid provenance is allowed only as
preregistered manifest provenance, not as evidence.

## Milestone 43 Use

Milestone 43 audits why v3 protocol families remain blocked under the same
fixed M34 criteria. It does not change M42 decisions. It does not retune
thresholds. It does not run stress tests. It does not execute v4. V4 design is
planned-only and requires a later execution milestone.
