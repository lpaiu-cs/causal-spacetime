# V4 To V5 Remediation Iteration Audit

Milestone 46 audits remediation-iteration risk because v2, v3, and v4 have all
failed fixed carry-forward criteria. Another remediation cycle is not
automatically justified.

Each proposed v5 change is classified as one of the following:

- semantic correction,
- diagnostic-completeness correction,
- protocol-design correction,
- profile-resolution correction,
- stability-design correction,
- coverage-design correction,
- overfitting-risk change.

Designs that only chase failed thresholds without a protocol-level rationale
are flagged. A proposal that changes thresholds, bypasses the stopped
stress-test boundary, removes hard cases after seeing outcomes, or searches
open-endedly for a passing family is treated as a remediation-iteration risk.

Report-only counterfactuals are allowed only as diagnostic summaries. They do
not change M45 decisions and cannot justify threshold retuning.

Repeated remediation is an overfitting risk unless justified by root-cause
analysis. M46 therefore separates justified upstream remediation from changes
that would amount to de facto tuning.
