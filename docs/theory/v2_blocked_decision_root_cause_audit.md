# V2 Blocked-Decision Root-Cause Audit

Milestone 39 audits why v2 families were blocked. It follows Milestone 38,
where the fixed Milestone 34 cross-family robustness criteria were applied to
the diagnostic-complete v2 bundle. No v2 family became carry-forward or
provisional, and `failed_controls_v2` remained a failed-control family.

Milestone 39 does not change those decisions. Blocked v2 families remain
reported. The purpose is to separate root-cause categories before any future
design work is executed.

Structural blocking and measured blocking are distinct. Structural blocking
means that a fixed criterion cannot be satisfied by the generated output
structure, such as a family with `manifest_count` below the fixed minimum.
Measured blocking means that a finite diagnostic fails a fixed threshold, such
as held-out violation or latent-order disagreement. Diagnostic blocking means a
required metric is missing or incomplete. Control-family blocking means a
failed or report-only control is intentionally not eligible.

Manifest-count failure is therefore different from heldout-violation failure.
The v2 bundle achieved diagnostic completeness, so missing metrics are not the
main explanation for the structured-family blocks. The audit asks which blocks
come from one-manifest family design and which remain after that structural
count issue is considered.

Structural counterfactuals are report-only. They do not alter Milestone 38
decisions and do not justify threshold changes. Milestone 39 does not retune
thresholds. Milestone 39 does not run stress tests.

Milestone 40 audits response-profile protocol invariance before v3 execution.
This is a pre-execution design correction, not threshold retuning. It does not
change M34, M38, or M39 decisions. It does not generate v3 manifests. It does
not run stress tests. It does not fit new representation models. A single
response profile must not mix measurement protocols. This patch preserves the
M39 structural-versus-measured blocking report.
