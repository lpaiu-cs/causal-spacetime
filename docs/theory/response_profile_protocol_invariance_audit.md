# Response-Profile Protocol Invariance Audit

Milestone 40 audits response-profile protocol invariance before v3 execution.
It is inserted before v3 execution because the M39 design was planned-only and
there is still time to correct profile semantics before production v3 manifests
exist.

This is a pre-execution design correction, not threshold retuning. It does not
change M34, M38, or M39 decisions. It does not generate v3 manifests. It does
not run stress tests. It does not fit new representation models.

## Why The Audit Is Needed

Pairwise response-profile dissimilarity assumes comparable columns. A response
profile used for pairwise response-profile dissimilarity may vary reference
chains inside one fixed measurement protocol. If emission, gate, echo rule,
spectrum type, subsampling, normalization, missing policy, tie policy, or
margin policy varies, those variants must form separate profile families or be
explicitly marked exploratory/report-only.

Measurement-protocol variation and reference-chain variation therefore have
different roles. Reference-chain variation can populate one protocol-invariant
multi-reference response profile. Measurement-protocol variation belongs in
cross-protocol robustness diagnostics.

## Decision Preservation

The audit does not change v2 carry-forward decisions. It also does not change
the M39 blocked-decision analysis. Existing blocked v2 families remain
reported, and no blocked family is reinterpreted as eligible.

## Mixed Contexts

Mixed contexts may be useful for exploratory diagnosis, but they must be
explicitly marked report-only. They are not admissible production inputs for
pairwise response-profile dissimilarity.
## Milestone 41 Extension

M41 executes the protocol-invariant patched v3 manifest-generation run after the
Milestone 40 protocol-invariance audit. It does not change M34, M38, M39, or
M40 decisions. It does not evaluate carry-forward decisions, does not run stress
tests, does not retune thresholds, and does not infer metric geometry.
