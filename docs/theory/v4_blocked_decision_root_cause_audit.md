# V4 Blocked-Decision Root-Cause Audit

Milestone 46 audits why v4 protocol families remain blocked. It follows from
Milestone 45 because no v4 protocol family received a carry_forward or
provisional decision under the fixed M34 criteria.

Milestone 45 was a fixed-criteria decision layer over the diagnostic-complete
v4 protocol bundle. Milestone 46 does not change M45 decisions. It decomposes
measured v4 failures into criterion-level root causes so that any future
remediation remains tied to upstream design issues rather than threshold
changes.

The root-cause categories include held-out failure, stricter-pass failure,
destructive-null gap failure, latent-order instability, coverage failure,
generalization-gap failure, restart instability, symmetry-control failure,
precondition failure, and control-family blocking.

The M45 result keeps stress-test design stopped. No stress-test design is
allowed unless a later fixed-criteria evaluation produces a carry_forward or
provisional family.

M46 does not retune thresholds. It does not run stress tests. It does not
execute v5. It does not infer metric geometry.
