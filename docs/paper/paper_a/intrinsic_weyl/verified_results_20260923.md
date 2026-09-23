# Verified coverage result and incomplete high-count diagnostic

Checked on 2026-09-23 using CPU only. These are results for the frozen local
4D vacuum plane-wave populations and the dimensionless, internally framed
ten-component Weyl tensor. They do not establish arbitrary-spacetime accuracy.

## Completed 512-case coverage campaign

The fixed population completed across its documented interruption: **509/512**
cases satisfy both the existing continuum-bias condition and nominal full-tidal
containment. The exact 95% two-sided binomial interval is
**[0.9829726248831611, 0.9987900221742452]**. Its lower endpoint exceeds 0.95,
so the frozen coverage gate is **SUPPORTED** (required 497/512).

All three failures remain in the denominator: cases 275 and 303 have a zero
measurement clock; case 370 returns a fit whose ball misses the flat truth.
This is a coverage result, not a statement that every estimate is accurate.
Earlier failed pilots keep their original verdicts.

The final checker used the original MAIN inverse at source stamp
`f001de6af42ba9dbc411291178cbbd3a3e695294`, before any worktree merge. It checked
512 assigned outcomes, 1534 artifact hashes, 15 frozen sources, the 913 preserved
pre-interruption files, population/seeds, input schemas, and saved bias and
confidence calculations. A pure inverse replay of unfavorable case 370 matched
every fit field exactly; maximum numerical difference was zero. Truth arrays
were kept closed, and no generator, oracle evaluator or GPU library was loaded.
This does not regenerate all orders or independently rederive all oracle truth.

Evidence in the main repository's outputs directory:

- `intrinsic_weyl_nested_coverage_v1_continued/summary.json`, SHA-256
  `238ac2400deda3d6553e7da628bc0c27fdb277aa6bb221d924a2d21d6ed86d27`.
- `intrinsic_weyl_nested_coverage_v1_continued.final_audit.json`.

## High-count results available so far

Seven of eight fixed geometries have all three prefixes: 21/24 assigned
results. The last geometry (A=0.75, w=-0.5) was interrupted during generation;
it is retained as incomplete. No completed full-diagnostic verdict is available.

At N=4194304, all five completed nonzero-curvature geometries have smaller
point errors than their fresh N=262144 prefix. The signed full-tidal errors
were independently recomputed from E and B, including reflection parity.

| A | w | Relative point error, N=262144 | Relative point error, N=4194304 | Error plus nominal radius, divided by truth norm |
|---:|---:|---:|---:|---:|
| 0.25 | 0 | 41.85% | 16.23% | 46.61% |
| 0.75 | 0 | 9.12% | 1.64% | 10.94% |
| 0.25 | +0.5 | 56.20% | 14.76% | 45.34% |
| 0.75 | +0.5 | 23.63% | 5.92% | 16.69% |
| 0.25 | -0.5 | 62.49% | 8.42% | 30.70% |

The final column is a descriptive triangle-bound calculation, not a newly
introduced pass/fail gate. It is below 25% in two completed geometries; the
three weak-curvature geometries do not meet that sufficient bound. The original
diagnostic has no powered accuracy gate, and these nested streams are not
independent replicates. A small point error alone does not certify precision.

All 21 available pure fits match their saved evaluation wrappers. Their input
hashes, declared import boundary, errors, covariances and nominal radii check.
All 23 available input prefixes preserve the original selection and exact
nesting, with valid landmark/event causal transitivity. All 132 checkpoint
files, 12 bound sources and 25 old-pilot files were verified unchanged.

Evidence: `outputs/intrinsic_weyl_high_count_v1.partial_audit_21.json`, SHA-256
`cfb82aced2cd2dc38826c8216e81af1071c84121db819e32af31707e2c1004ee`.
Its sibling `.py` is the runnable post-result audit. This verifies saved fits,
not 21 new inverse runs. The established truth normal-coordinate map is reused.

The remaining case follows [the fixed continuation](high_count_continuation.md).
After it completes, audit all assignments and measured costs before deciding
whether the sample-price extrapolation supports a separately preregistered,
independent powered accuracy campaign. The full research goal is not complete.
