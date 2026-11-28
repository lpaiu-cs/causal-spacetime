# Preserving the fixed 512-case study across an interruption

2026-09-23. This records an execution recovery before evaluating any unfinished
case or generating any remaining data. It changes no estimator, population,
seed, sample count, failure rule, confidence construction or statistical gate
in [the frozen protocol](nested_coverage_protocol.md).

## Observed interruption and immutable checkpoint

The former execution handle is missing. A process inspection found no coverage
runner or Python multiprocessing worker. The source directory has no final
summary, 226 complete result records (IDs 0--225), and four complete input/truth
pairs without fit results (226--229). Its last observed artifact modification
was 2026-09-22 14:31:46 UTC; that is not a verified interruption time or cause.
The original run remains **INCOMPLETE**, with no coverage verdict.

All 15 frozen source hashes and all 678 artifacts referenced by completed rows
match. All eight unfinished-input files are readable, including their compressed
members. The frozen Python 3.14.0, NumPy 2.5.3 and SciPy 1.18.1 versions match;
both recorded BLAS/OpenMP environment variables remain unset. The manifest hash
is `df14768daefa8f92a9ef8b9fae070cccb996d565d6228679cf30c5e54fc28480`.

## Exact continuation rule

Preserve the original directory. Copy its 913 files byte-for-byte into a new
continuation directory and record every source-file hash before dispatch.
Keep all 226 result records without re-evaluation or exclusion. Reuse the four
existing order archives without regeneration, and evaluate them with the frozen
evaluator. Generate only IDs 230--511, using their original assigned seeds and
population entries. Use four fresh Windows workers and the unchanged thread
environment. The isolated coordinator imports scientific modules from the
original, clean frozen-source checkout; the newer covariance optimization is
not used in this dataset.

Before dispatch, one fresh worker must reproduce every saved pure-fit and
confidence field of completed case 225 exactly, without importing the rotating
generator or oracle evaluator. Failure stops the continuation before new data.
The coordinator refuses changed checkpoint contents or existing result files.
Its committed source version and hash are recorded separately from the original
scientific source stamp. All original and copied checkpoint hashes are checked
again before a final combined summary can be written.

The full denominator remains 512, including every generated, selection, fit,
bias or noncoverage failure. At least 497 valid-and-covered cases are still
required. There is no early stop, replacement seed, enlarged sample, changed
radius, or new calibration population. The original interrupted directory never
receives a passing summary. A completed continuation must be described as the
same fixed population completed across an interruption, with the interruption
and reuse disclosed, and independently audited before any claim.

For the four reused inputs, reported case seconds cover continuation evaluation
only; their earlier generation time was not saved. Calendar time across the
interruption must not be reported as computation time. The dataset still has
the original maximum of 136,314,880 distinct sampled events plus 1024 tips.

The coordinator's `--check` command performs a read-only checkpoint, source,
runtime and population audit. It is not a scientific coverage result.

## Terminal evidence located after the continuation started

At 2026-09-22 23:27:57 UTC, a bounded inspection of the local execution record
found the original command's completion event at **14:31:54.833 UTC** that day.
The execution item is `exec-b150d07c-6c93-4465-9a5c-16f177c6d4e9`, process handle
`55336`. It records status `failed`, exit code **-1**, duration **18330.908958
seconds**, and empty stderr. All 226 stdout JSON records match the preserved
case-result files field-for-field; no traceback or case-failure record appears
in that command output.

This confirms terminal failure, beyond the later missing-handle observation.
The generic exit code and empty stderr do **not** establish its cause. The
original run remains incomplete; no scientific verdict or continuation rule
changes on this evidence.

The focused receipt is
`outputs/intrinsic_weyl_nested_coverage_v1_interruption_evidence.json`, SHA-256
`220a1ab6c904b172d2784c3eb9e10f4f990dd24a420a8186e7cf1b513a51d741`.
It includes the command, event identity and stdout hash. The original UTF-8
event record, including its newline, has SHA-256
`d03ea1a7b8d5944b1a0bce27bd623586d85e0d5ca9dbf2f74ad11f48c1402aa2`.
