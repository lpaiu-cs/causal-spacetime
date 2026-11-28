# Complete the interrupted fixed high-count diagnostic

2026-09-23. This continuation preserves the original diagnostic, population,
three measurement levels, estimator, masks, seeds, tolerances and failure rules.
It adds no accuracy gate. The user has reserved the GPU; use CPU only.

The recorded coordinator PID 620 is absent (Windows OpenProcess error 87).
Its last status says generation of case 7 at 3670016 events. Cases 0--6 have
all 21 results; case 7 has only its 262144 and 1048576 prefix archives. There
is no summary or case-7 generation result. The cause is unknown. Observation:
`outputs/intrinsic_weyl_high_count_v1.interruption_evidence.json`, SHA-256
`5c68aef9d1f47c0aeb8e1fee88991932e768cf4047af333921e29a3b98dfee91`.

Before running, commit the coordinator and this protocol. Verify the recorded
132 original files, all 12 bound sources, 25 old-pilot inputs, runtime and
original reachable source commit. Copy the 132 files into a new destination;
do not change the original or copied checkpoint files, including its stale
status. Write a separate coordinator status and continuation receipt.

Before any regenerated measurement, replay saved high-count case 6 at 262144
events in the original single-thread worker. Compare all finite-fit,
inference, continuum-fit and evaluation fields with the saved full wrapper
using the existing 1e-8 absolute/relative comparison and exact discrete flags.
Preserve its pure fit and comparison. A mismatch stops before regeneration.

Regenerate only case 7 using the original `default_rng([2026103007,1])`,
4096-event chunks, old selection and generator. Reconstruct the two saved
prefixes bit for bit, and stop on any mismatch. Never replace those archives.
Only the missing 4194304-event archive is written. This re-executes a fixed
stream; it is not a new seed, replacement case or independent replicate.
Retain ordinary guarded generation/fit failures and all three assignments.

Use the unchanged frozen worker, which writes the pure order-only fit before
importing truth evaluators. Keep OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1,
the other three BLAS overrides unset, and CUDA_VISIBLE_DEVICES=-1. No GPU
library or process is used. Run one generator and one sequential worker.

Recheck all preserved/source/input hashes at completion and combine the 21
unchanged results with the three assigned case-7 outcomes. Keep unique
archived events separate from computational sample operations: at least
3670016 previous case-7 draws/classifications must be charged in addition
to this replay. Its exact interrupted wall/CPU cost is unavailable, not zero.
Independently audit the completed continuation before calling the diagnostic
complete. A further interruption remains incomplete.
