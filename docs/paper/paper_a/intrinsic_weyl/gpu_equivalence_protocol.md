# GPU equivalence before independent weak-curvature accuracy data

2026-09-23. The user released the GPU restriction. The completed high-count
diagnostic is unchanged. This is a resource engineering check, not new
evidence of scientific accuracy. Commit the check and optional GPU modules
before reconstructing any old measurement stream.

Keep NumPy 2.5.3, SciPy 1.18.1, CPU physical-volume sampling and its 4096-event
chunks unchanged. The optional dependency is CuPy 14.2.0 with the installed
CUDA 13 toolkit. All device arithmetic is float64. Copy the four guarded
relation-query method bodies verbatim to the GPU array backend and compare
their ASTs. Combine eight sampling chunks per 32768-event GPU query batch.
Preserve event order, mathematical expressions and ambiguity guards.

Reconstruct all eight saved high-count 262144 prefixes using the original
seeds 2026103000+i and old selected coordinates. Require exact packed-bit
equality to every saved after/before prefix; no tolerated bit mismatch.
At lengths 1e-9, 1e-5, .1, .75, 1, require CPU/GPU action differences below
both conservative error guards. Construct null-margin queries and require
both to reject; check both resolved signs at 100 times the larger guard.

The estimator-side accelerator imports no generator or truth evaluator.
Retain the CPU grouping of causal-membership signatures. Transfer only the
weighted matrix products to float64 GPU GEMM. Inputs, products and all
nonnegative partial sums are integers below 2**53, hence exactly representable.
Test against direct ungrouped membership products and the old grouped helper
for repeated, distinct and single-pattern cases.

In a separate process, replay saved weak rotating case 4 at N=262144.
Allow np.load to open only that order archive. Require exact equality of
counts, masks, covariance, IDs and N between CPU and accelerated readers.
Then require every field of the full pure fit to equal the stored pure fit,
with no generator/truth module imported before return. Keep CPU BLAS and
OpenMP threads at one and the other three old overrides unset.

Time three warmed CPU/GPU packed-query repetitions on the first case-2 batch,
including GPU transfers, and record both complete clock-reader times.
Timing is descriptive under current load; it changes no scientific gate.
Bind the clean reachable commit, GPU/source hashes, existing frozen CPU source
and pilot hashes, runtime/device and all outcomes. A mismatch stops adoption
before independent accuracy data. It is not repaired by loosening a gate.

The later accuracy study will separately fix the population (including new
weak selection geometries and flat controls), confidence/error allocations,
variance caps, bias assumptions, power-based event counts and independent
seeds before any of its data. This engineering check cannot certify those
statistical assumptions or declare the weak-curvature goal complete.
