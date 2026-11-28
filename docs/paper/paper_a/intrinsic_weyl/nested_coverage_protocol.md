# Frozen independent coverage check for nested Weyl reconstruction

2026-09-22. Commit this protocol and implementation before drawing any of
its new random orders. The [eight-pair pilot](nested_pilot_results.md)
reduced actual tensor error in every pair, but one of its narrowed nominal
balls missed truth. No population coverage claim is made from that pilot.
This campaign tests the unchanged candidate with all failures included.

## Frozen estimator and population

Use the nested union of the 64- and 128-maximum connected selectors, preserving
all baseline landmarks and the same frame vectors. Keep the 192-interior cap,
selection count threshold two, full shared-membership covariance, eigenvalue
cutoff, four initialization starts, tolerances and optimizer budget unchanged.
Fit the fifth-degree vacuum jet plus the fixed curvature-square world-function
and volume-clock terms. All ten Weyl and 24 gradient components remain free.
Do not add profile/type priors, supply coordinates, or use a flat comparator.

The estimator receives only its six-field order archive. The generator and
truth evaluator remain separate; truth is read after the pure fit returns.
The nominal 95% ten-component confidence ball, residual-reflection quotient
and **2% relative bias allowance** are unchanged. No radius inflation, changed
clock mask or omitted failed case is permitted during the run.

Draw **512 iid cases** from the product mixture

```text
A uniform on {0,.25,.75}; w independent uniform on {-.5,0,.5}.
```

Use `default_rng(2026100950)`, draw all 512 amplitudes first, then all 512
rotation rates. Preserve the ordered population in the manifest before cases
run. Case i uses seed **2026101000+i**, with independent streams `[seed,0]`
for selection and `[seed,1]` for measurement. None overlaps earlier pilots
or the failed fifth-degree coverage campaign.

Each case receives **4096 selection and 262144 measurement events**, generated
in chunks of 4096. The fixed maximum is **136,314,880 sampled events plus 1024
supplied tips**. Generate only the nested arm; there is no paired baseline fit
in this campaign. All 512 assigned cases remain in the denominator, including
selection, empty-clock, rank, boundary and optimization failures. Record the
events actually generated if a selection failure prevents measurement.

## Primary gate and prospective power

For each returned fit, refit exact continuum volume clocks on the SAME
selected landmarks, clock mask and empirical covariance. Require 32/64
quadrature refinement below 1e-10. Define success as the conjunction of

1. successful generation and full-rank, converged, interior fitting;
2. relative full-tidal continuum bias <= .02, with numerical flat tolerance
   1e-10; and
3. truth contained in the nominal full-tidal ball.

Compute an exact **95% two-sided Clopper-Pearson interval** for the 512 iid
success indicators. Support the reliability gate only if its lower endpoint
is at least .95. This requires **at least 497 successes**, allowing at most
15 failures; 496 successes must fail the mechanical gate check.

The exact power at true end-to-end success probability **.98 is
.9442856165**, exceeding the repository's 90% planning target. At .975 the
power is only .7834228654, and at .99 it is .9999191612. The .98 alternative
is a planning assumption, not an estimate established by eight cases. This
uses more cases than the previous 256-case design, whose power at .98 was
only .5947069136. A correct 95% nominal construction alone does not establish
that the stronger .98 planning assumption holds.

No early stop, replacement, event extension, altered threshold or outcome-
dependent exclusion. Complete all cases even after the support gate becomes
unreachable. A process/storage interruption that prevents completion leaves
the run incomplete, not passed. Report stratum outcomes descriptively; the
primary claim concerns this mixture, not each stratum or all vacuum metrics.

## Execution, checks and limits

Use four fresh Windows worker processes with the pilot's unchanged inherited
BLAS/OpenMP environment. Record `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS`
(null means unset), without changing them. The observed host has 16 logical
processors. Computation remains substantial and its throughput is not assumed
from a single replay.

A pre-data engineering trial with four BLAS threads selected equivalent raw
axis signs and matched physical Weyl/covariance quantities to roundoff, but
its continuum nuisance gradient differed by 2.44e-8, exceeding the full-jet
1e-8 replay check. That resource change is **not adopted**; the guard is not
relaxed. Its one 74.14-second replay also did not improve elapsed time. These
old-data engineering trials do not re-score any scientific outcome.

Before new data, run the original confidence/generator/jet checks, the nested
exact-flat rank and inclusion check, and explicit 497/496 gate controls. Replay
old nested-pilot case 5 through one spawned worker in the unchanged environment:
require exact equality
of all six archive fields, and agreement of canonical physical coefficients,
covariance, coordinates, full covariant gradient and evaluation to **1e-8
absolute and relative numerical tolerance**, with identical rank/success
flags. This comparison is not an accuracy tolerance or a changed verdict rule.

The required spawned-worker replay passed before any new data: all six input
fields were exact, and canonical finite/continuum coefficients, covariance,
radius, coordinates and covariant gradients had zero numerical difference.
Evaluation and rank/success flags also matched. Both recorded thread variables
were unset. The audit is retained at
`outputs/intrinsic_weyl_nested_worker_default_replay_v1/audit.json`; its
coverage-runner source SHA256 (LF-normalized) is
`de6d2e83d79064d4c23040e5fd58ef03e78f8788636addc6dfdd77af35fc1146`.

Record clean reachable generating SHA, source hashes, population, counts,
runtime, all raw inputs, returned fits and failures. After completion verify
hashes, iid assignments, the independent gate calculation and a pure estimator
replay. Preserve noncoverage and bias failures even when the other check passes.

This validates reliability at N=262144 if supported. It is **not** a 25%
accuracy campaign, a calibrated high-count price, a uniform remainder bound,
or a claim for arbitrary-frequency, singular or globally complicated metrics.
The later accuracy gate still requires confidence containment and 90% power
sizing before its data. The regular-convex-vacuum and dimensionless-scale
conditions remain unchanged.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_jet_coverage.py --nested
```
