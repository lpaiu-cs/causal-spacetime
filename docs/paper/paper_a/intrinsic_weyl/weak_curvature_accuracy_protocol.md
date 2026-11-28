# Independent weak-curvature precision test

2026-09-23. Freeze this protocol and its runner before any new selection or
measurement event. The objective is to resolve the weak-wave precision failure
identified in the completed diagnostic: its three A=.25 cases had nominal
error-plus-radius ratios .307--.466 at N=4194304. Keep the full order-only,
internally framed ten-component target and the fifth-degree vacuum jet with
24 free gradient components, fixed curvature-square terms, nested selector,
clock mask, full covariance, four starts, cutoffs and optimizer unchanged.
No supplied coordinates, wave-profile prior or paired flat input is allowed
in the inverse. This does not promise uniform recovery as amplitude tends
to zero or over arbitrary frequencies, singularities or global geometries.

## Assigned cases and independent data

All eight cases are co-required. Retain all three difficult old weak geometries,
add three independently selected weak geometries, and retain both flat controls.
The existing pilot contributes only selection order/IDs and withheld geometry;
none of its measurement events, or the high-count measurement events, is reused.

| Case | A | w | Selection | Variance cap V | Final measurement N |
|---:|---:|---:|---|---:|---:|
| 0 | 0 | 0 | old pilot case 0 | 4000 | 25,108,480 |
| 1 | 0 | +.5 | old pilot case 1 | 4000 | 25,108,480 |
| 2 | .25 | 0 | old pilot case 2 | 5500 | 49,311,744 |
| 3 | .25 | +.5 | old pilot case 4 | 5700 | 51,101,696 |
| 4 | .25 | -.5 | old pilot case 6 | 2900 | 26,001,408 |
| 5 | .25 | -.5 | fresh selection | 8000 | 71,720,960 |
| 6 | .25 | 0 | fresh selection | 8000 | 71,720,960 |
| 7 | .25 | +.5 | fresh selection | 8000 | 71,720,960 |

Case i uses seed 2026112000+i. For cases 5--7 use an independent stream
default_rng([seed,0]) for 4096 physical-volume selection events, append the
same supplied outer tips, permute interior labels, form the original guarded
CPU selection order and apply selection(order, design="nested") unchanged.
Every case uses default_rng([seed,1]) for its entirely new measurements.
Keep CPU sampling chunks at 4096; group up to eight chunks for GPU relations.
The total is 391,794,688 measurement events plus 12,288 new selection events
and six supplied tips. Retained selections incur no new selection draws.

## Simultaneous confidence and interval accuracy gate

Set family confidence to **95%, central/two-sided**, with Bonferroni error
.05/8 per case. Each ten-component ball has nominal coverage .99375.
Use the existing covariance and reflection quotient. Let delta=.02 and
eta=1e-10. The radius is

    r = (sqrt(chi2_10(.99375)*lambda_max) + delta*||estimate|| + eta)/(1-delta).

The eta term makes the existing numerical bias slack explicit in the ball.
It is fixed before data. No confidence level, bias allowance or radius changes
are permitted after outcomes. The family statement uses all eight balls.

For each case require all of:

1. Completed assigned generation and the full-rank, interior, converged pure
   estimate, with a valid continuum refit using the SAME clock mask and
   measured covariance and a 32/64 quadrature relative change below 1e-10.
2. Oracle continuum bias <= delta*||truth||+eta. This is a model-scope check;
   it is not a data-derived bound for an unknown arbitrary profile.
3. The adjusted confidence ball contains truth.
4. **The whole ball is within the accuracy target:** error+r <= .25*||truth||
   for the six weak cases; error+r <= .125 for each flat control. Relative
   error at zero remains undefined. Errors use the original full-tidal norm
   and residual-reflection quotient.
5. The observed V=N*lambda_max is at most the case's fixed planning cap.

Support the precision gate only if every requirement holds in every assigned
case. Missing, selection, empty-clock, guard, fit, rank, bias, variance-cap and
accuracy failures remain in the result. Complete later cases after ordinary
failures; no replacement, redraw, reseeding, sample extension, early favorable
stop or post-result gate change. Infrastructure interruption means incomplete.

## Power, not break-even

Set target **90% joint power**, allocating .1/8 failure probability per case.
Let q_c=sqrt(chi2_10(1-.05/8)) and q_p=sqrt(chi2_10(1-.1/8)).
Their squared sum is 94.18944697204289. Under the stated Gaussian stochastic
error approximation, variance cap and bias bound, a sufficient count is

    weak: N >= V*(q_c+q_p)^2 / (.205*||truth||-2*eta)^2,
    flat: N >= V*(q_c+q_p)^2 / (.98*.125-2*eta)^2.

Here .205=(1-delta)*.25-2*delta, and ||truth||=.5 for A=.25.
Bias affects both the center error and the confidence radius. Round UP to
the next multiple of 4096. The table uses these power counts, not the smaller
counts that merely make the expected interval touch the boundary.

The retained weak caps exceed their completed high-count V values by about
25--29%; the flat caps exceed their two V values by at least 32%. The new
selection cap 8000 exceeds the largest returned weak-stratum V in the old
512-case campaign (6622.5671) by about 21%. These are conservative planning
assumptions, not proved population envelopes. Their observed checks cannot
retroactively certify unconditional power or the Gaussian approximation.
The completed 512-case coverage test supports the earlier nominal procedure
on its mixture at N=262144; it does not itself validate these new high counts.

## Execution, provenance and audit

Adopt GPU acceleration only after the separate unchanged-bit engineering
gate passes on all eight old prefixes and the complete pure-fit replay.
Preserve the first packing-API engineering failure. Bind the successful
receipt, its reachable source commit, current source hashes, all retained
input/truth hashes, exact case table, random seeds and runtime before data.

The successful engineering receipt is
outputs/intrinsic_weyl_gpu_equivalence_v2.json, SHA-256
baa7a69e68f2becb98eb41e776618a491618375fdde64483f5386408112a97c9.
Before any new scientific sample, also replay old high-count case 4 at N=262144
through this runner's actual fresh worker and its enforced archive-only loader.
Compare every finite fit, inference, continuum fit and evaluation field to the
old wrapper at the existing 1e-8 absolute/relative tolerance, with exact discrete
flags. Preserve its pure fit and full comparison. A mismatch stops before data.

Keep OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1, and unset the other three
old BLAS overrides. Use one generator and one sequential fresh estimator
process. GPU query expressions and guards are the audited float64 port;
CPU rejection sampling is unchanged. Only exact integer-weighted covariance
products move to GPU. The mathematical inverse and CPU solver are unchanged.

Use six-field uncompressed NPZ archives to avoid compression cost; preserve
all order bits. The estimator may open only its order archive until the pure
fit is saved. Import truth evaluators only afterward. Store the original 95%
evaluation wrapper as a reference and the new simultaneous-precision fields
separately; do not relabel an old coverage verdict.

Record sampled/classified events, direction predicates, wall/CPU time, peak
working set and archive bytes. Require at least 30 GiB free storage before
starting; predicted packed data are below 20 GiB even at the selector cap.
On guarded generation failure retain the classified prefix and the failure,
not an uninitialized full archive. Original source/data must stay unchanged.

After all eight outcomes, independently recompute the population, power
counts, confidence level, tensor errors, bias, covariance radii, variance caps
and interval-containment verdict. Verify all source/input/output hashes,
selection rules, causal transitivity, and a separate post-completion pure
replay of case 5 (if it failed, the first returned fit in case order).
Only then decide whether the weak-curvature precision objective has
been achieved on the stated tested regime. General weak-field universality,
new physical laws and superiority to published alternatives are not implied.
