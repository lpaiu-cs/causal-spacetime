# Fresh curvature-square pilot: bias repair, no precision gain

The [frozen eight-case comparison](square_pilot_protocol.md) completed from
clean commit **373ba57ea0fbce937c46c94c3b558b96a84dbb05**. All corrected and
reference fits returned successfully. All eight corrected cases meet the
unchanged 2% continuum-bias condition, as do these eight fifth-degree
references. This small development pilot has no coverage or accuracy verdict.

The correction reduces the continuum bias in five of six curved cases but
reduces actual finite-order error in only **two of six**. The median ratio
of corrected to reference largest variance is **1.000865**, with range
0.98556--1.01164. It is a useful bias correction, **not a demonstrated
precision lever**. Do not launch a large accuracy campaign solely on this
result or treat the previous coverage campaign as validation of this model.

| Case | A | w | Fifth-degree bias | Corrected bias | Fifth-degree finite error | Corrected finite error |
|---|---:|---:|---:|---:|---:|---:|
| 2 | .25 | 0 | 0.3205% | 0.0037% | 95.6884% | 93.0855% |
| 3 | .75 | 0 | 1.0807% | 0.0423% | 41.9117% | 42.6227% |
| 4 | .25 | .5 | 1.3033% | 1.2540% | 141.6054% | 145.2675% |
| 5 | .75 | .5 | 1.5937% | 1.2152% | 56.4262% | 55.3093% |
| 6 | .25 | -.5 | 1.1593% | 1.2662% | 86.7724% | 86.9110% |
| 7 | .75 | -.5 | 1.7601% | 1.2658% | 34.2889% | 35.9126% |

All errors use the full ten-component tidal norm modulo the residual frame
reflection. They are point errors, not interval-containment successes.
The two flat cases have zero continuum bias and absolute finite errors
0.38577 and 0.31179; relative error at zero is undefined. All eight nominal
confidence balls contain the truth, but their radii are 1.036--1.236, much
larger than a relative 25% accuracy target for either nonzero amplitude.

## Remaining bottleneck

In the eight corrected fits, **81.88--84.77%** of total nominal tidal
variance lies in the five magnetic components. The eigenvector of largest
variance has **99.36--99.97%** of its squared norm in that magnetic subspace.
These figures use the full tidal metric, including electric and magnetic
off-diagonal multiplicities; they are not raw coordinate diagonal sums.
The effective trace dimension, `trace(Sigma)/lambda_max(Sigma)`, is 4.12--5.24.
This is a descriptive covariance audit of these cases, not a population bound.

The known curvature-square terms do not repair that information shortage.
Adding further curvature derivatives would add nuisance parameters, and
improved truncation order alone would not establish useful finite recovery.
The next candidate must demonstrate a reduction in magnetic/full-tensor
uncertainty at a fixed event budget, or justify the actual larger event
budget with independent uncertainty and bias validation. The earlier wider-
center selector failed that test; no untested selector is adopted here.

The corrected model remains a fifth-degree jet with a partial sixth-degree
correction. Second covariant curvature derivatives are omitted. The rotating
cases' remaining 1.22--1.27% bias is consistent with a remaining truncation
term, but these numbers do not prove its mechanism or a uniform bound.
No arbitrary-frequency, singularity or global-topology claim follows.

## Provenance and verification

- All 32 finite/continuum fits across the two methods have Weyl rank ten,
  nuisance rank 277, successful optimizers and no boundary flags.
- A separate pure-estimator replay of corrected case 5 matched **every saved
  fit field exactly**, with no generator, profile or truth module imported.
- Fifteen source hashes were checked against the reachable generating commit;
  all 32 archived input/full-fit hashes match. Each input contains only the
  allowed six order fields. The saved evaluation and confidence fields were
  recomputed from the full fits and compared independently.
- The constant-wave physical scaling and full complex-step Jacobian checks
  passed before data. Additional generic-tensor rotation/reflection and
  world-term boost checks have discrepancies at most 5.43e-20; quadratic
  curvature scaling also holds. These check implementation covariance, not
  statistical coverage or a generic-spacetime remainder bound.

The sampled event count is **2,129,920**, plus 16 supplied tips. Sum of per-case
wall times under concurrent machine load is 671.40 seconds. No case, event,
seed or initialization was added after observing its results.

[Compact results and audit](square_pilot_summary.json) retain every outcome.
Raw archives and full fits are in `outputs/intrinsic_weyl_square_pilot_v1`.
The raw summary SHA256 is
`02f761dc339f694127543b9986e917bbc2b1fdedc8a6c2368feb19826ab79d95`.
The generating SHA and source hashes refer to the isolated worktree; the
original fifth-degree coverage process continues on its unchanged checkout.
