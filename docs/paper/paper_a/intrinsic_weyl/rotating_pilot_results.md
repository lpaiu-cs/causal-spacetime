# Rotating-profile pilot: completed, bias transfer fails

2026-09-22. The [prospective protocol](rotating_pilot_protocol.md) and
generator were committed before random data at
`8490f2de0ae8cff5b422d33ae0688ac7aa1a9555`. All eight fixed cases completed,
with 64 interior landmarks and projected curvature rank ten. No case,
seed or measurement clock was replaced. The estimator received only the
same six order-archive fields as the earlier connected pipeline.

**All four nonconstant curved cases exceeded the previous 2% model-bias
allowance.** Their continuum relative errors were 2.70--4.74%. The present
quartic model therefore does not inherit that allowance from constant waves.
This is a model-bias finding, not a formal failed accuracy verdict: the pilot
had no accuracy gate and the low-count finite-order errors remain large.

## Complete outcomes

The norm is the combined electric/magnetic Frobenius norm in the prescribed
internal frame, modulo z reflection. Continuum fits use the same selected
landmarks, clock mask and empirical covariance weights as the finite fit.

| Case | A | w | Finite relative error | Continuum relative error | Absolute confidence radius | 2% bias condition |
|---|---:|---:|---:|---:|---:|---|
| 0 | 0 | 0 | undefined | numerical zero | 1.31575 | met |
| 1 | 0 | .5 | undefined | numerical zero | .96555 | met |
| 2 | .25 | 0 | 72.49% | .2064% | .91041 | met |
| 3 | .75 | 0 | 27.05% | 1.2441% | .90595 | met |
| 4 | .25 | .5 | 86.43% | 2.7006% | 1.23251 | exceeded |
| 5 | .75 | .5 | 24.95% | 4.7356% | 1.07538 | exceeded |
| 6 | .25 | -.5 | 133.81% | 3.7357% | .94355 | exceeded |
| 7 | .75 | -.5 | 37.79% | 2.8441% | 1.00676 | exceeded |

Flat absolute finite errors were .96471777 and .30604632. Their continuum
errors were numerical zero. All nominal confidence balls contain their
targets, but they are broad and four violate the assumption used to construct
them. This is **not** a transferred coverage guarantee. Case 5's point error
below 25% is **not** an accuracy success: its radius alone exceeds the
corresponding .375 target radius.

The constant-profile coverage result in [coverage_results.md](coverage_results.md)
stands unchanged. It was explicitly limited to its frozen iid mixture.

## Verification and provenance

The actual-diamond sampler and causal predicate use the independently
checked [rotating generator](rotating_generator.md), not a constant-profile
approximation. Before data, its propagator, geodesic inverse map and path
action agreed with independent calculations to approximately `1e-15`.
Every case's interval-volume quadrature refinement was below `3e-15`.
These are numerical checks in the stated regular patch, not arbitrary-domain
causal certification or directed-rounding error bounds.

Raw archives and full fits are in `outputs/intrinsic_weyl_rotating_pilot_v1`.
`rotating_pilot_summary.json` preserves the complete summary normalized to LF.
Original summary SHA-256:
`0372a0e600dc8ab79e01fdb877473d65ead3a95ed14ff37f6132a3c04906f3b1`.
Committed LF summary SHA-256:
`e20cb249ae36a91b69d304302353604c7989df1028da0998441c0a5b335dac8d`.

A read-only check matched all **24 artifact hashes** and **10 source hashes**
against the generating commit and current files, verified the fixed eight
cases/seeds and archive field boundary, and recomputed saved confidence
calculations. Case 7 was replayed through the pure estimator without a
generator or truth input; coefficients, covariance and inferred coordinates
matched exactly (maximum difference **0.0**). This check did not regenerate
every random order or repeat every continuum fit.

The run used Python 3.14.0, NumPy 2.5.3 and SciPy 1.18.1. Sum of per-case
runtimes was 239.57 seconds for 2,129,920 sampled events plus tips.

## Next lever: information in the landmark design, with a bias budget

The 254 completed constant-profile coverage fits have median
`trace(G_B Sigma_BB)/trace(G_E Sigma_EE)=5.0444` (range 3.5329--7.5884).
Thus magnetic uncertainty dominates the current tidal error budget.
Increasing sample count treats that inefficiency expensively and cannot
remove the varying-profile continuum bias.

The existing three-slice identity shows why central spatial spread is a
candidate lever. For the magnetic part and a central point `(0,z)`, the
quartic response against `(t,r)` contains `2t r dot (z cross Bz)`, quadratic
in the central radius. Wider central landmarks can increase that response,
but shorter or fewer connecting intervals can also increase clock noise;
the net gain must be measured, not assumed from radius alone.

A post-result feasibility check on the already-used coverage cases 0--5
held their outer landmark IDs fixed and expanded only the central candidate
radius cap from .20 to .30, retaining at least four clocks to each side.
Each had at least 185 eligible expanded central candidates. Ranking by
`minimum_side_degree * radius_squared` gave 22 candidates with median
order-derived radius .249--.273 (previously .173--.191); their minimum
side degree was 4--6. This reuses development data and measures candidate
availability only. It does **not** establish lower variance, full rank after
joint fitting, or controlled bias for the new geometry.

The next bounded comparison should freeze the old and expanded selectors
before fresh data and hold event count constant. It must measure all-component
variance and error, preserve the same outer frame, and audit continuum bias
for both constant and varying waves. Do not enlarge the old 2% allowance
or reinterpret the present cases as successes. A regularity/window condition
or an explicitly new local model is still needed before a powered accuracy
claim beyond constant profiles.
