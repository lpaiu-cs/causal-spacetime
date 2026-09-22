# Completed independent coverage validation

2026-09-22. The frozen [protocol](coverage_protocol.md) ran to completion
from `e057d101bf98b7d04ea2f233e576b2ad58ac264f`. All 256 cases remain in
the denominator. **254/256 were valid and covered**. The exact two-sided
95% Clopper--Pearson interval is **[0.9720650858, 0.9990524657]**; its lower
endpoint exceeds the preregistered 0.95 threshold. The coverage gate is
**SUPPORTED for the frozen iid constant-wave mixture**.

This establishes a population-specific coverage/availability result. It does
not certify a 25% recovery error, each amplitude separately, arbitrary
landmark shapes, or varying profiles. All 254 completed fits satisfied the
explicit 2% continuum-bias condition; the largest relative bias was 1.55712%.

## Precision remains inadequate

Statistics below are descriptive, conditional on a returned fit. Errors
compare all ten Weyl components in the prescribed internal frame modulo
the remaining reflection. The true tidal norm is `2A`.

| A | Cases | Valid and covered | Median absolute error | Median relative error | Median radius | Maximum variance eigenvalue |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 98 | 98 | 0.46415 | undefined | 1.02389 | 0.090551 |
| 0.25 | 67 | 66 | 0.44770 | 89.54% | 1.01574 | 0.085039 |
| 0.75 | 91 | 90 | 0.48469 | 32.31% | 1.04428 | 0.108281 |

Even the smallest returned radius exceeds the 25% target radius in both
curved populations. Coverage is possible here partly because the balls are
wide. It must not be presented as useful-precision recovery. The earlier
development point errors and verdicts are unchanged.

The provisional later-study variance envelope `lambda <= .08` at this
measurement count was exceeded in **8** completed cases: 5 flat, 1 weak,
2 strong. Its maximum was .10828055. Therefore the preceding price based
on that envelope is not an established sufficient price for this population.
Do not silently enlarge the envelope or launch the expensive accuracy study
under its old assumptions. A new design must explicitly handle variation
in landmark information, bias, failure probability and powered precision.

## The two retained failures

Both cases failed before producing a fit because a clock selected from an
independent selection order had zero events in the measurement order.
Post-completion continuum calculations diagnose, but do not repair, them.

| Case | A | Seed | Landmark pair | Selection interval count | Expected measurement count | Exact conditional probability of zero |
|---|---:|---:|---|---:|---:|---:|
| 75 | .75 | 2026092375 | 15,36 | 2 | 1.65643239 | .19081753 |
| 151 | .25 | 2026092451 | 1,36 | 2 | .29782263 | .74243289 |

The last column is `(1-V_ij/V_pq)^262144`, using the actual selected
landmarks. These are sparse real intervals; neither zero requires a causal
predicate malfunction. The independent selection count does not guarantee
a well-populated measurement clock. Neither case was dropped or re-seeded.

## Provenance and independent read-only check

Raw order/truth archives and full fit arrays remain in
`outputs/intrinsic_weyl_coverage_v1`. The committed
`coverage_summary.json` preserves the complete summary with line endings
normalized to LF, including all case dispositions and artifact hashes.
The original raw summary's SHA-256 is
`1939a131288790615e7db91cd0f09798ebf01f6693e266ba768294119194d756`.
The committed LF copy's SHA-256 is
`ae84f3987c6017092f0c6255b944e1004d9c03d3fc8d733938860bbecdf520f3`.

`check_intrinsic_weyl_coverage.py` verified all 256 archive schemas, seeds,
population assignments and saved confidence calculations; all **766 artifact
hashes** and **8 source hashes** matched. The source hashes match both the
reachable generating commit and the files after completion. It replayed
case 255 through the pure order-only estimator without importing a generator
or opening a truth archive: maximum difference in coefficients, covariance
and inferred coordinates was **0.0**. It did not regenerate every order or
repeat every oracle bias calculation.

The recorded run used Python 3.14.0, NumPy 2.5.3 and SciPy 1.18.1. Sum of
per-case runtimes was 2906.47 seconds. The event count was 68,157,440 plus
tips. No extra measurement events were added after observing failures.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/check_intrinsic_weyl_coverage.py outputs/intrinsic_weyl_coverage_v1
```

Next work is actual precision improvement and finite-order validation on
varying waves. A new exact generator alone will not fulfill either task.
