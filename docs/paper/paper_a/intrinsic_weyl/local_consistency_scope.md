# What fitting the first curvature derivative can change

2026-09-22. This is a **conditional local error derivation**, not a new
geometrical identity, a proof of this selector/optimizer's consistency,
a minimax lower bound, or a finite-sample accuracy verdict. It explains
the possible substantive benefit of the enlarged inverse model and makes
the missing regularity and conditioning requirements explicit.

## Assumptions that cannot be dropped

Consider shrinking convex normal vacuum neighborhoods about a fixed regular
point, with physical tip separation h. Work in a bounded orthonormal frame.
The metric is sufficiently smooth for a uniform sixth-degree remainder;
curvature and its first two covariant derivatives remain bounded on these
neighborhoods. The interval-volume expansion also has a uniform remainder.
This excludes caustics, singular points, nonvacuum terms and an unbounded-
frequency family whose derivative bounds increase as h decreases.

Condition on independently selected landmarks with fixed bounded normalized
shape. Retained timelike intervals have volume fractions bounded away from
zero, and the two vectors defining the spatial frame remain nondegenerate.
After fixing the coordinate gauge, the joint geometry/curvature/gradient
Jacobian has uniformly full column rank in a neighborhood of the limiting
flat configuration. The limiting clock weight is positive definite on the
retained observation space, with bounded condition number. Require the fit
to stay on the local least-squares branch near the true geometry.

These last requirements are substantive. A numerical rank on a finite pilot
does not prove a uniform condition bound. The current connected selector
and four-start optimizer have not been shown to satisfy these assumptions
with probability tending to one. A shrinking covariance eigenvalue cutoff
cannot silently replace the positive-definite observation-space assumption.

## Derivation in normalized coordinates

Put physical normal coordinates x=hX and write the dimensionless target as
`c=h^2 C`, with nuisance derivative coefficients `d=h^3 nabla C`.
The retained normalized squared-time model has the form

```text
F_ij(X,c,d) = (X_i-X_j)^T eta (X_i-X_j)
             - c(X_i,X_j,X_i,X_j)/3
             - d_e(X_i,X_j,X_i,X_j)(X_i+X_j)^e/12.
```

The next world-function terms contain second curvature derivatives and
curvature products. Their normalized remainder is `O(h^4)` under the stated
bounds. This follows from the established expansion, including the order
convention, in [Brewin, section 9 and equations 11.21b--d](https://arxiv.org/html/0903.2087).

For a vacuum Alexandrov interval, the leading relative volume correction
is fourth order in its size; see [Wang, section VI, equation 79](https://arxiv.org/html/1904.01034).
Only its **order** is used here, not a numerical coefficient or a claim
that the curved diamond volume equals the flat volume. Consequently,
replacing normalized proper-time squares by `-sqrt(V_ij/V_outer)` also
adds `O(h^4)`, including the outer normalization error.

For N independent uniform physical-volume measurement events, each interval
frequency is an average of bounded indicators. With finitely many intervals
and probabilities bounded away from zero, their correlated square-root
clocks differ from the continuum clocks by `O_p(N^-1/2)`. Correlation does
not change that rate; treating overlapping intervals as independent would
still give an incorrect covariance. The estimator uses the shared-event
covariance. Its normalized weight must converge as assumed above.

Let theta=(X,c,d), and let J be the limiting gauge-fixed model Jacobian.
The local least-squares normal equation has derivative `J^T W J`, invertible
under the rank/weight assumptions. Expanding this equation around the true
geometry gives

```text
theta_hat - theta = O(h^4) + O_p(N^-1/2),
c_hat - h^2 C     = O(h^4) + O_p(N^-1/2).
```

This is the local-branch conclusion of the nonsingular normal equation,
not a guarantee that the implemented optimizer finds that branch. If an
independent physical clock supplies h, division by h^2 yields

```text
C_hat - C = O(h^2) + O_p(N^-1/2 h^-2).
```

Without that clock, only c is identifiable. For a fixed nonzero physical
tensor, the relative error of c has the same scaling, with constants
depending on `1/||C||`; no relative flat-case guarantee follows.

If the measurement density is rho in four dimensions, `N asymp rho h^4`:

```text
fifth-degree fit: O(h^2) + O_p(rho^-1/2 h^-4).
quartic fit:      O(h)   + O_p(rho^-1/2 h^-4), generically when nabla C != 0.
```

Thus a sufficient consistency regime for either local construction is
`h -> 0` and `rho h^8 -> infinity`. Balancing the two displayed terms gives
`h asymp rho^-1/12` and error `O_p(rho^-1/6)` for the fifth-degree model,
versus `h asymp rho^-1/10` and `O_p(rho^-1/10)` for the generic quartic model.
These are **conditional achievable-rate calculations**, with unknown
constants. They do not prove optimality, a 2% bias bound, or lower actual
cost on finite data. In particular, the extra nuisance parameters can
increase the finite variance, as the completed pilot actually observed.
For special constant profiles the quartic leading bias can already be
quadratic, so the generic exponent comparison does not apply there.

## Fixed deterministic check of the predicted bias orders

Before executing this diagnostic, fix the existing jet pilot case 7, its
landmarks and its clock mask/covariance weights. Use exactly
`s=[1,.5,.25,.125]`, with generator amplitudes `.75*s^2` and rotation rates
`-.5*s`. Fixed bounded dimensionless landmark coordinates approach a fixed
shape. Use exact continuum interval volumes and fit both unchanged models;
no true coordinates or coefficients enter either inverse. The separate
evaluator alone uses generator data to measure the errors.

Report every scale, ranks, bias and 32/64 quadrature discrepancy. This tests
whether halving s approaches bias ratios 2 and 4 respectively. It is a
deterministic oracle audit, not new finite-order data or a power-sized gate.
An unfavorable slope must be retained and investigated, not fitted away by
changing the ladder or landmark configuration. The independent 256-case
coverage campaign remains unchanged and separate.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_jet_scaling.py
```

## Completed deterministic check

The frozen four-scale diagnostic completed from
`6924a08ce91e9694337878bcc6bbdeb684847eb0`, without changing the ladder.

| Scale s | Quartic relative continuum bias | Fifth-degree relative continuum bias |
|---:|---:|---:|
| 1 | 3.60218% | 1.83314% |
| .5 | 1.71023% | .46700% |
| .25 | .87475% | .11693% |
| .125 | .44735% | .02921% |

Successive halving ratios are 2.1063, 1.9551, 1.9554 for the quartic model
and 3.9254, 3.9937, 4.0037 for the fifth-degree model. The latter's local
log2 slopes are 1.9728, 1.9977, 2.0013, consistent with the predicted
quadratic bias in this fixed configuration. All eight fits completed with
full coordinate/Weyl ranks; each fifth-degree nuisance rank was 277.

Historical source and input hashes were checked. The scale-one fifth-degree
fit was exactly equal to the previous pilot's complete continuum fit.
[The saved summary](jet_scaling_summary.json) records all scales, ranks,
ratios and hashes. The full raw output remains at
`outputs/intrinsic_weyl_jet_scaling_v1.json`, SHA-256
`42a18f38611424565a548eb8cd35c5d209548952c16db409957290487b8806bd`.
The LF summary hash is
`7957e827b0ee1f8a01d239dad85669b4222a17d1cb8cd4bbafec876eb36a7855`.

This supports the order-of-bias mechanism on one deterministic shape. It
does not validate the conditional density rate, the random selector,
confidence coverage, a universal 2% allowance, or finite-sample accuracy.
