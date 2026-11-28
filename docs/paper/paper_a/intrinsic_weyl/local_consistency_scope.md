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

## Finite-count scope of the square-root approximation

The total event count alone does not justify the clock approximation: the
expected count of each retained interval matters. Conditional on fixed
landmarks, one clock has `K ~ Binomial(N,p)` and `Y=-sqrt(K/N)`. Directly,

```text
E[Y^2] = p,
Var(Y) = p - (sum_k sqrt(k/N) BinomialPMF(k;N,p))^2.
```

Compare this with the known-p first-order value `v_delta=(1-p)/(4N)`.
The implemented diagonal instead plugs in K/N and uses N-1 from the
unbiased indicator covariance; neither is an exact transformed variance.
For the stipulated probability `p=3/262144`, numerical binomial sums give:

| N | Expected K | P(K=0) | Var(Y)/v_delta | Var(Y given K>0)/v_delta |
|---:|---:|---:|---:|---:|
| 262144 | 3 | .0497862 | 1.3599625 | .8445185 |
| 4194304 | 48 | 1.42477e-21 | 1.0080554 | 1.0080554 |

The sums used SciPy 1.18.1 binomial PMFs between the 1e-14 lower and upper
tail quantiles. Omitted probability was below 1e-14 in these two rows;
centered-moment calculations independently agreed with the displayed second-
moment identity to absolute 3e-19. No random samples were drawn.

These are illustrative marginal calculations at stated p, not an estimate
of p from an observed count of three or a diagnosis of pilot case 6. The
inverse rejects zero clocks, and the frozen reliability gate counts those
failures; conditioning just on one positive clock is not conditioning on a
successful full fit. Marginal moment correction would not establish joint
covariance or ten-component confidence coverage. The calculation changes
no clock mask, radius, failure accounting or ongoing coverage protocol.

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

## A finite-sample route for an independently fixed affine inverse

2026-09-23. This is a conditional statistical construction, not an implemented
replacement for the nonlinear inverse or a new coverage/accuracy verdict.
The ongoing 512-case campaign and all earlier failures remain unchanged.

Condition on the independently selected landmarks and a separate pilot fit.
For each subsequent measurement event, let I be the vector of interval
membership indicators. Its entries may be arbitrarily correlated. Freeze a
matrix A and offset b using the pilot only, and put

```text
Z = b + A I,       mu = E[Z | pilot],       estimate = mean_i Z_i.
```

Use ten orthonormal tidal components, so Euclidean distance is the existing
electric-plus-magnetic Frobenius distance. A row's range length is bounded by
`R_j = sum_k |A_jk|`, regardless of which indicator patterns are realizable.
Known constant indicators can be removed before taking this bound. Neither
the observed range nor an observed minimum interval count supplies a valid
population range bound. Conditional on the pilot, different sampled events
are iid; different clocks or components need not be independent.

This permits a finite-sample confidence construction. With m independent
calibration events, let s_j be the usual sample standard deviation of Z_j.
For K fixed cases and d=10 components, set a calibration failure budget eta.
Theorem 10(3) of [Maurer and Pontil (2009)](https://www.cs.mcgill.ca/~colt2009/papers/012.pdf)
gives simultaneous variance upper bounds, with probability at least 1-eta:

```text
u_j = min(R_j^2/4,
          (s_j + R_j sqrt(2 log(K d / eta)/(m-1)))^2).
```

Here the first term is the universal variance bound for an interval of length
R_j. Calibration must be independent of both pilot fitting and the final
measurement sample. For a sampling error budget gamma, define

```text
L_gamma = log(2 K d / gamma),
a_j = R_j L_gamma / (3 N),
q_j(N,gamma) = a_j + sqrt(a_j^2 + 2 u_j L_gamma/N),
r(N,gamma) = sqrt(sum_j q_j(N,gamma)^2).
```

Inverting the two-sided scalar Bernstein inequality and taking a union bound
over components and cases yields `||estimate-mu|| <= r(N,gamma)` except on
probability at most eta+gamma. The scalar inequality is stated in section 3.1
of [Martinez-Taboada and Ramdas (2024), v1](https://arxiv.org/html/2409.06060v1#S3.SS1).
This use does not require their more elaborate vector confidence sequence.
N may be selected from calibration results, before opening the independent
final sample; conditioning on calibration preserves the fixed-N argument.
It does not authorize optional stopping on the final sample.

The bound is for mu, not yet for the physical tensor theta. If a separately
justified bound `||mu-theta|| <= delta ||theta||`, delta<1, holds in the
chosen internal frame (or one fixed representative of its reflection orbit),
an observable radius for theta is

```text
R_theta = (r(N,gamma) + delta ||estimate||)/(1-delta).
```

The denominator follows from the triangle inequality; it is not optional.
For simultaneous 95% coverage and 90% planning power, one example allocates
eta=.01, gamma=.04 for coverage, and gamma_power=.09 for the noise event used
in the power proof. On the latter event,

```text
||estimate-theta|| + R_theta
 <= (r(N,.04) + r(N,.09) + 2 delta ||theta||)/(1-delta).
```

Thus the existing 25% containment target has a sufficient condition
`r(N,.04)+r(N,.09) <= .205 ||theta||` when delta=.02. This is a
power inequality, not just a break-even width. A future protocol would have
to fix calibration size, the count-selection rule, all error allocations,
and a maximum affordable N before data; exceeding the cap is a failure of
feasibility, not permission to drop difficult cases. No numeric sample price
is established here. End-to-end coverage and power claims must also account
for pilot, bias and cap failures wherever they prevent the required output;
the displayed noise calculation alone does not do so. Conditioning on a
favorable calibration or a successful pilot cannot establish unconditional
90% power. Flat targets require an absolute
bias/accuracy treatment; the relative formula cannot silently supply one.

The possible computational benefit is also specific. After A is frozen,
the final estimate needs only the interval counts, since
`mean_i Z_i = b + A mean_i I_i`. Calibration can process ten-dimensional
event scores. The large-count phase need not reconstruct the full clock-by-
clock covariance or repeat the nonlinear fit. This is an algebraic work
reduction, not a measured runtime or memory saving. Pilot fitting, calibration,
relation queries and any retained order archives must still be charged.

### Bounded diagnostic fixed before its calculation

Before considering fresh events or adopting this construction, use **all eight**
existing nested-pilot cases 0--7, including the missed ball in case 6, at
`outputs/intrinsic_weyl_nested_pilot_v1`. Draw no new events. The pure stage
receives each saved order archive and the previously returned `fit` subobject,
not its truth/evaluation fields. Keep its selected mask, fifth-degree model,
curvature-square terms, whitening cutoff and gauge.

At that pilot solution theta0, write f0 for the negative squared-time clock
model, J for its full geometry/curvature/gradient Jacobian and W for its saved-
data whitening map. Let E_C select the ten curvature parameters. Freeze

```text
A_C = E_C (W J)^+ W diag(1/(2 f0)),
b_C = C0 - A_C f0^2.
```

Require f0 to be strictly negative and its square to be a probability.
Canonicalize the pilot frame and apply the existing tidal norm transform to
A_C and b_C. Verify the local left-inverse identity to absolute 1e-7 and
replay the saved clock-linearized covariance to relative 1e-8. Preserve
rank, domain or replay failures; do not remove troublesome clocks.

Freeze and hash the resulting map before the evaluator opens the saved
truth archive. The evaluator then computes the same exact interval-volume
fractions at quadrature orders 32 and 64, requiring their existing 1e-10
relative agreement. Report the affine conditional mean b+A p, its full tidal
bias, every component range and its empirical score standard deviation, and
all eight outcomes. Label the latter variances descriptive: the same pilot
sample chose A, so they are not the independent calibration bounds above.

Compare mean bias with the **unchanged** 2% nonzero-target allowance and 1e-10
flat-zero tolerance for this diagnostic. A violation must be retained and
prevents claiming that the old bias condition transfers. This mean bias
includes pilot linearization error; the old continuum nonlinear-fit bias
alone does not bound it. The diagnostic cannot certify power, future coverage,
or recovery, and must not re-score the ongoing nonlinear campaign. Quadrature
agreement is a numerical check, not a rigorous bias remainder bound. Numerical
rounding error also needs a budget before a finite-sample implementation can
claim the mathematical guarantee.

### Completed eight-case affine diagnostic

The preceding calculation plan was committed as `081f0ae` before execution.
The implementation and its synthetic derivative/covariance check were
committed as `aea7eb61b286d415d8af6a2e648499ffd490b874`; all eight maps were
saved and hashed before any truth archive was opened or oracle module
imported. No new events were generated. The full result and two read-only
audits are preserved in [affine_diagnostic_result.json](affine_diagnostic_result.json).
The raw diagnostic summary SHA256 is
`4b2b7335611bf004292963b17c665d9357176c1e20c952f34ab50d511f14da13`.

All eight maps passed the domain, rank and covariance-replay requirements.
The largest local left-inverse discrepancy was 2.85e-13 and the largest
relative covariance discrepancy was 5.51e-14. Exact-volume quadrature
differences were below 2.89e-15. The unchanged bias condition held in **five
of eight** cases:

| Case | Amplitude | Rotation | Absolute mean bias | Relative mean bias | Old bias condition |
|---|---:|---:|---:|---:|---|
| 0 | 0 | 0 | 0.00569195 | undefined | FAIL |
| 1 | 0 | +0.5 | 0.00709783 | undefined | FAIL |
| 2 | 0.25 | 0 | 0.00517963 | 1.03593% | met |
| 3 | 0.75 | 0 | 0.00600480 | 0.40032% | met |
| 4 | 0.25 | +0.5 | 0.00672979 | 1.34596% | met |
| 5 | 0.75 | +0.5 | 0.01666881 | 1.11125% | met |
| 6 | 0.25 | -0.5 | 0.01418166 | 2.83633% | FAIL |
| 7 | 0.75 | -0.5 | 0.02312156 | 1.54144% | met |

These are conditional means with exact probabilities, so increasing only the
future measurement count cannot remove these discrepancies with the pilot
map held fixed. The flat cases violate the fixed 1e-10 tolerance. Case 6
violates the 2% allowance even without future sampling noise. This finding
rejects transfer of the old bias condition to this particular affine
construction at the existing pilot size; it does not rule out all one-step
estimators, better pilots, or separately justified bias bounds.

A post-result comparison with the previously saved nonlinear continuum fits
locates the additional error. The old nonlinear mean satisfies the same
bias condition in all eight cases. In case 6 its bias is 0.00640067 (1.28013%),
while its full-tidal distance from the affine mean is 0.01038240. These vector
errors do not add as scalars; their combined bias is the 0.01418166 above.
Both flat nonlinear means are zero, whereas their affine means are not.
This comparison used saved results only, performed no new fits, and does
not re-score either pipeline's original finite-order outcome.

An independent direct event-score replay of case 6 reconstructed all 262144
score vectors from the packed order bits and the frozen map. Its empirical
covariance agreed with the delta-method reconstruction to 2.10e-15 relative
error and its mean to 3.51e-14 absolute error. This rules out a covariance
conversion/counting mismatch in that replay; it is not an independent
calibration sample. The observed score ranges are not population bounds.

There is also a separate precision obstacle: across the eight cases the
valid triangle-inequality range bounds are 1188--2297 times the same-pilot
empirical component standard deviations. A valid bound can be this loose;
this ratio is not a sample-price estimate or proof that all concentration
methods are unaffordable. No calibrated finite-sample radius or power price
has been obtained. The diagnostic therefore does **not** justify replacing
the nonlinear inverse, claiming a measured cost saving, or running the
large-count accuracy study with this map. No threshold or clock was changed
to rescue its failed cases.

The statistical ingredients are established prior work. One-step estimation
from an initial estimator is described, with its regularity assumptions, by
[Brouste, Dutang and Noutsa Mieniedou (2021)](https://journal.r-project.org/articles/RJ-2021-044/).
First-order nuisance cancellation and sample splitting also have extensive
precedent in [Chernozhukov et al.](https://arxiv.org/abs/1608.00060). Neither
source establishes the geometric remainder bound needed here. The possible
research contribution is a validated order-only Weyl instrument and its
cost, not rediscovery of those statistical methods.
