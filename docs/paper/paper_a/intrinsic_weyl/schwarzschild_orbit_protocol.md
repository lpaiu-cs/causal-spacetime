# Independent electric and magnetic parts: a circular-orbit Schwarzschild test of the unchanged order-only Weyl inverse

2026-11-21. Freeze this protocol, the generator, the runner and the checking
script before any test selection or measurement event. The engineering replay
of section 4 uses a different interval and seed and is not test data.

## 1. Question

Every test of the fifth-degree order-only inverse so far had one of two
special tidal structures. The plane waves are of Petrov type N, with electric
and magnetic parts of equal norm, `|E| = |B|`. The
[Schwarzschild test](schwarzschild_results.md) put the tips on a radial
geodesic at rest at the apex, so that its observer saw a purely electric field,
`B = 0`. The Weyl note therefore says that "a field with independent electric
and magnetic parts ... was not tested".

This protocol asks one question about that gap: **given only the order of an
interval whose tips move on a circular Schwarzschild geodesic, does the same
frozen inverse (a) detect the field, (b) detect its magnetic part, while
detecting none at the same event for an observer at rest there, and (c) resolve
that the electric part is larger than the magnetic one?** Claim (c) concerns
the invariant `|E|^2 - |B|^2`, which is the same for every observer at the apex;
it says something about the moving observer's split only together with (b).

The design uses one event and two observers. At the apex the Weyl tensor is one
geometric object; its split into `E` and `B` depends on the observer, and the
estimand of the inverse is that split in the frame of the tips (Sect. 2 of the
note: time axis along the tip direction). The observer on the circular geodesic
sees `|E| != |B| != 0`. The observer at rest at the same event sees the field of
the earlier test, `B = 0`. Nothing on the inverse side changes.

Scope, stated now so that no sentence later outgrows it: the magnetic part here
is that of a static type-D field seen by a moving observer. It is a genuine
magnetic part of the Weyl tensor in the frame the inverse reconstructs, and `E`
and `B` are not tied together as in a plane wave, but this is not the
gravitomagnetic field of a rotating source. It is one geometry, one event, one
interval size and two cases per observer. It is not an accuracy test.

## 2. Geometry, observers and truth

Isotropic coordinates as in the static test, `M = 1`:

    ds^2 = -A^2 dt^2 + B^4 dX.dX,  A = (1 - M/2rho)/(1 + M/2rho),  B = 1 + M/2rho.

The apex is the event `(t, X) = (0, rho_a, 0, 0)` with areal radius `r_a = 6M`,
`rho_a = 4.949490`. Three classes of interval share it, all with tip proper
time `T = 3M` and all built by one generator class:

| Class | Tips' geodesic through the apex | Tips `P, Q = (t, x, y, z)` |
|---|---|---|
| orbit | circular, in the `xy` plane, speed `v = sqrt(M/(r_a - 2M)) = 1/2` along `+y` relative to the static observer | `(-/+2.121320, 4.898022, -/+0.711919, 0)` |
| rest | radial, momentarily at rest at the apex (the static test's geometry) | `(-/+1.838721, 4.917860, 0, 0)` |
| flat | `M = 0`, straight, the same velocity `v = 1/2` along `+y` | `(-/+1.732051, 6, -/+0.866025, 0)` |

For the orbit the coordinate half-time is `(T/2)/sqrt(1 - 3M/r_a) = 1.5 sqrt(2)`
and the azimuth of the tips `+/- Omega t = +/-0.144338` with
`Omega = sqrt(M/r_a^3)`. Tip proper time is the unit, as in every earlier test.

**Truth.** In the static orthonormal frame at the apex (MTW sign, so that
`E_ij = C_0i0j`), `T^2 E = (T^2 M/r_a^3) diag(-2, 1, 1)` and `B = 0`, with
`T^2 M / r_a^3 = 1/24`. The comoving frame of the orbit is the static frame
boosted by `v` along the direction of motion. With `B_ij = eps_ikl C_kl0j / 2`
(the note's definition) and `K = gamma^2 v^2 = 1/3`, in axes (radial, direction
of motion, `z`):

    orbit:  T^2 E = (1/24) diag(-(2 + 3K), 1, 1 + 3K) = diag(-1/8, 1/24, 1/12),
            T^2 B = (3 gamma^2 v / 24)(e_x e_z + e_z e_x),  B_xz = 1/12;
    rest:   T^2 E = diag(-1/12, 1/24, 1/24),  B = 0;
    flat:   0.

| | `|T^2 E|` | `|T^2 B|` | full norm `s` | `T^4(|E|^2 - |B|^2)` | `E.B` | eigenvalues of `T^2 E` | eigenvalues of `T^2 B` |
|---|---:|---:|---:|---:|---:|---|---|
| orbit | `sqrt(14)/24 = 0.15590` | `sqrt(8)/24 = 0.11785` | `sqrt(22)/24 = 0.19543` | `6/576 = 0.010417` | 0 | `-0.1250, 0.0417, 0.0833` | `-0.0833, 0, 0.0833` |
| rest | `0.10206` | 0 | `0.10206` | `6/576` | 0 | `-0.0833, 0.0417, 0.0417` | 0 |

`|B|/|E| = sqrt(4/7) = 0.756` for the orbit. `|E|^2 - |B|^2` and `E.B` are
invariants of the Weyl tensor, the same for both observers at the apex.

The truth was confirmed before any data in five ways, each reaching it from the
metric by a different route (generator selfcheck, Sect. 3): the closed form above; the
coordinate Riemann tensor from nested finite differences of the metric,
contracted with the comoving tetrad (`2e-10`); the relative motion of particles
released at the apex with the tips' four-velocity (electric part, `1.5e-5`
relative); the world function itself, i.e. the estimand of the inverse, between
geodesically shot events in the generator's normal coordinates, which returns
all ten components including the sign of `B` (`1.2e-8` relative after two-scale
extrapolation); and the proper volume of a small moving interval, whose excess
over `pi T^4/24` follows `T^4 (17|E|^2 + 3|B|^2)/10080` to `2.4e-4` at `T = M/2`.
The volume law weights `|E|^2` and `|B|^2` differently, while their difference
is the boost invariant `6 (T^2 M/r_a^3)^2` known from the static frame, so this
last check determines both norms without using the boost formula.

**Variation across the interval.** Both intervals reach areal radii from about
`4.82M` to `7.25M` (boundary found by bisection on 2 x 6144 rays: orbit
`4.819M`-`7.246M`, rest `4.817M`-`7.248M`), over which `M/r^3` ranges from
`1.93` to `0.57` times its apex value; in the comoving frame the field also
turns with the orbit. The inverse models a constant Weyl tensor, its
24 vacuum first derivatives and the curvature-square terms; second and higher
derivatives are not in the model. Exact interval volumes between landmarks are
not implemented for this geometry, so the continuum bias is not measured
separately: it is part of the reported error, as in the static test.

## 3. What is new and what is unchanged

New, generator side only: `experiments/positive_control/schwarzschild_orbit.py`.
It subclasses the frozen generator `schwarzschild_diamond.py`, which is not
modified; the metric, the optical-distance solver, the geodesic integrator, the
Chebyshev surrogates of the squared optical distance and the guarded relation
decisions are inherited unchanged. Every place where the old class assumed
radial tips at rest is replaced:

- **Tips.** The geodesic through the apex with four-velocity
  `u = gamma (e_t + v e_y)` in the static orthonormal frame, shot to proper time
  `T/2`; the past tip is the image of the future one under
  `(t, x, y, z) -> (-t, x, -y, z)`, an isometry that fixes the apex and reverses
  `u` (checked against backward integration: identical to the bit).
- **Region and surrogate domain.** With foci `P, Q` (the tips' spatial points)
  the spatial section is the optical ellipsoid `d(P, X) + d(X, Q) <= t_Q - t_P`.
  It lies in the coordinate ball about `c = (P + Q)/2` of radius
  `t_h / n(|c| + t_h)`; that ball replaces the old one about the tips' common
  point, bounds the surrogates and their guard test points, and is checked to
  contain the boundary found by the direct solver.
- **Sampling.** The spatial density is `A B^6 h(X)` with
  `h = t_Q - d(X, Q) - t_P - d(P, X)`, by rejection from a proposal proportional
  to `1 - |X - c|/R`, whose envelope `h <= 2 t_h (1 - |X - c|/R)` follows from
  `|X - P| + |X - Q| >= 2|X - c|`; `t` is then uniform on
  `[t_P + d(P, X), t_Q - d(X, Q)]`. A violated envelope raises an error; the
  runner records that case as failed and continues.
- **Quadrature.** The old quadrature was axially symmetric about the radial
  line. Moments of the measure now split into one term per focus, because the
  time range above `X` is bounded by a function of `d(P, X)` below and of
  `d(X, Q)` above; each term is smooth along rays from its own focus. Ray
  lengths come from bisection with the direct solver.
- **Frame.** Normal coordinates at the apex in the comoving frame
  `(u, e_x, gamma (v e_t + e_y), e_z)`, the static ones boosted along `y`.
- **Truth.** The static-frame Weyl tensor boosted by `v` (Sect. 2).

Selfcheck (deterministic, run by the runner before data; values from the
uncommitted tree, about four minutes):

- at `v = 0` the class is the old generator: tips, bounding ball, a sample of
  4096 events from a common seed, normal coordinates and the tidal tensor are
  equal to the bit;
- orbit tips on the circular geodesic (`2e-15` in radius, time and azimuth),
  at `(-/+1/2, 0, 0, 0)` in normal coordinates (`5e-15`);
- the five confirmations of the truth listed in Sect. 2;
- relation bits from the surrogate equal those of the direct solver for the
  orbit and the moving flat class, every sampled event lies inside the tips'
  interval, and the order is antisymmetric and transitive;
- on a grid of 16 x 24 directions from either focus, sampled at 64 radii,
  every ray crosses the boundary once, inside the bounding ball, and the
  largest rejection ratio on 20,000 probe points is `0.76`;
- six moments of `2e6` orbit events, among them `<t^2>` and the tilt `<t y>`
  that only a moving interval has, agree with the quadrature (`|z| <= 0.48`);
- the orbit volume is the same on 24x32x24 and 32x48x32 nodes per focus
  (tolerance `1e-12`, observed identical), and at `v = 0` the new quadrature
  reproduces the old axially symmetric one to `3e-15` (tolerance `1e-12`);
- `normal_coordinates`, whose frame the runner uses to express the truth,
  returns the vectors along which events were shot with the comoving tetrad, to
  `3e-12`: the frame of the truth checks and the frame of the truth coincide;
- the moving flat class reproduces Lorentz-boosted tips, normal coordinates and
  the volume `pi T^4/24` (`2e-15`);
- the parent class still passes its own selfcheck.

The runner adds one check of its own before data: for each of the three
classes, the relation bits evaluated on the GPU, which are the ones archived,
equal those of the direct solver on 262,144 pairs.

Unchanged: `selection(order, design="nested")`, the interval clocks, the exact
shared-membership covariance, `fit_order_archive(..., model_order=5,
square_correction=True)`, four starts and all optimizer settings, the
`confidence_ball` construction, and the blind boundary: the pure fit is saved
before any generator or truth module is imported, and only the assigned order
archive may be opened. No coordinates, mass, radius, velocity or paired flat
sample enter the inverse.

## 4. Engineering replay, not test data

One orbit case on a larger interval at the same apex, `T = 4M`
(`T^2 E = (2/27) diag(-3, 1, 2)`, `T^2 B_xz = 4/27`, `s = 0.347`), seed
`2026112290`, at `N = 2^20, 2^22, 2^24`, to check plumbing, run time, memory and
the checking script. It ran from the uncommitted tree before this freeze. The
pre-data selfchecks took about five minutes; selection and `2^24` events took
642 s while another session kept about nine of the sixteen CPU threads busy;
the three fits took 69, 82 and 121 s with a peak of 12.2 GiB. The selection
kept 129 interior landmarks (131 with the tips) and 3811 clocks; no measurement
margin needed the direct solver.

| N | `Q` | `Q_B` | 95% interval of `I` (truth `0.0329`) | relative error | nominal 95% radius / `s` | ball contains truth | `Q` about truth | `N lambda_max` | fit `chi2` z |
|---:|---:|---:|---|---:|---:|---|---:|---:|---:|
| 1,048,576 | 69.1 | 10.4 | `[-0.0009, 0.1021]` | 37.7% | 68.6% | yes | 10.6 | 2965 | 1.88 |
| 4,194,304 | 317.0 | 79.7 | `[-0.0014, 0.0563]` | 14.3% | 35.5% | yes | 7.1 | 2974 | 1.03 |
| 16,777,216 | 1156.2 | 297.2 | `[0.0073, 0.0350]` | 7.5% | 18.7% | yes | 9.2 | 2969 | 2.11 |

At `2^24` the estimated eigenvalues are `(-0.204, 0.068, 0.137)` for `E`
against `(-0.222, 0.074, 0.148)`, and `(-0.149, -0.000, 0.149)` for `B` against
`(-0.148, 0, 0.148)`. The `S^-1`-weighted readings along the true pattern are
`0.938 +/- 0.028` (full), `0.916 +/- 0.032` (electric half) and
`1.004 +/- 0.058` (magnetic half): the magnetic part is read on target and the
electric part about 8% low, 2.6 nominal sigma. That low electric reading is why
`I_hat = 0.021` is 64% of the truth: `I` is a difference of two squared norms,
so a bias in one of them is amplified. The static test's engineering replay, at
the same `T = 4M`, read `0.94 +/- 0.05` of its purely electric tensor. The
unmodelled curvature variation scales as `T^2`, so an electric bias of about
`(3/4)^2 x 8% = 4.5%` would be expected at the test's `T = 3M` if it has the
same origin. This is a reading at one selection, not a measurement of the
continuum bias, which cannot be computed for this geometry; no correction is
applied.

The checking script's functions were exercised on this replay (not as a
verdict): the truth from areal coordinates matched the generator's to `2e-13`;
the run's truth vectors and those from the O(3) alignment differed by at most
0.25% of the norm; all recorded statistics were reproduced; 998 of 998 order
bits checked with the Clairaut solver, among them the 500 pairs nearest the
light cone, matched, none inside its error bound; the moments of the first
`2^20` events against the quadrature with Clairaut bounds and the areal volume
element had `|z| <= 1.48`, the Clairaut bounds differing from the generator's by
at most `6e-9` and putting the generator's boundary nodes on the null boundary
to `1.2e-8`; nested prefixes and a CPU replay of the first 4096 events were
bit-exact; and a fresh single-threaded refit of the `2^20` row was identical
field for field.

The committed runner must reproduce this replay field for field, from a clean
worktree, before any test event is drawn; a mismatch stops the run.

## 5. Assigned cases and independent data

All six cases are assigned at every level.

| Case | Class | Seed |
|---:|---|---:|
| 0 | orbit | 2026112200 |
| 1 | orbit, independent selection and events | 2026112201 |
| 2 | rest | 2026112202 |
| 3 | rest, independent | 2026112203 |
| 4 | flat, moving tips | 2026112204 |
| 5 | flat, moving tips, independent | 2026112205 |

Each case draws 4096 proper-volume selection events from
`default_rng([seed, 0])`, appends the two tips, permutes interior labels and
applies the nested selection. Measurement events come from
`default_rng([seed, 1])` in chunks of 4096. The levels are nested prefixes of
one stream per case, the four of the earlier tests:

    N = 1,048,576;  4,194,304;  16,777,216;  67,108,864.

At most `6 * 67,108,864 = 402,653,184` measurement events and 24,576 selection
events are drawn. Results at different levels of one case are not independent.

## 6. Frozen gates

For a completed fit let `c` be the ten fitted Weyl coefficients and `S` their
nominal covariance, both after `confidence_ball`'s axis-sign convention, with
`c_E, c_B` the electric and magnetic halves and `G` the Gram matrix of the five
symmetric basis tensors, so that `|E|^2 = c_E^T G c_E`. A fit is valid only if
`confidence_ball` accepts it (optimizer success, interior solution, full
coordinate, nuisance and projected curvature rank, positive definite
covariance); an invalid fit is a failure of its row. The statistics below use
the fit alone and are saved before the truth is evaluated; a failure in the
truth step is recorded but does not fail the row.

Statistics:

- `Q = c^T S^-1 c`, the detection statistic of every earlier test;
- `Q_B = c_B^T (S_BB)^-1 c_B`, the same statistic for the magnetic half alone,
  with the marginal `5 x 5` block of `S`;
- `I = |E|^2 - |B|^2`, estimated without the noise bias of the squared norms,
  `I_hat = c^T A c - tr(A S)` with `A = diag(G, -G)`, with standard error
  `sqrt(4 c^T A S A c + 2 tr(A S A S))` (the Gaussian variance of a quadratic
  form with the estimate in place of the mean, which overstates it on average)
  and the interval `I_hat +/- 1.959964 se`.

Thresholds: `q*_10 = chi2_10(1 - 1e-3) = 29.5883` and
`q*_5 = chi2_5(1 - 1e-3) = 20.5150`. All three statistics are invariant under
rotations and reflections of the spatial frame, so none needs the reflection
convention.

**G1, field (cases 0, 1, 4, 5).** For each level: DETECTED if `Q > q*_10` in
both orbit cases and in neither flat case; MIXED if exactly one orbit case
exceeds and no flat case does; NOT_DETECTED if neither does and no flat case
does; CONTROL_FAILED if a flat case exceeds; INCOMPLETE if any of the four
failed.

**G2, magnetic part (all six cases).** For each level: MAGNETIC_DETECTED if
`Q_B > q*_5` in both orbit cases and in none of the four controls (two rest,
two flat); MIXED if exactly one orbit case exceeds and no control does;
NOT_DETECTED if neither does and no control does; CONTROL_FAILED if any control
exceeds; INCOMPLETE if any of the six failed. A rest control that exceeds means
that the inverse reads a magnetic part where the tips' observer sees none; the
level then supports no statement about the magnetic part, whatever the orbit
cases show. A rest control that does not exceed is a non-detection, not a proof
of absence: at `2^26` it detects a spurious magnetic part of norm `|T^2 B|`
0.026 with probability 0.5 and 0.036 with probability 0.9 (calibration
covariances, random magnetic patterns; 22% and 30% of the orbit's 0.118), at
`2^24` 0.052 and 0.071. The rest interval also has a weaker, axially symmetric
electric field, so it does not control for leakage peculiar to the orbit
interval.

**G3, unequal norms (cases 0, 1).** For each level: UNEQUAL if the 95% interval
of `I` lies above zero in both orbit cases; REVERSED if it lies below zero in
either (the inverse would then resolve `|B| > |E|`, the wrong ordering); MIXED
if exactly one lies above and none below; NOT_RESOLVED if neither; INCOMPLETE if
either failed.

Precedence, for all three gates: bad news first. CONTROL_FAILED (G1, G2) or
REVERSED (G3) is decided on the valid fits and reported even if another row of
the gate failed; only then does a failed row make the gate INCOMPLETE; then the
remaining outcomes. Exceedance is strict: a statistic equal to its threshold
does not exceed it, and an interval bound equal to zero does not lie above or
below it.

G3 alone says nothing about the moving observer, because `I` is the same for
every observer at the apex; an inverse that returned the resting observer's
tensor for the orbit interval would pass it. It is read only together with G2.

**Primary claim, fixed now.** The headline of this test is the conjunction at
the top level, `N = 67,108,864`: G1 DETECTED, G2 MAGNETIC_DETECTED and G3 UNEQUAL.
SUPPORTED if all three hold; INCOMPLETE if any of the three is INCOMPLETE;
NOT_SUPPORTED otherwise. The lower levels and the other outcomes are reported in
full as the detection curve, and none of them can replace the primary claim.

House rule. G1 and G2 are existence gates and G3 is an ordering gate, decided on
the house-default 95% two-sided interval with the method above fixed now; its
power in Sect. 8 is computed for that same construction. None is an accuracy
gate; no sentence of the form "recovers `X` to within `t`" may be quoted from
this run, and every error reported below is quoted with its interval or
radius. Sect. 8 says why no accuracy gate is included and what one would
cost.

Selection, empty-clock, guard, rank and fit failures stay in the result. No
replacement, redraw, reseeding, sample extension, level change or threshold
change after outcomes. An infrastructure interruption (a crashed process, a
driver or allocation error, a power loss; never an outcome) means incomplete. A
rerun after an engineering repair is allowed only for such an interruption, must
be disclosed with the original stop, and counts only if it reproduces every row
the stopped run completed: order archives bit for bit and fits field for field
after the axis-sign gauge is aligned. Otherwise the stopped run's INCOMPLETE
stands.

`Q` and `Q_B` test departure from flat and from `B = 0`; they do not test
agreement with the truth, and no detection is called a recovery of this field.
Agreement is reported by one descriptive label fixed now, CONSISTENT at a level
if the nominal 95% ball of `confidence_ball` contains the truth in both orbit
cases, INCONSISTENT if not, UNAVAILABLE if an orbit fit or its truth evaluation
is missing, together with the quantities of Sect. 7. At
`2^26` those balls have radii near 18% of `s`, so CONSISTENT is not an accuracy
statement; INCONSISTENT is reported as a disagreement with the truth.

## 7. Descriptive quantities, not gates

Reported for every completed fit, with no pass/fail meaning:

- the full-tidal error modulo the residual reflection, with the nominal 95%
  radius of `confidence_ball`, never the point error alone. The truth is the
  apex tensor of Sect. 2 expressed in the frame that the inverse's own rule
  builds from the true normal coordinates of its first two landmarks;
- frame-free: the eigenvalues of the estimated `E` and `B`, their norms, `I_hat`
  with its interval, and `E.B` with its standard error (`E.B` is a pseudoscalar:
  its sign flips with the unresolved reflection; the truth is 0);
- for the orbit and rest cases, the `S^-1`-weighted reading of the estimate
  along the true pattern, for the full tensor, its electric half and (orbit) its
  magnetic half, each with its nominal sigma, using the reflection candidate
  nearer to the estimate; the error of each half; and the statistic of the
  estimate about the truth, the smaller over the two reflection candidates;
- the nominal signal-to-noise ratio along each true pattern, taking the less
  favourable reflection, and `N_1 = N/SNR^2`;
- `N * lambda_max`, margins resolved by the direct solver, run time, peak
  memory and archive sizes;
- the boost check: `I_hat` of the rest cases and of the orbit cases estimate the
  same invariant `6/576`.

## 8. Power and predictions, before data

All numbers assume the nominal Gaussian model. The inverse fixes its spatial
frame by two landmarks whose directions are random, so the truth of Sect. 2 is
applied to existing covariances in random orientations. Three sources:

- (a) the 63 completed fifth-degree covariances of the published flat
  calibration at `N = 2^24`, 400 orientations each, half of them reflected;
- (b) the eight covariances of the static test at `2^24` and `2^26` (four
  cases);
- (c) the engineering replay at `2^24`, rescaled from `T = 4M` by `(3/4)^4` for
  noncentralities and `(3/4)^2` for the standard error of `I`.

| Rate, zero model bias | (a) 5% - median - 95% | (b) | (c) |
|---|---|---|---:|
| noncentrality of `Q` per `10^6` events | 18.5 - 24.8 - 31.0 | 17.7 - 28.3 | 24.7 |
| noncentrality of `Q_B` per `10^6` events | 3.6 - 5.0 - 6.6 | 3.9 - 6.0 | 5.6 |
| `se(I_hat) x sqrt(N)` | 14.8 - 16.7 - 19.6 | 15.4 - 19.5 | 16.3 |

The flat calibration also gives the null behaviour of `Q_B`, which no earlier
stage gated: over its replicates the mean of `Q_B / 5` was 1.053 at `2^20` and
0.976 at `2^24` (exact 95% intervals `[0.91, 1.24]` and `[0.84, 1.15]`), with no
value above `q*_5` among 127. The static test's eight curved fits had `Q_B`
between 3.3 and 9.8 (mean 5.5): at that strength the electric field leaked no
resolvable magnetic part.

Probabilities from (a), two independent orbit cases, zero model bias, before
the control factor below:

| N | `E[Q]`, orbit | `E[Q_B]`, orbit | `I / se` | G1 DETECTED | G2 MAGNETIC_DETECTED | G3 UNEQUAL |
|---:|---|---|---:|---|---|---:|
| 1,048,576 | 36 (29-43) | 10 (9-12) | 0.64 | 0.47 (MIXED 0.43) | 0.003 | 0.01 |
| 4,194,304 | 114 (88-140) | 26 (20-33) | 1.28 | 1.00 | 0.46 (MIXED 0.44) | 0.06 |
| 16,777,216 | 426 (321-530) | 89 (65-117) | 2.56 | 1.00 | 1.00 | 0.52 |
| 67,108,864 | 1675 (1253-2088) | 340 (244-451) | 5.12 | 1.00 | 1.00 | 0.996 |

Parentheses are 5%-95% over selections and orientations. The control factor is
the probability that no control exceeds: 0.998 for G1 (two flat cases at the
nominal `1e-3`) and 0.996 for G2 (four controls), 0.994 if the magnetic block's
variance scale is at the calibration's 1.05. At `2^22` the probability that one
orbit case has `Q_B > q*_5` ranges from 0.43 to 0.85 over selections (5%-95%).

G3 is the gate that model bias can defeat, because `I` is a difference of
squared norms. At `2^26` its probability is 0.996 without bias (0.98 for the 5%
least favourable selections); 0.96 (0.87) if the electric part is read 4.5% low,
as the engineering replay suggests for `T = 3M`; 0.81 (0.64) if it is read 8% low,
as at `T = 4M`; 0.80 with opposing 5% biases on both norms (`E` low, `B` high)
and 0.09 with 10%; about 1 for biases in the other direction. At `2^24` it is
0.52, 0.28 and 0.14 for the first three cases. A NOT_RESOLVED or MIXED at `2^26`
could come from a bias of that size, from an unfavourable selection, or from an
inverse that does not resolve the ordering; this run cannot tell these apart,
because it does not measure the bias.

The repository's 90% power convention is met for G1 from `2^22`, for G2 from
`2^24`, and for G3 at `2^26` without bias and with the engineering-suggested
bias; it is not met for G3 if the electric bias at `T = 3M` is as large as at
`T = 4M`. The primary claim (Sect. 6) inherits these: with the control factors
its probability is 0.99 without bias, 0.95 with the engineering-suggested
electric bias (0.86 for the 5% least favourable selections) and 0.80 (0.63)
with the `T = 4M` bias.

**Why there is no accuracy gate.** A gate that the `S^-1`-weighted magnetic
reading lies in `[1 - t, 1 + t]` on its 95% interval in both orbit cases would
be the natural one. At `2^26` that reading has a nominal sigma of 0.047-0.065
(5%-95% over selections, median 0.055), so 90% joint power needs `t >= 0.21` for
the median selection without bias, `t >= 0.25` for the 5% least favourable one,
and `t >= 0.28` with a 5% bias. At such a tolerance the gate could not separate
`|B|` from `|E|`, which differ by 24% of `|E|`; that ordering is what G3 tests
directly, with more power. A tighter tolerance would need several times `2^26`
events per case, more than one fit can hold in memory, and the continuum bias an
accuracy gate must be sized against cannot be computed for this geometry. So no
accuracy gate is set, and none of the descriptive quantities of Sect. 7 may be
quoted as an accuracy result.

**Predictions, so that the run can contradict them.**

- G1: MIXED or DETECTED at `2^20` (about even), DETECTED from `2^22`.
- G2: NOT_DETECTED at `2^20`; MIXED or MAGNETIC_DETECTED at `2^22` (about
  even); MAGNETIC_DETECTED at `2^24` and `2^26`.
- G3: NOT_RESOLVED at `2^20` and `2^22`; any of NOT_RESOLVED, MIXED or UNEQUAL
  at `2^24`; UNEQUAL at `2^26`.
- Orbit, top level: nominal 95% radius 16-20% of `s`; point error of order
  5-15% of `s`; the estimated `E` eigenvalues recognizably `(-3, 1, 2)/24` and
  the `B` eigenvalues `(-2, 0, 2)/24`; magnetic reading within two of its sigma of
  1; electric reading between 0.90 and 1.00; `I_hat` below `6/576` if the
  electric bias seen in the engineering replay is real.
- Rest (not gated): `Q` about 19, 45, 150, 580 at the four levels from the
  calibration covariances (the static test, same geometry, realized 32 and 19,
  42 and 42, 140 and 158, 462 and 534); `Q_B` distributed as `chi2_5` at every
  level; `I_hat` estimating the same `6/576` with `se` about `7.2/sqrt(N)`
  (0.0009 at `2^26`).
- Flat: `Q` and `Q_B` distributed as `chi2_10` and `chi2_5`.
- Run time about two to three hours; peak memory about 50 GiB at the top level.

## 9. Execution and provenance

Commit this protocol, the generator, the runner and the checking script; run
once from the clean worktree. The manifest records the commit and the
LF-normalized hashes of the bound sources; the runner refuses a dirty worktree
and an existing destination, and after the last fit refuses to write the summary
if a bound source changed meanwhile.
Fits run in fresh single-threaded worker processes (`OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`); the relation queries and the exact integer second moments
use the GPU as before. **No other job may use the GPU while the stage runs**:
on 2026-11-20 a concurrent analysis on the same device made one calibration
replicate fail to allocate page-locked memory. Archives stay under `outputs/`;
only summaries and the results note are committed.

The checking script `check_intrinsic_weyl_schwarzschild_orbit.py` is committed
with this protocol. It was written before the run and exercised on the
engineering replay and on two six-case plumbing runs of this runner at
`N = 2^16, 2^17` with disjoint seeds (2026112300 and 2026112400), its generator
selfcheck stubbed, none of them test data. In the first, too few selection
events were used and every selection failed; the run then stopped at its end
because the protocol was edited while it ran (the source guard of the runner),
the summary assembled from its saved rows with the runner's own functions
reported every gate INCOMPLETE, and the script reproduced that and passed as
designed. The second completed all twelve rows; it
showed that the Clairaut solver can differ from the generator by up to `4e-7`
in coordinate time at a few quadrature nodes, exactly within the error bound
that solver reports (a refined generator solution agreed with the generator to
`5e-16`), so the script judges solver disagreements against that bound, not
against a fixed `1e-7`. The script restates the three gates with the same
precedence, the primary claim and the consistency label, and it tolerates the
failures this protocol keeps: a case without archives passes its checks only
if the run recorded that failure. Its refit comparison aligns the axis-sign
gauge before comparing fields, because the calibration audit showed that raw
fields can differ by such signs alone; that rule is fixed here, before the
outcome, as are the two rows it refits: case 0 at `N = 2^22` and case 1 at
`N = 2^24`. Corrections made to the script after the run will be visible in the
history.

## 10. Scope

Whatever the verdicts, this run establishes at most, for one static type-D
field at one event and one interval size, with two cases per observer and
nominal Gaussian uncertainty: detection of the field by the frozen inverse,
detection of the magnetic part that the moving observer sees together with a
non-detection, at a stated sensitivity, for the observer at rest, and the
ordering `|E| > |B|`, which concerns the moving observer's split only together
with the magnetic detection. It does not establish
accuracy at this strength, the size of the model bias, behaviour at other radii,
speeds or interval sizes, the gravitomagnetic field of a rotating source,
general vacuum geometry, or sample efficiency against any other method.
