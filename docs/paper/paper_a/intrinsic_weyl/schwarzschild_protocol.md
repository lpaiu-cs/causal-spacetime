# Schwarzschild detection test for the unchanged order-only Weyl inverse

2026-11-12. Freeze this protocol, the generator and the runner before any test
selection or measurement event. The engineering replay of section 3 uses a
different interval and seed and is not test data.

## 1. Question

Every quantitative test of the fifth-degree order-only inverse so far used an
exact vacuum plane wave: Petrov type N, electric and magnetic parts of equal
norm, curvature depending on one null coordinate. The
[benchmark of 2026-11-11](bombelli_benchmark_results.md) is of that class too,
so the record has to say "plane waves only".

This protocol asks one question outside that class: **does the same frozen
inverse detect the tidal field of a Schwarzschild exterior — static, Petrov
type D, purely electric for the observer through the tips, and varying by a
factor of more than three across the interval — at a comparable dimensionless
strength and comparable event counts?** Nothing on the inverse side changes.

It is one geometry, one position and one interval. It is not an accuracy test,
and the model bias cannot be measured separately here (section 2).

## 2. Geometry and unit map

Isotropic coordinates, `M = 1`:

    ds^2 = -A^2 dt^2 + B^4 dX.dX,  A = (1 - M/2rho)/(1 + M/2rho),  B = 1 + M/2rho.

The tips `P`, `Q` lie on the radial timelike geodesic that is momentarily at
rest at areal radius `r_a = 6M` (the apex, `rho_a = 4.949490`), at proper times
`-T/2` and `+T/2` from the apex with `T = 3M`:

    P, Q = (t, rho) = (-/+ 1.838721, 4.917860),  on the x axis.

As in the plane-wave tests the tip proper time is the unit. The apex is the
midpoint of the tip geodesic, and the static orthonormal frame there is the rest
frame of that geodesic. In it

    T^2 E = (T^2 M / r_a^3) diag(-2, 1, 1) = diag(-1/12, 1/24, 1/24),  B = 0,

with the radial axis first. The full tidal norm is `sqrt(6)/24 = 0.10206`; the
benchmark wave had `0.0987`. `r_a` and `T` were chosen for that match and to
keep the interval well outside the photon sphere; no Schwarzschild fit existed
when they were chosen.

The interval reaches areal radii `4.816M` to `7.248M` on the slice `t = 0`, so
`M/r^3` ranges from `1.93` to `0.57` times its apex value. The inverse models a
constant Weyl tensor, its 24 vacuum first derivatives and the curvature-square
term. Second and higher derivatives are **not** in the model. For the plane
waves their effect was measured by a continuum refit on exact interval volumes.
No exact interval volumes are implemented for this geometry, so here the model
bias is not measured separately: it is part of the reported error.

The spacetime is static, so `q` is causally after `p` exactly when
`t_q - t_p` is at least the optical distance of their spatial points, the
Fermat minimum of the integral of `n = B^2/A`. Events are drawn uniformly in
proper volume, `sqrt(-g) = A B^6`.

The flat controls use the same generator class with `M = 0` (tips at
`t = -/+ 1.5`), so that selection, sampling, relation and packing code is shared
with the curved cases. They cannot reveal an error that vanishes with `M`; the
deterministic checks of section 3 are what guard the curved measure and order.

## 3. What is new and what is unchanged

New, generator side only: `experiments/positive_control/schwarzschild_diamond.py`.
Its deterministic selfcheck must pass before data:

- the derivatives of the index, by complex step;
- radial pairs against the exact tortoise-coordinate difference (`3e-15`);
- the optical-distance solver (Newton iteration on ten polynomial deviation
  modes from the chord, 32 Gauss nodes) against a 16-mode, 64-node solution
  (`2e-15`), with Bouguer's invariant `n rho sin(psi)` constant along the solved
  path to `4e-9`; the time error is second order in that path error;
- symmetry and the triangle inequality of the optical distance;
- an independent integration of the ray equation with shooting (`1e-12`, the
  accuracy of that integrator);
- the geodesic integrator: norm, energy and angular momentum are conserved on
  radial and on generic non-radial geodesics;
- the tidal tensor from the relative motion of particles released at rest
  (`1e-5`), and the tips at `(-/+ 1/2, 0, 0, 0)` in normal coordinates;
- the sampler: every event is inside the tips' interval by the direct solver,
  and four moments of `2e6` events agree with a Gauss quadrature of the
  proper-volume measure that uses the direct solver only (`|z| < 4.5`);
- the proper volume of a small interval (`T = M/2`) exceeds the flat value
  `pi T^4/24` by the curvature-square law `17 T^4 E_ab E^ab / 10080` of
  [square_correction.md](square_correction.md) to `1e-3`. This checks metric,
  solver, volume element and tip proper time together, and it tests the electric
  coefficient of that law apart from the magnetic one, which a plane wave
  cannot;
- relation bits from the surrogate equal those of the direct solver, and the
  order is antisymmetric and transitive;
- the `M = 0` class reproduces Minkowski distances, tips, normal coordinates
  and volume.

The runner adds one check of its own before data: for both classes the relation
bits evaluated on the GPU, which are the ones archived, equal those of the
direct solver on 262,144 pairs.

A development check outside the selfcheck: the older solver of
`s1_schwarzschild_cost.py`, which works in areal coordinates with the Clairaut
constant and shares no code with the generator, agreed within its own error
bound on 1500 of 1500 random pairs of the interval (largest bound `1.9e-9`,
median difference `7e-14`).

Relation decisions. For each landmark the squared optical distance to a point is
a smooth function of two variables; a degree-24 Chebyshev series of it is
evaluated on the GPU. Its guard is eight times its largest error on 512 test
points plus `1e-12`. A pair inside the guard is decided by the direct solver. If
that margin is within `1e-13` the generator stops and the case is a recorded
failure; nothing is guessed. About `2e-3` such stops are expected in the whole
run. Events closer than `1e-10` in coordinate time to the boundary of the
interval are not drawn, which removes a fraction of order `1e-10` of the
measure.

Unchanged: `selection(order, design="nested")`, the interval clocks, the exact
shared-membership covariance, `fit_order_archive(..., model_order=5,
square_correction=True)`, four starts and all optimizer settings, the
`confidence_ball` construction, and the blind boundary (the pure fit is saved
before any generator or truth module is imported; only the assigned order
archive may be opened). No coordinates, mass, radius or paired flat sample enter
the inverse.

Engineering replay, not test data: one case on a larger interval at the same
apex, `T = 4M` (tidal norm `0.1814`, areal radii `4.45M`–`7.67M`), seed
`2026111290`, at `N = 2^20, 2^22, 2^24`, to check plumbing, run time and peak
memory. A first replay ran from the uncommitted tree before this freeze
(124 landmarks, 3381 clocks, ten minutes in all):

| N | Q | relative error | nominal 95% radius / norm | contains truth | N*lambda_max | fit chi2, z | peak memory |
|---:|---:|---:|---:|---|---:|---:|---:|
| 1,048,576 | 9.2 | 82.6% | 121.1% | yes | 2598 | -0.01 | 4.2 GiB |
| 4,194,304 | 86.5 | 35.8% | 61.9% | yes | 2610 | 0.42 | 7.5 GiB |
| 16,777,216 | 393.7 | 19.1% | 32.0% | yes | 2619 | 1.53 | 12.6 GiB |

With this selection the true tensor has a noncentrality of `2.6e-5` per event
(27, 109 and 434 at the three levels). Along the true tensor the fit reads
`0.39 +/- 0.19`, `0.85 +/- 0.10` and `0.94 +/- 0.05` of it; the levels are nested
and therefore correlated. The fitted landmark coordinates agree with the true
normal coordinates to an rms of `9e-4`, `4e-4` and `2e-4` of the tip time, and
the fitted cloud is the mirror image of the true one, as the reflection
ambiguity allows.

The lowest level is far below its expectation of `10 + 27`. To separate chance
from a low-count effect, the first `2^24` engineering events were cut into 16
disjoint blocks of `2^20` and each was fitted with the same landmarks. `Q` ranged
from 9.2 to 66.0 with mean 37.0, against 37.2 expected; the nested prefix is the
lowest of the sixteen. Twelve of sixteen exceeded `q*`. The reading along the
truth averaged `0.93 +/- 0.06`. The statistic of the estimate about the truth
averaged 13.4 where 10 is nominal, so at `2^20` events the nominal covariance
may understate the scatter by about 15% in standard deviation; sixteen blocks
of one selection do not settle that.

The committed runner must reproduce this replay field for field, from a clean
worktree, before any test event is drawn; a mismatch stops the run.

## 4. Assigned cases and independent data

All four cases are co-required at every level.

| Case | Geometry | Seed |
|---:|---|---:|
| 0 | Schwarzschild, `M = 1`, `r_a = 6M`, `T = 3M` | 2026111200 |
| 1 | the same, independent selection and events | 2026111201 |
| 2 | flat control, `M = 0`, same `T` | 2026111202 |
| 3 | flat control, independent | 2026111203 |

Each case draws 4096 proper-volume selection events from
`default_rng([seed, 0])`, appends the two tips, permutes interior labels and
applies the nested selection. Measurement events come from
`default_rng([seed, 1])` in chunks of 4096. The levels are nested prefixes of
one stream per case, the same four as in the benchmark:

    N = 1,048,576;  4,194,304;  16,777,216;  67,108,864.

At most `4 * 67,108,864 = 268,435,456` measurement events and 16,384 selection
events are drawn. Results at different levels of one case are not independent.

## 5. Frozen existence gate

This is an existence gate, not an accuracy claim; the house rule on accuracy
intervals is therefore not the governing rule here, and no accuracy sentence
may be quoted from this run without its interval.

For a completed fit let `c` be the ten fitted Weyl coefficients and `S` their
nominal covariance. The statistic is

    Q = c^T S^-1 c,

compared with `q* = chi2_10(1 - 1e-3) = 29.5883`. A fit is valid only if
`confidence_ball` accepts it (optimizer success, interior solution, full
coordinate, nuisance and projected curvature rank).

For each level:

- **DETECTED** if `Q > q*` in both curved cases and `Q <= q*` in both controls;
- **NOT_DETECTED** if `Q <= q*` in all four;
- **MIXED** if exactly one curved case exceeds `q*` and no control does;
- **CONTROL_FAILED** if either control exceeds `q*`; the level then supports no
  detection statement, whatever the curved cases show;
- **INCOMPLETE** if any assigned generation, selection or fit failed.

Selection, empty-clock, guard, rank and fit failures stay in the result. No
replacement, redraw, reseeding, sample extension, level change or threshold
change after outcomes. An infrastructure interruption means incomplete; a
rerun after an engineering repair must be disclosed with the original stop.

`Q` tests departure from the flat model, not agreement with the Schwarzschild
tensor. A detection is therefore described as a recovery of this field only
together with the comparison of section 6; a detected signal that disagrees
with the truth is reported as that.

## 6. Descriptive quantities, not gates

Reported for every completed fit, with no pass/fail meaning:

- the full-tidal error modulo the residual reflection, together with the
  nominal 95% radius, never the point error alone. The truth is the apex tensor
  of section 2, expressed in the frame that the inverse's own rule builds from
  the true normal coordinates of its first two landmarks;
- frame-free: the eigenvalues of the estimated electric part (truth
  `-1/12, 1/24, 1/24`) and the norm of the estimated magnetic part (truth `0`);
- `N * lambda_max`;
- for the curved cases, the nominal signal-to-noise ratio along the tidal
  pattern of the truth, taking the less favourable reflection, and the implied
  `N_1 = N / SNR^2`;
- margins resolved by the direct solver, run time, peak memory, archive sizes.

Recorded before data so the run can contradict it:

- The covariances of the sixteen completed benchmark fits, applied to a type-D
  electric pattern of norm `0.10206` in a random orientation, give a
  noncentrality of `6.1–11.4` per `1e6` events (5%–95% over orientations and
  the four selections; medians `7.4–9.8`). If the inverse read the apex tensor
  without bias, `Q ~ 10 + noncentrality`: about `16–22`, `36–58`, `110–200`
  and `420–780` at the four levels against `q* = 29.6`. Expected verdicts: no
  detection at the lowest level; probable detection at the second, where each
  curved case should exceed `q*` with probability roughly `0.6–0.97`, so a
  mixed verdict there would not be surprising; detection at the upper two.
- Per unit norm this field should be easier to detect than the benchmark wave,
  which realized `4.2–4.6` per `1e6` events, because the magnetic directions
  carry most of the variance and this field has no magnetic part.
- At the top level the nominal 95% radius should be about 30% of the norm, as
  in the benchmark, with a point error of roughly 10–25%; this run therefore
  cannot support an accuracy sentence at this strength.
- The eigenvalue pattern of the estimated electric part should be recognizably
  `(-2, 1, 1)` at the top level, and the estimated magnetic norm should lie
  inside the nominal radius.
- The unmodelled variation of the curvature across the interval is expected to
  cost several per cent of the norm; the engineering replay, on a larger
  interval, bounds it loosely. No correction is applied.
- The engineering replay gives `7.9e-4` per event and unit squared norm, inside
  the range `5.9e-4`–`1.09e-3` of the benchmark-based figure above, and reads
  `0.94 +/- 0.05` of the apex tensor on an interval whose unmodelled
  second-order variation is `(4/3)^2 = 1.8` times larger than here. For the
  test interval this suggests a noncentrality near `8e-6` per event and a
  reading within a few per cent of the apex tensor: `Q ~ 19, 44, 150, 560`.
- In terms of verdicts, with noncentral chi-square probabilities over the
  range above: NOT_DETECTED at the first level with probability about
  `0.75–0.9` (otherwise MIXED); at the second level DETECTED with about
  `0.5–0.8` and MIXED with about `0.2–0.4`; DETECTED at the upper two. The
  engineering blocks show that a single case can fall far below its
  expectation at a noncentrality near 30.

## 7. Execution and provenance

Commit this protocol, the generator and the runner; run once from the clean
worktree; the manifest records the commit and LF-normalized source hashes, and
the runner refuses a dirty worktree and an existing destination. Fits run in
fresh single-threaded worker processes (`OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`); the relation queries and the exact integer second moments
use the GPU as in the earlier campaigns. Archives stay under `outputs/`; only
summaries and this protocol's results note are committed.

The checking script `check_intrinsic_weyl_schwarzschild.py` is committed with
this protocol. It was written before the run, and its functions were exercised
on the engineering replay only: recorded values reproduced at the three
levels, 1997 of 1997 sampled order bits confirmed by the Clairaut solver, and
the first `2^20` engineering events consistent with the areal volume element
(`|z| <= 2.4` on four moments; `0.7` on all `2^24`). Corrections made to it
after the run will be visible in the history.

## 8. Scope

Whatever the verdict, this run establishes at most detection of one static
type-D field at one radius on one interval by one frozen inverse, with nominal
Gaussian uncertainty. It does not establish accuracy at this strength, the size
of the model bias, behaviour at other radii or interval sizes, behaviour for a
field with both electric and magnetic parts that is not a plane wave (a
rotating source, for example), general vacuum geometry, or sample efficiency
against any other method.
