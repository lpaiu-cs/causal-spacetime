# Ibrahim–Bombelli benchmark for the unchanged order-only Weyl inverse

2026-11-11. Freeze this protocol and its runner before any benchmark selection
or measurement event. Engineering checks on the audited constant-profile
generator are not benchmark data and are listed separately in section 3.

## 1. Question

[Ibrahim and Bombelli, arXiv:2611.09056](https://arxiv.org/abs/2611.09056)
estimate that detecting one Weyl component of a plane gravitational wave by
discrete geodesic deviation needs `N_total ~ 1e8–1e9` sprinkled elements, with
an ideal element-to-element distance estimator and endpoint separations known
without error, for amplitude `h0 = 1e-2` on a region of height `lambda/2`. No
recovery from a sprinkled causal set is performed there.

On 2026-11-11 this repository recorded only an **extrapolation** to that
condition from the 2026-09-23 weak-curvature campaign
([results](weak_curvature_accuracy_results.md)): about `3e6–7e6` events in the
full diamond for a signal-to-noise ratio of one on a single component. That
figure used the worst-direction variance `N*lambda_max ~ 6700` measured at tidal
norm `.5` on constant and slowly rotating profiles. The estimator has never
been run on this wave.

This protocol replaces the extrapolation by a measurement: **does the frozen
order-only inverse detect the benchmark wave, and at what event count?**

## 2. Geometry and unit map

Use the chart, tips and conventions of `rotating_plane_wave`:

    ds^2 = 2 du dv + a(u) (x^2 - y^2) du^2 + dx^2 + dy^2,
    tips P = (0, 1/2, 0, 0), Q = (1, -1/2, 0, 0), tip proper time T = sqrt(2).

Every profile `a(u)` is an exact vacuum solution. Along the central geodesic
`du/dtau = 1/sqrt(2)`, so the transverse deviation obeys
`d^2 xi/dtau^2 = (1/2) a(u) diag(1, -1) xi`, and the dimensionless tidal
eigenvalue in outer-tip units is `T^2 |E| = |a|`. The full tidal norm used by
the estimator is `2|a|` at the midpoint.

The benchmark wave has `d^2 xi^x/dt^2 = -(1/2) omega^2 h0 cos(omega (t - z)) xi^x`
(their Eq. 10). Their region has height `H = lambda/2`, so `omega T = pi`, and
equal separations at its top and bottom place a crest at mid-height. Hence

    a(u) = A0 sin(pi u),  A0 = pi^2 h0 / 2 = 0.0493480220...  (h0 = 0.01).

Three stated differences from their setup, none tuned after data:

- **Exact vacuum instead of the linearized metric.** Their metric satisfies
  the vacuum equations to first order in `h0`; its residual Ricci tensor is of
  relative order `h0`. The inverse assumes vacuum, so the exact plane wave with
  the same tidal profile on the axis is used. Their outlook names exact pp-waves
  as the intended next step.
- **Polarization orientation.** `a(u) diag(1, -1)` is their wave rotated by a
  quarter turn about the propagation axis. The inverse imposes no orientation.
- **Region.** They count elements in a thin slab around two geodesic tubes; this
  inverse uses the whole causal diamond between the tips. Event counts are
  therefore also reported as a density per `H^4`.

The truth at the midpoint is the constant-wave tensor of amplitude `A0`; the
half-cycle variation across the window is **not** in the inverse's fifth-degree
model. Its effect is measured by the continuum refit, not assumed small.

## 3. What is new and what is unchanged

New, generator side only: `experiments/positive_control/profile_plane_wave.py`.
It solves the two Jacobi equations by a power series in 60-digit decimal
arithmetic and stores Chebyshev series. Its deterministic selfcheck must pass
before data: the constant profile reproduces `RotatingWave(.25, 0)` in samples
(same random stream), order, interval volumes and normal map; the half-cycle
profile agrees with independent ODE solutions for the series, the action, the
causal relation and the geodesic round trip; CPU and GPU relation bits agree.
Margins inside the floating guard are decided in decimal arithmetic and
counted, never guessed.

Unchanged: `selection(order, design="nested")`, the interval clocks, the exact
shared-membership covariance, `fit_order_archive(..., model_order=5,
square_correction=True)`, four starts and all optimizer settings, the
`confidence_ball` construction, and the blind boundary (the pure fit is saved
before any generator or truth module is imported; only the assigned order
archive may be opened). No coordinates, profile, amplitude or paired flat sample
enter the inverse.

Engineering replay, not benchmark data: one constant-profile case `A = .25`
with seed `2026111190` at `N = 2^20, 2^22, 2^24`, to check plumbing, run time
and peak memory against the known behaviour of that case. A first replay ran
from the uncommitted tree before this freeze (127 landmarks, 3944 clocks):

| N | relative error | nominal 95% radius / norm | contains truth | continuum bias | N*lambda_max | peak memory |
|---:|---:|---:|---|---:|---:|---:|
| 1,048,576 | 32.7% | 46.2% | yes | 0.005% | 2648 | 5.4 GiB |
| 4,194,304 | 15.9% | 24.2% | yes | 0.005% | 2673 | 10.0 GiB |
| 16,777,216 | 6.7% | 13.1% | yes | 0.005% | 2668 | 16.6 GiB |

This matches the 2026-09-23 behaviour of that geometry (weak cases at
`N = 4,194,304`: 8.4–16.2% error; constant-profile continuum bias 0.003–0.004%).
The committed runner must reproduce this replay field for field, from a clean
worktree, before any benchmark event is drawn; a mismatch stops the run.

## 4. Assigned cases and independent data

All four cases are co-required at every level.

| Case | Geometry | Seed |
|---:|---|---:|
| 0 | benchmark wave, `A0 = pi^2 h0/2` | 2026111100 |
| 1 | benchmark wave, independent selection and events | 2026111101 |
| 2 | flat control | 2026111102 |
| 3 | flat control, independent | 2026111103 |

Each case draws 4096 physical-volume selection events from
`default_rng([seed, 0])`, appends the two supplied tips, permutes interior
labels and applies the nested selection. Measurement events come from
`default_rng([seed, 1])` in chunks of 4096. The levels are nested prefixes of
one stream per case:

    N = 1,048,576;  4,194,304;  16,777,216;  67,108,864.

At most `4 * 67,108,864 = 268,435,456` measurement events and 16,384 selection
events are drawn. Results at different levels of one case are not independent.
The levels were fixed from the engineering replay alone; its memory slope
predicts about 43 GiB at the top level on a machine with 125 GiB.

## 5. Frozen existence gate

This is an existence gate, not an accuracy claim; the house rule on accuracy
intervals is therefore not the governing rule here, and no accuracy sentence
may be quoted from this run without its interval.

For a completed fit let `c` be the ten fitted Weyl coefficients and `S` their
nominal covariance. The statistic is

    Q = c^T S^-1 c,

compared with `q* = chi2_10(1 - 1e-3) = 29.5883`. `Q` is invariant under the
axis-sign conventions. A fit is valid only if `confidence_ball` accepts it
(optimizer success, interior solution, full coordinate, nuisance and projected
curvature rank).

For each level:

- **DETECTED** if `Q > q*` in both wave cases and `Q <= q*` in both controls;
- **NOT_DETECTED** if `Q <= q*` in all four;
- **MIXED** if exactly one wave case exceeds `q*` and no control does;
- **CONTROL_FAILED** if either control exceeds `q*`; the level then supports no
  detection statement, whatever the wave cases show;
- **INCOMPLETE** if any assigned generation, selection or fit failed.

Selection, empty-clock, guard, rank and fit failures stay in the result. No
replacement, redraw, reseeding, sample extension, level change or threshold
change after outcomes. An infrastructure interruption means incomplete; a
rerun after an engineering repair must be disclosed with the original stop.

## 6. Descriptive quantities, not gates

Reported for every completed fit, with no pass/fail meaning:

- the full-tidal error modulo the residual reflection, together with the
  nominal 95% radius, never the point error alone;
- the continuum refit on exact interval volumes with the same clocks and
  covariance, which measures the model bias from the unmodelled half cycle;
- `N * lambda_max`;
- for the wave cases, the nominal signal-to-noise ratio along three fixed
  tidal patterns of the truth (electric part, magnetic part, both), taking the
  less favourable of the two reflection candidates, and the implied
  `N_1 = N / SNR^2`. The electric figure is the analogue of the Ibrahim–Bombelli
  condition `sigma < |signal|` for one targeted component; unlike theirs it is a
  nominal estimator standard deviation, with no distance oracle;
- margins decided in decimal arithmetic, run time, peak memory, archive sizes.

Recorded before data so the run can contradict it:

- The extrapolation of section 1 was `N_1 ~ 3e6–7e6`. It used the worst
  direction and a single-component amplitude, so the electric-pattern figure
  may be smaller.
- The engineering replay gives `N_1 = 10,800` for the electric pattern at
  `A = .25`. Scaling by `(0.25/A0)^2 = 25.7` predicts `N_1 ~ 2.8e5` here, larger
  by `1/b^2` if the inverse reads a fraction `b` of the crest curvature.
- Because the landmark bands sample about `0.8–0.9` of the crest curvature, the
  continuum refit is expected to read low by roughly 10–20%; this expectation is
  not a correction and is not applied.
- With the engineering selection's variance and `b = 0.8–0.9`, the expected
  statistic is `Q ~ 15, 33, 100, 370` at the four levels against `q* = 29.6`:
  no detection at the lowest level, a marginal second level, detection at the
  upper two. The 2026-09-23 selections varied by up to a factor of three in
  `N*lambda_max`, which could move the detection onset up by one level.

## 7. Execution and provenance

Commit this protocol, the generator and the runner; run once from the clean
worktree; the manifest records the commit and LF-normalized source hashes, and
the runner refuses a dirty worktree and an existing destination. Fits run in
fresh single-threaded worker processes (`OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`); the relation queries and the exact integer second moments
use the GPU as in the 2026-09-23 campaign. Archives stay under `outputs/`;
only summaries and this protocol's results note are committed.

## 8. Scope

Whatever the verdict, this run establishes at most detection of one prescribed
wave on one window by one frozen inverse, with nominal Gaussian uncertainty.
It does not establish accuracy at this amplitude, behaviour for other
wavelength-to-window ratios, phases or polarizations, general vacuum geometry,
sample efficiency against any other method, or any statement about the
geodesic-deviation estimator itself, which is not implemented here. The two
methods use different inputs and different regions, so their event counts are
compared only with those differences stated.
