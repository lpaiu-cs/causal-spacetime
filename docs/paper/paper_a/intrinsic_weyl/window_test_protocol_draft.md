# Finite-N test of the order-only five-window correction (DRAFT)

2026-11-21. **Draft for decisions, not frozen.** Nothing here is binding until
the freeze commit named in Sect. 10, which must precede every test selection
and measurement event. Section 11 lists the decisions still open; the numbers
in Sects. 6-7 use the proposed choices.

## 1. Question

On the Ibrahim–Bombelli half-cycle wave the fifth-degree order-only inverse
reads the crest tensor about 12% low in the continuum, because the jet omits
the second covariant derivatives of the curvature
([sixth-degree record](second_derivative_results.md)). Two deterministic
continuum studies found a cheaper remedy than the sixth-degree jet:

- [Multi-window study](window_study_results.md): windows displaced along the
  observer's worldline (proper-time offsets 0, ±1/2, ±1 of the tip
  separation) measure the time-time second derivative; a five-point
  second-order correction of the centre window with a minimum-norm operator,
  which assumes no wave type, leaves 0.54% and 0.60% of the truth norm at 5.7
  to 6.1 times the events of one window.
- [Frame link](frame_link_results.md): windows that share the landmarks of
  their inner neighbour are aligned by a Procrustes rotation of the shared
  landmarks' inferred coordinates to 0.03°, and the corrected results equal
  those of oracle alignment, so the correction is order-only end to end.

Both were continuum calculations with nominal Gaussian errors. This test asks:
**given only the causal order of five windows along one observer's worldline
in the half-cycle wave, does the order-only five-window correction recover
the crest tensor's amplitude to within ±t, on the 95% two-sided interval of
its mean over independent replicates, while the single-window inverse on the
same data reads it low by more than t?**

Scope, fixed now: the amplitude is the component of the estimate along the
true tensor (Sect. 4), one number. Full-tensor accuracy is reported
descriptively and is not claimed. One wave, one interval size, one amplitude.

## 2. Geometry, windows and truth

- **Wave.** The generator-side linearly polarized plane wave of the benchmark,
  `a(u) = A0 sin(pi u)`, `A0 = pi^2 h0 / 2`, `h0 = 0.01`
  (`profile_plane_wave.py`, unchanged).
- **Windows.** Five windows along the geodesic `x = y = 0` of the Brinkmann
  chart, with tips at proper times `delta -/+ 1/2` for
  `delta in {-1, -1/2, 0, 1/2, 1}`. The window at `delta` is the standard
  window (tips `(0, 1/2, 0, 0)`, `(1, -1/2, 0, 0)`) of the shifted profile
  `a(u + delta)`; an event's chart coordinates in two windows differ by
  `(delta1 - delta2, delta2 - delta1, 0, 0)`. Both facts were checked in the
  continuum studies (window isometry of actions to `3e-16`, squared intervals
  of shared events to `3e-12`).
- **Truth.** The crest tensor at the centre of window 0 (amplitude `A0`),
  expressed in window 0's inferred frame by the O(3) alignment of window 0's
  inferred landmark coordinates to their true normal coordinates, including
  the parity of that alignment (evaluation only, as `aligned_tensor_error`).
- **Continuum expectations** (benchmark landmark layouts, cases 0/1, run
  covariances): single-window reading `0.8829 / 0.8880`; five-point second
  order with the minimum-norm operator `0.9959 / 0.9957`, with the plane-wave
  operator `0.9940 / 0.9944`; three-point first order `0.9697 / 0.9713`.

## 3. What is new and what is unchanged

New, in a runner to be written (`intrinsic_weyl_window_test.py`):

- **Landmarks of a window off the centre.** The unchanged nested selection on
  the window's own 4096 selection events, plus the shared landmarks: the
  landmarks that the inner neighbour's own selection chose (for `±1/2` the
  centre window's, for `±1` those of `±1/2`), kept if they lie in the window's
  interval and if at least two of the window's selection events lie between
  each tip and the landmark. This rule is order-only; the continuum study's
  rule used flat proper times and is replaced.
- **Clocks.** Pairs with at least two of the window's selection events, the
  unchanged rule, applied to the window's landmarks including the shared
  ones; every landmark must keep both tip clocks.
- **Archive.** The window's selection order is extended by the candidate
  shared events, so that the landmarks and clocks can be re-derived from the
  archive alone, as `read_order_clocks` does now.
- **Correction (order-only, after the five pure fits).**
  1. Links: the orthogonal Procrustes map of the shared landmarks' inferred
     spatial coordinates, **over O(3)**: the order does not fix the
     handedness of an embedding, so two windows may have converged to mirror
     embeddings; an improper link maps the Weyl coefficients with the parity
     the tensor transformation gives. The `±1` windows are chained through
     `±1/2`.
  2. Alias map: Weyl rows of `J5^+ Jd` at window 0's inferred coordinates and
     zero curvature, whitened with window 0's clock covariance.
  3. Minimum-norm operator from the alias map (time-time slice, regularity
     norm of the multi-window study).
  4. Five-point second-order stencil, `h = 1/2`:
     `c0 - B D2[c] + B^2 D4[c]`.

Unchanged: generation of selection and measurement events and of relation
bits, the exact integer second moments, the fifth-degree jet with
curvature-square terms, four starts and all optimizer settings. Blind
boundary: the five pure fits and the corrected estimate are saved before any
generator or truth module is imported in the analysing process.

## 4. Statistics

For replicate `r` let `c_r` be the corrected estimate and `u_r` the window-0
estimate alone, both in window 0's inferred frame, and `c` the truth there.
With the full-tidal inner product `<x, y>`:

- corrected reading `a_r = <c, c_r> / <c, c>`;
- single-window reading `s_r = <c, u_r> / <c, c>`.

Over the completed replicates, the means `a`, `s` and their 95% two-sided
Student-t intervals with `R - 1` degrees of freedom from the sample standard
deviations (house default coverage; the method is fixed now). The sample
standard deviation is used instead of the nominal covariance, so the
intervals do not rely on the calibration of the augmented windows.

## 5. Gates

- **G1, accuracy (primary).** ACCURATE if the interval of `a` lies inside
  `[1 - t, 1 + t]`; INACCURATE if it lies entirely outside; otherwise
  UNRESOLVED.
- **G2, control.** BIASED_LOW if the interval of `s` lies entirely below
  `1 - t`; otherwise NOT_BIASED. If the single window is not biased beyond the
  tolerance, an ACCURATE G1 cannot be credited to the correction.
- **Primary claim.** SUPPORTED if G1 is ACCURATE and G2 is BIASED_LOW;
  INCOMPLETE if fewer than 90% of the replicates complete; otherwise
  NOT_SUPPORTED.

House rule: G1 is an accuracy gate, decided on the 95% two-sided interval and
sized for 90% power in Sect. 7 with that same construction. G2 is a sign gate.
A replicate fails if any of its five windows fails (selection, empty clock,
rank, fit, or a link with fewer than eight shared landmarks); failures stay in
the result and are never replaced. Interruption, rerun and precedence rules
as in the orbit protocol.

Descriptive, no verdict: the plane-wave and three-point readings; the mean
error vector and its full-tidal norm with a bootstrap interval; electric and
magnetic readings; per-window fits; shared counts, link parities and
Procrustes residuals; the angle between the order-only and the oracle link
(after the blind step).

## 6. Replicates, seeds, events

Proposed: `R = 64` replicates of five windows, `N = 2^25` measurement events
per window, seeds `2026112600 + r`; window `w` of replicate `r` draws its
selection from `default_rng([seed, w, 0])` and its measurement from
`default_rng([seed, w, 1])`. Total `64 x 5 x 2^25 = 1.07e10` measurement
events and `1,310,720` selection events. Windows are independent sprinkles;
in one causal set they would share events where they overlap, which changes
the variance but not the expected reading, and the interval uses the
replicates' own spread.

## 7. Power and cost (before data)

Inputs (`intrinsic_weyl_window_test_power.py`, commit `fd4e06e`, output
`window_test_power.json`): the one-window reading has standard deviation
0.113 at `2^24` (rms over 63 calibration selections x 100 orientations;
5-95%: 0.099-0.132); the corrected reading's variance is 1.29 times that
(minimum-norm, median over 400 orientations, larger of the two layouts), so
0.128 at `2^24` and 0.091 at `2^25`, assuming variance `~ 1/N`.

Replicates for 90% power of G1 at `t = 0.05`, minimum-norm operator:

| events per window | mean reading 0.996 | 0.986 | 0.976 |
|---|---|---|---|
| `2^24`, `kappa = 1` / 1.2 | 96 / 112 | 160 / 192 | over 256 |
| `2^25` | 48 / 56 | 80 / 96 | 160 / 160 |
| `2^26` | 28 / 32 | 40 / 48 | 80 / 80 |

`kappa` is a variance understatement factor; 0.986 and 0.976 allow a finite-N
bias of one and two percent below the continuum reading. At the proposed
design (`2^25`, `R = 64`) the power of G1 is 0.98 (`kappa = 1`) and 0.94
(`kappa = 1.2`) at the continuum reading, 0.87 and 0.81 if the reading is one
percent lower, and 0.61 and 0.53 if two percent lower. A finite-N bias of two
percent would therefore most likely give UNRESOLVED, not a false ACCURATE.
The three-point correction cannot pass at `t = 0.05` (continuum reading
0.970). G2 has power near 1 (single-window reading 0.885).

Cost model, scaled from one calibration replicate at `2^24` and to be checked
by the pilot: about 16 GPU hours and 15 CPU hours of fits at the proposed
design (fits run in parallel with generation), peak memory about 37 GB per
window while its events are generated.
`t = 0.06` would need about 40 replicates at `2^25` (10 GPU hours); `t = 0.04`
about 96 at `2^25` (24 GPU hours).

## 8. Engineering pilot (before the freeze, not test data)

Two replicates with disjoint seeds (`2026112700`, `2026112701`) at `2^22` and
`2^24` events per window, from the uncommitted runner:

- plumbing, archives, the blind boundary and the checking script;
- shared counts, link parities and residuals, fit success at finite `N`;
- run time and peak memory per window, to replace the cost model;
- the nominal standard deviation of each replicate's corrected reading
  against the design value.

The committed runner must reproduce the pilot field for field before any
test event.

## 9. Predictions (before data)

- `a` between 0.98 and 1.01; `s` between 0.86 and 0.90.
- Per-replicate standard deviation at `2^25`: corrected 0.08-0.10, single
  0.07-0.09.
- G1 ACCURATE, G2 BIASED_LOW, primary claim SUPPORTED.
- Shared landmarks 30-45 per link; links improper in about half of the
  replicates; Procrustes relative residual about `1e-3`.

## 10. Execution and provenance

Freeze this protocol, the runner and the checking script in one commit after
the pilot; run once from the clean worktree with the GPU exclusive. Archives
are kept for two designated replicates for the audit and hashed for the
others. Summaries and a results note are committed.

## 11. Open decisions (draft only)

1. **Tolerance `t`.** Proposed 0.05. 0.06 is cheaper (about 40 replicates at
   `2^25`); 0.04 costs about 1.6e10 events.
2. **Events per window and replicates.** Proposed `2^25` and 64. `2^24` needs
   about twice the replicates and fits for the same GPU time and is more
   exposed to finite-N bias; `2^26` halves the fits but needs about 74 GB per
   window during generation.
3. **Primary operator.** Proposed minimum-norm (no wave type assumed).
   Plane-wave with the estimated principal null direction is slightly
   cheaper (variance ratio 1.21 instead of 1.29).
4. **Amplitude.** Proposed the benchmark `h0 = 0.01`. Twice that would need
   about a quarter of the events, but would no longer be the Ibrahim–Bombelli
   condition.
5. **Sprinkling.** Proposed independent windows. One sprinkle over the union
   of the five windows is closer to one causal set but needs a new sampler and
   gives windows variable event counts.

## 12. Scope

Whatever the verdict, this run would establish at most, for one plane wave at
one amplitude and interval size: that the order-only five-window correction
removes the fifth-degree bias of the amplitude along the true tensor at
finite `N` to within the tolerance, with the single window as control. It
would not establish full-tensor accuracy, the behaviour for other waves or
fields, accuracy of one replicate (each replicate's own error is about ten
times larger), or the cost of a single-sprinkle implementation.
