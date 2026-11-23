# Finite-N test of the five-window correction at spacing T/3, where the inner side windows and the links carry the reading

2026-11-23. **Frozen** in the commit that renames this file from
`window_link_test_protocol_draft.md` to `window_link_test_protocol.md`. That
commit freezes this protocol together with the runner
`experiments/positive_control/intrinsic_weyl_window_link_test.py` and the
checking script `experiments/positive_control/check_intrinsic_weyl_window_link_test.py`,
both last changed at commit `6b4e1a4`. It follows the engineering pilot of
Sect. 8 and precedes every test selection and measurement event.

How the draft (commits `a6f7ae0`, `7e6921d`, `d83287a`, `871668c`) became this
protocol:

- On 2026-11-21 the user took the proposed choices of Sect. 12.1: spacing
  `T/3`, `t = 0.05`, 96 replicates, F2 as proposed, G3 secondary, four fit
  workers, the seeds of Sect. 6.
- An independent design review and an independent code review, run before
  the freeze at the user's request, led to the revisions of Sect. 12.2.
- At the freeze on 2026-11-23 the user confirmed those revisions, including
  G3 in the primary claim and the secondary gates G4 and G5.

The engineering pilot ran on 2026-11-21 with the earlier runner, and the
checking script passed on it (Sect. 8). The revised runner changes no reading
of a completed replicate. The pilot replay at this commit is therefore the
frozen runner's pilot, and it must reproduce the earlier one (Sect. 10). All
design numbers come from the design aid
`experiments/positive_control/intrinsic_weyl_window_link_test_power.py`
(commit `6b4e1a4`), whose recorded run is
[`window_link_test_power.json`](window_link_test_power.json). It reproduces
every number of the aid's first recorded run (commit `a58965c`) and adds G4.
Neither run drew a measurement event; their only draws are the same selection
events of a dry run (Sect. 6).

## 1. Why a second test

The [frozen five-window test](window_test_protocol.md) (commit `1ca7af2`,
[record](window_test_results.md), commit `02c21e2`) is SUPPORTED: on the
Ibrahim–Bombelli half-cycle wave the order-only five-window correction read the
crest amplitude along the true tensor as 0.9976 [0.9772, 1.0180] over 64
replicates, while the single window read 0.8870 [0.8689, 0.9051]. That record
stands as frozen and is not re-scored here.

A later independent review found the test nearly blind to its side windows and
links. The post hoc split of the recorded reading
(`intrinsic_weyl_window_test_sensitivity.py`, commit `c87fdab`, run from a
clean worktree; record
[`window_test_sensitivity.json`](window_test_sensitivity.json) and the revised
[record](window_test_results.md), commit `d4121a2`, merged into this branch at
`0be9ce2`) shows why:

- With windows at proper-time offsets 0, ±T/2, ±T (spacing `h = T/2`, `T` the
  tip separation) on a wave of wavelength `2T`, the ±T/2 windows sit on nodes
  and read about 0, and the ±T windows carry a weight of only 0.0067 along the
  truth. The corrected reading is about `1.139 x` the centre window's own
  estimate `- 0.012`, where 1.139 is the component of `I + 10B + 96B^2` along
  the truth, computed from the centre window's own inferred geometry.
- Deleting the side windows gives 1.0086 [0.9883, 1.0290]; identity links give
  1.0098 [0.9895, 1.0302]; the 64 Haar-random side-link draws have means
  1.0062–1.0114. Every such variant would have passed. The accuracy of the order-only links rests only on descriptive
  comparisons (parity 256/256, angle errors at most 0.054°).

The cause is a coincidence of `lambda = 2T` and `h = T/2`, not of the phase.
For a profile `f(delta) = cos(k delta + phi)` the five-point stencil sees the
pair sums `f(h) + f(-h) = 2 cos(kh) f(0)` and `f(2h) + f(-2h) = 2 cos(2kh) f(0)`.
At `kh = pi/2` the inner pair sum vanishes **for every phase**, and the centre
coefficient `30/(12 h^2) = 10` happens to be close to `k^2 = pi^2 = 9.87`, so
the centre window's own operator does the work. A phase shift only splits the
inner pair into equal and opposite parts, which identity or random links
leave unchanged in the mean. A longer wavelength moves the nodes but shrinks
the single-window bias as about `0.5 (T/lambda)^2` (5.5% at `lambda = 3T`),
too little for a control at `t = 0.05`. **The spacing is the remedy.**

This test asks: **with the side windows at 0, ±T/3, ±2T/3 on the same wave,
where the side windows read ±0.44 of the crest and the centre term alone
overshoots to about 1.20, does the order-only five-window correction recover
the crest amplitude along the true tensor to within ±0.05 on the 95%
two-sided interval of its mean over independent replicates, while (control)
the single window on the same data reads low by more than 0.05,
(falsification) the centre term alone and the correction with identity,
handedness-flipped or Haar-random side links all read above 1.05, and (link
signs) each side window's linked reading lies on its predicted side of zero?**

Scope, fixed now: the amplitude along the true tensor, one number. One wave,
one amplitude, one interval size. Full-tensor accuracy is descriptive.
"Order-only" refers to the links, the estimates and the operator `B`. The
correction also takes the spacing `h` and the tip separation `T` as inputs.
Those are proper-time intervals that the observer reads off their own clock,
as in the frozen test.

## 2. Choice of spacing (design aid, Part 1)

Continuum fits on the benchmark landmark layouts (cases 0 / 1, run clock
covariances), oracle frames, readings relative to the centre window's truth.
"Ratio" is the variance of the corrected reading over that of one window
(median over 400 orientations).

| `h/T` | inner / outer window reads | corrected | centre term alone | flipped links | Haar links (mean) | ratio |
|---|---|---|---|---|---|---|
| 1/4 | 0.62 / 0.00 | 1.0049 / 1.0042 | 1.556 / 1.522 | 1.365 / 1.351 | 1.55 / 1.52 | 3.30 / 3.04 |
| 5/16 | 0.49 / −0.34 | 1.0036 / 1.0029 | 1.261 / 1.245 | 1.177 / 1.170 | 1.26 / 1.25 | 2.07 / 1.96 |
| **1/3** | **0.44 / −0.44** | **1.0030 / 1.0024** | **1.205 / 1.193** | **1.140 / 1.134** | **1.20 / 1.19** | **1.87 / 1.79** |
| 3/8 | 0.34 / −0.62 | 1.0017 / 1.0012 | 1.125 / 1.118 | 1.086 / 1.083 | 1.12 / 1.12 | 1.62 / 1.56 |
| 7/16 | 0.17 / −0.82 | 0.9992 / 0.9988 | 1.052 / 1.049 | 1.036 / 1.034 | 1.05 / 1.05 | 1.41 / 1.37 |
| 1/2 (frozen) | 0.00 / −0.88 | 0.9959 / 0.9957 | 1.008 / 1.007 | 1.004 / 1.003 | 1.01 / 1.01 | 1.29 / 1.27 |

The `h = 1/2` row reproduces the frozen design aid (readings to `4e-10`,
variance ratio to `2e-15`), and the design aid's fits reproduce the
multi-window study's to `6e-8` relative.

**Proposed `h = T/3`.** Both pairs read far from zero and with opposite signs
(`cos(pi/3) = 1/2`, `cos(2 pi/3) = -1/2`). The centre term alone and every
wrong-link variant read at least 0.08 above 1.05, while the variance ratio
rises from 1.29 to 1.87, which sets the price (Sect. 7). At `h = T/3`
(case 0 / 1):

- Weights along the truth: centre 1.369 / 1.347, each inner window
  −0.208 / −0.195, each outer window +0.023 / +0.022.
- Contributions: centre +1.205 / +1.193, inner pair −0.181 / −0.171, outer
  pair −0.020 / −0.019.
- Corrected reading 1.0030 / 1.0024 (full-tensor error 0.42% / 0.44%). With
  the plane-wave operator, 1.0012 / 1.0011.
- With flipped links the side windows still contribute about −0.06 instead
  of −0.20, because the operator `B` mixes electric and magnetic parts; the
  variant reads 1.14, not the centre term's 1.20.
- **The outer pair is not load-bearing for the amplitude.** Deleting it gives
  1.0234 / 1.0214, which would pass. The same holds at every spacing in the
  table: the outer pair contributes between 0 (at `T/4`, where the outer
  windows sit on nodes) and −0.020. The five-point stencil's outer
  coefficients are small. Only gate G3 below tests the outer windows and
  their chained links, through the sign of their linked readings.

Smaller spacings make the side windows heavier but cost more variance (ratio
3.0–3.3 at `T/4`, where the outer pair sits on nodes again; about 16 at
`T/6` by the first-order estimate). Dyadic alternatives (`5/16`, `3/8`) are
listed in Sect. 12.1. `h = 1/3` is not exactly representable; offsets are the
floats `k * (1/3)`, and the chart move between neighbours adds their float
difference (Sect. 3).

## 3. What changes and what does not

The runner and checking script are reviewed copies of the frozen ones
(`intrinsic_weyl_window_test.py`, `check_intrinsic_weyl_window_test.py`)
with only the changes below. The frozen files are not modified.

Changes:

1. **Windows.** `WINDOWS = (-2/3, -1/3, 0, 1/3, 2/3)` as the floats
   `k * (1/3)`, `k = -2..2`; window `w` is `build_window("sinusoid", A0, 2,
   delta_w)`. Inner neighbours as before (`INNER = {0: 1, 1: 2, 3: 2, 4: 3}`).
   An event's chart coordinates in windows `delta1` and `delta2` satisfy
   `x(delta1) - x(delta2) = (delta2 - delta1, delta1 - delta2, 0, 0)` in
   float arithmetic (the window at offset `delta` carries the profile
   `a(u + delta)`). The frozen runner's chart-move selfcheck (equal orders and
   exact squared intervals to `1e-10`) is run on the new pairs. Window labels
   use the index or `f"{delta:+.4f}"` (the frozen `+.1f` labels would read
   `+0.3`, `+0.7`).
2. **Stencil.** Five-point second order at `h = 1/3`:
   `W_0 = I + 22.5 B + 486 B^2`, `W_(+-h) = -12 B - 324 B^2`,
   `W_(+-2h) = 0.75 B + 81 B^2`, by the formula of `correction_specs`, which
   only offers `h` in {1/4, 1/2, 1}; the design aid's `stencil_weights`
   reproduces `correction_specs` bit for bit at 1/4 and 1/2. The checker
   recovers `B = -h^2 (W_h + 4 W_(2h))` and checks the structure. Saved
   alongside, descriptive: the single window, five-point second order with
   the plane-wave operator, three-point first order (windows 0, ±1/3) with
   both operators.
3. **Evaluation** (after the corrections are saved) adds the counterfactual
   readings of Sect. 4 and the linked readings of the side windows. It first
   replays the primary reading from the saved weights, links and fits. A
   difference above `1e-9`, or a failure of the exact Weyl maps of the
   counterfactuals, is a code defect and stops the run (Sect. 5).
4. **Gates** of Sect. 5: F1 and F2 with their verdicts named for what they
   compute; G3 part of the primary claim; secondary G4 (side share) and G5
   (link parity against the oracle).
5. **Failures and interruptions** (Sect. 5): a `MemoryError` anywhere and a
   worker's non-zero exit stop the run instead of failing a replicate. The
   frozen runner recorded both as replicate failures, contrary to its
   protocol's interruption rule; that discrepancy is noted in the frozen
   record and is not carried over.
6. Seeds, replicates, destinations, four fit workers (Sect. 6).
7. Selfchecks added: the stencil against `correction_specs` (bit for bit at
   1/4, 1/2, 1) and its structure at 1/3; the chart move, the GPU relation
   bits and the links on the new offsets; the new gates. `--check-cpu` runs
   every selfcheck that needs no GPU, the generator's and the jet inverse's
   included.
8. **Checking script** (Sect. 11): every interval the record quotes, the link
   parities against the oracle, the correction keys and shifts, the run's
   constants against this protocol, its source digests against the stamped
   commit.

Unchanged from the frozen protocol, Sects. 2–6: the wave and amplitude
(`a(u) = A0 sin(pi u)`, `A0 = pi^2 h0 / 2`, `h0 = 0.01`), tips, 4096
selection events per window and their streams, the nested selection, the
order-only shared-landmark rule (at least two own selection events between
each tip and the candidate), the clock rule, the archive, measurement
generation and relation bits, the bounded exact integer second moments, the
fifth-degree inverse with curvature-square terms and its settings, the blind
boundary, links over O(3) by orthogonal Procrustes chained outward, at least
eight shared landmarks per link, the alias map and minimum-norm operator from
the centre window, the truth, the readings, the statistics, strictness, and
the list of replicate failures and the rerun rule (restated in Sect. 5).

## 4. Readings

For replicate `r`: `c_s` the pure fit of window `s` in its own inferred frame,
`R_s` its order-only link to the centre window (O(3), chained), `M(R)` the Weyl
map of `R` (the magnetic part with the pseudotensor sign), `W_s` the saved
stencil weights, and `c` the truth in the centre window's inferred frame.
With the full-tidal inner product and `<x> = <c, x> / <c, c>`:

- corrected `a_r = <sum_s W_s M(R_s) c_s>`; single `s_r = <c_0>`;
- **F1 variant**, centre term alone: `a0_r = <W_0 c_0>`;
- **F2 variants**, only the side links replaced, everything else as in `a_r`:
  (i) identity links, `M(R_s) -> I` (each side estimate in its own inferred
  frame, unrelated to the centre's); (ii) flipped handedness,
  `R_s -> -R_s` (the order-only link composed with spatial inversion, which
  flips the magnetic part); (iii) Haar-random links, `R_s` replaced by
  independent Haar-random elements of O(3): for replicate `r`,
  `default_rng([2026113050, r])`, 64 draws, each draw the four side windows in
  the order of their offsets, each element the QR of a standard normal 3x3
  matrix with the sign fix (the post hoc script's construction, seeded per
  replicate so that replicates evaluated in parallel draw the same matrices);
- **G3 readings**, each side window's linked reading `rho_(r,s) = <M(R_s) c_s>`;
- **G4 reading**, the side windows' share
  `d_r = a_r - a0_r = <sum_(s != 0) W_s M(R_s) c_s>`;
- **G5 comparison**, for each side link the oracle relative frame
  `O_s = A_s^T A_0`, where `A_w` is the O(3) alignment of window `w`'s inferred
  coordinates to its true normal coordinates over the interior landmarks
  (computed after the corrections are saved). The link agrees in parity if
  `det(R_s O_s^T) > 0`; its angle error is the rotation angle of `R_s O_s^T`.

What the variants can show:

- **Variants (i) and (iii)** put every side estimate in a frame unrelated to
  the centre's. The Haar mean of a spin-2 map is zero, so their side
  contributions average out, and they are expected to read within about 0.01
  of the centre term alone. They show that the side windows count only
  through frames related to the centre's, and add little to F1 beyond that.
- **Variant (ii)** is the only one that keeps each link's rotation and changes
  only its handedness. It models just one form of handedness error
  (Sect. 7).
- **The side share `d`** is measured far more precisely than `a`. In the
  frozen run re-weighted at `T/3`, its per-replicate standard deviation is
  0.017 against 0.097 for `a`, and it is uncorrelated with `a0`
  (correlation 0.01).

Over the completed replicates: means and their 95% two-sided Student-t
intervals with `R - 1` degrees of freedom from the sample standard deviations
(house default; method fixed now, as in the frozen test).

## 5. Gates

- **G1, accuracy (primary).** ACCURATE if the interval of `a` lies inside
  `[0.95, 1.05]`; INACCURATE if entirely below 0.95 or entirely above 1.05;
  otherwise UNRESOLVED.
- **G2, control.** BIASED_LOW if the interval of `s` lies entirely below 0.95;
  otherwise NOT_BIASED.
- **F1, the centre term alone misses.** CENTRE_TERM_MISSES if the interval of
  `a0` lies entirely above 1.05; otherwise NOT_SHOWN. F1 uses the centre
  window only. That the side windows carry the result follows from F1
  together with G1 ACCURATE, not from F1 alone.
- **F2, wrong links miss.** WRONG_LINKS_MISS if the intervals of variant (i),
  of variant (ii) and of each of the 64 Haar draws (iii) all lie entirely
  above 1.05; otherwise NOT_SHOWN. As Sect. 4 explains, (i) and (iii) nearly
  restate F1, and (ii) is the independent part.
- **G3, link sign per side window.** For each side window, LINKED if the
  interval of `rho_s` lies entirely on the predicted side of zero (positive
  for ±T/3, negative for ±2T/3); otherwise NOT_LINKED. It is an aggregate
  sign test, and the only gate on the outer windows and their chained links.
- **G4, side share (secondary).** AS_PREDICTED if the interval of `d` lies
  inside `(-0.22, -0.17)`; NOT_AS_PREDICTED if it lies entirely below −0.22
  or entirely above −0.17; otherwise UNRESOLVED.
  - The band is the continuum value of the two layouts (−0.2018 / −0.1904),
    widened by 0.02 on each side.
  - The expected finite-N shifts fall well inside that widening. The operator
    `B` reads about 3% low at finite `N` (along the truth `0.01231 +- 0.00021`
    in the frozen run, against 0.0127 in the continuum), which moves `d` by
    about +0.006. The side windows' finite-N bias of about 1% moves it by
    about 0.002.
  - G4 is secondary because its band rests on only two continuum layouts.
- **G5, link handedness against the oracle (secondary).** ALL_AGREE if every
  side link of every complete replicate agrees in parity with the oracle
  relative frame; DISAGREE if any link disagrees (the count is reported);
  UNAVAILABLE if no link disagrees but some complete replicate has no oracle
  frames. Angle errors are descriptive.
- **Primary claim.** INCOMPLETE if fewer than `ceil(0.9 R) = 87` replicates
  complete. Otherwise SUPPORTED if G1 is ACCURATE, G2 BIASED_LOW, F1
  CENTRE_TERM_MISSES, F2 WRONG_LINKS_MISS and G3 LINKED for all four side
  windows; otherwise NOT_SUPPORTED.
- **What the record says.**
  - It states that the side windows and the links carry the result only if
    the claim is SUPPORTED.
  - If G1 and G2 pass but F1, F2 or G3 does not, it names the gate (for G3,
    also the window) and says that the test did not show the side windows or
    the links to carry the result.
  - G4 and G5 never change the primary claim. The record reports them
    whatever they show. When G4 is not AS_PREDICTED, or G5 is not ALL_AGREE,
    it says so next to the primary claim, with the interval of `d` or the
    count of disagreeing links.
  - All gates are computed on the completed replicates and reported whatever
    the primary claim.

Strictness as frozen: an interval bound equal to 0.95 or 1.05 (or to 0 for
G3, or to a band edge for G4) does not lie inside, below or above it. F1 and
F2 name the predicted side (above): removing the side windows' negative
contributions raises the reading, so a variant below 0.95 would contradict
the mechanism and does not count.

House rule:

- G1 and G4 are accuracy gates. Each is decided on its 95% two-sided interval
  and sized for power in Sect. 7 with that construction. G1 sets the number
  of replicates; G4's power at that number is given there.
- F1 and F2 are claims that a variant misses the band. They use the same
  interval, and their power is given in Sect. 7.
- G2 and G3 are sign gates, and G5 is a count.

**Failures.** A replicate fails if any of the following happens:

- one of its windows fails selection (the nested selection or the tip-clock
  requirement);
- a clock has zero measurement count;
- the inverse's validity rule rejects a fit (optimizer failure, rank
  deficiency, boundary excess, non-positive full-tidal covariance);
- a link has fewer than eight shared landmarks;
- the truth cannot be evaluated for the centre window.

A window that fails in generation stops its replicate. The runner records as
such a failure every `ValueError`, `ArithmeticError` or `LinAlgError` raised
in a window's selection, generation or fit, or in the replicate's correction
or truth evaluation. That is the frozen classification, and the listed causes
raise these types. Failures stay in the result and are never replaced: no
redraw, reseeding, sample extension, level change or tolerance change after
outcomes.

**Interruptions** are never outcomes. The runner stops the whole run, writes
no summary and no verdict, and leaves status STOPPED with the reason, on any
of these:

- a `MemoryError` anywhere, CuPy's out-of-memory error included;
- a fit or correction worker that exits with a non-zero code (an uncaught
  exception, a crash);
- a driver error;
- a failure of the exact counterfactual maps;
- a replay of the primary reading that differs by more than `1e-9`.

A power loss stops the run too. The generation loop stops at its next window;
fits already running finish, and queued ones are cancelled. A stopped run is
INCOMPLETE.

A rerun after an engineering repair is allowed only after such a stop. It
goes to a new destination (Sect. 6) and must be disclosed together with the
original stop. It counts only if it reproduces every replicate the stopped
run completed: order archives bit for bit where kept, selection records bit
for bit, fits and corrections field for field. Otherwise the stopped run's
INCOMPLETE stands.

## 6. Replicates, seeds, events

- `R = 96` replicates of five windows, `N = 2^25` measurement events per
  window, seeds `s_r = 2026112800 + r`. Total `96 x 5 x 2^25 = 1.61e10`
  measurement events and 1,966,080 selection events.
- Engineering pilot (Sect. 8): seeds 2026113000, 2026113001.
- **Destinations.**
  - The test: `outputs/intrinsic_weyl_window_link_test_v1`.
  - A rerun after a stop (Sect. 5): `_v2`, then `_v3`, each disclosed with
    the stop before it.
  - The pilot replay at the freeze commit:
    `outputs/intrinsic_weyl_window_link_test_engineering_v2`, compared with
    the draft pilot `outputs/intrinsic_weyl_window_link_test_engineering_v1`.
  - The runner refuses an existing destination.
- Seeds already drawn by the design aid's selection dry run, excluded from
  the test: 2026112900–2026112903 (selection events only, `h = 1/3`) and
  2026112700 (a replay of the frozen pilot's selection at `h = 1/2`). Both
  recorded runs of the aid drew these same events.
- Archives kept for replicates 0 and 48, hashed for the others; selection
  records, truth landmarks, pure fits and generation records kept for all.
- No data of the frozen test enter this test. The frozen run's saved fits
  entered only the design aid (Part 2, post hoc, for the variance and the
  centre-term reading); the frozen record's verdicts are unchanged.

**Dry run at `h = 1/3`** (four design replicates, 20 windows; it replays the
frozen runner's selection steps, which it reproduces bit for bit at `h = 1/2`
on the frozen pilot's seed): side windows 160–186 landmarks including tips,
31–53 shared per link (median 43), 4,986–7,481 clocks (median 6,147); centre
windows 120–133 landmarks, 3,485–4,195 clocks. The frozen run's side windows
had 23–42 shared landmarks (median 39) and a median of 6,209 clocks, so a
window costs about what it did and no cap on shared landmarks is needed.

## 7. Power and cost (before data)

**Standard deviation of one replicate's reading at `2^25`,** the largest of
the available estimates (design aid, Part 4):

| reading | calibration single rms (0.0797) x continuum ratio | frozen run re-weighted (Part 2) | design value |
|---|---|---|---|
| single | 0.0797 (realized in the frozen run: 0.0726) | — | 0.0797 |
| corrected | 0.109 | 0.097 | **0.109** |
| centre term alone | 0.107 | 0.095 | 0.107 |
| Haar links (one draw) | 0.117 | 0.106 | 0.117 |

The design value takes the calibration's single-window rms (0.0797, about 10%
above the frozen run's realized 0.0726) times the continuum variance ratio;
with the realized 0.0726 instead, the first column would read 0.099 / 0.097 /
0.109.
The re-weighted frozen run uses the recorded replicates' own noise with the
new spacing's weights, each from that replicate's operator `B` (along the
truth `0.01231 +- 0.00021`); its centre-term-alone reading, exact for those
centre windows, is 1.1990. Per-replicate standard deviations are resampled
from the calibration's selection-by-orientation spread; `kappa` inflates the
variance.

**G1, replicates for 90% power at `t = 0.05`, `2^25`** (continuum reading
1.0027, the mean of the two layouts):

| mean reading | `kappa = 1` | `kappa = 1.2` |
|---|---|---|
| continuum + 0.01 (1.0127) | 96 | 112 |
| continuum (1.0027) | 72 | **80** |
| continuum − 0.01 (0.9927) | 80 | 88 |
| continuum ± 0.02 | 192 / 120 | 256 / 160 |

**Sizing rule, fixed now:** the smallest `R` with 90% power in each of three
scenarios: `kappa = 1.2` at the continuum reading, and `kappa = 1` with the
mean 0.01 above or below it. This gives **`R = 96`**. At `R = 96` the power of
G1 is 0.986 / 0.960 (`kappa` 1 / 1.2) at the continuum reading, 0.917 / 0.856
at +0.01, 0.966 / 0.934 at −0.01, 0.68 / 0.60 at +0.02 and 0.83 / 0.76 at
−0.02. In the
`kappa = 1.2` continuum scenario the price is: break-even (50% power) 48, 80%
power about 64, 90% power 80, so `R = 96` buys the ±0.01 allowance.

Disclosure: an earlier development run of the aid sized the test at
`R = 88`, under a narrower rule (90% power at `kappa = 1.2` with the mean
0.01 below the continuum).

- It ran at commit `c87fdab` with the aid not yet committed, on the same
  dry-run seeds. Its output, `outputs/intrinsic_weyl_window_link_test_design/power_debug.json`,
  is not kept in the repository.
- The rule was widened to the three scenarios above before the first
  recorded run.
- The change moved the design to more replicates, and neither run drew a
  measurement event.

**Other gates at `R = 96`, `2^25`, `kappa` up to 1.2:** G2, F1, F2 (each
variant and each Haar draw), and G3 for each side window all have simulated
power 1.0000 in 20,000 draws.

- Margins at `kappa = 1.2`, for the expected lower end of the interval above
  1.05: about ten standard errors of the mean for the centre term alone
  (1.199), nine for the Haar variant (1.198), five for the flipped-link
  variant (1.137).
- The single window's (0.885) expected upper end lies about five below 0.95.
- G3 now belongs to the primary claim at no cost: the linked readings lie
  about 0.44 from zero, against an interval half-width of about 0.02.

The joint power of the primary claim is therefore that of G1 to the quoted
precision.

**G4 (secondary) at `R = 96`, `2^25`** (design aid, Part 4; drawn after every
other power, so that the earlier numbers are unchanged):

- **Standard deviation of `d`.** The design value is the larger of two
  estimates, **0.0225**:
  - the calibration single rms times the side share's continuum variance
    ratio (0.0798 / 0.0677 for the two layouts), which gives 0.0225;
  - the frozen run re-weighted, which gives 0.0172.

  The pilot's nominal 0.017 (Sect. 8) agrees with the smaller one.
- **Power under G1's sizing rule,** applied to each layout's continuum value
  (−0.2018 / −0.1904):
  - 1.000 at `kappa = 1.2` at the continuum values;
  - at `kappa = 1`, with the mean 0.01 below / above: 0.942 / 1.000 for
    case 0, 1.000 / 0.995 for case 1.

  So the `R` that G1 sets gives G4 at least 94% power under the same rule.
- **Outside that rule,** with the mean shifted 0.01 toward the nearer edge
  (case 0 below, case 1 above):
  - at `kappa = 1.2`, 0.900 and 0.984;
  - at a stress standard deviation of 0.03, 0.756 and 0.921, against
    1.000 at the continuum values.

  G4 is secondary; these are its weak corners.

**What the gates can see of the links.**

- **Fine accuracy is out of reach.** The amplitude is second order in a small
  link error. Rotating every side window's link by a fixed angle about a
  random axis moves the corrected reading by +0.006 (10°), +0.024 (20°),
  +0.050 (30°), +0.11 (45°) and +0.25 (90°) at case 0 (mean of 256 draws;
  case 1 alike). So G1 and F2 certify the links only to roughly 20–30° in
  aggregate.
- **G3 is an aggregate sign test.** It fails if a side window's linked
  reading has the wrong sign on average. A window whose links are all flipped
  (`R -> -R`) reads about ∓0.015, on the wrong side of zero (continuum, both
  layouts: a window reads the electric part lower than the magnetic, so the
  flipped projection takes the opposite sign). Under a rotation error about a
  random axis, the mean linked reading changes sign near 72° (the zero of the
  spin-2 character).
- **G3 does not certify each link's handedness.** A real parity error of the
  Procrustes step is a reflection `R sigma` through some plane. Applied to the
  true tensor, a reflection through a plane that contains the propagation
  direction and a polarization axis of this wave (`sigma_x`, `sigma_y` in the
  generator frame) leaves its projection onto the truth at exactly 1.
  `R -> -R`, which is `sigma_z` composed with a half-turn about the
  propagation direction, sends that projection to 0. (Checked with the
  repository's `weyl_rotation`; the half-turn alone leaves it at 1.)
- **Errors confined to a minority of replicates barely move the means.** With
  10% of replicates carrying flipped links, `a` rises by only about
  0.012–0.014 (frozen run re-weighted / continuum), and an inner window's
  linked reading falls only to about 0.39.
- **What does catch them.** G5 compares each link's parity with the oracle
  directly. G4 also responds: each flipped replicate moves `d` by about
  +0.12 to +0.14, so `d` leaves its band once roughly 15–25% of the
  replicates are flipped, depending on the layout. Angle errors against the
  oracle stay descriptive.

**Cost** (model from the frozen run's measured times; it gives 8.53 h for the
frozen run, which took 8.62 h). Generation 96 s per window at `2^25`,
independent of the landmark count in the frozen run; side-window fits
`311 s x (clocks/6209)^1.63`. At `R = 96`, `2^25`:

- **1.61e10 measurement events**, 1.97e6 selection events;
- **GPU busy 12.8 h**, the GPU reserved for the whole run;
- wall time **about 12.8 h** by this model (up to about 16 h by the
  engineering pilot's timings, Sect. 8) with four fit workers (generation-bound; peak
  memory at most 5.1 GB per fit worker in the frozen run), up to 16.3 h with
  three if every side window had the dry run's largest clock count;
- 49 CPU hours of fits.

For comparison, `R = 80` costs 1.34e10 events and 10.7 h; `R = 112` costs
1.88e10 events and 14.9 h.

## 8. Engineering pilot (before the freeze, not test data)

Two replicates, seeds 2026113000 and 2026113001, at `2^22` and `2^24` events
per window (nested prefixes of one sprinkle), engineering mode of the draft
runner, run on 2026-11-21 from the clean worktree at commit `7e6921d`
(`outputs/intrinsic_weyl_window_link_test_engineering_v1`; runner SHA-256
under the LF rule
`1323216bea8c243c6f19a9872ba73844a8d3922180a8193566428f5c08774068`), four
fit workers. Other applications held the GPU at about 24% utilization
throughout; there was no other compute job.

- Completed in 27.4 minutes; all four replicate-levels complete, no failure,
  every fit valid. Every selfcheck passed at start, including the GPU
  relation bits on the new windows and the bounded second moments.
- Shared landmarks per link 39–48. Procrustes relative residual 0.14–0.22% at
  `2^24`, 0.28–0.40% at `2^22`.
- Mirror embeddings again: at every level three windows of replicate 0 (the
  centre among them) and two of replicate 1 converged to the mirror embedding
  of the truth. All 16 links had the oracle's parity, with angle errors
  0.017–0.091°, and the primary readings with oracle links equal the
  order-only ones to `7.6e-5`.
- Nominal standard deviation at `2^24`: primary 0.132 and 0.140, single
  window 0.097 and 0.103 (variance ratio about 1.85; the design used 1.79–1.87
  in the continuum and 0.154 for the primary at `2^24`).
- Generation 50–76 s per window to `2^24`, peak 1.6 GB; fit workers 86–119 s
  (centre) and 193–422 s (side windows, 5,560–6,888 clocks), one at 765 s,
  peak at most 2.9 GB; the correction and evaluation 6–7 s. Generation took
  up to 1.3 times and the fits 1.5–2 times as long as in the frozen pilot
  (three fit workers, no background load).
  Scaled from the frozen run (8.62 h for 64 replicates), the test needs about
  13 h; if its fits run as slowly as these, up to about 16 h.
- The checking script passed (receipt
  `outputs/intrinsic_weyl_window_link_test_engineering_v1.audit.json`,
  SHA-256 `8adb8ef7ded369b503786d00172d2224ef2841fa99d69cee6d8627d4bab3c32c`):
  all ten windows' landmarks and shared rows re-derived; links, estimates,
  readings and every counterfactual reading (the 64 Haar draws from their
  seed) reproduced at all four replicate-levels; the `h = 1/3` stencil
  structure to `4.4e-16`; the minimum-norm operator recomputed to `1.3e-14`;
  the four centre-window clock covariances (3,872 and 4,199 clocks)
  recomputed with the inverse's CPU float64 second moment, bit-identical to
  the saved ones; nested prefixes and the CPU replay of the first 4096
  measurement events of every window bit for bit; 5,997 order bits checked
  with the ODE integrator, all agree; no excluded candidate inside its
  window; shared landmarks are the moved landmarks of the neighbour, with
  exact squared intervals equal in both charts to `3.5e-13`; refits of window
  3 of replicate 0 at `2^24` and window 0 of replicate 1 at `2^22` reproduce
  every field.

For the record only (engineering data, two replicates, nothing is inferred
from them): primary readings 0.935 and 0.953 at `2^22`, 1.043 and 0.735 at
`2^24`; single window 0.869 and 0.770, 0.917 and 0.679. The side share `d`
was −0.207 and −0.100 at `2^22`, −0.200 and −0.188 at `2^24` (continuum
−0.196). Its nominal standard deviations, from the side windows' nominal
covariances with `B` held fixed, were 0.047–0.048 and 0.023–0.024, about
0.017 when scaled to `2^25`. The side windows' linked readings were +0.33 to
+0.59 (±T/3) and −0.37 to −0.53 (±2T/3).

The pilot ran the earlier runner (`7e6921d`). The revised runner (`6b4e1a4`)
changes only the failure and interruption paths, the gates and the
selfchecks; no field of a completed replicate's record changes. The checks
after the reviews (2026-11-22, CPU only, no events):

- `--check-cpu` passes.
- An interruption test runs the real dispatcher in engineering mode with
  every GPU step stubbed. A worker that crashes stops the run before the
  second replicate: 5 of 40 windows generated, status STOPPED with the
  reason, no summary. A `MemoryError` inside a fit propagates, while a
  `ValueError` is still recorded as the replicate's failure.
- On synthetic 96-replicate levels (eight cases, including INCOMPLETE, a
  failed G3, G5 disagreement and levels with one or two complete replicates),
  the revised checking script agrees item by item with the runner's gates. A
  tampered interval in the record is caught.
- On the pilot's four replicate-levels, the checking script reproduces every
  link's parity and angle against the oracle, and the primary reading with
  oracle links. The pilot's constants and source digests match this protocol
  and its stamp. Its requirement that the bound sources at its own checkout
  equal the stamped ones fails there for the runner, as it must; the pilot
  replay at the freeze commit meets it.

Done before the pilot, CPU only, no events: the runner's selfchecks that need
no GPU pass (gates, selection, links with mixed parity on the new offsets,
stencil, chart moves between the new neighbours to `6e-12`). A plumbing test
relabelled the frozen run's replicates 0 and 31 (whose centre clocks are kept)
to the new offsets and ran the draft correction and evaluation on them: the
primary reading replays exactly, and the draft checking script reproduces the
links, estimates, readings and every counterfactual reading including the 64
Haar draws, with the five-point structure at `h = 1/3` to `4e-16`; the
runner's and the checker's verdicts agree. Its readings mean nothing, since
those side windows were drawn at the old offsets; as expected, the relabelled
outer windows (old ±T, signal −0.88) read −0.87 to −0.97 and pass G3, and the
relabelled inner windows (old nodes) read about 0 and fail it.

## 9. Predictions (before data)

- `a` between 0.98 and 1.02 (continuum 1.0027; the frozen run's finite-N
  offset was +0.002). `s` between 0.865 and 0.905.
- Centre term alone `a0` between 1.17 and 1.23 (continuum 1.199; frozen
  centre windows 1.199).
- Identity links and each Haar draw between 1.16 and 1.24. Flipped links
  between 1.11 and 1.17 (continuum 1.137).
- The identity prediction rests on independent windows having unrelated
  inferred frames. Re-weighting the frozen run at `T/3` with identity links
  gives 1.202 against its centre term's 1.199. Aligned frames would have
  moved it by about −0.04 through the ±T windows (signal −0.88, weight
  +0.023 each).
- Linked readings: ±T/3 windows about +0.44 (0.40–0.48), ±2T/3 about −0.44.
- Side share `d` between −0.21 and −0.18 (continuum −0.196; frozen run
  re-weighted −0.195).
- Per-replicate standard deviation at `2^25`: corrected 0.09–0.11, single
  0.07–0.08, centre term alone 0.09–0.11, side share 0.015–0.025.
- Shared landmarks 30–55 per link; links improper in about half of the
  replicates; Procrustes relative residual about `1e-3`; every link with the
  oracle's parity and an angle error below 0.1°.
- Verdicts: G1 ACCURATE, G2 BIASED_LOW, F1 CENTRE_TERM_MISSES, F2
  WRONG_LINKS_MISS, G3 LINKED for all four side windows; primary claim
  SUPPORTED. G4 AS_PREDICTED, G5 ALL_AGREE.

## 10. Execution and provenance

As frozen, with this protocol's destinations (Sect. 6):

1. **The freeze.** Protocol, runner and checking script are committed
   together after the pilot.
2. **The pilot replay.** The pilot is replayed from the clean worktree at the
   freeze commit into `outputs/intrinsic_weyl_window_link_test_engineering_v2`
   and compared with `--compare` against
   `outputs/intrinsic_weyl_window_link_test_engineering_v1` before the test
   starts. This is the frozen runner's own pilot: it must reproduce the order
   bits exactly and every fit, correction and evaluation field to `1e-8`. The
   checking script then runs on it, at the freeze commit.
3. **The test.** One run from the clean worktree with no other GPU job, into
   `outputs/intrinsic_weyl_window_link_test_v1`. The runner binds the hashes
   of its sources and of this protocol, and refuses a dirty worktree or an
   existing destination.
4. **The audit.** The checking script runs with the refits
   `rep00:3:33554432` and `rep48:0:33554432`, from a checkout whose bound
   sources equal the stamped ones.
5. **The record.** Summaries, the audit receipt and a results note are
   committed.

Merge commits only, so the `code_version` stamps stay reachable. The pilot
replay (about half an hour of GPU) and the test (about 13 hours, up to about
16, with the GPU reserved) each wait for the user's go-ahead, since other
work may need the GPU.

## 11. Checking script

**Inherited from the frozen checker,** with the new spacing: links by polar
decomposition, Weyl maps through the 3x3 electric and magnetic tensors,
estimates, readings, intervals, verdicts, selection records, archives of two
replicates, order bits by an ODE integrator, minimum-norm operator, clock
covariance and refits.

**Added for this test:**

- The `h = 1/3` stencil structure, recovering `B = -h^2 (W_h + 4 W_(2h))`.
  Every correction's key and shifts are checked against Sect. 3, which pins
  `h` for all replicates. The structure check alone cannot tell a wrong `h`
  apart from a rescaled `B`.
- Every counterfactual reading and every linked side reading, recomputed
  from the saved fits, links and weights (the 64 Haar draws from their seed),
  and compared with both the evaluation records and the summary rows.
- Each side link's parity and angle against the oracle relative frames, and
  the primary reading with oracle links, recomputed from the truth landmarks
  by polar decomposition.
- **Every gate, restated.** The verdicts of G1–G5 and the primary claim, and
  every interval the record quotes: the corrected, single, centre-term,
  identity, flipped and side-share readings, the Haar summary (draws missing
  the band, lowest lower bound, range of means) and the four link-sign
  intervals. Each is compared item by item to `1e-9`. Fewer than two
  complete replicates give no interval and no crash.
- **Protocol, not the run.** The run's seeds, levels, windows, spacing,
  tolerance, coverage, minimum count, kept archives, primary correction and
  Haar seed are checked against this protocol, not taken from the run. It
  must have started from a clean worktree.
- **Provenance.** Its source digests are checked against the blobs of the
  stamped commit. The bound sources at the checking script's own checkout
  must equal the stamped ones, so that the refits run the stamped code.

**Not independent of the runner:**

- the frozen inverse's nested selection and fit (refits);
- the generator's window geometry (truth side);
- the alias map, whitening and minimum-norm operator of
  `intrinsic_weyl_window_study`.

The last of these recomputes `B` only for the two replicates whose clocks are
kept. For the others the checking script verifies only the five-point
structure of the saved weights. `B` sets both the centre term (1.20) and the
side weights, and the reading moves by about 10 per unit change of `B` along
the truth. The frozen test's record notes the same limit.

## 12. Decisions

### 12.1 Taken by the user on 2026-11-21 (the proposed choices)

1. **Spacing.** Proposed `T/3`. `3/8` is dyadic and cheaper (ratio 1.62;
   roughly 85 replicates by scaling the variance ratio, not computed
   directly) but its flipped-link variant reads only 1.086, 0.036 above the
   band, and its centre term 1.12. `5/16` is dyadic and more decisive (centre
   term 1.26) but costs ratio 2.07 (roughly 105 replicates).
2. **Tolerance.** Proposed `t = 0.05`, as the frozen test. `0.06` needs 64
   replicates by the same sizing rule (8.5 h) and still separates the single
   window (0.885) and every variant (flipped links 1.137, 0.077 above 1.06);
   `0.04` needs 192 (25.6 h).
3. **Replicates.** Proposed 96 by the sizing rule of Sect. 7. 80 if sized as
   the frozen test was (`kappa = 1.2` at the continuum reading only).
4. **F2 composition.** Proposed identity, flipped and 64 Haar draws, all
   required. A minimal form keeps only the Haar draws.
5. **G3.** Proposed secondary, and taken; revised in 12.2.
6. **Fit workers.** Proposed four (generation-bound, about 12.8 h); three as
   frozen (up to 16.3 h). Neither changes any number.
7. **Seeds.** As in Sect. 6.

### 12.2 Revisions after the independent reviews (2026-11-22)

A design and statistics review and a code review read the draft, the runner,
the checking script, the design aid and the pilot. They found no error in the
design numbers and one blocking defect.

**Proposed after the reviews, confirmed by the user at the freeze
(2026-11-23):**

1. **G3 in the primary claim.** It is the only gate on the outer windows and
   their chained links, and its power is 1.0000, so it costs nothing. As a
   secondary gate it would have let the claim stand with a window on the
   wrong side of zero.
2. **G4, side share (secondary).** It measures the side windows' share with
   an interval about a fifth the width of `a`'s.
3. **G5, link parity against the oracle (secondary).** G1 and G3 cannot
   certify each link's handedness (Sect. 7).

**Corrected (defects and overstatements):**

- **Interruptions** (blocking). The runner recorded a `MemoryError` or a
  worker's non-zero exit as a replicate failure, while the protocol treats
  them as interruptions. A shared GPU could have removed replicates silently,
  plausibly those with the most clocks. They now stop the run (Sect. 5). The
  frozen runner has the same discrepancy, recorded in its post hoc note.
- **F1 and F2** now report CENTRE_TERM_MISSES and WRONG_LINKS_MISS. A variant
  that misses the band does not by itself make anything load-bearing, and
  F2's identity and Haar variants nearly restate F1 (Sect. 4).
- **Handedness claims** are narrowed to aggregate statements (Sect. 7).
- **The checking script** now compares every interval the record quotes,
  checks the run's constants, source digests and correction shifts against
  the protocol and the stamp, recomputes the links' oracle parities, and
  discloses the operator code it shares with the runner (Sect. 11).
- **Provenance.** The Sect. 1 numbers now cite their record (`d4121a2`,
  merged here), and its Haar range is corrected. The earlier development run
  of the design aid is disclosed (Sect. 7).
- **Wording.** The outer pair's contribution at `T/4`, a mislabelled table
  column, the chart-move sign, the destinations, and `h` and `T` as clock
  inputs.

## 13. Scope

Whatever the verdict, this run would establish at most the following, for
one plane wave at one amplitude and interval size:

- that the order-only five-window correction at spacing `T/3` removes the
  fifth-degree amplitude bias along the true tensor at finite `N` to within
  ±0.05, with the single window as control;
- that this result depends on the inner side windows (the centre term alone
  overshoots);
- that it depends on the order-only links being related to the centre's
  frame, right in aggregate handedness and roughly right in orientation
  (identity, flipped and random links fail, and every side window's linked
  reading has its predicted sign).

It would not establish any of the following:

- **Each link's handedness.** G5 compares every link with the oracle, but
  the amplitude gates certify handedness only in aggregate, and some
  reflections leave the reading unchanged.
- **The links' fine accuracy.** Sect. 7: they are certified only to about
  20–30°.
- **A load-bearing role for the outer pair in the amplitude.** It contributes
  −0.02, and G3 tests only the sign of its linked readings.
- **Full-tensor accuracy**, or the behaviour for other waves or fields.
- **The accuracy of one replicate.** Each replicate's own corrected reading
  has a standard deviation of about 0.1, twice the tolerance.
- **The cost of a single-sprinkle implementation.**
- **A correction free of clock inputs.** The correction takes `h` and `T`
  from the observer's clock.
