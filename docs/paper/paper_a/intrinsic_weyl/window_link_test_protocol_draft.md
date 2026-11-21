# Finite-N test of the five-window correction with load-bearing side windows and links (DRAFT)

2026-11-21. **Draft, not frozen.** Nothing here binds until a freeze commit
that renames this file to `window_link_test_protocol.md` together with the
runner and the checking script and precedes every test selection and
measurement event (Sect. 10). The user took the proposed choices of Sect. 12
on 2026-11-21 (spacing `T/3`, `t = 0.05`, 96 replicates, F2 as proposed, G3
secondary, four fit workers, the seeds of Sect. 6). The runner
`experiments/positive_control/intrinsic_weyl_window_link_test.py` and the
checking script `check_intrinsic_weyl_window_link_test.py` are written and
committed as drafts; the engineering pilot (Sect. 8) comes next and waits for
the user's go-ahead. All design numbers come from the design aid
`experiments/positive_control/intrinsic_weyl_window_link_test_power.py`
(commit `a58965c`), whose recorded run is
[`window_link_test_power.json`](window_link_test_power.json). That run drew no
measurement event; its only draws are selection events of a dry run (Sect. 6).

## 1. Why a second test

The [frozen five-window test](window_test_protocol.md) (commit `1ca7af2`,
[record](window_test_results.md), commit `02c21e2`) is SUPPORTED: on the
Ibrahim–Bombelli half-cycle wave the order-only five-window correction read the
crest amplitude along the true tensor as 0.9976 [0.9772, 1.0180] over 64
replicates, while the single window read 0.8870 [0.8689, 0.9051]. That record
stands as frozen and is not re-scored here.

A later independent review found the test nearly blind to its side windows and
links. The post hoc split of the recorded reading
(`intrinsic_weyl_window_test_sensitivity.py`, commit `c87fdab`) shows why:

- With windows at proper-time offsets 0, ±T/2, ±T (spacing `h = T/2`, `T` the
  tip separation) on a wave of wavelength `2T`, the ±T/2 windows sit on nodes
  and read about 0, and the ±T windows carry a weight of only 0.0067 along the
  truth. The corrected reading is about `1.139 x` the centre window's own
  estimate `- 0.012`, where 1.139 is the component of `I + 10B + 96B^2` along
  the truth, computed from the centre window's own inferred geometry.
- Deleting the side windows gives 1.0086 [0.9883, 1.0290]; identity or
  Haar-random side links give 1.0086–1.0098. Every such variant would have
  passed. The accuracy of the order-only links rests only on descriptive
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
the single window on the same data reads low by more than 0.05, and
(falsification) the centre term alone and the correction with identity,
handedness-flipped or Haar-random side links all read above 1.05?**

Scope, fixed now: the amplitude along the true tensor, one number. One wave,
one amplitude, one interval size. Full-tensor accuracy is descriptive.

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
  1.0234 / 1.0214, which would pass. This holds at every spacing in the table
  (outer pair −0.012 to −0.020): the five-point stencil's outer coefficients
  are small. Only gate G3 below tests the outer windows and their chained
  links, through the sign of their linked readings.

Smaller spacings make the side windows heavier but cost more variance (ratio
3.0–3.3 at `T/4`, where the outer pair sits on nodes again; about 16 at
`T/6` by the first-order estimate). Dyadic alternatives (`5/16`, `3/8`) are
listed in Sect. 12. `h = 1/3` is not exactly representable; offsets are the
floats `k * (1/3)`, and the chart move between neighbours adds their float
difference (Sect. 3).

## 3. What changes and what does not

The runner and checking script are to be written after the decisions of
Sect. 12, as reviewed copies of the frozen ones
(`intrinsic_weyl_window_test.py`, `check_intrinsic_weyl_window_test.py`)
with only the changes below. The frozen files are not modified.

Changes:

1. **Windows.** `WINDOWS = (-2/3, -1/3, 0, 1/3, 2/3)` as the floats
   `k * (1/3)`, `k = -2..2`; window `w` is `build_window("sinusoid", A0, 2,
   delta_w)`. Inner neighbours as before (`INNER = {0: 1, 1: 2, 3: 2, 4: 3}`);
   an event's chart coordinates in windows `delta1`, `delta2` differ by
   `(delta1 - delta2, delta2 - delta1, 0, 0)` in float arithmetic. The frozen
   runner's chart-move selfcheck (equal orders and exact squared intervals to
   `1e-10`) is run on the new pairs. Window labels use the index or
   `f"{delta:+.4f}"` (the frozen `+.1f` labels would read `+0.3`, `+0.7`).
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
   replays the primary reading from the saved weights, links and fits; a
   difference above `1e-9` is an `ArithmeticError` and fails the replicate.
4. Seeds, replicates, destination, four fit workers (Sect. 6).
5. Selfchecks added: the stencil against `correction_specs` (bit for bit at
   1/4, 1/2, 1) and its structure at 1/3; the chart move, the GPU relation
   bits and the links on the new offsets; the new gates (`--check-cpu` runs
   every selfcheck that needs no GPU).

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
the failure, interruption and rerun rules.

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
- **G3 readings**, each side window's linked reading `rho_(r,s) = <M(R_s) c_s>`.

Over the completed replicates: means and their 95% two-sided Student-t
intervals with `R - 1` degrees of freedom from the sample standard deviations
(house default; method fixed now, as in the frozen test).

## 5. Gates

- **G1, accuracy (primary).** ACCURATE if the interval of `a` lies inside
  `[0.95, 1.05]`; INACCURATE if entirely below 0.95 or entirely above 1.05;
  otherwise UNRESOLVED.
- **G2, control.** BIASED_LOW if the interval of `s` lies entirely below 0.95.
- **F1, side windows load-bearing.** LOAD_BEARING if the interval of `a0`
  lies entirely above 1.05; otherwise NOT_SHOWN.
- **F2, links load-bearing.** LOAD_BEARING if the intervals of variant (i),
  of variant (ii) and of each of the 64 Haar draws (iii) all lie entirely
  above 1.05; otherwise NOT_SHOWN.
- **G3, link sign per side window (secondary).** For each side window,
  LINKED if the interval of `rho_s` lies entirely on the predicted side of
  zero (positive for ±T/3, negative for ±2T/3); otherwise NOT_LINKED.
  Reported per window; not part of the primary claim.
- **Primary claim.** INCOMPLETE if fewer than `ceil(0.9 R)` replicates
  complete; otherwise SUPPORTED if G1 is ACCURATE, G2 BIASED_LOW, F1 and F2
  LOAD_BEARING; otherwise NOT_SUPPORTED. If G1 and G2 pass but F1 or F2 does
  not, the record says that the test did not show the side windows or the
  links to carry the result. All gates are computed on the completed
  replicates and reported whatever the primary claim.

Strictness as frozen: an interval bound equal to 0.95 or 1.05 (or to 0 for
G3) does not lie inside, below or above it. F1 and F2 name the predicted side
(above): removing the side windows' negative contributions raises the reading,
so a variant below 0.95 would contradict the mechanism and does not count.

House rule: G1 is an accuracy gate, decided on its 95% two-sided interval and
sized for 90% power in Sect. 7 with that construction. F1 and F2 are claims
that a variant misses the band; they use the same interval and their power is
given in Sect. 7. G2 and G3 are sign gates.

Failures as frozen: a replicate fails if any window fails selection, a clock
has zero count, a fit fails the inverse's validity rule, a link has fewer than
eight shared landmarks, or the truth cannot be evaluated; failures stay in the
result, no redraw, reseeding or extension.

## 6. Replicates, seeds, events

- `R = 96` replicates of five windows, `N = 2^25` measurement events per
  window, seeds `s_r = 2026112800 + r`. Total `96 x 5 x 2^25 = 1.61e10`
  measurement events and 1,966,080 selection events.
- Engineering pilot (Sect. 8): seeds 2026113000, 2026113001.
- Seeds already drawn by the design aid's selection dry run, excluded from
  the test: 2026112900–2026112903 (selection events only, `h = 1/3`) and
  2026112700 (a replay of the frozen pilot's selection at `h = 1/2`).
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

| reading | realized single x continuum ratio | frozen run re-weighted (Part 2) | design value |
|---|---|---|---|
| single | 0.0726 (frozen run) | — | 0.0797 (calibration rms) |
| corrected | 0.109 | 0.097 | **0.109** |
| centre term alone | 0.107 | 0.095 | 0.107 |
| Haar links (one draw) | 0.117 | 0.106 | 0.117 |

The design value takes the calibration's single-window rms (0.0797, about 10%
above the frozen run's realized 0.0726) times the continuum variance ratio.
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

**Other gates at `R = 96`, `2^25`, `kappa` up to 1.2:** G2, F1, F2 (each
variant and each Haar draw), and G3 for each side window all have simulated
power 1.0000 in 20,000 draws. Margins at `kappa = 1.2`: the expected lower
end of the interval lies about ten standard errors of the mean above 1.05 for
the centre term alone (1.199) and nine for the Haar variant (1.198), five for
the flipped-link variant (1.137); the single window's (0.885) expected upper
end lies about five below 0.95. The joint power of the primary claim is
therefore that of G1 to the quoted precision.

**What the gates can see of the links.** The amplitude is second order in a
small link error: rotating every side window's link by a fixed angle about a
random axis moves the corrected reading by +0.006 (10°), +0.024 (20°), +0.050
(30°), +0.11 (45°), +0.25 (90°) at case 0 (mean of 256 draws; case 1
alike). So G1 and F2 certify the links
only to roughly 20–30° in aggregate. Their reported agreement with the oracle
(about 0.05°) remains a descriptive comparison. G3 catches a link with the
wrong handedness: it leaves the side window's linked reading at about ∓0.015,
on the wrong side of zero (continuum, both layouts; a window reads the
electric part lower than the magnetic, so the flipped projection takes the
opposite sign). For a rotation error about a random axis the
mean linked reading changes sign near 72° (the zero of the spin-2 character),
so G3 also catches gross orientation errors.

**Cost** (model from the frozen run's measured times; it gives 8.53 h for the
frozen run, which took 8.62 h). Generation 96 s per window at `2^25`,
independent of the landmark count in the frozen run; side-window fits
`311 s x (clocks/6209)^1.63`. At `R = 96`, `2^25`:

- **1.61e10 measurement events**, 1.97e6 selection events;
- **GPU busy 12.8 h**, the GPU reserved for the whole run;
- wall time **about 12.8 h** with four fit workers (generation-bound; peak
  memory at most 5.1 GB per fit worker in the frozen run), up to 16.3 h with
  three if every side window had the dry run's largest clock count;
- 49 CPU hours of fits.

For comparison, `R = 80` costs 1.34e10 events and 10.7 h; `R = 112` costs
1.88e10 events and 14.9 h.

## 8. Engineering pilot (before the freeze, not test data)

Two replicates, seeds 2026113000 and 2026113001, at `2^22` and `2^24` events
per window (nested prefixes), engineering mode of the draft runner. It
checks the runner end to end at the new spacing: selection and links (shared
counts, parity against the oracle, angle errors), the new stencil's
structure, the counterfactual readings, the checking script, timings and
memory. Nothing is inferred from its readings.

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
  centre windows 1.199). Identity links and each Haar draw between 1.16 and
  1.24. Flipped links between 1.11 and 1.17 (continuum 1.137). The identity
  prediction rests on independent windows having unrelated inferred frames:
  re-weighting the frozen run at `T/3` with identity links gives 1.202
  against its centre term's 1.199, where aligned frames would have moved it
  by about −0.04 through the ±T windows (signal −0.88, weight +0.023 each).
- Linked readings: ±T/3 windows about +0.44 (0.40–0.48), ±2T/3 about −0.44.
- Per-replicate standard deviation at `2^25`: corrected 0.09–0.11, single
  0.07–0.08, centre term alone 0.09–0.11.
- Shared landmarks 30–55 per link; links improper in about half of the
  replicates; Procrustes relative residual about `1e-3`; angle error against
  the oracle below 0.1°.
- G1 ACCURATE, G2 BIASED_LOW, F1 and F2 LOAD_BEARING, G3 LINKED for all four
  side windows; primary claim SUPPORTED.

## 10. Execution and provenance

As frozen: protocol, runner and checking script committed together (the
freeze) after the pilot; the pilot replayed from the clean worktree at the
freeze commit and compared before the test starts; one run from the clean
worktree with no other GPU job; the runner binds the hashes of its sources and
of this protocol and refuses a dirty worktree or an existing destination;
checking script with the refits `rep00:3:33554432` and `rep48:0:33554432`;
summaries, audit receipt and a results note committed. Merge commits only, so
the `code_version` stamps stay reachable.
The pilot (tens of minutes of GPU) and the test (about 13 hours with the GPU
reserved) each wait for the user's go-ahead, since other work may need the
GPU.

## 11. Checking script

The frozen checker's independent reconstruction (links by polar
decomposition, Weyl maps through the 3x3 electric and magnetic tensors,
estimates, readings, intervals, verdicts, selection records, archives of two
replicates, order bits by an ODE integrator, minimum-norm operator, clock
covariance, refits) with the new spacing, plus: the `h = 1/3` stencil
structure (recovering `B = -h^2 (W_h + 4 W_(2h))`), the readings of every
counterfactual and of the linked side windows from the saved fits, links and
weights, the 64 Haar draws from their seed, and the verdicts of F1, F2 and G3.

## 12. Decisions (the proposed choices, taken by the user on 2026-11-21)

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
5. **G3.** Proposed secondary. It could be made a co-requirement of the
   primary claim; its power is 1.0000 either way.
6. **Fit workers.** Proposed four (generation-bound, about 12.8 h); three as
   frozen (up to 16.3 h). Neither changes any number.
7. **Seeds.** As in Sect. 6.

## 13. Scope

Whatever the verdict, this run would establish at most, for one plane wave at
one amplitude and interval size: that the order-only five-window correction at
spacing `T/3` removes the fifth-degree amplitude bias along the true tensor at
finite `N` to within ±0.05, with the single window as control, and that this
result depends on the side windows (the centre term alone overshoots) and on
the order-only links being right in handedness and roughly in orientation
(identity, flipped and random links fail). It would not establish the links'
fine accuracy (Sect. 7: certified only to about 20–30°), a load-bearing role
for the outer pair in the amplitude (it contributes −0.02; G3 tests only the
sign of its linked readings), full-tensor accuracy, the behaviour for other
waves or fields, the accuracy of one replicate (each replicate's own corrected
reading has a standard deviation of about 0.1, twice the tolerance), or the
cost of a single-sprinkle implementation.
