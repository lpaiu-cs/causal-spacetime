# P12: does the instrument read CURVED geometry? — power-first preregistration

Status: **DESIGN v1.1 (2026-07-28), for in-session review. Nothing
below has been run.** Dates are local (UTC+9); commit timestamps
carry their +09:00 offset.

v1.1 corrections, all pre-freeze and pre-data, from design review
(each independently reproduced before being applied):
(i) **the patch widens to `|x| <= 1.5` and box-inside-patch joins
eligibility** — under v1.0's constants the frozen protocol completed
only 12% / 22% of samples against a pin demanding 99.9%, because six
disjoint boxes of area ~0.18 do not fit a patch of UV-area 2; with
the wider patch and the eligibility rule, completion measured
100% / 100%. The geometric budget did not transfer from P11 along
with the machinery, and the packing demonstration this design
demands of its own Stage B addendum was the thing Stage A itself
needed. (ii) **the causal box must lie inside the patch**: 37% of
band pairs had boxes leaving the sprinkled region, where there are no
events to be found, so their chains were truncated by an
`N`-independent geometric artifact -- mostly cancelling in `Delta`
but corrupting the constant, the flat arm, and (fatally) Stage B's
volume counts. P11's unit diamond is UV-convex so this was free
there; this patch is not. (iii) the null convention is pinned and the
flat arm's formula corrected (a spurious factor 2), with the
flat arm's expected level frozen as a prediction and the
convention identity asserted by the degeneracy test *(applied in
v1.1a: the v1.1 commit claimed this correction in its log while
the body still carried the old formula -- a patch whose match
silently failed. The narrative-versus-artifact defect class this
programme hunts, committed against its own correction log; every
substitution is now assert-guarded)*. (iv) the Stage B
addendum gains two requirements against a second-order signal being
swamped by an amplified bias. (v) wording and cross-reference pins.

Lineage. P11 established, over its tested ladder, that continuum-unit
proper time, spacelike distance, and coordinates are reconstructible
from order-plus-count with accuracy that improves as density grows.
Every one of those results is in FLAT Minkowski. P12 asks the first
question beyond flatness: **does the same instrument read the
geometry when the geometry is curved — and can the curvature itself
be recovered from order plus count?**

Why this is the right next lever, and why it is cheap: in 1+1D every
metric is conformally flat, so a curved sprinkling has the SAME
causal-order structure as a flat one read in null coordinates. What
changes is the sprinkling DENSITY (the volume element). Curvature
therefore enters exactly where the programme's thesis says geometry
lives — in the counting measure — and the P11 machinery (estimators,
gates, power design, certificates) transfers without redesign. If the
instrument reads curvature, "order + number = geometry" holds beyond
the flat case; if it cannot, the boundary is located precisely.

---

## 1. Power design (first, by house rule)

### 1.1 Ladder and primary contrast

Ladder: expected sprinkling count `N in {600, 1200, 2400}` over a
FIXED curved patch, so density scales by 4 end to end (as in P11) and
the curvature radius is held constant — the discreteness scale moves,
the geometry does not.

Per-sample statistic, Stage A:
`y = log10( median over K pairs of |tau_hat - tau_curved| /
tau_curved )`, with `tau_curved` the exact geodesic proper time of
Section 3. Rung statistic: the MEAN of `y` (the CLT-exact aggregator
of P11 1.2). Primary contrast
`Delta = mean_y(2400) - mean_y(600)`.

### 1.2 Effect size, sample size, verdict table

The chain estimator's error is governed by the same longest-chain law
(Section 7), so at fixed continuum separation

    Delta*_A = -(1/3) log10(4) = -0.2007 dex,

derived, identical to P11's — the ladder and the estimator are the
same; only the ambient geometry changes. Sample size, verdict table,
equivalence margin, and the ex-ante flat-availability declaration are
P11 Sections 1.2 and 1.4 **applied verbatim**, including the
calibrated Bonett variance bound and its published-vs-consumed rule
(P11's revision list, item (1) of the deferred-P2 batch: pilots
publish the calibrated bound they feed to power -- 9.6 is the
separate wall-time rule). Stage P-12 pilots BOTH endpoint rungs, 200 samples each,
cross-rung statistics forbidden.

### 1.3 The outcome this design must be able to report

A curvature bias is physically possible and must be representable:
the chain law is asymptotic in `interval size / curvature radius`, so
a flat-normalized estimator can carry a bias that does NOT shrink
with density. That would show as FLAT-WITHIN-MARGIN or as a
positive-`Delta` DEGRADES, and it is an informative finding, not a
failure — it would locate the boundary of the instrument's reach
exactly where theory says the boundary is. The frozen band of
Section 4 keeps `tau / ell ~ 0.3`: curvature is unmistakably present,
and the law is still expected to apply.

Verdict-word caveat carried from P11 Stage C, where `n_eq = 168`
exceeded the cap: if the pilot declares the flat verdict unavailable
ex ante, a bias floor reads UNRESOLVED rather than FLAT-WITHIN-
MARGIN. The information then lives in the labelled per-rung medians
and the slope, not in the verdict word, and this record says so in
advance so that an UNRESOLVED is not later mistaken for an
uninformative campaign.

---

## 2. The ensemble: a TRUE Poisson sprinkling in 2D de Sitter

Patch: the flat slicing of `dS_2`,

    ds^2 = (ell^2 / eta^2) ( -d eta^2 + dx^2 ),   eta in [eta_0, eta_1],

with Ricci scalar `R = 2 / ell^2` (constant, positive). Null
coordinates `U = eta + x`, `V = eta - x` give
`ds^2 = -(4 ell^2 / (U+V)^2) dU dV`: the conformal factor is the only
thing that differs from Minkowski.

- **Causal order**: `p < q` iff `U_p < U_q` and `V_p < V_q` — the
  SAME rule as P11, because the conformal factor cannot change the
  causal structure in 2D. Every P11 order-side routine transfers
  unchanged.
- **Sprinkling**: a genuine inhomogeneous Poisson process of
  intensity `rho * sqrt(-g)`, i.e. proportional to
  `Omega^2 = 4 ell^2 / (U+V)^2`, realized by thinning a uniform
  proposal on the coordinate patch. This is a REAL Poisson
  sprinkling, not the conditioned rank grid of P11 Section 3 — so
  P12 also retires that ensemble caveat rather than inheriting it.
- **Frozen patch constants**: `ell = 1`, `eta in [-2.0, -1.0]`,
  `|x| <= 1.5` (v1.1: widened from 0.5 so six disjoint supports
  fit; the cost is stated rather than hidden -- the bottom rung's
  interval count drops to `m ~ 18`, comparable to P11's `m ~ 24`,
  so the ladder stays in the same pre-asymptotic regime). Proper volume of the patch is finite and computed
  exactly by the frozen quadrature of Section 3; `rho` is set so the
  expected count equals the rung's `N`.
- **Density calibration**: as in P11, the instrument is given the
  count; explicitly, `rho = (realized event count) / (frozen patch
  proper volume)`, the volume computed once by the Section 3
  quadrature and frozen with the constants. The
  claim is reconstruction from order plus density calibration — and
  Section 5 is explicit that recovering `R` from that pair is the
  strong form of "order + number = geometry", not a weaker one.

## 3. Truth (exact, frozen, and cross-checked before any run)

For two points in the flat slicing, the de Sitter invariant

    Z = [ eta_p^2 + eta_q^2 - (x_p - x_q)^2 ] / ( 2 eta_p eta_q )

gives the timelike geodesic proper time

    tau_curved = ell * arccosh(Z),    Z >= 1.

(Short-separation check, done at design time:
`Z = 1 + (d eta^2 - d x^2) / (2 eta^2) + ...` reproduces
`tau = (ell/|eta|) sqrt(d eta^2 - d x^2)`, the local metric.) Before
Stage P-12 may run, a regression must confirm this closed form
against numerical integration of the geodesic equation to `1e-9`
relative on a frozen grid of separations — the same "verify the
assumption before freezing" rule that killed Stage C v1.

## 4. Stage A: does the chain estimator read curved proper time?

- Estimator: **unchanged from P11 Stage A** —
  `tau_hat = (L - 2) / sqrt(2 rho)` through the shared definition
  `estimate_tau_from_longest_chain_1p1`, with `L` the closed-interval
  longest chain from the shared `longest_chain_length`. Nothing about
  it knows the geometry is curved; that is the point.
- Pairs: `K = 6` per sample from the precomputed pool of related
  pairs satisfying BOTH frozen eligibility conditions: `tau_curved`
  in the band `[0.28, 0.34]` (so `tau / ell ~ 0.3`), and **the
  pair's causal box lies inside the patch** (v1.1 -- a geometric
  condition on coordinates, decided before any measurement, so
  completeness still never conditions on a measured value). Drawn
  without replacement with the P11 greedy disjoint-support rule
  (supports are the pairs' `(U, V)` boxes), 200-draw cap on overlap
  rejections, first-`n`-complete fill, skip cap 20.
- Gate: the 1.4 verdict table on `Delta`. Labelled, never gating:
  the slope of mean `y` on `log10 N` against `-1/3`; the median
  relative error per rung; and **the flat-comparison arm** — the
  same samples scored against the FLAT proper time
  `sqrt(dU dV)` a naive analyst would use, whose error should be
  large and NOT improving, quantifying how much of the result is
  curvature actually being read. Its expected level is fixed by the
  patch geometry alone, so it is FROZEN here as a prediction rather
  than read off afterwards: **relative error `0.42` to `0.49`**
  (design-check measurement at the v1.1 constants), large and flat
  in `N`.

**Null convention, pinned (v1.1).** This document uses FULL-null
coordinates `U = eta + x`, `V = eta - x`, in which the flat proper
time is `sqrt(dU dV)`. P11's shared routines use the HALF-null
convention `u = U/2`, `v = V/2`, in which the same number reads
`2 sqrt(du dv)`. Every call into a P11 routine therefore passes
`U/2, V/2`, and the `ell -> infinity` degeneracy test of Section 8
asserts that identity numerically. v1.0 mixed the two (it quoted the
half-null formula against full-null coordinates), which would have
inflated the flat arm by a factor 2 and diluted the curvature
diagnostic; the assertion exists so the clash cannot recur
silently.

## 5. Stage B: recover the curvature itself (declared; frozen by addendum)

Declared now, estimator frozen by addendum after Stage A reports and
before any Stage B data — the staged-freezing precedent of P11
Sections 10 and 11.

The idea: the chain estimator measures proper time; the interval
COUNT measures volume; and in a curved spacetime the volume of a
causal diamond of given proper time differs from the flat value by a
curvature term. Comparing the two therefore yields `R` without any
metric input beyond the density calibration. In 2D the small-diamond
expansion has the shape `V(tau) = (tau^2/2)(1 - c R tau^2 + ...)`,
but the coefficient convention is exactly the kind of unverified
constant that killed Stage C v1, so the addendum must:

1. compute the dS_2 diamond volume relation EXACTLY (numerically, to
   a pinned tolerance) rather than rely on an expansion, and freeze
   the inversion `R_hat = F(|I|, tau_hat, rho)` it implies;
2. derive `Delta*_B` from the propagated rates of its inputs, not
   borrow `-0.2007`;
3. demonstrate on quarantined design seeds that six disjoint pairs
   per sample are realizable at the bottom rung;
4. state the expected `|R_hat - R| / R` per rung BEFORE running
   (v1.1): curvature enters the diamond volume only at relative
   order `c R tau^2 ~ 0.002-0.008` here, and the inversion amplifies
   that by two to three orders of magnitude, so an uncorrected chain
   bias (`~0.89 m^(-1/3)`, i.e. 19-34% on this ladder) can leave
   `|R_hat - R| / R >> 1` at every rung while `Delta_B` still reads
   `~ -0.2007` from bias decay alone -- a gate that passes while
   recovering nothing. The budget is what makes that outcome
   visible;
5. use the **flat-twin difference** as the primary construction
   (v1.1): subtract the same chain-versus-count discrepancy measured
   on a FLAT ensemble matched in band and in `m`, so the universal
   bias constant cancels at first order and the curvature signal is
   what remains. The Stage A flat arm is thus promoted from a
   labelled check to Stage B's load-bearing control.

Gate: the 1.4 table applied to `Delta_B` for
`y_B = log10( |R_hat - R| / R )`. A pass means the instrument
recovers the SCALAR CURVATURE from order plus count, with accuracy
improving as density grows.

## 6. Protocol and seed windows

Everything procedural is P11's, verbatim: verification pin
(`>= 1998 of 2000` complete per rung, plus measured wall times),
first-`n`-complete fill with reserve slots and published skip
identities, clean-worktree preflight, stamp-equality within a chain,
cross-stage gate on P11 Stage A's frozen IMPROVES (P12's instrument
IS that estimator, so its flat certification is the prerequisite),
and one-significant-figure quoting of every measured timing.

| block | base | fills | window span |
|---|---|---|---|
| verification-12 (non-experimental) | 1000000 | 2000 per rung | 1000000-1005999 |
| Stage P-12 pilot, N=600 | 1020000 | 200 of 220 | 1020000-1063999 |
| Stage P-12 pilot, N=2400 | 1064000 | 200 of 220 | 1064000-1107999 |
| Stage A-12, N=600 | 1120000 | <= 60 of 80 | 1120000-1135999 |
| Stage A-12, N=1200 | 1136000 | <= 60 of 80 | 1136000-1151999 |
| Stage A-12, N=2400 | 1152000 | <= 60 of 80 | 1152000-1167999 |
| Stage B-12 blocks | 1200000-1999999 | — | frozen with the Stage B addendum |

All spans sit above every range the programme has used (documented
maxima: P11's experimental windows to 705999, its design checks at
900000+). **Design-check space for P12 is `2000000-5999999`**
(v1.1), disjoint from every experimental window above; the v1.1
completeness and truncation measurements used seeds in `2000000+`
and produced no experimental data.

## 7. Theory appendix: why the chain law should transfer, and where it stops

Myrheim's proposal and the Brightwell-Gregory theorem read the
longest chain between two events of a sprinkling as the geodesic
proper time, with the flat-space constant, in the regime where the
interval is small compared with the curvature radius: the interval
looks locally Minkowskian, the chain statistic is dominated by the
local structure, and the leading behaviour `E L = 2 sqrt(m)` with
`m` the interval count carries over. The correction is controlled by
`tau / ell` and by the discreteness `m^(-1/3)`; this design holds the
first fixed at about `0.3` while sweeping the second by a factor of
4, which is exactly the separation of effects the gate needs. The
same references bound where the transfer must fail: at
`tau / ell = O(1)` the diamond is not locally flat and no
flat-normalized reading is expected to survive — a lever for a later
ladder in `tau / ell` rather than in density.

## 8. Implementation plan (after this design is approved)

New module `experiments/positive_control/p12_curved.py`: the
inhomogeneous sprinkler (thinning, frozen patch constants), the exact
truth of Section 3, the pair pool with BOTH eligibility conditions,
and stage runners verify/pilot/A. Estimators come from the shared definitions used by
P11 — imported, never re-implemented. Tests before any run: the
truth formula against numerical geodesic integration; the sprinkler's
realized density profile against `Omega^2` (chi-square, frozen
tolerance); the causal-order rule against the flat P11 routine on a
degenerate patch (`ell -> infinity` must reproduce Minkowski,
including the assertion `sqrt(dU dV) == 2 sqrt(du dv)` at
`u = U/2, v = V/2`, which pins the Section 4 convention); pool
band membership; window privacy and freshness; verdict-table
precedence. Stage P-12 runs only from a clean commit containing all
of it, and only after the verification pin passes.

---

## 9. Stage records (results)

### 9.1 Verification-12 and Stage P-12 (2026-07-28)

- **Verification-12**: 2000/2000 complete at every rung against the
  1998 pin — the v1.1 completeness fix holds at pin resolution, where
  v1.0's constants completed 12-22%. Measured wall times ~0.009 /
  ~0.03 / ~0.1 s per sample (one significant figure, per P11 9.6).
- **Stage P-12**: 200 samples per endpoint rung, zero skips. The
  calibration fired on both rungs (coverage 0.948 / 0.865 at nominal
  -> `z = 1.693 / 2.522`; the top rung's 0.865 is the largest
  under-coverage the programme has recorded, and the calibration
  absorbing it is exactly why P11 v1.9 replaced the nominal bound).
  `n_sup = 12`, `n_eq = 127 > cap` — flat verdict declared
  unavailable ex ante, so this campaign's fourth verdict is
  UNRESOLVED (the caveat Section 1.3 wrote down in advance).
  `n_per_rung = 12`; projected Stage A ~2 s. FEASIBLE.

### 9.2 Stage A-12 (2026-07-28): gate **IMPROVES**, mechanism checks say **bias floor**

12 complete samples per rung, zero skips, no selection caveat.

    Delta = -0.104 dex,  95% CI [-0.177, -0.041]  ->  IMPROVES

The gate passes as frozen. **The labelled checks do not support
reading this as the theory rate**, and the record says so plainly:

| labelled check | reading | expectation |
|---|---|---|
| mean `y` by rung | -0.567 / -0.712 / -0.671 | monotone |
| middle rung between endpoints | **NO** — the top rung is worse than the middle | yes |
| median relative error | 0.282 / 0.218 / 0.230 | falling |
| slope vs `log10 N` | **-0.173, CI [-0.295, -0.060]** | contains `-1/3` |
| `Delta` against derived `Delta*_A` | -0.104, about HALF of -0.2007 | at or beyond |
| mean interval count `m` | 18.8 / 38.7 / 75.8 | (as designed) |

Three of these fail together and in the same direction: the error
falls from the bottom rung to the middle and then STOPS, the slope
excludes the predicted `-1/3`, and the contrast is half the derived
size. That is the signature Section 1.3 described in advance —
**a curvature bias floor**: the discreteness error keeps shrinking
with density, but the flat-normalized estimator carries a bias set by
`tau / ell`, not by `m`, and once discreteness falls below it the
total stops improving. The gate's IMPROVES is real but is measuring
the approach to that floor, not convergence to the truth.

**Prediction miss, recorded rather than absorbed.** The flat-comparison
arm was frozen at `0.42-0.49` (Section 4, from the design check); it
read **0.401 / 0.390 / 0.381** — outside the frozen band and slightly
FALLING where the prediction said flat. The likely cause is post hoc
and labelled as such: the design check measured all band-eligible
pairs, while the campaign scores only pairs that survive disjoint
packing, which is not a uniform sample of the band. The arm still does
its job — the flat reading is far worse than the curved one at every
rung, so curvature is unmistakably being read — but the frozen number
was wrong, and freezing it is what made that visible.

**What this establishes, and what it does not.** Establishes: over
this ladder, the P11 chain estimator applied to a genuinely curved
sprinkling reads CURVED proper time to ~20-28% and improves with
density; the flat reading is much worse at every rung, so the
instrument is reading geometry rather than a flat template. Does not
establish: that the improvement follows the longest-chain rate here —
the slope excludes `-1/3` and the sequence is non-monotone at the top
rung, both consistent with a bias floor at `tau / ell ~ 0.3`. Claims
bind to this estimator, this patch, and this band.

**The lever this hands over**, already named in Section 7: the next
ladder should sweep `tau / ell` at fixed density rather than density
at fixed `tau / ell`. If the floor is curvature bias, the error at
fixed `m` must grow with `tau / ell`, and the exponent should return
as `tau / ell` shrinks. Stage B (curvature recovery) is unaffected in
design — its flat-twin differencing (Section 5, requirement 5) is
aimed precisely at a universal bias of this kind — but it now has a
measured reason to exist rather than a hypothetical one.

Artifacts: `docs/prereg/frozen/p12/` (verification, pilot, Stage A),
stamp recorded in each.
