# P12: does the instrument read CURVED geometry? — power-first preregistration

Status: **STAGE A HAS RUN and returned IMPROVES (Section 9). STAGE B
IS FROZEN BY THE SECTION 10 ADDENDUM (2026-07-31) and has NOT run.**
Sections 1 through 8 are design v1.1 as frozen before Stage A; Section
9 is Stage A's record and is not re-read; Section 10 is the Stage B
addendum that Section 5 declared and deferred. Dates are local
(UTC+9); commit timestamps carry their +09:00 offset.

*(The header previously read "Nothing below has been run" while
Section 9 carried a completed campaign — the narrative-versus-artifact
inconsistency this programme hunts, committed against its own record.
Corrected here with the addendum rather than silently.)*

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

### 9.2 Stage A-12 (2026-07-28): gate **IMPROVES**; the campaign's
apparent plateau was small-sample fluctuation

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

Three of these depart together and in the same direction: the error
falls from the bottom rung to the middle and then stops, the slope
excludes the predicted `-1/3`, and the contrast is half the derived
size.

**First reading, WITHDRAWN on review (2026-07-28).** This record
originally called that pattern the signature of a curvature bias
floor and wrote that "the gate's IMPROVES is measuring the approach
to that floor, not convergence to the truth." Two objections, both
correct, and the second decisive:

1. The three checks are not three observations. All three are
   functions of the SAME three rung means at `n_per_rung = 12` (the
   floor value of Section 1.2), so they are one correlated
   observation, and the entire signature rests on a single 12-sample
   rung. The middle-to-top reversal is under one rung standard error.
2. A non-decaying floor is not available theoretically: at fixed
   geometry the chain law's convergence `L / sqrt(2 rho) -> tau_geo`
   holds in curved spacetime too (Section 7), so curvature can inflate
   a PRE-ASYMPTOTIC plateau but cannot stop convergence. The
   withdrawn sentence asserted a mechanism the data could not
   distinguish.

**The adjudicating diagnostic (post hoc, quarantined, labelled).**
Run at 200 samples per rung in the design-check seed space
(`2000000+`, disjoint from every experimental window), with an extra
`N = 4800` rung, through the runner's own code path:

| `N` | mean `y` | median relative error | mean `m` |
|---|---|---|---|
| 600 | -0.5673 +/- 0.0084 | 0.278 | 18.7 |
| 1200 | -0.6578 +/- 0.0094 | 0.228 | 37.5 |
| 2400 | -0.7605 +/- 0.0084 | 0.179 | 76.2 |
| 4800 | -0.8386 +/- 0.0077 | 0.149 | 152.0 |

**Monotone across all four rungs; slope `-0.304` over the four, and
`-0.321` over the campaign's three — both consistent with `-1/3`;
`Delta(600 -> 2400) = -0.193`, i.e. the derived `-0.2007`.** The
plateau was small-sample fluctuation at `n = 12`. No bias floor is
detected out to `m ~ 152` and `tau / ell ~ 0.3`.

**What the campaign therefore establishes**: over this ladder the
unchanged P11 chain estimator reads CURVED proper time on a genuinely
curved sprinkling, at 15-28% relative error, improving with density,
with the flat template far worse at every rung. The frozen gate says
IMPROVES and the diagnostic says the improvement runs at the
longest-chain rate. **What it does not establish**: anything about
`tau / ell` regimes other than ~0.3, and nothing about a floor —
neither its presence nor its absence beyond `m ~ 152`.

**The flat arm, rebuilt — and the frozen prediction vindicated.**
Review found that the arm as first implemented computed
`|tau_flat - tau_curved| / tau_curved`, which never touches
`tau_hat`: a property of the PATCH, not of the instrument, and
therefore incapable of supporting the claim it was cited for. The arm
is now the estimator's error against the flat template,
`|tau_hat - tau_flat| / tau_flat`, and the old quantity is retained
under an honest name as what it always was — how much curvature the
band carries at all.

| rung | error vs CURVED truth | error vs FLAT template | geometric gap |
|---|---|---|---|
| 600 | 0.282 | 0.464 | 0.401 |
| 1200 | 0.218 | 0.430 | 0.390 |
| 2400 | 0.230 | 0.442 | 0.381 |

The estimator sits about twice as close to the curved truth as to the
flat template at every rung — the claim the arm was cited for, now
resting on a quantity that involves the estimator.

**But the "flat and vindicated" reading of it is withdrawn (same
round, second finding).** Once the arm is `|tau_hat - tau_flat| /
tau_flat` it is DENSITY-DEPENDENT by definition: as `tau_hat`
converges to `tau_curved` the arm tends to
`|tau_curved - tau_flat| / tau_flat`, not to a constant. Here
`tau_flat = 1.39 tau_curved` on the band, so that limit is
`0.39 / 1.39 ~ 0.28`, and the arm must FALL toward it. The
quarantined 200-per-rung diagnostic confirms it does —
`0.489 / 0.457 / 0.410 / 0.391` across `N = 600 / 1200 / 2400 / 4800`,
monotone — so the apparent flatness of the campaign's three
12-sample medians was noise, and the agreement with the frozen band
`0.42-0.49` was coincidence of the ladder's middle, not confirmation.

The frozen prediction is therefore RETIRED for this quantity rather
than scored: Section 4's rationale ("its expected level is fixed by
the patch geometry alone") is true of the geometric gap and false of
the estimator-based arm, so the prediction described one quantity
while the claim needed the other. What the arm does support, and
more strongly than flatness would, is that the SEPARATION widens with
density: curved-truth error `0.278 -> 0.149` while the flat-template
error stays above `0.39`, so the reading tracks the curved metric
better and better while never approaching the flat one.

**Two further corrections from the same round**, both silent-drift
risks rather than wrong numbers: the estimator now CALLS
`estimate_tau_from_longest_chain_1p1` instead of re-deriving its
endpoint correction inline (a later recalibration would otherwise
have split P11 and P12 while this module claimed they shared an
instrument), and the truth regression now integrates the GEODESIC
equation — `x` is cyclic, so `k = Omega^2 dx/dtau` is conserved and
`Delta x` and `tau` follow by quadrature, shot on `k` — asserting
`1e-9` as Section 3 requires, where the first version integrated a
straight coordinate path at `2e-3`. The closed form passes at `1e-9`.
The gate is unchanged to the last digit by all three fixes
(`Delta = -0.104`, CI `[-0.177, -0.041]`), which is the evidence that
they were publication and hygiene defects rather than measurement
ones.

**The lever this hands over**, already named in Section 7: sweep
`tau / ell` at fixed density. That experiment no longer rests on a
measured floor — there is none to rest on — but on the sharper
question the diagnostic leaves open: how large `tau / ell` may grow
before the flat-normalized reading degrades, and whether the
exponent survives there.

Artifacts: `docs/prereg/frozen/p12/` (verification, pilot, Stage A),
stamp recorded in each.

---

## 10. Stage B addendum (frozen 2026-07-31, before any Stage B data)

Section 5 declared Stage B and deferred its estimator, listing five
things the addendum must do before any Stage B datum exists. This
section does all five, plus one thing Section 5 did not ask for and
that P13's record argues is necessary (10.7). Every number quoted here
is measured on quarantined design-check seeds and is **disclosed
pre-freeze**, per the v1.1 and P13 precedents, so that a later pass
cannot be mistaken for a surprise.

Design-check artifacts, both committed with this addendum and both
containing no experimental seeds:
`docs/prereg/frozen/p12/p12_stage_b_volume_check.json` (deterministic
quadrature, no sprinkling at all) and
`.../p12_stage_b_ensemble_check.json` (P12's design-check space
`2000000-5999999`, Section 6).

### 10.1 Item 1: the diamond volume relation, exactly

Section 5 refused the small-diamond expansion because its coefficient
is "exactly the kind of unverified constant that killed Stage C v1".
The exact relation is available in closed form. In the flat slicing
with `Omega^2 = ell^2/eta^2` and null coordinates `U = eta + x`,
`V = eta - x`, the diamond of a timelike pair is the coordinate
rectangle `[U_p, U_q] x [V_p, V_q]`, and

    Vol = 2 ell^2 int int dU dV / (U+V)^2
        = 2 ell^2 ln[ (U_q+V_p)(U_p+V_q) / ((U_p+V_p)(U_q+V_q)) ]
        = 2 ell^2 ln( (Z+1)/2 )
        = 4 ell^2 ln cosh( tau / (2 ell) ),

using `(eta_p+eta_q)^2 - (Delta x)^2 = 2 eta_p eta_q (Z+1)` with `Z`
the de Sitter invariant of Section 3, and `Z = cosh(tau/ell)`.
**Position drops out entirely**, as maximal symmetry requires — the
volume is a function of proper time alone.

Verified before freezing, per Section 3's rule:

| check | result |
|---|---|
| closed form vs brute-force quadrature of `Omega^2` over the diamond, 5 positions | worst relative error `4.1e-9` (the midpoint rule's own grid error) |
| same `tau` at three different patch positions | volume spread `1.5e-15` |
| inversion recovers `R` on exact inputs, four `tau/ell` values | relative error `<= 1.0e-12` |

The flat limit is `Vol -> tau^2/2 - tau^4/(48 ell^2)`, so Section 5's
`(tau^2/2)(1 - c R tau^2)` has **`c = 1/48`** in this convention. That
number is recorded to retire the ambiguity, and is not used: the
campaign uses the exact form.

### 10.2 Item 1 continued: the inversion, frozen

Write `g = Vol / tau^2`. With `s = tau / (2 ell)`,

    g = ln cosh(s) / s^2  =:  G(s),

strictly decreasing from `G(0+) = 1/2` (flat) toward 0, hence
invertible on `(0, 1/2)`. The instrument's per-pair observable is
order-plus-count only:

    g_hat = m_open / ( rho_hat * tau_hat^2 )

with `m_open` the open interval cardinality (order data), `rho_hat` the
Section 2 density calibration, and `tau_hat` the **unchanged** P11
chain estimator called through
`estimate_tau_from_longest_chain_1p1`. The frozen inversion is

    s_hat = G^-1( Q_hat / 2 ),      Q_hat = g_bar_curved / g_bar_flat
    R_hat = 8 s_hat^2 / tau_hat_bar^2

with the ratio `Q_hat` taken against the flat twin of 10.6 at matched
`m`. `R_hat = F(|I|, tau_hat, rho)` as Section 5 required, with the
twin supplying the normalization.

**Boundary convention, frozen.** The positive-curvature domain is
`Q_hat/2 < 1/2`. A measurement outside it has recovered no curvature,
so the boundary-constrained estimate is `R_hat = 0` and the relative
error is exactly 1. This is the edge of the parameter space, not a
clamp over noise, and the **boundary-hit fraction is published per
rung** — 28.2% / 27.0% / 14.8% in design check, so it is a large and
declining share and must never be silent.

### 10.3 Item 1 continued: why a RATIO and not a difference

Section 5 item 5 specified subtracting the flat arm's discrepancy. The
addendum freezes a **ratio** instead, because the bias it must remove
is multiplicative and the design check shows how large it is:

| rung (`m`) | `g_bar` curved | `g_bar` flat | continuum flat value |
|---|---|---|---|
| 19.0 | 0.97247 | 1.04001 | 0.5 |
| 38.1 | 0.80787 | 0.86685 | 0.5 |
| 75.4 | 0.70700 | 0.77497 | 0.5 |

Each arm's `g_bar` is biased by **54% to 108%** against its own
continuum value (the curved arm's is `G(0.75) = 0.45914`, the flat
arm's is 0.5) — entirely because
`tau_hat` underestimates `tau` (`tau_hat/tau` = 0.705 / 0.760 / 0.803),
and `tau_hat` enters squared. A difference would leave that bias in the
units of the answer. The ratio cancels it, and what survives is the
signal:

| rung | `1 - Q_hat` | continuum truth |
|---|---|---|
| 19.0 | +0.06494 | 0.08172 |
| 38.1 | +0.06804 | 0.08172 |
| 75.4 | +0.08771 | 0.08172 |

A 55%-biased pair of arms yielding an 8.8% ratio signal against an
8.2% truth is the whole case for the ratio, and it is measured rather
than argued.

### 10.4 Item 2: `Delta*_B`, derived and NOT borrowed

`g_hat` inherits two independent dispersions whose EXPONENTS are
theory: a Poisson `1/m` from the count, and a BDJ `m^(-2/3)` from the
chain (`sd(L) ~ m^(1/6)` against `E L ~ 2 sqrt(m)`, doubled by the
square). So

    rel_sd(g)^2 = a/m + b m^(-2/3).

The two COEFFICIENTS are calibrated on design-check seeds rather than
asserted, because `m` and `L` rise together inside a box and their
ratio therefore disperses LESS than independence would give. Fitted:
**`a = 1.6019`, `b = 0.4338`**, against an uncorrelated analytic
`b = 0.8132` — the shortfall is that correlation, named rather than
absorbed. The fit reproduces the measurements it was built from
(measured 0.3825 / 0.2755 / 0.2225 versus model 0.3809 / 0.2836 /
0.2134).

Propagating over the frozen ladder (`m` 19.02 -> 75.42):

    Delta*_B = log10[ rel_sd(g; m_top) / rel_sd(g; m_bot) ]
             = -0.2516 dex.

**Steeper than `-0.2007`, and that is the evidence it was not
borrowed**: the count contributes a Poisson component alongside the
chain's, so the ratio's dispersion decays faster than the chain's
alone. The equivalence margin follows P11's convention,
`delta_eq_B = |Delta*_B| / 3 = 0.0839` dex.

**A note against a plausible future error.** The `1.645` in the
equivalence sample-size slot is `Phi^-1(0.95)`, i.e. **already the
two-sided quantile** for `beta = 0.10`, because ROBUST-style
equivalence needs BOTH interval bounds inside the margin and so has
power `2 Phi(z) - 1`. P13's review found exactly this slot "upgraded"
with a one-sided quantile; the note exists so it is not upgraded again.

### 10.5 Item 3: six disjoint pairs at the bottom rung

Measured on design-check seeds, 400 samples per rung, both arms:

| rung | mean `m` | completion, curved | completion, twin |
|---|---|---|---|
| 600 | 19.0 | **1.000** | 1.000 |
| 1200 | 38.1 | 1.000 | 0.998 |
| 2400 | 75.4 | 1.000 | 1.000 |

The bottom rung — the one Section 5 item 3 singled out, where the pool
is thinnest — completes 400 of 400. The packing budget is not the
binding constraint at this operating point, which is a consequence of
10.8's move: P13's `tau/ell = 1.5` patch is large (`X = 17`) precisely
because six boxes had to stack along `x` there.

### 10.6 Item 5: the flat twin, promoted and matched

Section 5 item 5 promotes the Stage A flat arm from a labelled check to
Stage B's load-bearing control, and here it is inside the estimator
rather than beside it: `g_bar_flat` is the denominator of `Q_hat`. Its
construction is P13's twin verbatim — same coordinate rectangle, a
uniform intensity calibrated so realized mean `m` matches, a band on
`tau_flat = sqrt(dU dV)` centred on the curved rung's mean box side,
same eligibility, `K`, rejection cap and fill rule.

Two requirements follow, both frozen:

1. **The `m`-matching gate applies to BOTH arms**, each against its own
   grand mean, at `+/- 5%`. This is P13's review C5 correction carried
   forward before it can recur: a twin rung drifting in discreteness
   would break the bias cancellation the ratio depends on, in the arm
   the estimator divides by. The cross-arm `m` LEVEL offset is
   published as a labelled field, not gated.
2. **`g_bar_flat` against 0.5 is published per rung** as a labelled
   check. Its ratio to 0.5 IS the measured `(tau/tau_hat)^2` bias
   (2.08 / 1.73 / 1.55 in design check), so the cancellation the
   estimator relies on is auditable from the artifact rather than
   asserted here.

`n_twin = n_per_rung`. The twin's per-pair dispersion is measured
equal to the curved arm's (0.381 / 0.295 / 0.220 against 0.383 / 0.276
/ 0.222), so equal `n` equalizes the two contributions to `se(Q_hat)`.
The twin is NOT separately power-sized as P13's was, and the reason is
that P13's twin supplied a GATED CONTRAST while this one supplies a
NORMALIZER — its precision enters the estimator, where equal `n` is
the right rule, not a verdict, where equivalence sizing would be.

### 10.7 Item 4: the budget, and the co-requirement it forces

Section 5 item 4 predicted the failure mode in advance: an uncorrected
chain bias could leave `|R_hat - R| / R >> 1` at every rung while
`Delta_B` still improves from bias decay alone — "a gate that passes
while recovering nothing". The design check confirms the first half.

| rung | `m` | per-sample median `\|R_hat-R\|/R` | boundary hits | pooled `\|R_hat-R\|/R` |
|---|---|---|---|---|
| 600 | 19.0 | 2.56 | 28.2% | 0.554 |
| 1200 | 38.1 | 1.56 | 27.0% | 0.406 |
| 2400 | 75.4 | 1.10 | 14.8% | 0.683 |

**No single sample recovers `R`**, at any rung: the per-sample median
relative error exceeds 1 everywhere. And the pooled dimensionful
recovery is **not monotone** — 0.554 / 0.406 / 0.683 — because it is a
product of two errors with opposite signs. Splitting it shows why:

| rung | `R tau^2` recovered | truth | relative error | `tau_hat/tau` |
|---|---|---|---|---|
| 600 | 3.4725 | 4.5 | 0.228 | 0.705 |
| 1200 | 3.6580 | 4.5 | 0.187 | 0.760 |
| 2400 | 4.8816 | 4.5 | **0.085** | 0.803 |

`R tau^2 = 8 s^2` needs **no length at all** — it comes from `Q_hat`
alone — and it recovers monotonically to 8.5%. Dividing by
`tau_hat_bar^2` to reach `R` in continuum units then multiplies the
error by `(tau/tau_hat)^2`, which is 1.55 at the top rung, and the
partial cancellation at the middle rung is arithmetic coincidence
rather than trend.

**So the addendum adds a co-requirement Section 5 did not ask for.**
The budget alone makes the failure visible; it does not stop the
verdict word from being read as recovery, which is precisely what P13's
Section 11 record cost. Frozen:

> **A pass requires BOTH.** (i) the Section 1.4 table returns IMPROVES
> on `Delta_B` for `y_B = log10(|R_hat - R| / R)`, and (ii) the top
> rung's DIMENSIONLESS recovery satisfies
> `|R tau^2 hat - R tau^2| / (R tau^2) <= 0.25`, reported with its
> bootstrap CI. If (i) passes and (ii) fails, the record reads
> **RATE-ONLY** and states in the same sentence that no curvature was
> recovered — the gate measured a decaying bias and nothing else.

`0.25` is a design choice, stated as one. Design check measured 0.085
at the top rung, so (ii) is expected to pass and is disclosed here so
that its passing is not read as a discovery. The threshold is loose
enough to be purchasable at the frozen `n` and tight enough to exclude
"recovering nothing", which sits near 1.

### 10.8 The operating point moves to `tau/ell = 1.5`

Section 5 was written when Stage A's band was the only certified one,
so it assumed `tau/ell ~ 0.3`. The volume check measured what that
costs. The inversion's amplification — relative error in `R_hat` per
relative error in `g` — is

| `tau/ell` | `1 - 2G` (the signal) | amplification |
|---|---|---|
| 0.30 | 0.00373 | **268.8x** |
| 0.60 | 0.01465 | 68.9x |
| 1.00 | 0.03908 | 26.2x |
| 1.50 | 0.08172 | **12.9x** |

At `tau/ell = 0.3` the signal is 0.37% and every input error is
amplified 269-fold, which is Section 5 item 4's "two to three orders of
magnitude" made exact. **P13 has since certified the same estimator out
to `tau/ell = 1.5`** (CURVATURE-ROBUST, its Section 13, with the
caveats that record carries), so Stage B takes that rung: the signal is
22x larger and the amplification 21x smaller, a 480-fold reduction in
the sample count a given recovery costs.

Frozen constants, P13's 1.50 rung verbatim: `ell = 1`,
`eta in [-7.0, -1.0]`, `|x| <= 17.0`, band `tau_curved in [1.35, 1.65]`
(`+/- 10%`), twin band centre `4.2262`. The density ladder is fixed
geometry with `rho` scaling 4x end to end as Stage A's did:
`rho in {19.025, 38.05, 76.1}` giving `m ~ 19 / 38 / 75`, and the twin
`rho_twin in {2.1277, 4.2554, 8.5108}` (P13's `8.5108` at the top,
halved and quartered down the ladder).

**Scope this moves, stated plainly.** Stage B no longer measures at
Stage A's band, so it does not inherit Stage A's IMPROVES as evidence
about its own operating point — it inherits P13's. Both are recorded
gates with reachable stamps, and the cross-stage gate of 10.10 requires
both.

### 10.9 What Stage B can and cannot claim

A pass under 10.7 supports: **the scalar curvature of the ambient
geometry is recoverable from order plus the count, over this ladder, to
the accuracy the record reports.** That is the strong form of "order +
number = geometry" for a curvature invariant, and it is the first
quantity in this programme that is not a length or a time.

It does NOT support the general claim that curvature is recoverable
from causal order in general. P13's Section 13.3 gives the reason and
it applies here verbatim: **1+1D metrics are conformally flat**, so the
curvature reaches the order only through the conformal factor — that
is, through exactly the volume this estimator reads. Stage B measures
that channel accurately; it does not show there is no other channel,
and `d >= 3` is where that question lives. The record must carry this
sentence next to any pass.

### 10.10 Gates, power, and seed windows

**Cross-stage gate.** Stage B runs only after P12 Stage A's frozen
IMPROVES and P13 Stage A-13C's frozen CURVATURE-ROBUST, both with
reachable stamps — the first because the estimator is that one, the
second because the operating point is that one.

**Verification-B.** 2000 single-stream samples per rung, completeness
pin `>= 1998`, measured wall times for the feasibility projection,
quoted to one significant figure. Must pass before Stage P-B may run.

**Stage P-B.** Both endpoint rungs and both twin endpoint rungs, 200
samples each, cross-rung statistics forbidden, calibrated Bonett bounds
published as consumed. Then P11's frozen formulas with `Delta*_B` and
`delta_eq_B` of 10.4:

    n_sup = ceil( S^2_90 * (1.960 + 1.282)^2 / Delta*_B^2 )
    n_eq  = ceil( S^2_90 * (1.960 + 1.645)^2 / delta_eq_B^2 )
    n_per_rung = clamp( max(n_sup, n_eq), 12, 1300 )

with the FLAT-WITHIN-MARGIN row available only if `n_eq <= 1300`.
Indicative, uncalibrated, from design check: `S^2 = 0.4606`,
`n_sup = 77`, `n_eq = 852` — so **both verdicts look purchasable at
the cap**, which is the affordability statement P13's Section 1.3
requires a design to make about itself before running. The pilot's
calibrated bounds decide.

**Seed windows (frozen). Section 6's `1200000-1999999` reservation is
SUPERSEDED**, for two reasons stated rather than worked around: six
blocks at this addendum's sample sizes need about 1.6M seeds and that
reservation holds 800k, and `2000000+` is P12's own design-check space
which this addendum's checks have now consumed. Stage B takes a fresh
decade above every range the programme has used (P13 v3 ends at
21503999; P13's Section 12.4 reserves `22000000+` for design checks,
left intact):

| block | base | slots | span |
|---|---|---|---|
| verification-B, 600 | 30000000 | 2000, consecutive | 30000000-30001999 |
| verification-B, 1200 | 30002000 | 2000, consecutive | 30002000-30003999 |
| verification-B, 2400 | 30004000 | 2000, consecutive | 30004000-30005999 |
| pilot, curved 600 | 30100000 | 320 | 30100000-30163999 |
| pilot, curved 2400 | 30200000 | 320 | 30200000-30263999 |
| pilot, twin 600 | 30300000 | 320 | 30300000-30363999 |
| pilot, twin 2400 | 30400000 | 320 | 30400000-30463999 |
| Stage B, curved 600 | 31000000 | 1320 | 31000000-31263999 |
| Stage B, curved 1200 | 31400000 | 1320 | 31400000-31663999 |
| Stage B, curved 2400 | 31800000 | 1320 | 31800000-32063999 |
| Stage B, twin 600 | 32200000 | 1320 | 32200000-32463999 |
| Stage B, twin 1200 | 32600000 | 1320 | 32600000-32863999 |
| Stage B, twin 2400 | 33000000 | 1320 | 33000000-33263999 |

Slots are `cap + 20` at stride 200. Design-check space from here on is
`40000000+`. Everything else procedural is P12 Section 6's, verbatim:
first-`n`-complete fill with published skip identities, skip cap 20,
clean-worktree preflight, stamp equality within a chain, and
one-significant-figure quoting of every measured timing.

### 10.11 Implementation plan (after this addendum merges)

Extend `experiments/positive_control/p12_curved.py` with a Stage B
module that imports rather than reimplements: P13's per-rung sprinkler
and eligibility for the 1.50-rung geometry, the shared chain estimator,
P11's power and verdict machinery. The only new code is `G`, its
bisection inverse with the boundary convention, `g_hat`, and the
twin-ratio aggregation.

Tests before any run, each written against the case it is meant to
catch:

1. the closed-form volume against quadrature at the frozen tolerance,
   and its position independence at fixed `tau`;
2. `G` strictly decreasing on the tested range, and `G^-1(G(s)) = s`
   to `1e-10`;
3. the boundary convention: `Q_hat/2 >= 1/2` yields `R_hat = 0` and
   relative error exactly 1, and the boundary-hit count is reported;
4. the `m`-gate covering BOTH arms, failing on a drifting twin beside a
   clean curved arm (P13's C5 regression, ported);
5. `Delta*_B` recomputed from the frozen `a`, `b` and the frozen `m`
   ladder, asserted equal to the `-0.2516` this section freezes, so the
   number cannot drift from its derivation;
6. the equivalence coefficient's two-sided identity, `2 Phi(z) - 1 =
   0.90` at the frozen constant (P13's C2 regression, ported);
7. window privacy and freshness against every documented spent range,
   asserted against the ranges themselves rather than against a floor;
8. the co-requirement of 10.7 returning RATE-ONLY when (i) passes and
   (ii) fails.

Stage P-B runs only from a clean commit containing all of it, and only
after Verification-B's pin passes.
