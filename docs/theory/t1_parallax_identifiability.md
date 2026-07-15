# T1: Parallax identifiability and stability of bracket-width echo profiles

Status: **THEORY DRAFT v0.3 — statements and proof programs; nothing frozen**
(v0.3 2026-07-16 KST: G1 closed, Theorem 1 upgraded to `[PROVED]`;
v0.2 2026-07-16 KST / 2026-07-15 UTC; v0.1 2026-07-15 — revised after PR
reviews, see Revision notes).

Each claim below carries a proof-status tag:

- `[PROVED]` — argument written out here and believed complete.
- `[PROVABLE]` — route is standard; details to be written and checked.
- `[CONJECTURED]` — plausible, not yet proved; must not be cited as a result.

Nothing in this document may be cited in a manuscript above its tag.

## 1. Goal

State and prove conditions under which the order-level bracket-width echo
profiles of PC-V1 Section 5 determine the targets' spatial order up to
reflection (and overall scale), and a stability bound under Poisson
sprinkling with expected error decreasing in sprinkling density. This
complements the embedding-uniqueness theorem (uniqueness *given* an
embedding) with a finite, order-intrinsic reconstruction criterion.

## 2. Setup and definitions

Work in the 1+1D causal diamond `M = { (t,x) : |x| + |t| < 1 }` with the
Minkowski *causal* order (null-inclusive):

```
(t,x) < (t',x')   iff   t' - t > 0  and  (t'-t)^2 >= (x'-x)^2.
```

This is J (causal precedence), not the strict chronological order I —
deliberately: it is the convention `causal_matrix_minkowski()` implements
(`dt > 0 and dt^2 >= dx^2`), and Lemma 2's rank-gap identity is exact
*because* null-aligned ticks count as predecessor/successor. Under the
strict order a tick exactly on the light cone would be neither, and the
identity fails on endpoint-aligned configurations — which the
deterministic grid of Model D can realize (review supplied the
counterexample). Theory and instrument must share this convention. Fix:

- A Poisson sprinkling `S` of density `rho` on `M`.
- `R` observer chains `C_1..C_R` at fixed spatial coordinates `x0_r`.
  Ticks are the chain's elements, indexed by rank; rank is the only clock.
  Two *clock models* for how tick times arise:
  - **Model D (the instrument).** What PC-V1 actually constructs:
    `build_positive_control_scene()` calls
    `make_stationary_observer_chain_1p1()`, which appends an exact
    vertical worldline whose tick times are a deterministic uniform grid
    (`np.linspace`, `ticks_per_chain = 96` by default), *independent of
    the sprinkling and of rho*. Tick spacing is `delta = T/(K-1)` for `K`
    ticks over time span `T`; write `lambda = 1/delta` for the tick rate.
  - **Model P (stochastic idealization).** Tick times form a simple
    stationary point process of intensity `lambda` on the worldline
    (Poisson being the special case used for concentration bounds).
    No code constructs such chains today (Section 6, G2).
- Targets `e_1..e_n` in `S` at coordinates `(t_j, x_j)`.

Observable (PC-V1 Section 5, `find_radar_ticks_from_order`):

```
p[j,r] = max { i : tick_i of C_r precedes e_j }
s[j,r] = min { i : e_j precedes tick_i of C_r }
W[j,r] = s[j,r] - p[j,r]          (integer tick-rank units)
```

`W[j,r]` is defined only when both sides exist (`reachable[j,r]`).

Parallax centering (PC-V1 Section 6, amendment D1):

```
P[j,r] = W[j,r] - mean_{r'} W[j,r']
```

## 3. Continuum lemmas

### Lemma 1 (radar bracket width) `[PROVED]`

For an inertial observer at spatial coordinate `x0` and an event
`e = (t_e, x_e)`, with `dx = x_e - x0`, the past radar point is
`(t_e - |dx|, x0)` and the future radar point is `(t_e + |dx|, x0)`, so the
proper-time width of the continuum radar bracket is `2 |dx|`, independent
of `t_e`.

*Proof.* `(t, x0) < (t_e, x_e)` iff `t_e - t >= |dx|` and `t < t_e`
(null-inclusive order), so the set of preceding observer times is the
closed ray `t <= t_e - |dx|` and its supremum `t_e - |dx|` is attained;
symmetrically the succeeding times are `t >= t_e + |dx|`. For an
inertial worldline coordinate time is proper time. ∎

### Lemma 2 (rank calibration)

The engine is a deterministic rank-gap identity, not an overshoot
argument (v0.1 added two renewal overshoots and double-counted the
constant; review caught it).

**Identity `[PROVED]`.** Let the ticks of a chain be simple (no
coincident times), and let no tick be *spacetime-coincident* with the
target — automatic when `|dx| > 0`, since ticks live at `x0 != x_e`.
Under the null-inclusive order of Section 2, every tick time `tau` then
falls in exactly one of three classes relative to the target:
predecessor (`tau <= t_e - |dx|`), spacelike
(`t_e - |dx| < tau < t_e + |dx|`, the *open* radar interval), or
successor (`tau >= t_e + |dx|`) — a partition, with no endpoint
orphans. (The coincidence hypothesis is not decorative: at `|dx| = 0` a
tick exactly at `t_e` has `dt = 0` and is unrelated under `dt > 0` —
neither predecessor, spacelike, nor successor — and then `W = N + 2`,
residual exactly `+1`, outside the band below. Measure zero for
sprinkled targets, but realizable by construction; the verification
harness pins it and its band predicate rejects it.) Let `N` count the
spacelike ticks and assume at least one predecessor and one successor
exist (reachability). The predecessor of
maximal rank `k` is followed in rank order by the `N` spacelike ticks
(`k+1..k+N`) and then by the first successor at rank `k+N+1`, so

```
W = N + 1        (exactly, realization by realization).
```

The convention is load-bearing: under the *strict* (chronological)
order, a tick lying exactly on the light cone of the target is neither
predecessor nor successor, the three classes no longer partition, and
`W` can exceed `N + 1` — e.g. a target `(0, 0.5)` against uniform ticks
`-0.75, -0.5, ..., 0.75` gives `N = 3` but strict `W = 6`. The
instrument's `causal_matrix_minkowski()` is null-inclusive and yields
`W = 4 = N + 1` on that configuration, which is why the theory adopts
its convention rather than the reverse.

**Model D `[PROVED]`.** For the uniform grid with spacing `delta`, the
number of grid points in the *open* interval of length `L = 2|dx|` is
`L/delta + theta` with `theta in [-1, 1)` depending on the interval's
position relative to the grid (its *phase*; `theta = -1` exactly when
both endpoints are grid-aligned, the reviewer-example case). Hence

```
W = 2|dx|/delta + 1 + theta,      |W - (2 lambda |dx| + 1)| <= 1.
```

`W` is deterministic given the target, affine in `|dx|` up to a bounded
quantization term -- but `theta` depends on `t_e` and on `|dx|` through
the phase, so it is *observer-dependent* and does not cancel across
chains. This is not hypothetical: in the actual fixed-rank PC-V1 grid,
moving a target at `x = 0.10` from `t = 0` to `t = 0.003` shifts the
centered six-chain profile entries by about one third of a rank.

**Model P `[PROVED]`.** For any simple *stationary* tick point process
of intensity `lambda` on an (effectively) infinite worldline, Campbell's
formula gives `E[N] = 2 lambda |dx|` exactly, so with the identity above

```
E[ W ] = 2 lambda |dx| + 1        (exactly; no O(1), no moment condition).
```

Stationarity is what Model D lacks (its grid has a fixed phase), which
is precisely where the `theta` term comes from; averaging `theta` over a
uniform random phase recovers the stationary answer. For the *Poisson*
case additionally `N ~ Poisson(2 lambda |dx|)`, which is the variance
input Theorem 2 needs -- stationarity alone fixes the mean, not the
fluctuations.

Boundary caveat (both models): near the diamond boundary the bracket
interval may leave the chain's tick support; those targets fail
reachability and are excluded, matching the code's `reachable` mask.

### Lemma 3 (what centering does and does not remove)

Two claims of different strength, previously conflated (review caught
the conflation).

**3a — gauge invariance `[PROVED]`.** If some mechanism adds to each
target an arbitrary per-target scalar `c_j` common across observers,
i.e. `W'[j,r] = W[j,r] + c_j`, then `P'[j,r] = P[j,r]` identically:
centering removes exactly the shared-scalar gauge. Immediate from
linearity of the row mean. This operationalizes the underdetermination
principle: a single shared scalar across observers is not distance
structure. (The `+1` of Lemma 2 is such a common mode and cancels.)

**3b — expectation profile, per clock model.** Under Model P with equal
rates `lambda` across chains (protocol invariance),

```
E P[j,r] = 2 lambda ( |x_j - x0_r| - mean_{r'} |x_j - x0_r'| )
```

*exactly*, by Lemma 2 (Model P) and linearity `[PROVED conditional on
the Model P hypotheses]`. Under Model D — the instrument — the same
display holds only up to a *centered quantization residue*
`theta[j,r] - mean_{r'} theta[j,r']`, bounded by 2 ranks but
observer-dependent and `t_j`-dependent: equal rates alone do not make
the phase terms common across observers, so centering does not cancel
them. "No temporal quantity survives" is therefore true only up to this
O(1) residue; the numeric example in Lemma 2 (profile shift of ~1/3
rank under a pure time translation of the target) is that residue in
the flesh. Any consumer of `P` on Model D data must budget for it.

## 4. Identifiability

### Theorem 1 (spatial order from exact profiles) `[PROVED]`

Let `R >= 2` observers have distinct positions `x0_1 < ... < x0_R`, all
targets lie in the open hull `(x0_1, x0_R)`, and profiles equal their
continuum expectations (Lemma 2, common `lambda`). Then the map

```
x  |->  ( E W[.,1], ..., E W[.,R] )
```

restricted to `(x0_1, x0_R)` is injective, and the difference
`E W[j,1] - E W[j,R] = 2 lambda ( |x_j - x0_1| - |x_j - x0_R| )`
is strictly increasing in `x_j` on the hull. Consequently the target
order along the slice is decodable from profile differences alone; from
*unlabeled* profiles (or from the dissimilarity `D` of PC-V1 Section 6)
the order is determined up to global reversal, and positions up to a
positive affine map (scale `1/(2 lambda)` unknown by construction).

*Proof.* On `[x0_1, x0_R]` the function
`f(x) = |x - x0_1| - |x - x0_R|` equals `2x - x0_1 - x0_R`: strictly
increasing with slope 2 (this is where "targets inside the hull" is
needed — outside the hull `f` is constant and order is NOT identifiable
from these two observers). Strict monotonicity of the flanking
difference gives injectivity of the labeled map and order decoding by
sorting that difference. The unlabeled step — order from `D` alone,
without observer labels — is Lemma 4(c) below. Reversal ambiguity is
the isometry group of the line acting on unlabeled configurations, and
the positive-affine position ambiguity is the unknown rank-per-length
scale `2 lambda` together with the unknown origin (Corollary). ∎

Sharpness `[PROVED by example]`: with `R = 1`, `|x - x0_1|` folds left
and right of the observer — two targets symmetric about `x0_1` have
identical profiles. One observer is never enough; the roadmap's "parallax"
is exactly the second observer breaking the fold.

Model D refinement `[PROVED]`: on the instrument's deterministic grid
the flanking difference carries a quantization error of less than 2
ranks (Lemma 2, Model D: each width carries `theta in [-1, 1)`, so
`theta_1 - theta_R in (-2, 2)`), while its slope is `2 lambda` per unit
position for each of the two flanking observers. Comparing two targets
`i, j`, the signal is `4 lambda (x_i - x_j)` against a total error
below 4 ranks, so order decoding is guaranteed once the spatial gap
exceeds one tick spacing, `|x_i - x_j| > delta`; below that the
quantization phase can invert a comparison. Identifiability on the
instrument is resolution-limited, not noise-limited.

### Lemma 4 (unlabeled decoding: seriation from `D` alone) `[PROVED]`

This closes G1. Data model: the decoder receives only the dissimilarity

```
D(i,j) = RMS_r ( P[i,r] - P[j,r] )
```

over the columns common to both targets (all `R` in the exact model) —
the exact-model counterpart of `profile_dissimilarity_matrix()`. `D` is
invariant under any single column permutation applied to all targets at
once, so "decodable from `D`" implies "decodable from unlabeled
profiles" a fortiori: `D` is the weaker data. Throughout, `R >= 2`
observers at distinct positions, targets at distinct positions in the
open hull, profiles equal to their Lemma 2 expectations.

**(a) Gap geometry.** Write `phi(x)` for the centered profile vector of
a target at position `x`: `phi_r(x) = 2 lambda ( |x - x0_r| - h(x) )`
with `h(x)` the mean of `|x - x0_r|` over `r`. `phi` is piecewise
linear in `x` with kinks exactly at the interior observer positions.
On the gap `G_k = (x0_k, x0_{k+1})` its direction (per unit `x`, in
units of `2 lambda`) is the vector `u_k` with entries `+2 b_k / R` on
the `a_k` observers to the left of `G_k` (moving right increases their
distances) and `-2 a_k / R` on the `b_k = R - a_k` observers to the
right. For gaps `k <= l` the entries multiply in three blocks (left of
`G_k`; between; right of `G_l`):

```
u_k . u_l = (4/R^2) [ a_k b_k b_l  -  (a_l - a_k) a_k b_l  +  a_k a_l b_l ]
          = (4/R^2) a_k b_l ( b_k + a_k )
          = 4 a_k b_l / R   >  0        on the hull (a_k >= 1, b_l >= 1).
```

The hull condition is exactly the positivity condition: outside
`[x0_1, x0_R]` one of the counts vanishes and with it the direction
(`|u_k|^2 = 4 a_k b_k / R = 0`) — the fold again.

**(b) Strict Robinson property.** For `x < y < z` in the hull,
`phi(y) - phi(x)` and `phi(z) - phi(y)` are nonnegative combinations of
the gap directions (with the traversed lengths as coefficients, not all
zero), so their inner product is strictly positive by (a), and with
`D = |.| / sqrt(R)`:

```
D(x,z)^2  >  D(x,y)^2 + D(y,z)^2      (strictly),
```

in particular `D(x,z) > max( D(x,y), D(y,z) )`, and `D(x,y) > 0` for
`x != y`. So `D` is a *strictly Robinson* dissimilarity in the true
spatial order — the seriation structure: every row increases strictly
moving away from the diagonal.

**(c) Decoding.** Finite target set. (i) The pair maximizing `D` is
unique and equals the extreme pair `{x_min, x_max}`: any other pair is
dominated through a strict (b)-chain. (ii) Anchored at either extreme
`a`, the map `x -> D(a, x)` is strictly monotone in position by (b), so
sorting the anchor row of `D` recovers the spatial order; anchoring at
the other extreme yields exactly the reversal. No observer labels, no
coordinates, no access to the profiles themselves — `D` alone, and
`R >= 2` suffices. ∎

**(d) What `R >= 3` actually buys.** Not decodability — (c) holds at
`R = 2`, where the expectation in v0.2's G1 that `R >= 3` would be the
clean hypothesis turned out pessimistic. Two real benefits remain.
*Falsifiability:* centering confines profiles to the sum-zero
hyperplane, dimension `R - 1`; at `R = 2` the geometry is automatically
one-dimensional and (b) tests nothing, while at `R >= 3` the strict
quadrance superadditivity is a nontrivial signature that fails for
configurations not lying on a monotone curve. *Redundancy:* the
instrument's `D` is computed over common *reachable* columns
(`min_common_columns = 4` in PC-V1), so extra observers keep pairs
comparable when reachability drops columns.

**(e) Model D version `[PROVED, constants not optimized]`.**
*Shared-centering hypothesis:* both targets are centered over the same
observer set that the exact comparison uses — e.g. every chain reaches
every target. PC-V1 enforces this by construction: scene validity
requires `min_bracketing_chains = 6 = R`, so on the instrument the
hypothesis holds for every admitted target, and the harness asserts it
rather than assuming it. Under the hypothesis: each measured width
carries `theta in [-1, 1)` (Lemma 2, Model D); centering moves each
entry by less than 2, the difference of two targets' centered errors is
below 4 per entry, and the RMS is a norm, so

```
| D_measured(i,j) - D_exact(i,j) |  <  4.
```

The hypothesis is not decorative: `profile_dissimilarity_matrix()`
centers each target over its *own* reachable columns, and two targets
with different reachable sets pick up a row-mean shift of order
`lambda` times the spread of the dropped `|dx|` values — a
coordinate-scale error, not a rank-scale one, even at zero quantization
error. Lemma 4d's redundancy remark is about keeping pairs
*comparable* when columns drop; this quantitative bound additionally
needs the centering masks to agree.

Hence every strict comparison in (c) survives measurement whenever its
exact-model margin exceeds 8 — this margin-level statement is the
primary Model-D claim, and the one the harness pins. A position-level
sufficient condition follows via (b): since all `u_k . u_l >= 4/R` on
the hull (`a_k b_l >= 1`), any chord obeys `D(y,z) >= (4 lambda / R) g`,
and the anchor-row gap obeys
`D(a,z) - D(a,y) > D(y,z)^2 / (2 D_max) >= (4 lambda g / R)^2 / (2 D_max)`,
so gaps `g > R sqrt(D_max) / lambda` always decode correctly. That
constant is loose twice over (the uniform `a_k b_l >= 1` and the
chord-to-margin step) and quadratically worse than the labeled flanking
decoder's `g > delta` (Theorem 1, Model D refinement): the anchor
decoder compares long chords, where quantization has more room to hide.
On the actual PC-V1 scene the measured perturbation is about 1 rank,
far below the worst-case 4 (Section 7).

### Corollary (what "up to reflection and scale" means here)

Reflection = global reversal of the recovered order (line isometry);
scale = the unknown `2 lambda` rank-per-length factor. Both are gauge:
no order-intrinsic observable in this setup can fix them, so Theorem 1's
conclusion is the strongest possible of its type.

## 5. Stability, per clock model

The error law depends on where the randomness lives, and the two clock
models put it in different places.

### Model D (the instrument): resolution-limited, not noise-limited

On the deterministic grid, `W[j,r]` is a *function* of the target
coordinates — zero sampling noise. The entire position uncertainty is
the quantization of Lemma 2 (Model D): at most 2 ranks, i.e. a
positional error `<= delta = 1/lambda`, decreasing like `1/K` in the
tick count `K`, **not** like any `rho^{-1/2}` — and unaffected by the
bulk sprinkling density, because `ticks_per_chain` does not scale with
`n_events`. The sprinkling enters only through *which targets exist and
where*, not through the measured widths. `[PROVED, given Lemma 2
Model D]`

### Theorem 2 (Model P: order recovery with high probability)
`[CONJECTURED, route below]`

Ticks Poisson of rate `lambda`. For targets `i, j` in the hull with
spatial separation `g = |x_i - x_j|`, the flanking profile-difference
estimator recovers their relative order except with probability at most

```
2 exp( - c * lambda g^2 / L ),
```

where `L` bounds the bracket widths involved; a union bound over pairs
gives full order recovery w.h.p. once
`min gap >~ sqrt( L * log(n) / lambda )`.

*Route.* By the Lemma 2 identity, `W - 1` is a Poisson count over the
bracket interval — no overshoot terms — so Bennett/Chernoff gives
fluctuations `~ sqrt(lambda L)` against a signal of `4 lambda g` (slope
2 per observer, two flanking observers). To write: independence
bookkeeping (the same ticks enter several brackets — interval
disjointness for well-separated targets, or negative association), and
the boundary caveat.

### The rho^{-1/2} target is conditional on future instrumentation

The roadmap's `error ~ 1/sqrt(rho * area)` law presumes tick statistics
coupled to the sprinkling density (`lambda ~ rho * ell`, e.g. chains
harvested from the sprinkling). PC-V1's chains are not that (Model D),
so **that law is currently untestable in this codebase**: varying
`n_events` in the existing generator does not vary the clock. Testing it
requires a genuine density-coupled tick protocol — a harvested-chain
constructor with its own audit — listed in Section 7 as future
instrumentation, not as a v0.2 verification target.

## 6. Known gaps (the honest list)

- **G1 — unlabeled decoding.** **CLOSED (v0.3, Lemma 4).** Order from
  `D` alone is proved for `R >= 2` — the expected `R >= 3` hypothesis
  turned out unnecessary in the exact model. The mechanism is the strict
  Robinson (seriation) structure of `D`: gap-direction inner products
  `4 a_k b_l / R > 0` on the hull make `D` strictly increasing away from
  the diagonal, so the max pair identifies the extremes and one anchor
  row sorts the rest. What `R >= 3` actually buys is falsifiability of
  the Robinson signature and reachability redundancy (Lemma 4(d)).
- **G2 — the stochastic clock model has no instrument.** (Rewritten in
  v0.2: the v0.1 text claimed the code harvests maximal-ish paths from
  the sprinkling with longest-path tick statistics — it does not.
  `build_positive_control_scene()` appends *exact* vertical worldlines
  via `make_stationary_observer_chain_1p1()` with deterministic
  `np.linspace` ticks, `ticks_per_chain = 96` fixed regardless of
  `n_events`.) The instrument is Model D through and through; Model P
  and every `rho`-scaling statement hang on a density-coupled tick
  protocol (harvested chains or Poisson-thinned clocks) that would have
  to be built and audited first. Until then, Model-P results are theory
  about a future instrument.
- **G3 — dependence between brackets.** Overlapping brackets share ticks;
  the union bound in Theorem 2 needs the independence bookkeeping made
  precise. `[PROVABLE, standard but fiddly]`
- **G4 — 1+1D only.** The roadmap's ">= 3 non-collinear observers"
  phrasing anticipates 2+1D, where the hull condition and the reflection
  group change. Out of scope for v0.1; the 1+1D statements are the ones
  the existing instrument exercises.

## 7. Numerical verification plan

With existing generators (Model D — the right knob is `ticks_per_chain`,
not `n_events`):

1. **Slope/quantization check** (Lemma 2, Model D): regress measured
   `W[j,r]` on `|x_j - x0_r|`; confirm slope `2 lambda`, intercept
   `1 + theta`, and that every residual is within the proved +/-1 band.
2. **Fold demonstration** (Theorem 1 sharpness): single-observer profiles
   for mirrored target pairs must coincide (exactly, on Model D, up to
   quantization); two observers must separate them.
3. **Resolution scaling** (Section 5, Model D): positional RMSE and
   order-recovery error vs `ticks_per_chain` on a log-log axis; target
   slope `-1` in `K` (quantization law). Varying `n_events` at fixed
   `K` must leave the width errors unchanged — a direct falsifier for
   any accidental density dependence.
4. **Centered-residue exhibit** (Lemma 3b): reproduce the reviewer-class
   example — pure time translations of a target moving the centered
   profile by O(1) ranks — and check the residue never exceeds the
   proved bound.
5. Reuse `measure_bracket_echo_profiles` / `PositiveControlScene`
   unchanged; the harness is analysis-only and carries no gate — it
   verifies theory, it does not certify the instrument.
6. **Unlabeled decoding check** (Lemma 4): exact-model `D` on random
   hull configurations at `R = 2, 3, 6` — gap-direction inner products
   match `4 a_k b_l / R`, strict quadrance superadditivity on every
   triple, anchor decoder recovers the order up to reversal (including
   at `R = 2`). Pipeline: `|D_measured - D_exact| < 4` on the PC-V1
   scene, and every anchor comparison with exact-model margin above 8
   is ordered correctly by the measured decoder (the Model-D claim at
   its proved strength — full-order recovery is *not* asserted, since
   sprinkled targets may sit closer than the decodable separation).
   The shared-centering hypothesis of Lemma 4e is asserted, not
   assumed: the check fails unless every chain reaches every target,
   because with unequal reachable sets the per-target centering shifts
   row means at coordinate scale and the `< 4` bound does not apply.

Future instrumentation (required before any `rho^{-1/2}` claim): a
density-coupled tick protocol (harvested chains or Poisson-thinned
clocks) with its own audit; only then does Theorem 2's scaling become
testable. Out of scope for the existing-generator harness.

### Execution outcome (2026-07-16 KST / 2026-07-15 UTC)

The harness (`experiments/theory/t1_verification.py`, regression tests in
`tests/test_t1_verification.py`) ran all checks below. All passed:

1. Quantization band: 400 controlled targets and 240 full-pipeline
   measurements (`build_positive_control_scene` seed 0), residuals in
   `[-0.93, 0.99]` — inside the proved `[-1, 1)` band, zero violations.
2. Fold: mirrored targets identical on the centered observer (`W = 14`
   both), separated by 27 ranks on an offset observer (predicted 27.14).
3. Resolution: RMSE log-log slope `-1.017` against target `-1` over
   `K = 12..384`, every pointwise error within the proved `delta/2`
   bound; widths bit-identical across bulk sizes 0/300/3000 (the density
   falsifier — Model D's clock does not know `rho`).
4. Centered residue: the review's counterexample reproduced *exactly*
   (`t = 0` profile `[77/3, 35/3, -7/3, -49/3, -49/3, -7/3]`, matching
   the numbers quoted in the PR #4 thread); max residue 0.79 < 2.
5. Edge and convention pinning: the tick-coincident orphan measures
   `W = 2`, residual exactly `+1`, and the band predicate rejects it
   (Revision note 5); Lemma 2's null-aligned example `(0, 0.5)` measures
   `W = 4`, residual exactly `-1` (`theta = -1`, the closed band edge,
   in band), with the strict-order counterfactual measuring `W = 6` —
   and a meta-regression swaps the strict relation in at runtime and
   asserts the check then fails (Revision note 6).
6. Density falsifier, builder level: scenes built by
   `build_positive_control_scene()` at `n_events` = 300/900/2700 place
   bit-identical chain worldlines and measure bit-identical widths for
   one fixed appended target set; the hand-built-order variant (item 3)
   is kept as a unit invariant of the order machinery (Revision note 6).
7. Unlabeled decoding (Lemma 4, v0.3): exact model — gap-direction
   inner products match `4 a_k b_l / R`, strict quadrance
   superadditivity holds on every triple, and the anchor decoder
   recovers the order up to reversal at `R = 2, 3, 6` (the `R = 2`
   sufficiency of Lemma 4c pinned explicitly). Pipeline (seed-0 scene,
   40 targets) — the shared-centering hypothesis holds and is asserted
   (all 40 targets reachable on all 6 chains, per scene validity);
   `max |D_measured - D_exact| = 0.998 < 4`, and all 515 anchor
   comparisons with exact margin above 8 sort correctly, zero
   violations. Diagnostic, not asserted: the measured argmax pair *is*
   the true extreme pair, and only 4 of 780 target pairs invert in the
   full measured order, all with true separations below one tick
   spacing (max 0.0053 against `delta = 0.0147`) — the resolution
   limit, exactly where Lemma 4e permits it.

The `[PROVED]` Model-D statements of Lemmas 1-3 are therefore also
verified against the instrument, and the band/fold/density assertions are
pinned in CI as exact (non-statistical) regressions.

## Revision notes (after PR reviews; notes 1-6 are v0.1 -> v0.2, note 7 is v0.3)

1. G2 rewritten: the v0.1 description of the observer chains was wrong
   about the code — PC-V1 appends deterministic uniform-grid worldlines
   (`ticks_per_chain` fixed, independent of `rho`); it does not harvest
   paths from the sprinkling. Consequently the `rho^{-1/2}` law was
   untestable as planned; Sections 5 and 7 now separate the
   deterministic instrument (resolution law `1/K`) from the stochastic
   idealization (conditional on future instrumentation).
2. Lemma 2's overshoot argument double-counted the Poisson constant.
   Replaced by the reviewer's direct rank-gap identity `W = N + 1`,
   which is exact realization-by-realization; the Poisson expectation is
   `2 lambda |dx| + 1` exactly, and the statement strengthens to any
   simple stationary tick process via Campbell's formula (the renewal
   hypothesis and its unproven uniformity are gone).
3. Lemma 3 split: the gauge-invariance algebra (proved) is now separate
   from the expectation display (conditional on the clock model), and
   the observer-dependent quantization residue that centering does NOT
   remove is stated with its bound, matching the numeric counterexample
   raised in review.
4. Order convention aligned with the instrument (second review round):
   Section 2 now uses the null-inclusive causal order that
   `causal_matrix_minkowski()` implements. Under v0.1's strict
   (chronological) order the rank-gap identity fails on ticks exactly
   on the target's light cone — probability zero under Model P, but
   realizable on Model D's deterministic grid, per the review's
   `(0, 0.5)` counterexample. The identity's partition argument, the
   open-interval count, and `theta`'s range (`[-1, 1)`) are restated
   under the shared convention.
5. Coincidence hypothesis added to the identity (harness review round):
   at `|dx| = 0` a tick exactly at the target's time has `dt = 0` and is
   unrelated even under the null-inclusive order, breaking the partition
   with `W = N + 2` and residual exactly `+1`. The hypothesis "no tick
   spacetime-coincident with the target" (automatic for `|dx| > 0`)
   excludes it; the verification harness pins the case, and its band
   predicate keeps the `[-1, 1)` upper edge genuinely open — a
   symmetric float tolerance would have admitted the true `+1`.
6. Convention regression and builder-level density falsifier (second
   harness review round). The reviewer swapped `causal_matrix_1p1` for
   the strict relation at runtime and every harness check stayed green:
   all sampled targets sit in general position, where the two
   conventions agree, so the convention note 4 calls load-bearing was
   asserted by the document but not pinned by the harness. The harness
   now measures Lemma 2's own null-aligned example (`W = 4` inclusive
   vs `6` strict, residual `-1` vs `+1`) and a meta-regression performs
   the swap and demands failure. Separately, the density falsifier only
   appended bulk events to hand-built chains — an invariance the
   pairwise causal matrix guarantees no matter what the scene builder
   does — so a builder-level check now runs
   `build_positive_control_scene()` at three densities and asserts
   bit-identical chains and widths for a fixed appended target set; the
   direct-order form is kept as a unit invariant.
7. G1 closed (v0.3): Lemma 4 proves unlabeled decoding from `D` alone
   via the strict Robinson structure, with `R >= 2` sufficient — v0.2's
   expectation that `R >= 3` would be the clean hypothesis was
   pessimistic; `R >= 3`'s real roles (signature falsifiability,
   reachability redundancy) are stated as Lemma 4(d). Theorem 1 and its
   Model D refinement upgrade `[PROVABLE] -> [PROVED]` with the labeled
   and quantized arguments written out. The Model-D decoding bound is
   stated at margin level (exact margin > 8 survives measurement) with
   an explicitly loose position-level corollary; the harness pins the
   margin-level form. Review (PR #6 round) caught a missing hypothesis
   in the Model-D bound: `profile_dissimilarity_matrix()` centers each
   target over its own reachable columns, so with unequal reachable
   sets the row-mean shift is coordinate-scale, not rank-scale, and
   `< 4` fails structurally. The shared-centering hypothesis is now
   stated in Lemma 4e (PC-V1 satisfies it by scene validity,
   `min_bracketing_chains = R`) and asserted by the harness instead of
   assumed.

## 8. Relation to the frozen program

T1 makes no claim about the P3-P7 gates, G, or any frozen artifact. If
Theorems 1-2 close, they upgrade the *interpretation* of PC-V1's
reconstruction from "the fitter recovers it" to "any consistent decoder
must recover it, with this error law" — the finite, order-intrinsic
criterion the external review asked for. Any manuscript use waits until
every cited tag reads `[PROVED]`.
