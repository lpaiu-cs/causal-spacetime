# T1: Parallax identifiability and stability of bracket-width echo profiles

Status: **THEORY DRAFT v0.5 — statements and proof programs; nothing frozen**
(v0.5 2026-07-16 KST: G2 instrumented, `rho^{-1/2}` shown
protocol-dependent; v0.4 2026-07-16 KST: G3 closed, Theorem 2 upgraded
to `[PROVED]`; v0.3 2026-07-16 KST: G1 closed, Theorem 1 upgraded to
`[PROVED]`; v0.2 2026-07-16 KST / 2026-07-15 UTC; v0.1 2026-07-15 —
revised after PR reviews, see Revision notes).

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
    Realized as theory-track instrumentation in v0.5 by
    `make_poisson_clock_chain_1p1()` (Section 6, G2); the frozen PC-V1
    scene still uses only Model D.
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
is strictly increasing in `x_j` on the hull. Consequently, from
*labeled* profiles the target order along the slice is decodable from
the flanking difference alone, and positions follow up to a positive
affine map (the flanking difference is affine in `x_j` with positive
slope `4 lambda`; slope and origin unknown by construction). From
*unlabeled* profiles (or from the dissimilarity `D` of PC-V1
Section 6) the order is determined up to global reversal — and nothing
metric: `D` alone provably does not determine spacings, even up to a
positive affine map (Lemma 4f).

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

**(f) Sharpness: `D` determines order, not spacings `[PROVED by
example]`.** With unknown observer positions, nothing metric survives
in `D`. Two hypothesis-satisfying configurations (`R = 6`,
`lambda = 1`) — observers `[0, .1, .21, .34, .44, .65]` with targets
`[.02, .28, .53]`, versus observers
`[0, .0292803932..., .0752379602..., .3002515990..., .3119385925...,
.65]` with targets `[.02, .2411728090..., .5114951311...]` — produce
the *same* dissimilarity matrix (entrywise to float precision, max
difference `~1e-15`), while their affine invariants
`(x_2 - x_1)/(x_3 - x_1)` differ: `0.5098...` versus `0.45`. No
positive affine map relates the two target sets, so no decoder reading
`D` alone can output spacings. The mechanism is (a) itself: the
per-gap speeds `4 a_k b_l / R` depend on the unknown observer layout,
which an adversary can re-tune to absorb spacing information while
preserving every pairwise `D`. Order up to reversal (c) is therefore
the exact content of the `D`-based claim — Theorem 1's positive-affine
clause belongs to the *labeled* flanking decoder only. For unlabeled
profiles (more data than `D`) spacing recovery is neither claimed nor
refuted here. (Review supplied the counterexample; the harness
reproduces it to float precision.)

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
no order-intrinsic observable in this setup can fix them. For the
*labeled* decoder the conclusion — order plus positions up to a
positive affine map — is the strongest possible of its type. For the
`D`-only decoder the strongest possible conclusion is order up to
reversal: Lemma 4f exhibits two configurations with identical `D` and
affinely inequivalent target sets, so spacings are genuinely not in
`D`.

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

### Theorem 2 (Model P: order recovery with high probability) `[PROVED]`

Hypotheses: each chain an independent simple Poisson process of rate
`lambda` on its worldline, tick support containing every radar interval
below with at least one tick on each side (reachability; boundary
caveat at the end). Targets `i, j` in the open hull with `x_j < x_i`,
gap `g = x_i - x_j`, *arbitrary* times. Estimator: the flanking
difference `D_m = W[m,1] - W[m,R]`, order declared by sorting `D`.
Write `I_{m,r}` for target `m`'s open radar interval on chain `r`
(length `2 d_{m,r}`, `d_{m,r} = |x_m - x0_r|`) and
`L = 2 max(d_{i,1}, d_{j,1}, d_{i,R}, d_{j,R})`.

**Claim 1 (pairwise, ties counted as errors).**

```
P( D_i <= D_j )  <=  exp( - 2 lambda g^2 / (L + g/3) ).
```

**Claim 2 (same-slice refinement).** If additionally `t_i = t_j`, then
`D_i - D_j >= 0` on *every* realization — the estimator never strictly
inverts a same-slice pair — and

```
P( D_i = D_j )  =  exp( - 4 lambda g )        exactly.
```

**Claim 3 (full order).** For `n` targets with minimum adjacent gap
`g_min`, and with the *global* length bound
`L_all = 2 max_{m, r in {1,R}} d_{m,r}` (over all `n` targets and both
flanking chains — the pairwise `L` above is local to one pair):
`P(recovered order != true order) <=
(n-1) exp( - 2 lambda g_min^2 / (L_all + g_min/3) )`, so order
recovery holds w.h.p. once `g_min >~ sqrt( L_all log n / lambda )` —
improving to `g_min >~ log n / (4 lambda)` on a common slice, where
the per-pair term is the exact tie probability of Claim 2.

*Proof.*

**Step 1 — the G3 bookkeeping: shared regions cancel.** By the
rank-gap identity (Lemma 2; applicable a.s., since Poisson ticks are
simple and a.s. avoid both the target's light cone and its position,
so the null-alignment and coincidence hypotheses hold with probability
1), `W[m,r] = N_r(I_{m,r}) + 1` with `N_r` the chain-`r` counting
measure. On one chain,

```
N_r(I_i) - N_r(I_j) = N_r(I_i \ I_j) - N_r(I_j \ I_i):
```

the count over the shared region `I_i ∩ I_j` cancels *exactly*,
realization by realization. What remains are counts over disjoint
regions of one Poisson process — independent — and the two flanking
chains are independent of each other. So with
`A_r = I_{i,r} \ I_{j,r}`, `C_r = I_{j,r} \ I_{i,r}`,

```
Delta := D_i - D_j = N_1(A_1) - N_1(C_1) - N_R(A_R) + N_R(C_R)
```

is a signed sum of four *mutually independent* Poisson counts.
Overlapping brackets leave no residual dependence in a pairwise
comparison — this is the entire content of the worry G3 recorded.

**Step 2 — mean and variance.** `|I_{i,r}| - |I_{j,r}| = |A_r| - |C_r|`,
and inside the hull `d_{i,1} - d_{j,1} = g`, `d_{i,R} - d_{j,R} = -g`,
so `E Delta = 2 lambda g - (-2 lambda g) = 4 lambda g`. Independence
gives `Var Delta = lambda (|I_{i,1} Δ I_{j,1}| + |I_{i,R} Δ I_{j,R}|)
<= 4 lambda L` — only the symmetric differences enter; the shared
regions dropped out of the variance too.

**Step 3 — Chernoff.** For `theta in [0, 3)`, the downward-tail factors
obey `E exp(-theta (N(A) - lambda|A|)) = exp(lambda|A|(e^{-theta} - 1 +
theta)) <= exp(lambda|A| theta^2/2)` and the upward ones
`E exp(theta (N(C) - lambda|C|)) = exp(lambda|C|(e^{theta} - 1 -
theta)) <= exp(lambda|C| theta^2 / (2(1 - theta/3)))`. The standard
Bernstein optimization then yields, with `s = 4 lambda g` and
`v = Var Delta`,

```
P(Delta <= 0) = P(Delta - E Delta <= -s)
             <= exp( - s^2 / (2 (v + s/3)) )
             <= exp( - 2 lambda g^2 / (L + g/3) ).
```

(Constants not optimized.)

**Step 4 — same slice.** With `t_i = t_j` the intervals are concentric:
`I_{j,1} ⊂ I_{i,1}` and `I_{i,R} ⊂ I_{j,R}`, so `C_1 = A_R = ∅` and
`Delta = N_1(A_1) + N_R(C_R) >= 0` pathwise — a sum of two independent
`Poisson(2 lambda g)` variables (`|A_1| = |C_R| = 2g`). `Delta = 0`
iff both annuli are empty: probability `exp(-2 lambda g)^2 =
exp(-4 lambda g)`, exactly.

**Step 5 — union bound; where cross-pair dependence lives and why it
is harmless.** Different pairs *do* share ticks (`D_i` enters every
comparison involving target `i`), but Claim 3 needs no independence:
sorting real numbers recovers the true order as soon as `D` is
strictly increasing across the `n - 1` *adjacent* true pairs
(transitivity of `<`), and the union bound over those events uses
subadditivity of probability only. Each adjacent pair `k` has its own
gap `g_k >= g_min` and its own local length bound `L_k <= L_all`, and
the Claim 1 exponent `2 lambda g^2 / (L + g/3)` is increasing in `g`
and decreasing in `L`, so every per-pair bound is dominated by the
single displayed term with `(g_min, L_all)`:
`g_k^2/(L_k + g_k/3) >= g_min^2/(L_k + g_min/3) >=
g_min^2/(L_all + g_min/3)`. ∎

Boundary caveat `[PROVED, crude]`: on a finite tick support,
"`W` exists" is itself an event. Each width needs one tick on each
side of its radar interval; if every interval sits at distance `>= m`
from the support's ends, each one-sided failure is an empty Poisson
segment of length `>= m`, probability `<= exp(-lambda m)`, so a union
bound adds `2 x (number of widths involved) x exp(-lambda m)` to each
bound above (`8 exp(-lambda m)` for the pairwise Claim 1, `4n
exp(-lambda m)` for Claim 3). Targets failing reachability are
excluded by the instrument's `reachable` mask rather than imputed,
matching the code.

### The rho^{-1/2} law is protocol-dependent (instrumented in v0.5)

The roadmap's `error ~ 1/sqrt(rho * area)` law presumes tick statistics
coupled to the sprinkling density as `lambda ~ rho * ell`. PC-V1's
frozen chains are not that (Model D), but the coupling is now
*instrumented* (`density_coupled_clocks.py`, audited) in two protocols,
and the law turns out to hold for exactly one of them
(`experiments/theory/t1_g2_density_scaling.py`):

- **Poisson-thinned clock** (`lambda = rho * ell`, exact worldline):
  realizes Model P with the presumed coupling. Here
  `W - 1 ~ Poisson(2 lambda d)` by the Lemma 2 identity, so the
  distance estimator `d_hat = (W - 1)/(2 lambda)` is unbiased with
  `sd = sqrt( d / (2 lambda) ) ~ rho^{-1/2}` — a `[PROVED]` corollary
  of Theorem 2's setup, and `1/sqrt(rho * area)` with `area ~ d * ell`.
  Measured: rate exponent `+0.995`, RMSE exponent `-0.463`, per-density
  RMSE tracking the exact prediction (grid-mean ratio within 15%,
  intra-chain correlation accounted for).
- **Harvested chain** (longest causal chain of *sprinkled events* in a
  spatial tube — the order-intrinsic clock): its rate is a measurement,
  not a choice, and it comes out `lambda ~ rho^{0.49..0.55}` — the
  discreteness scale `sqrt(rho)`, **not** `rho`. The `rho^{-1/2}` law
  therefore does *not* transfer: measured error exponent `-0.317` with
  a discreteness-scaled tube (`w ~ 3/sqrt(rho)`), flattening to
  `-0.165` with a fixed-width tube (whose position wiggle adds an
  error floor). `[MEASURED]`

Mechanism note, open: with `lambda ~ sqrt(rho)`, a Poisson-rate guess
gives exponent `-1/4`; a maximal path is *more regular* than Poisson
(KPZ-like concentration would suggest `-1/3`). The measured `-0.32`
sits in that band, nearer `-1/3`, but distinguishing the fluctuation
class needs a dedicated study — recorded in G2 as the remaining open
question. Any density-scaling claim must name its protocol.

## 6. Known gaps (the honest list)

- **G1 — unlabeled decoding.** **CLOSED (v0.3, Lemma 4).** Order from
  `D` alone is proved for `R >= 2` — the expected `R >= 3` hypothesis
  turned out unnecessary in the exact model. The mechanism is the strict
  Robinson (seriation) structure of `D`: gap-direction inner products
  `4 a_k b_l / R > 0` on the hull make `D` strictly increasing away from
  the diagonal, so the max pair identifies the extremes and one anchor
  row sorts the rest. What `R >= 3` actually buys is falsifiability of
  the Robinson signature and reachability redundancy (Lemma 4(d)). The
  `D`-based conclusion is order-only: spacings are provably not in `D`
  (Lemma 4f counterexample), and Theorem 1's positive-affine clause
  belongs to the labeled flanking decoder.
- **G2 — the stochastic clock model has no instrument.**
  **INSTRUMENTED (v0.5); fluctuation class open.** (History: v0.1
  wrongly claimed the code harvests paths from the sprinkling; v0.2
  corrected that the frozen instrument is Model D through and through —
  `build_positive_control_scene()` appends exact `np.linspace`
  worldlines with `ticks_per_chain = 96` fixed regardless of
  `n_events`.) v0.5 builds the missing protocols as theory-track
  instrumentation (`density_coupled_clocks.py`, audited: chain
  property, containment, simplicity, determinism): a Poisson-thinned
  clock realizing Model P with `lambda = rho * ell`, and a
  harvested-chain constructor whose ticks are sprinkled events.
  Outcome (Section 5): the `rho^{-1/2}` law holds for the thinned
  protocol as a proved corollary and is measured at exponent `-0.46`;
  harvested chains couple at the discreteness scale
  (`lambda ~ sqrt(rho)`, measured `0.49..0.55`) and give a distinctly
  shallower error law (measured `-0.32`), so the roadmap's law is
  protocol-dependent. *Remaining open:* the harvested chain's
  fluctuation class (`-1/4` Poisson-rate guess vs `-1/3` KPZ-like;
  measurement sits between, nearer `-1/3`). The frozen PC-V1
  instrument is unchanged; any confirmatory use of these protocols
  would need its own prereg freeze.
- **G3 — dependence between brackets.** **CLOSED (v0.4, Theorem 2).**
  The bookkeeping turned out cleaner than "standard but fiddly": in a
  pairwise flanking comparison the shared interval region cancels
  *exactly* (`N(I_i) - N(I_j) = N(I_i \ I_j) - N(I_j \ I_i)`), leaving
  four mutually independent Poisson counts — no negative-association
  argument needed; and the full-order union bound needs subadditivity
  only, over the `n - 1` adjacent pairs. Bonus from the same
  cancellation: same-slice pairs are *pathwise* monotone (the estimator
  can tie but never strictly invert), with tie probability exactly
  `exp(-4 lambda g)` — exponentially better in the gap than the
  conjectured `exp(-c lambda g^2 / L)` form, which remains the correct
  rate for arbitrary-time pairs.
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
   Sharpness (Lemma 4f): the two same-`D` affinely-inequivalent
   configurations are pinned as a deterministic regression — `D`
   carries order, not spacings.
7. **Theorem 2 simulation check** (Model P): direct seeded Monte Carlo
   of the *stated stochastic model* — Poisson tick chains simulated in
   the harness, widths via the rank-gap identity, cross-checked
   exactly against `find_radar_ticks_from_order()` on the same draws.
   Asserted: same-slice pairs never strictly invert (pathwise claim —
   a single inversion falsifies Step 4); the tie rate, mean and
   variance match `exp(-4 lambda g)` / `4 lambda g` / `4 lambda g`
   within exact binomial/MC tolerances; arbitrary-time error rates are
   dominated by the Claim 1 bound; full-order failure rates are
   dominated by the Claim 3 union bound at parameters where that bound
   is nontrivial (`< 1`). Every simulated width's reachability flag is
   collected and asserted (zero unreachable draws per section — a
   nonzero count fails the check loudly instead of letting synthetic
   widths into the statistics). This verifies the *theorem*, not the
   instrument: no code constructs Poisson chains as an instrument
   (G2), and no `rho`-scaling claim is tested here.

8. **G2 instrumentation + density-scaling characterization**
   (`experiments/theory/t1_g2_density_scaling.py`): audit the
   harvested-chain constructor (chain property, tube containment,
   simplicity, determinism — the run refuses to measure if any fails),
   then measure `lambda(rho)` and `RMSE(rho)` exponents for three arms
   (thinned `lambda = rho * ell`; harvested fixed tube; harvested
   discreteness-scaled tube). The thinned arm is asserted against its
   exact `[PROVED]` prediction (`sd = sqrt(d / 2 lambda)`), with
   tolerances derived from the intra-chain correlation (targets on one
   chain share its tick realization — Theorem 2 Step 1's shared
   regions); the harvested arms are characterization with sanity bands
   only, and the fluctuation-class question is left open, not settled
   by a fit.

The density-coupled protocols are theory-track instrumentation: no
frozen gate consumes them, and any confirmatory use would need its own
prereg freeze. The `rho^{-1/2}` claim itself is settled as
protocol-dependent (Section 5).

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
   violations. Sharpness (Lemma 4f) — the review's same-`D`
   counterexample reproduced: max entrywise difference `~1.4e-15`
   between the two configurations' dissimilarity matrices, affine
   invariants `0.5098` vs `0.45`. Diagnostic, not asserted: the measured argmax pair *is*
   the true extreme pair, and only 4 of 780 target pairs invert in the
   full measured order, all with true separations below one tick
   spacing (max 0.0053 against `delta = 0.0147`) — the resolution
   limit, exactly where Lemma 4e permits it.
8. Theorem 2 simulation (Model P, v0.4): rank-gap-identity widths equal
   `find_radar_ticks_from_order()` widths exactly on shared Poisson
   draws. Same-slice pair (`lambda = 40`, `g = 0.02`): **0 strict
   inversions in 4000 trials** (the pathwise claim), tie rate `0.046`
   against the exact `exp(-4 lambda g) = 0.041` (inside MC tolerance),
   mean and variance matching `4 lambda g`. Arbitrary-time pair with
   partial bracket overlap: error rate `0.22` dominated by the Claim 1
   bound `0.76`, variance matching Step 2's symmetric-difference
   bookkeeping. Full order (`n = 6`, `lambda = 300`): failure rate
   `0.000` dominated by the nontrivial union bound `0.478` (computed
   with the global `L_all`, as Claim 3 requires). Zero unreachable
   draws in every section, asserted rather than assumed.
9. G2 density scaling (v0.5): harvested-chain audit clean
   (causal-chain property, containment, simplicity, determinism).
   Thinned protocol (`lambda = rho * ell`, `ell = 0.1`,
   `rho = 500..32000`, 16 seeds): rate exponent `+0.995`, RMSE
   exponent `-0.463`, per-density RMSE tracking the exact
   `sqrt(d / 2 lambda)` prediction (grid-mean ratio inside 15%; every
   per-density ratio inside `[0.6, 1.6]`). Harvested chains: rate
   exponents `+0.547` (fixed tube) and `+0.494` (scaled tube) — the
   discreteness scale; RMSE exponents `-0.165` (fixed tube: wiggle
   floor visible) and `-0.317` (scaled tube). Zero unreachable
   measurements in all arms. Full table in
   `outputs/theory/t1_g2_density_scaling.json`.

The `[PROVED]` Model-D statements of Lemmas 1-3 are therefore also
verified against the instrument, and the band/fold/density assertions are
pinned in CI as exact (non-statistical) regressions.

## Revision notes (after PR reviews; notes 1-6 are v0.1 -> v0.2, notes 7-8 are v0.3, note 9 is v0.4, note 10 is v0.5)

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
8. Theorem 1's `D`-clause narrowed to order-only (PR #6 review round
   2). The v0.2 statement let "positions up to a positive affine map"
   attach to the unlabeled/`D` conclusion; the review exhibited two
   hypothesis-satisfying configurations with identical `D` and affinely
   inequivalent target sets, so that reading is false. The affine
   clause now belongs explicitly to the *labeled* flanking decoder
   (where it is proved), the `D`-only conclusion is order up to
   reversal, and the counterexample is recorded as Lemma 4f and
   reproduced to float precision by the harness. Spacing recovery from
   unlabeled profiles (more data than `D`) is left explicitly open —
   neither claimed nor refuted.
9. G3 closed (v0.4): Theorem 2 upgrades `[CONJECTURED] -> [PROVED]`
   with the independence bookkeeping written out. The feared
   dependence dissolves in two different ways at the two places it
   could enter: pairwise, the shared interval region cancels exactly
   inside `N(I_i) - N(I_j)`, leaving four mutually independent Poisson
   counts (no negative-association machinery); across pairs, the
   union bound needs only subadditivity over the `n - 1` adjacent
   comparisons. New same-slice refinement: concentric intervals make
   the flanking estimator pathwise monotone — ties possible
   (probability exactly `exp(-4 lambda g)`), strict inversions
   impossible — strengthening the conjectured `g^2/L` rate to a linear
   `g` rate on a common slice. Boundary caveat stated with a crude
   additive `4 exp(-lambda m)` term. The theorem is verified by direct
   seeded simulation of the stated model (Section 7, item 7); this
   does not touch G2 — no instrument realizes Model P.
10. G2 instrumented (v0.5): the density-coupled tick protocols the gap
   demanded now exist as theory-track instrumentation
   (`density_coupled_clocks.py` — Poisson-thinned clock and
   harvested-chain constructor, with a constructor audit the scaling
   run refuses to proceed without). The headline finding is that the
   roadmap's `rho^{-1/2}` law is *protocol-dependent*: it holds for
   the `lambda = rho * ell` thinned clock as a proved corollary of
   Theorem 2's setup (measured RMSE exponent `-0.463`, tracking the
   exact prediction), while the order-intrinsic harvested chain
   couples at the discreteness scale (`lambda ~ sqrt(rho)`, measured
   `0.49..0.55`) and yields a distinctly shallower error law
   (measured `-0.317`, discreteness-scaled tube). The harvested
   chain's fluctuation class (`-1/4` vs KPZ-like `-1/3`) is the
   remaining open question in G2. The frozen PC-V1 instrument is
   untouched.

## 8. Relation to the frozen program

T1 makes no claim about the P3-P7 gates, G, or any frozen artifact. If
Theorems 1-2 close, they upgrade the *interpretation* of PC-V1's
reconstruction from "the fitter recovers it" to "any consistent decoder
must recover it, with this error law" — the finite, order-intrinsic
criterion the external review asked for. Any manuscript use waits until
every cited tag reads `[PROVED]`.
