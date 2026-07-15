# T1: Parallax identifiability and stability of bracket-width echo profiles

Status: **THEORY DRAFT v0.1 — statements and proof programs; nothing frozen**
(2026-07-15).

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
Minkowski order `(t,x) < (t',x')` iff `t'-t > |x'-x|`. Fix:

- A Poisson sprinkling `S` of density `rho` on `M`.
- `R` observer chains `C_1..C_R`: for T1 v0.1 an observer is an inertial
  worldline at fixed spatial coordinate `x0_r` (the code's chains are
  sprinkled-point paths hugging such a worldline; the gap between the two
  models is addressed in Section 6, G2). Ticks are the chain's elements,
  indexed by rank; rank is the only clock.
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

*Proof.* `(t, x0) < (t_e, x_e)` iff `t_e - t > |dx|`, so the supremum of
preceding observer times is `t_e - |dx|`; symmetrically the infimum of
succeeding times is `t_e + |dx|`. For an inertial worldline coordinate
time is proper time. ∎

### Lemma 2 (rank calibration) `[PROVABLE]`

If ticks along `C_r` form a one-dimensional Poisson process of rate
`lambda_r` in proper time (or any renewal process with mean spacing
`1/lambda_r`), then conditional on reachability and away from the diamond
boundary,

```
E[ W[j,r] ] = 2 lambda_r |x_j - x0_r| + O(1),
```

with the `O(1)` a boundary/overshoot term bounded uniformly in `|dx|`.
In particular `E W` is affine in `|dx|`, hence strictly increasing.

*Route.* The bracket rank count is the number of ticks in the proper-time
interval of Lemma 1, plus the two overshoot terms of a renewal process
(inspection paradox gives the `O(1)`; for Poisson it is exactly `+1` in
expectation on each side). Boundary truncation only shrinks brackets and is
excluded by the reachability condition. Details to write: the exact
overshoot constant and uniformity of the `O(1)`.

### Lemma 3 (centering removes the common mode) `[PROVED]`

Suppose every chain shares tick rate `lambda` (protocol invariance) and
some mechanism adds to each target an arbitrary per-target scalar `c_j`
common across observers, i.e. `W'[j,r] = W[j,r] + c_j`. Then
`P'[j,r] = P[j,r]`: the centered profile is invariant, and in continuum
expectation

```
E P[j,r] = 2 lambda ( |x_j - x0_r| - mean_{r'} |x_j - x0_r'| ),
```

which depends on `x_j` only — no temporal quantity survives.

*Proof.* Immediate from linearity of the row mean and Lemma 2. ∎

This operationalizes the underdetermination principle: a single shared
scalar across observers is not distance structure, and centering removes
exactly that gauge.

## 4. Identifiability

### Theorem 1 (spatial order from exact profiles) `[PROVABLE]`

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

*Route.* On `[x0_1, x0_R]` the function
`f(x) = |x - x0_1| - |x - x0_R|` equals `2x - x0_1 - x0_R`: strictly
increasing with slope 2 (this is where "targets inside the hull" is
needed — outside the hull `f` is constant and order is NOT identifiable
from these two observers). Injectivity and order-decoding follow;
reversal ambiguity is the isometry group of the line acting on unlabeled
configurations. Details to write: the unlabeled-profile step (recovering
order from `D` rather than from labeled columns), which is where
`R >= 3` earns its keep as redundancy and as a check that the profile
geometry is one-dimensional.

Sharpness `[PROVED by example]`: with `R = 1`, `|x - x0_1|` folds left
and right of the observer — two targets symmetric about `x0_1` have
identical profiles. One observer is never enough; the roadmap's "parallax"
is exactly the second observer breaking the fold.

### Corollary (what "up to reflection and scale" means here)

Reflection = global reversal of the recovered order (line isometry);
scale = the unknown `2 lambda` rank-per-length factor. Both are gauge:
no order-intrinsic observable in this setup can fix them, so Theorem 1's
conclusion is the strongest possible of its type.

## 5. Stability under sprinkling `[CONJECTURED, route below]`

### Theorem 2 (order recovery with high probability)

Same geometry, ticks Poisson of rate `lambda` per unit proper time.
There exist constants such that for any two targets `i, j` in the hull
with spatial separation `g = |x_i - x_j|`, the profile-difference
estimator recovers their relative order except with probability at most

```
2 exp( - c * lambda^2 g^2 / (lambda * L) )  =  2 exp( - c * lambda g^2 / L ),
```

where `L` bounds the bracket widths involved; and a union bound over all
pairs gives full order recovery w.h.p. once

```
min gap  >~  sqrt( L * log(n) / lambda ).
```

With `lambda ~ rho * ell` for chains harvested from a sprinkling of
density `rho` (Section 6, G2), the positional standard error scales as
`1 / sqrt(rho * area)`, matching the roadmap's target bound.

*Route.* `W[j,r]` is a Poisson count over the Lemma 1 interval plus O(1)
overshoot; Poisson concentration (Bennett/Chernoff) gives fluctuations
`~ sqrt(lambda L)`. The signal separating the order of `i` and `j` in the
flanking-difference statistic is `4 lambda g` (slope 2 per observer, two
observers). Chernoff + union bound. To write: independence bookkeeping
(the same ticks enter several brackets — use the interval-disjointness of
brackets for well-separated targets, or a negative-association argument),
and the boundary-truncation caveat.

### What must NOT be claimed

Theorem 2 with its constants is a concentration statement about the
*idealized* tick process. It becomes a statement about the code's scenes
only after G2 (below) is closed. Until then, the numerical harness
(Section 7) tests the *scaling*, not the constants.

## 6. Known gaps (the honest list)

- **G1 — unlabeled decoding.** Theorem 1's last step (order from `D`
  alone, without observer labels) needs its own small argument; `R >= 3`
  expected to be the clean hypothesis. `[PROVABLE, to write]`
- **G2 — chains are not inertial worldlines.** The code's observer chains
  are maximal-ish paths in the sprinkled order near a worldline, not
  exact worldlines; their tick process is not exactly Poisson in proper
  time (path fluctuates in space, spacing follows longest-path statistics
  with the Myrheim-Meyer constant). Needed: either a coupling bounding
  the discrepancy, or restate Lemma 2 for harvested chains directly.
  This is the main modelling gap between theory and instrument.
- **G3 — dependence between brackets.** Overlapping brackets share ticks;
  the union bound in Theorem 2 needs the independence bookkeeping made
  precise. `[PROVABLE, standard but fiddly]`
- **G4 — 1+1D only.** The roadmap's ">= 3 non-collinear observers"
  phrasing anticipates 2+1D, where the hull condition and the reflection
  group change. Out of scope for v0.1; the 1+1D statements are the ones
  the existing instrument exercises.

## 7. Numerical verification plan (existing generators, no new physics)

1. **Monotonicity/slope check** (Lemma 2): scenes from the PC-V1
   generator over a density ladder; regress measured `W[j,r]` on
   `|x_j - x0_r|`; confirm affine fit, slope ratio `2 lambda`, and the
   `O(1)` intercept.
2. **Fold demonstration** (Theorem 1 sharpness): single-observer profiles
   for mirrored target pairs — must be statistically indistinguishable;
   two observers must separate them.
3. **Error scaling** (Theorem 2): order-recovery error rate and positional
   RMSE vs `rho` on a log-log axis; target slope `-1/2` in
   `sqrt(rho * area)`. Any systematic deviation flags G2 as material
   rather than technical.
4. Reuse `measure_bracket_echo_profiles` / `PositiveControlScene`
   unchanged; the harness is analysis-only and carries no gate — it
   verifies theory, it does not certify the instrument.

## 8. Relation to the frozen program

T1 makes no claim about the P3-P7 gates, G, or any frozen artifact. If
Theorems 1-2 close, they upgrade the *interpretation* of PC-V1's
reconstruction from "the fitter recovers it" to "any consistent decoder
must recover it, with this error law" — the finite, order-intrinsic
criterion the external review asked for. Any manuscript use waits until
every cited tag reads `[PROVED]`.
