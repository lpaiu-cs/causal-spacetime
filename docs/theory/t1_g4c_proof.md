# G4c: a proof program for the general-dimension rigidity law

Companion to `t1_parallax_identifiability.md` Section 6 (gap G4c). Same
proof-status tags, same rule: nothing here may be cited above its tag.

- `[PROVED]` — argument written out here and believed complete.
- `[PROVABLE]` — route is standard; details to be written and checked.
- `[CONJECTURED]` — plausible, not yet proved; must not be cited as a result.

Rounds 1-3 (`t1_g4c_predictions{,_round2,_round3}.json`) established the
law by measurement across spatial dimensions one through five. This
document asks what can be *derived*. The answer is: most of it. One
step resists, it is the step rigidity theory always gets stuck on, and
it is named in Section 7 rather than papered over.

Every object below is computed and checked numerically by
`experiments/theory/t1_g4c_proof_check.py`; the lemmas are verified as
lemmas, not merely through the conclusion they support.

## 1. Setup

Observers `p_1, ..., p_R` and targets `x_1, ..., x_n`, all in `R^d`, with
no target on an observer. Write the scene as
`theta = (x_1, ..., x_n, p_1, ..., p_R) in R^{d(n+R)}`.

    phi_r(x) = |x - p_r|,   Phi(x) = (phi_1(x), ..., phi_R(x)) in R^R

Let `M = I_R - (1/R) 1 1^T`, the orthogonal projector onto

    V = {v in R^R : sum_r v_r = 0},     m := dim V = R - 1

and let `w(x) = M Phi(x) in V`, `w_j := w(x_j)`. The instrument consumes

    D(j,k) = ||w_j - w_k|| / sqrt(R)

and nothing else: no observer labels, no observer positions.

An *infinitesimal flex* is a `delta theta` with `delta D(j,k) = 0` for
every pair. The *nullity* is the dimension of the flex space. The
*gauge* is the rigid motions of the whole scene, `d(d+1)/2` of them —
translations and rotations, and **not** scale, since `D` is homogeneous
of degree one in the scene and a global dilation multiplies it.

## 2. Lemma A — the reduction `[PROVED]`

> `delta theta` is a flex **iff** the induced `delta w_j` satisfies
> `<w_j - w_k, delta w_j - delta w_k> = 0` for every pair `j, k`.

*Proof.* `R * D(j,k)^2 = ||w_j - w_k||^2`. Differentiating,
`2 R D(j,k) delta D(j,k) = 2 <w_j - w_k, delta w_j - delta w_k>`. For
`D(j,k) != 0` the two vanish together. (`D(j,k) = 0` means two targets
share a centered profile exactly; excluded on the generic set, and the
squared observable is in any case the object the pipeline forms.) ∎

This is the whole idea. **The flex problem for the scene is the
infinitesimal rigidity problem for the point cloud `{w_j}` in
`V ≅ R^m`, with every pair constrained — pulled back along the map that
builds the cloud out of the scene.** Nothing about `d` has entered.

## 3. Lemma B — flexes of a complete-graph cloud `[PROVED]`

> Let `y_1, ..., y_n in R^m` have affine span of dimension `q`. Then
> `F = { (delta y_j) : <y_j - y_k, delta y_j - delta y_k> = 0 for all j,k }`
> has dimension `n(m - q) + q(q+1)/2`.

*Proof.* Translate so the affine span is a linear subspace `W`, `dim W = q`,
containing every `y_j`. Split `delta y_j = a_j + b_j` with `a_j in W`,
`b_j in W^perp`. Each `y_j - y_k` lies in `W`, so the constraint reads
`<y_j - y_k, a_j - a_k> = 0` and the `b_j` are entirely free: that is
`n(m - q)` dimensions.

For the `W` part, translate again so `y_1 = 0`; then `{y_j}` spans `W`
linearly. The pair `(1, j)` gives `<y_j, c_j> = 0` where `c_j := a_j - a_1`.
Expanding a general pair and substituting that,

    <y_j, c_k> + <y_k, c_j> = 0   for all j, k.                    (*)

Choose a basis `y_{j_1}, ..., y_{j_q}` of `W` and define a linear map `A`
on `W` by `A y_{j_i} = c_{j_i}`. This is well defined: if
`sum_i lambda_i y_{j_i} = 0` then for every `k`, by (*),
`<y_k, sum_i lambda_i c_{j_i}> = - <sum_i lambda_i y_{j_i}, c_k> = 0`,
and the `y_k` span `W`, so `sum_i lambda_i c_{j_i} = 0`. Applying (*)
again, `<y_j, A y_k> + <y_k, A y_j> = 0` on a spanning set, so
`A + A^T = 0`: `A` is skew. Hence `a_j = a_1 + A y_j`, an infinitesimal
isometry of `W`, and distinct `(A, a_1)` give distinct `(a_j)` because
the configuration spans. That is `q(q-1)/2 + q = q(q+1)/2` dimensions. ∎

## 4. Lemma C — the decomposition `[PROVED]`

Let `L : R^{d(n+R)} -> V^n` be the derivative of `theta |-> (w_1,...,w_n)`.
By Lemma A the flex space is `L^{-1}(F)`, so

> **`nullity = dim ker L + dim(Im L ∩ F)`.**

*Proof.* Linear algebra: `dim L^{-1}(F) = dim ker L + dim(Im L ∩ F)`. ∎

Equivalently, with `Psi(delta theta, f) := L delta theta - f` on
`R^{d(n+R)} × F`, the projection `(delta theta, f) |-> delta theta` is
injective on `ker Psi` (because `f = L delta theta` is determined), so

> **`nullity = dim ker Psi`.**                                      (**)

Verified exactly on 15 measured cells, including rigid ones — see
Section 8.

## 5. Genericity hypotheses

Two conditions, both checkable and both checked numerically.

**(G1)** `rank(M U(x)) = min(d, m)`, where `U(x)` has rows
`u_r(x)^T = ((x - p_r)/|x - p_r|)^T`, so `D_x w = M U(x)`. `[PROVABLE]`
(`U` has rank `min(d,R)` generically, and `M` removes only the `1`
direction, which costs rank only on a proper subvariety of observer
placements.)

**(G2)** The cloud `{w_j}` affinely spans `V` when `m <= d`, and more
generally has affine span `q = min(m, ...)` as measured. `[PROVABLE]`
for `m <= d`: by (G1) the image of `w` is then an open subset of `V`, and
`n >= m + 1` generic points of an open set affinely span.

## 6. Theorem 1 — the regimes `R <= d + 1` `[PROVED under (G1), (G2)]`

> For `m = R - 1 <= d` and generic scenes,
>
>     nullity = dR + n(d - R + 1) + R(R-1)/2.

*Proof.* Write `delta w_j = M U(x_j) delta x_j - M W_j delta P`, the two
terms being the target and observer contributions.

*(i) `ker L`.* By (G1) `M U(x_j) : R^d -> V` has rank `m`, hence is onto
`V`. The observer term lies in `V`. So for **every** `delta P` there is a
solution `delta x_j`, unique up to `ker(M U(x_j))`, of dimension `d - m`.
Therefore `dim ker L = dR + n(d - m)`.

*(ii) `L` is onto.* `dim Im L = d(n+R) - dR - n(d-m) = nm = dim V^n`.
Hence `Im L ∩ F = F`.

*(iii) `dim F`.* By (G2), `q = m`, so Lemma B gives `dim F = m(m+1)/2`.

Lemma C sums these: `dR + n(d-m) + m(m+1)/2`, which is the stated form
since `m = R - 1`. ∎

**This derives the counting law.** Rounds 1-3 recorded
`nullity = dR + d(d+1)/2` in the saturated regime and called it a
heuristic calibrated on one point. It is the case `m = d` of Theorem 1:
the second term vanishes and `m(m+1)/2` **is** `d(d+1)/2` — the
coincidence that made the saturated flex count look like the scene gauge
is just `m = d` there. It is the gauge of the *profile space*, not of the
scene, and the two agree only in that one regime.

## 7. Theorem 2 — the threshold, necessity and sufficiency

Let

    N(d, R) := ceil( [ dR + m(m+1)/2 - d(d+1)/2 ] / (m - d) ),   m = R - 1 > d.

### 7a. Necessity `[PROVED under (G2)]`

> If the scene is infinitesimally rigid, then `n >= N(d, R)`.

*Proof.* By (**), `nullity = dim ker Psi >= dim domain - dim target`,
that is

    nullity >= d(n + R) + dim F - mn.

Rigidity means `nullity = d(d+1)/2`. With `dim F = m(m+1)/2` this forces
`n(m - d) >= dR + m(m+1)/2 - d(d+1)/2`, i.e. `n >= N(d,R)`. ∎

Note what this does **not** need: no genericity of `Psi`, no curvature
argument, no claim about attainment. It is a hard lower bound on how
many targets any rigid configuration must carry.

### 7b. Sufficiency `[CONJECTURED]`

> For `R >= d + 2`, `d >= 2` and generic scenes with `n >= N(d,R)`,
> `nullity = d(d+1)/2`.

Equivalently: `Psi` attains its maximal possible rank generically. This
is the step rigidity theory always gets stuck on — the gap between
"the count permits rigidity" and "rigidity happens", the same gap that
separates Maxwell counting from Laman-type theorems. It is **not**
supplied by anything above.

Two things are known about it here, and they point in opposite
directions, which is why the tag is `[CONJECTURED]`.

**In favour.** The count is attained *exactly*, with no slack, in every
`(d, R)` cell ever run — the ten of rounds 1-3, and the five that round 4
committed in advance (marked ✎, `t1_g4c_predictions_round4.json`):

| `d` | `R` | `m - d` | `N(d,R)` | measured | |
|---|---|---|---|---|---|
| 2 | 4 | 1 | 11 | 11 | |
| 2 | 5 | 2 | 9 | 9 | |
| 2 | 6 | 3 | 8 | 8 | |
| 2 | 7 | 4 | 8 | 8 | ✎ |
| 2 | 8 | 5 | 9 | 9 | |
| 3 | 5 | 1 | 19 | 19 | |
| 3 | 6 | 2 | 14 | 14 | |
| 3 | 7 | 3 | 12 | 12 | ✎ |
| 3 | 8 | 4 | 12 | 12 | |
| 4 | 6 | 1 | 29 | 29 | |
| 4 | 7 | 2 | 20 | 20 | ✎ |
| 4 | 8 | 3 | 17 | 17 | |
| 5 | 7 | 1 | 41 | 41 | |
| 5 | 8 | 2 | 27 | 27 | ✎ |
| 6 | 8 | 1 | 55 | 55 | ✎ |

Not one cell needs a target more than the count allows. Two of the
preregistered cells — `(4,7)` and `(5,8)` — have non-exact divisions
(`39/2`, `53/2`), so their agreement tests something the exact ones
cannot: that the bound is attained at the *first integer above* it, not
merely somewhere above it.

**Against, and it is instructive.** Sufficiency is **false at `d = 1`**.
There `m = R - 1 > 1` for every `R >= 3`, the necessary condition is met
with room to spare, and the configuration is nonetheless never rigid at
any `(n, R)` — measured, and Lemma 4f of the main document. The reason is
exactly a failure of generic rank: in one dimension `phi_r(x) = |x - p_r|`
is piecewise linear, so `w` is affine on each cell between observers, the
profile "surface" is a polyline of exactly zero curvature, and per-cell
slopes trade off against spacings to produce flexes the count cannot see.

So the curvature of the profile surface is not decoration in the
argument — it is the entire content of the sufficiency step, and `d = 1`
is the standing proof that the count alone does not suffice.

### 7c. The closed form at `R = d + 2` `[PROVED, given 7a]`

At `R = d + 2` we have `m - d = 1` and the ceiling is vacuous:

    N(d, d+2) = d(d+2) + (d+1)(d+2)/2 - d(d+1)/2
              = d^2 + 2d + (d+1)
              = d^2 + 3d + 1.

Round 3 fitted `d^2 + 3d + 1` to the three thresholds `11, 19, 29` and
labelled it *speculative*, noting that three points determine a quadratic
so the fit carried no information. That tag was correct then and is now
superseded: the quadratic is the `m - d = 1` case of the count, where the
division is trivial. The out-of-sample confirmation at `d = 5`
(predicted 41, measured 41) was, in hindsight, a test of the count.

## 8. What is checked numerically, and how

`experiments/theory/t1_g4c_proof_check.py` verifies the *lemmas*, not
just the conclusion:

1. **Lemma A** — the flex space computed from `D` directly and the space
   computed as `L^{-1}(F)` agree as subspaces, not merely in dimension.
2. **Lemma B** — `dim F` measured against `n(m-q) + q(q+1)/2` with `q`
   the measured affine span.
3. **Lemma C** — `nullity = dim ker L + dim(Im L ∩ F)` on every cell,
   rigid and flexible alike.
4. **Theorem 1** — the closed form against measured nullity for every
   `(d, R, n)` with `R <= d + 1`.
5. **Theorem 2a** — that no measured rigid cell violates the lower bound.
6. **The `d = 1` counterexample to sufficiency** — the necessary
   condition holds while rigidity fails.

## 9. Status

| statement | tag |
|---|---|
| Lemma A (reduction to profile-cloud rigidity) | `[PROVED]` |
| Lemma B (complete-graph flex dimension) | `[PROVED]` |
| Lemma C (decomposition, and `nullity = dim ker Psi`) | `[PROVED]` |
| (G1), (G2) genericity | `[PROVABLE]` |
| Theorem 1 (`R <= d+1` closed form) | `[PROVED]` under (G1), (G2) |
| Theorem 2a (threshold necessity) | `[PROVED]` under (G2) |
| Theorem 2b (threshold sufficiency) | `[CONJECTURED]`, 15/15 cells |
| `R >= d + 2` as *the* rigidity threshold | `[MEASURED]`, `d = 1..6` |

What this changes relative to v1.2: the counting law and the threshold
formula had four and one confirmations respectively and **no
derivation**. They now have one. The regime law's *necessary* half is a
theorem in every dimension at once, which is what a general-`d` result
was supposed to buy.

What it does not change: the scope. Everything here is infinitesimal
rigidity in the exact model. Global uniqueness, noisy or
`delta`-quantized `D`, and `D` harvested through the instrument from a
measured causal set are all still outside, and no amount of further
dimensions would bring them in.

**The one open step is 7b**, and it is not a matter of more measurement.
`d = 1` shows the count is not self-sufficient, so a proof has to use
curvature of the profile surface — plausibly by exhibiting, at some
convenient configuration, a `Psi` of maximal rank and invoking lower
semicontinuity of rank, which is the standard route and the reason the
tag is `[CONJECTURED]` rather than `[PROVABLE]`: the convenient
configuration has to be constructed, and in `d >= 2` it must be one
where the profile surface genuinely curves.
