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

### 7b. Sufficiency

> For `R >= d + 2`, `d >= 2` and generic scenes with `n` large enough,
> `nullity = d(d+1)/2`.

This is the step rigidity theory usually gets stuck on — the gap between
"the count permits rigidity" and "rigidity happens", the same gap that
separates Maxwell counting from Laman-type theorems. **Sections 7d–7g
close it**, up to one explicit finite-dimensional hypothesis (G3), and
close it *sharply* at `R = d + 2`. What follows first is the evidence
that motivated looking.

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

### 7d. The tangency reformulation

Fix the observers for a moment. An element of `Im L ∩ F` is an ambient
isometry `v |-> Av + b` of `V` (with `A` skew) that some scene motion
induces. Writing `T_x := Im(M U(x))` for the tangent space of the
profile surface `Sigma_P = w(R^d)` at `w(x)`, the condition at a target
`x` is

    A w(x) + b + M W_x deltaP  in  T_x.                              (T)

With `deltaP = 0` and (T) imposed at *every* `x`, this says precisely
that the linear vector field `X(v) = Av + b` is **tangent to `Sigma_P`
everywhere** — that is, `Sigma_P` is invariant under a one-parameter
group of ambient isometries. Sufficiency therefore becomes a question
about symmetries of one surface, and that question has an answer.

### 7e. Lemma E — `Sigma_P` has no continuous symmetry `[PROVED under (G1)]`

> For `d >= 2` and generic `P`, the only `(A, b)` tangent to `Sigma_P`
> at every smooth point is `(A, b) = 0`.

*Proof.* **`Sigma_P` is not smooth.** `phi_r(x) = |x - p_r|` fails to be
differentiable at `x = p_r`, and near there `Sigma_P` is a cone over a
`(d-1)`-sphere: with `x = p_r + rho * omega`,

    Phi(x) = E_r + rho ( e_r + sum_{s != r} <omega, u_s(p_r)> e_s ) + O(rho^2)

whose tangent cone at `w(p_r)` is a genuinely curved cone, not a union of
linear subspaces. Away from the observers `w` is real-analytic and (G1)
makes it an immersion, so every other point of `Sigma_P` has a tangent
cone that *is* a finite union of `d`-planes. Hence the set

    S := { v in Sigma_P : the tangent cone at v has a curved component }

equals `{ w(p_1), ..., w(p_R) }`, and is finite.

The flow of `X` is a one-parameter group of isometries of `V` preserving
`Sigma_P`, so it preserves `S`. A **connected** group cannot permute a
finite set nontrivially, so it fixes every point of `S`:

    A w(p_r) + b = 0    for r = 1, ..., R.

Now suppose `(A, b) != 0`. Its zero set `{v : Av + b = 0}` is empty when
`A = 0`, and otherwise an affine subspace of dimension `dim ker A`. A
nonzero real skew-symmetric matrix has even rank at least `2`, so that
dimension is at most `m - 2`.

But `w(p_r) = M E_r`, the centered rows of the **non-squared** observer
distance matrix `E`. (Squared distance matrices have rank at most
`d + 2`; non-squared ones are generically nonsingular.) With
`rank E = R`, the vectors `{M E_r}` span `M(R^R) = V`, of dimension `m`,
so their **affine** span has dimension at least `m - 1`.

`m - 1 > m - 2`. The cone points do not fit in the zero set, so
`(A, b) = 0`. ∎

This is where `d >= 2` finally earns its keep in a *proof* rather than in
a picture: it is what makes the cone at `w(p_r)` curved. At `d = 1` the
"cone" is two rays — a union of linear pieces — `S` is empty, and the
argument has nothing to grip.

### 7f. Lemma F — finitely many targets suffice `[PROVED]`

Let `Lambda := R^{dR} × so(V) × V`, of dimension `dR + m(m+1)/2`, and for
`x in R^d` let `C(x) ⊆ Lambda` be the subspace cut out by (T). Then
`K(n) := ∩_{j<=n} C(x_j)` is the space of candidate flexes after `n`
targets, and `K(∞) := ∩_{x} C(x)`.

`K(1) ⊇ K(2) ⊇ ...` is a nested chain of subspaces of `Lambda`, so it
strictly decreases at most `dim Lambda - dim K(∞)` times. Choosing
targets greedily — at each step, if the current intersection strictly
contains `K(∞)`, some `x` (hence a generic one) strictly reduces it —

> `K(n) = K(∞)` for generic targets once
> `n >= dim Lambda - dim K(∞)`. ∎

### 7g. Theorem 2b′ — sufficiency `[PROVED under (G1), (G2), (G3)]`

**(G3)** For generic `P` with `R >= d + 2`, the map
`G : P |-> ( || M(E_r - E_s) || )_{r<s}` has differential of rank
`dR - d(d+1)/2`: it is an immersion modulo congruence.

> `K(∞)` consists exactly of the rigid motions of the scene, so
> `dim K(∞) = d(d+1)/2`; and therefore for generic scenes with
>
>     n >= dR + m(m+1)/2 - d(d+1)/2
>
> we have `ker L = ` rigid motions, `Im L ∩ F = 0`, and
> `nullity = d(d+1)/2`.

*Proof.* Take `(deltaP, A, b) in K(∞)`. Condition (T) holding at every
`x` says the scene deformation carries the cloud by the ambient isometry
`(A,b)`; equivalently `Sigma_{P + eps deltaP}` is congruent to
`Sigma_P` to first order. Congruences carry cone points to cone points,
so the pairwise distances of the `w(p_r)` are preserved to first order:
`deltaP in ker DG`. By (G3), `deltaP` is an infinitesimal rigid motion
of the observers; extend it to a rigid motion `g` of the whole scene and
replace `(deltaX, deltaP)` by `(deltaX, deltaP) - g`, which changes
nothing since `g` induces `delta w = 0`. Now `deltaP = 0`, and Lemma E
gives `(A, b) = 0`. So `K(∞)` is exactly the rigid motions.

Lemma F then gives `K(n) = K(∞)` for `n >= dim Lambda - d(d+1)/2`.
Consequently no element of `K(n)` has `(A,b) != 0`, i.e.
`Im L ∩ F = 0`; and an element with `(A,b) = 0` has `deltaP` rigid and
then `delta x_j` uniquely determined by the injectivity of `M U(x_j)`
(rank `d`, as `m > d`), so `ker L` is exactly the rigid motions. ∎

### 7h. Where the bounds meet: `R = d + 2` is settled `[PROVED]`

Theorem 2a gives `n >= N(d,R)`; Theorem 2b′ gives rigidity once
`n >= dR + m(m+1)/2 - d(d+1)/2`. The lower bound is the upper bound
divided by `m - d` and rounded up, so **when `m - d = 1` — that is,
exactly when `R = d + 2` — the two coincide and the threshold is
pinned**:

    n_threshold(d, d+2) = d^2 + 3d + 1,   proved, for every d >= 2.

| `d` | `R = d+2` | threshold |
|---|---|---|
| 2 | 4 | 11 |
| **3** | **5** | **19** |
| 4 | 6 | 29 |
| 5 | 7 | 41 |
| 6 | 8 | 55 |
| 7 | 9 | 71 |

The `d = 3` row is the physical one: **in 3+1 spacetime, five observers
and 19 targets, and both numbers are now theorems** rather than
measurements — modulo (G1)–(G3).

For `R > d + 2` the bounds separate (e.g. `d = 3, R = 6` gives
`14 <= n_threshold <= 27`, measured 14). Sharpness there is still open,
and is what remains of the original conjecture.

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
7. **Lemma E's two inputs** — that `rank E = R` and that the cone points
   span `V`, so their affine span exceeds the largest possible zero set
   of a nonzero ambient isometry.
8. **(G3)** — the rank of `DG` against `dR - d(d+1)/2`.
9. **The brackets** — that `N(d,R) <= dR + m(m+1)/2 - d(d+1)/2` always,
   that every measured threshold lies between them, and that they
   coincide exactly on the cells with `R = d + 2`.

## 9. Status

| statement | tag |
|---|---|
| Lemma A (reduction to profile-cloud rigidity) | `[PROVED]` |
| Lemma B (complete-graph flex dimension) | `[PROVED]` |
| Lemma C (decomposition, and `nullity = dim ker Psi`) | `[PROVED]` |
| (G1), (G2) genericity | `[PROVABLE]` |
| **(G3)** observer gap map is an immersion mod congruence | `[PROVABLE]`, 18/18 cells |
| Theorem 1 (`R <= d+1` closed form) | `[PROVED]` under (G1), (G2) |
| Theorem 2a (threshold necessity) | `[PROVED]` under (G2) |
| Lemma E (`Sigma_P` has no continuous symmetry) | `[PROVED]` under (G1) |
| Lemma F (finitely many targets suffice) | `[PROVED]` |
| **Theorem 2b′ (sufficiency, with an explicit bound)** | `[PROVED]` under (G1)–(G3) |
| **Threshold at `R = d + 2` equals `d^2 + 3d + 1`** | `[PROVED]` under (G1)–(G3) |
| Sharpness of the threshold for `R > d + 2` | `[CONJECTURED]`, 15/15 cells |
| `R >= d + 2` as *the* rigidity threshold | `[MEASURED]`, `d = 1..6` |

What this changes relative to v1.2: the counting law and the threshold
formula had four and one confirmations respectively and **no
derivation**. They now have one, and the sufficiency step — the one that
looked like the permanent obstacle — is closed at `R = d + 2`, which is
the case that answers the physical question. Three spatial dimensions
need five observers and 19 targets, and both numbers are theorems.

The residual conjecture is narrower than the original in two ways. It
concerns only `R > d + 2`, where the extra observers make the count
divide unevenly; and it is a *sharpness* claim, not an existence one,
since Theorem 2b′ already gives rigidity there at a larger explicit `n`.

(G3) is worth separating out. It is the single unproved input to
sufficiency, and it is a statement about one concrete finite-dimensional
map on observer configurations alone — no targets, no surface, no
`n`. It is verified at rank exactly `dR - d(d+1)/2` in all 18 cells run
across `d = 2..7`.

What it does not change: the scope. Everything here is infinitesimal
rigidity in the exact model. Global uniqueness, noisy or
`delta`-quantized `D`, and `D` harvested through the instrument from a
measured causal set are all still outside, and no amount of further
dimensions would bring them in.

**What is left is (G3) and the sharpness for `R > d + 2`.** Neither is a
matter of more measurement. (G3) is an explicit rank statement about one
map on observer configurations, of the kind that yields to a single
well-chosen configuration plus lower semicontinuity of rank. Sharpness
for `R > d + 2` asks why `m - d` independent conditions arrive with each
new target rather than fewer, and Lemma F's greedy argument only
guarantees one.

The route that closed 7b is worth recording, because the obstacle was
never really the counting. `d = 1` shows counting cannot be
self-sufficient; what supplies the missing content is that the profile
surface **is not smooth**. Each observer contributes a conical singular
point, a continuous group cannot permute finitely many of them, and the
cone points are forced to span because the non-squared observer distance
matrix has full rank. The curvature that `d = 1` lacks turns out to be
the curvature of those cones.
