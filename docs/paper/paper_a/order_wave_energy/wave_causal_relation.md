# The causal relations of the wave and of the background, written out

2026-11-26. **Proof note: the exact causal relation of a Rosen plane wave,
its form for the circular family, the vacuum condition, and the
background's conformal chart.**

- **Why this note exists.** [one_point_strip_bound.md](one_point_strip_bound.md)
  proves the one-point bound from the wave's exact causal relation, in the
  form `2 Delta v >= Delta y^T (mI + E)^-1 Delta y`. Its fresh-context AI
  referee named that imported relation the step it trusts least
  ([independent_audit/one_point_bound/review.md](independent_audit/one_point_bound/review.md)):
  every constant and every check assumes it.
  [quantized_order_bound_r2.md](quantized_order_bound_r2.md) Section 8
  lists three related imports:
  - the wave relation and the form of `E_lambda`;
  - the background relation (2);
  - the `a^4` density in the `(T, z)` chart.
- **What supported them before.**
  - [design.md](design.md) Section 1 states the Rosen Ricci formula.
  - [full_order_resolution_probe.md](full_order_resolution_probe.md)
    states the family and the matrix `M`.
  - The audit referee's report
    ([proof1_lower_bound/report.md](independent_audit/proof1_lower_bound/report.md)
    Section 3.1) sketches the relation in four lines.
  - Its script `check_a_exact_family.py` checks the Ricci tensor, the
    inverse, the null geodesics and a numerical minimisation.
  - No document wrote the "if and only if" with the time orientation,
    curves that pause in `u`, and the domain stated.
- **What this note does.**
  - Theorem 2 proves the causal relation of any Rosen plane wave.
    Lemma 1 fixes the time orientation it rests on.
  - Lemma 3 puts the circular family into the form the other notes use.
  - Lemma 4 derives the Ricci tensor and the vacuum condition.
  - Lemma 5 proves the background relation and density through a global
    conformal chart.
  - Section 8 lists deterministic checks. The first run failed three
    curvature checks: a finite-difference step that ignored the pole of
    `f` left too much truncation error near it. After the step was fixed,
    all 21 pass.
- **Authorship and review.** AI-written. Reviewed by an AI referee in a
  fresh context:
  [independent_audit/wave_relation/review.md](independent_audit/wave_relation/review.md).
  - **Mathematics:** CORRECT WITH MINOR FIXES. Every step was re-derived,
    and none was false.
  - **Checks and runs:** CORRECT WITH MINOR FIXES. The run-1 mechanism
    was misdescribed, two numbers were misquoted, and several checks
    cannot fail.

  Section 9 lists the corrections.
- **What it is not.** It is not human-reviewed. The results are standard
  in plane-wave geometry. The note's job is to write them in this
  programme's notation, with every hypothesis visible. It claims no
  novelty.

**Notation warning.** Here, as in the probe, the oracle note and the
referee reports, `c = cosh(2 beta)`, `s = sinh(2 beta)` and `t = s/c`.
In one_point_strip_bound.md `s` means `sqrt(q)` instead. `t` means the
same in both.

## 1. Setting

Let `I` be an open interval and `gamma : I -> R^(2x2)` a `C^2` map into
symmetric positive-definite matrices. The **Rosen plane wave** is
`S = I x R x R^2`, with coordinates `(u, v, x)` and

\[
g=-2\,du\,dv+\gamma_{ij}(u)\,dx^i dx^j .
\]

`Xi = d_u + d_v` has `g(Xi, Xi) = -2`, so it is timelike; it is declared
future-directed. A nonzero causal vector `w` is future-directed iff
`g(w, Xi) < 0`.

A **future-directed causal curve** is a locally absolutely continuous
curve `C : [0, 1] -> S` whose velocity is, at almost every parameter
value, zero or a future-directed causal vector. Piecewise-`C^1` curves are
included. Write `X <= Y` (`X` causally precedes `Y`) if `X = Y` or such a
curve runs from `X` to `Y`. This is the relation `J^+`.
- **The model's order.** The model's order is taken to be `J^+` on
  distinct sampled events. Every user document writes its relation with
  `>=`, which is `J^+`. This is a definition, listed in Section 7.
- **`J^+` and the chronological relation agree almost surely** on sampled
  events (written out after review, finding F1).
  - **Strict inequality gives a timelike curve.** If (a) of Theorem 2 holds
    strictly, the curve built in its proof has `eps > 0` and
    `g(C', C') = -2 eps / Delta u < 0`, so `Y` is in `I^+(X)`.
  - **The two relations differ only on a null set.** For fixed `X`, the
    pairs in `J^+` but not in `I^+` lie in two sets. One is the graph
    `2 Delta v = Delta x^T M^-1 Delta x` over `u > u_1`. The other is the
    slice `u = u_1`. Both have measure zero.
  - **Sampled pairs avoid that set.** The two events of a pair are drawn
    from a product of absolutely continuous laws, so almost surely no pair
    lies in it.

For `u_1 < u_2` in `I`, put

\[
M(u_1,u_2)=\int_{u_1}^{u_2}\gamma^{-1}(w)\,dw ,
\]

a symmetric positive-definite matrix.

## 2. Time orientation

**Lemma 1.** A vector `w = (w^u, w^v, w^x)` is future-directed causal iff
`w != 0`, `2 w^u w^v >= w_x^T gamma w_x` and `w^u + w^v > 0`. Then
`w^u >= 0` and `w^v >= 0`. If moreover `w^u = 0`, then `w_x = 0` and
`w^v > 0`.

*Proof.* `g(w, w) = -2 w^u w^v + w_x^T gamma w_x`, and
`g(w, Xi) = -(w^u + w^v)`. Causal means `g(w, w) <= 0`, which is the
inequality. Since `gamma > 0`, it gives `w^u w^v >= 0`, so `w^u` and
`w^v` do not have opposite signs. With `w^u + w^v > 0`, both are
non-negative. If `w^u = 0`, the inequality gives
`w_x^T gamma w_x <= 0`, so `w_x = 0`, and then `w^v > 0`. []

**Remark 1.1 (why this orientation; added after review, finding F2).**
`S` is connected and `Xi` is a global timelike field, so `S` has exactly
two time orientations: `Xi` and `-Xi`.
- **The model's tips select `Xi`.** Under `Xi`, Theorem 2(a) gives
  `p_- <= p_+` for `p_- = (-1/2, 0, 0)` and `p_+ = (1/2, h, 0)`, since
  `2h >= 0`.
- **The other choice is empty.** Under `-Xi` every relation reverses, so
  the diamond `{p_- <= X <= p_+}` would be empty.

The choice is part of the model's definition (Section 7).

## 3. The causal relation of a Rosen plane wave

**Theorem 2.** Let `X = (u_1, v_1, x_1)` and `Y = (u_2, v_2, x_2)` be
points of `S`. Write `Delta v = v_2 - v_1` and `Delta x = x_2 - x_1`.
Then `X <= Y` iff one of the following holds:

- (a) `u_1 < u_2` and `2 Delta v >= Delta x^T M(u_1,u_2)^-1 Delta x`;
- (b) `u_1 = u_2`, `Delta x = 0` and `Delta v >= 0`.

In particular, `X <= Y` never holds when `u_1 > u_2`.

*Proof of "only if".* Let `C(tau) = (u, v, x)(tau)` be a future-directed
causal curve from `X` to `Y`.
- **`u` does not decrease.** By Lemma 1, `u' >= 0` almost everywhere.
  `u` is absolutely continuous, so it is non-decreasing. Hence
  `u_2 >= u_1`, and `u(tau)` stays in `[u_1, u_2]`.
- **Case `u_1 = u_2`.** Then `u` is constant, so `u' = 0` almost
  everywhere. By Lemma 1, `x' = 0` and `v' >= 0` almost everywhere. So
  `Delta x = 0` and `Delta v >= 0`: this is (b).
- **Case `u_1 < u_2`.** Fix `p` in `R^2` and put

  \[
  F(\tau)=v(\tau)-p^Tx(\tau)+\tfrac12\,p^TM(u_1,u(\tau))\,p .
  \]

  `M(u_1, .)` is `C^1` on the compact interval `[u_1, u_2]`, hence
  Lipschitz there, so `F` is absolutely continuous. Almost everywhere,

  \[
  F'=v'-p^Tx'+\tfrac12\,u'\,p^T\gamma^{-1}(u)\,p .
  \]

  - Where `u' > 0`, Lemma 1 gives `v' >= x'^T gamma x' / (2u')`, and so
    \[
    F'\ \ge\ \frac{1}{2u'}\Big(x'^T\gamma x'-2u'\,p^Tx'+u'^2\,p^T\gamma^{-1}p\Big)
    =\frac{1}{2u'}\,\big|\gamma^{1/2}x'-u'\gamma^{-1/2}p\big|^2\ \ge 0 .
    \]
  - Where `u' = 0`, Lemma 1 gives `x' = 0` and `v' >= 0`, or the velocity
    is zero. Either way `F' = v' >= 0`.

  So `F(1) >= F(0)`, that is,
  `Delta v - p^T Delta x + (1/2) p^T M p >= 0` for every `p`, with
  `M = M(u_1, u_2)`. The choice `p = M^-1 Delta x` gives
  `Delta v >= (1/2) Delta x^T M^-1 Delta x`: this is (a).

*Proof of "if".*
- **Case (b).** If `Delta v = 0`, then `X = Y`. Otherwise
  `tau -> (u_1, v_1 + tau Delta v, x_1)` has velocity `Delta v d_v`. This
  is null with `g(., Xi) = -Delta v < 0`, so the curve is future-directed.
- **Case (a).** Put `p = M^-1 Delta x`, `Delta u = u_2 - u_1` and
  `eps = Delta v - (1/2) Delta x^T M^-1 Delta x >= 0`. Parametrise by
  `u` in `[u_1, u_2]`:

  \[
  x(u)=x_1+M(u_1,u)\,p,\qquad
  v(u)=v_1+\tfrac12\,p^TM(u_1,u)\,p+\varepsilon\,\frac{u-u_1}{\Delta u}.
  \]

  - Endpoints: `x(u_2) = x_2` and `v(u_2) = v_2`.
  - The velocity is `(1, v', x')`, with `x' = gamma^-1 p` and
    `v' = (1/2) p^T gamma^-1 p + eps / Delta u`.
  - Causal: `2 v' - x'^T gamma x' = 2 eps / Delta u >= 0`.
  - Future-directed: `1 + v' > 0`.

  So the curve is a `C^1` future-directed causal curve from `X` to
  `Y`. []

**Remark 2.1 (the boundary is reached by null geodesics).** When
`eps = 0`, the curve of the proof is a null geodesic with affine
parameter `u`. The Christoffel symbols of Lemma 4 make `u'' = 0`. The
transverse equation is `(gamma x')' = 0`, which holds since
`gamma x' = p`. The `v` equation is `v'' = -(1/2) x'^T gamma' x'`, which
holds since `(gamma^-1)' = -gamma^-1 gamma' gamma^-1`.
- **For `u > u_1`.** The boundary of `J^+(X)` in `S` is
  `2 Delta v = Delta x^T M^-1 Delta x`, and it is swept by these null
  geodesics from `X`.
- **At `u = u_1`.** The boundary also contains the ray
  `{u = u_1, x = x_1, v >= v_1}`, which is the null geodesic along `d_v`.
  (Added after review, finding F4.)

**Remark 2.2 (which spacetime).** Theorem 2 concerns `S`, the whole
domain of the Rosen chart. A causal curve between two points with `u` in
a subinterval stays in that subinterval. So the relation between points of
the model's window `|u| <= 1/2` is decided inside the window.

The wave can also be viewed as part of a larger plane-wave spacetime,
such as its Brinkmann form. The same conclusion then needs two
conditions:
- `u` extends to a function that does not decrease along future causal
  curves;
- the extension's slab `u^-1([u_1, u_2])` is exactly the slab of `S`.

For the Brinkmann form both hold, since the Rosen-to-Brinkmann map takes
`S` onto `{U in I}`. This note does not use or prove any extension.
(Corrected after review, finding F3: the first version omitted the slab
condition.)

## 4. The circular family

Fix `0 < q < pi^2`, `lambda > 0`, a phase `phi` and `beta > 0`. Put
`a(u) = cos(sqrt(q) u)`, `f = a^-2`, `c = cosh(2 beta)`, `s = sinh(2 beta)`,
`t = s/c`, and

\[
\gamma_\lambda=a^2R(\theta)\,\mathrm{diag}(e^{2\beta},e^{-2\beta})\,R(\theta)^T,\qquad
\theta=\tfrac12(u/\lambda+\phi),\qquad
T(\psi)=\begin{pmatrix}\cos\psi&\sin\psi\\ \sin\psi&-\cos\psi\end{pmatrix}.
\]

`R` is the rotation matrix. Take `I = (-pi/(2 sqrt q), pi/(2 sqrt q))`,
on which `a > 0`. It contains the window `[-1/2, 1/2]` exactly because
`q < pi^2`.

**Lemma 3.**
- (i) `gamma_lam^-1 = a^-2 (c I - s T(u/lambda + phi))` and
  `det gamma_lam = a^4`.
- (ii) `M(u_1,u_2) = c m I - s J`, where `m = m(u_1,u_2) = integral f`
  and `J = integral_(u_1)^(u_2) f(w) T(w/lambda + phi) dw`.
- (iii) In the chart `x = sqrt(c) y`,
  `Delta x^T M^-1 Delta x = Delta y^T (m I + E)^-1 Delta y`, with
  `E = E_lambda(u_1,u_2) = -t J`.
- (iv) `E` is symmetric and trace-free, with eigenvalues `+/- e`, where
  `e = t |G|` and `G = e^(i phi) integral_(u_1)^(u_2) f(w) e^(iw/lambda) dw`.
- (v) The volume density is `sqrt|det g| = a^2` in `(u, v, x)`, and
  `c a^2` in `(u, v, y)`. The tips `(u, v, 0)` keep `x = y = 0`.
- (vi) `T(psi + phi) = R(phi/2) T(psi) R(phi/2)^T`. So changing the phase
  conjugates every `E` by the same rotation, which is what Remark 4 of
  one_point_strip_bound.md uses.

*Proof.*
- (i) Write `diag(e^(-2beta), e^(2beta)) = c I - s diag(1, -1)`. For any
  angle, `R(theta) diag(1,-1) R(theta)^T = T(2 theta)` (multiply out).
  Inverting `gamma` inverts `a^2` and the diagonal factor, which gives
  (i). The determinant is `a^4 e^(2beta) e^(-2beta) = a^4`.
- (ii) Integrate (i).
- (iii) `c m I - s J = c (m I - t J)`, so
  `M^-1 = c^-1 (m I + E)^-1`. Then
  `Delta x^T M^-1 Delta x = c Delta y^T c^-1 (mI + E)^-1 Delta y`.
- (iv) Each `T(psi)` is symmetric with trace 0, so `J` is too, and
  `J = [[A, B], [B, -A]]` with `A + iB = G`. Then `J^2 = |G|^2 I`.
- (v) `|det g| = det gamma` because the `(u, v)` block has determinant
  `-1`. The map `x = sqrt(c) y` has Jacobian `c`.
- (vi) `T(psi) = R(psi/2) diag(1,-1) R(psi/2)^T`, and rotations compose
  by adding angles. []

## 5. Curvature and the vacuum condition

**Lemma 4.** For the Rosen metric of Section 1, the only non-zero
Christoffel symbols are

\[
\Gamma^v_{ij}=\tfrac12\gamma'_{ij},\qquad
\Gamma^i_{uj}=\Gamma^i_{ju}=\tfrac12\,(\gamma^{-1}\gamma')^i{}_j .
\]

The only non-zero Ricci component is

\[
R_{uu}=-\tfrac12\operatorname{tr}P'-\tfrac14\operatorname{tr}P^2,\qquad P=\gamma^{-1}\gamma'.
\]

For `gamma = a^2 Sigma` with `det Sigma = 1`, this becomes
`R_uu = -2 a''/a - (1/4) tr Q^2`, with `Q = Sigma^-1 Sigma'`. For the
circular family, `tr Q^2 = 2 s^2/lambda^2`, so

\[
R_{uu}=2q-\frac{s^2}{2\lambda^2}.
\]

So the family is vacuum exactly when `s = sinh(2 beta) = 2 lambda sqrt q`.
The background `gamma = a^2 I` has `R_uu = 2q`.

*Proof.*
- **Christoffel symbols.** The only non-constant components are
  `g_ij = gamma_ij(u)`, with inverse components `g^uv = -1` and
  `g^ij = gamma^ij`.
  - `Gamma^u_bc` contains only derivatives of `g_vc`, which are zero.
  - `Gamma^v_bc = (1/2) d_u g_bc`, which is non-zero only for
    transverse `b, c`.
  - `Gamma^i_bc` is non-zero only for `{b, c} = {u, j}`.
- **The trace.** `Gamma^a_ab` is non-zero only for `b = u`, where it
  equals `(1/2) tr P`.
- **Ricci components.** Use
  `R_bd = d_a Gamma^a_bd - d_d Gamma^a_ab + Gamma^a_ae Gamma^e_bd - Gamma^a_de Gamma^e_ab`.
  - These Christoffel symbols are zero:
    - every `Gamma^u_bc`;
    - every `Gamma^a_vb`;
    - every `Gamma^a_uu`;
    - every `Gamma^k_ij` with all indices transverse.
  - Nothing depends on `v` or `x`, so only `d_u` derivatives survive.
    `d_u Gamma^a_bd` is non-zero only for `a = v` with `b, d` transverse,
    or for `a = i` with `{b, d} = {u, j}`.
  - Every term in `R_ij`, `R_ui`, `R_uv` and `R_vv` contains one of these
    zero factors. Some vanish only because every `Gamma^u` is zero, for
    example `Gamma^k_ju Gamma^u_ki` and `d_u Gamma^u_ij` in `R_ij`. So
    these components vanish. (The list of zero factors was completed after
    review, finding F5.)
  - `R_uu = -d_u Gamma^a_au - Gamma^i_uj Gamma^j_ui`, which is the
    formula above.
- **Reduction for `gamma = a^2 Sigma`.** Here
  `P = 2 (a'/a) I + Q`, and `tr Q = (log det Sigma)' = 0`. So
  `tr P' = 4 (a'/a)'` and `tr P^2 = 8 (a'/a)^2 + tr Q^2`. Then use
  `(a'/a)' + (a'/a)^2 = a''/a`.
- **The circular family.** Write `Sigma = R D R^T`, with
  `D = diag(d, 1/d)`, `d = e^(2beta)` and `theta' = 1/(2 lambda)`.
  - `R' = R K theta'`, with `K = [[0, -1], [1, 0]]`.
  - So `Q = theta' R (D^-1 K D - K) R^T`.
  - `D^-1 K D - K = [[0, 1 - d^-2], [d^2 - 1, 0]]`.
  - Its square has trace `2 (d - 1/d)^2 = 8 s^2`.
  - So `tr Q^2 = 8 s^2 theta'^2 = 2 s^2/lambda^2`.
  - Finally, `a = cos(sqrt(q) u)` gives `a''/a = -q`. []

This agrees with [design.md](design.md) (1), where `Q = s^2 theta'^2`,
and with the audit referee's `2q - sinh^2(2 beta)/(2 lambda^2)`. Theorem 2
and Lemma 3 do not use the vacuum condition. They hold for every `beta`.

## 6. The background relation and its density

**Lemma 5.** Let `gamma = a^2 I` on `I`. Put
`m(u) = integral_(-1/2)^u a^-2`, `T = (m + v)/sqrt 2`,
`Z = (m - v)/sqrt 2` and `z = (Z, x)`.
- (i) `Phi : (u, v, x) -> (T, Z, x)` is a diffeomorphism of `S` onto
  `R^4`, and `g_0 = a^2 (-dT^2 + dZ^2 + |dx|^2)`.
- (ii) `X <=_0 Y` iff `T_Y - T_X >= |z_Y - z_X|`. This is (2) of
  quantized_order_bound_r2.md.
- (iii) The volume density is `a^2 du dv d^2x = a^4 dT dZ d^2x`.

*Proof.*
- (i) `m' = a^-2 > 0`, and `m(u) = (tan(sqrt(q) u) + tan(sqrt(q)/2))/sqrt q`
  tends to `+/- infinity` at the ends of `I`. So `m` maps `I` onto `R`,
  and `Phi` is a diffeomorphism onto `R^4`. Since
  `dT^2 - dZ^2 = 2 dm dv = 2 a^-2 du dv`, the metric
  `a^2 (-dT^2 + dZ^2 + |dx|^2)` equals `-2 du dv + a^2 |dx|^2`.
- (ii) `(S, g_0)` is conformal, through `Phi`, to all of Minkowski
  space. A conformal factor does not change which vectors are causal.
  Time orientations agree, since `dT(d_u + d_v) = (a^-2 + 1)/sqrt 2 > 0`.
  So the causal curves correspond and the relation is Minkowski's.
- (iii) The Jacobian of `(u, v) -> (T, Z)` has absolute value `m' = a^-2`.

**Cross-check against Theorem 2.** With `gamma = a^2 I`, `M = Delta m I`,
and `Delta T^2 - Delta Z^2 = 2 Delta m Delta v`.
- Theorem 2(a) gives `Delta v >= 0`, hence `Delta T > 0` and
  `Delta T^2 >= Delta Z^2 + |Delta x|^2`.
- Conversely, `Delta T >= |Delta z|` forces `Delta m >= 0` and
  `Delta v >= 0`. If `Delta m = 0`, it forces `Delta x = 0`, which is
  Theorem 2(b).

So the two descriptions coincide. []

**Remark 5.1 (what r2 also asserts; added after review, finding F7).**
- **The transverse coordinate.** quantized_order_bound_r2.md writes (2) in
  the wave's `y` chart. In the coupling, the reference model is placed in
  the wave's `y` coordinates. Its transverse coordinate is `y`, with the
  same tips (coordinate_oracle_barrier.md Section 1). Lemma 5 is stated in
  `x` for the metric `-2 du dv + a^2 |dx|^2`, and it applies with `x`
  renamed `y`. This identification is part of the model's definition
  (Section 7).
- **The density near the closed diamond.** r2 also asserts that the
  normalised density is positive, bounded and `C^1` there, with bounded
  `d rho / dT`.
  - The density is proportional to `a^4` with `a = cos(sqrt(q) u(m))`,
    and `m = (T + Z)/sqrt 2`.
  - So `d rho / dT` is proportional to `(1/sqrt 2) 4 a^3 a' du/dm`, which
    is `2 sqrt 2 a^5 a'`.
  - On the closed diamond `|u| <= 1/2`, so `a` lies between
    `cos(sqrt(q)/2) > 0` and 1. Every one of these quantities is therefore
    bounded.

## 7. What is now written, and what remains

- **Written here:**
  - the causal relation of the wave, for every phase, `lambda` and
    `beta`;
  - the form `mI + E` and the eigenvalues of `E`;
  - the volume densities;
  - the vacuum condition;
  - the background relation (2) and its `a^4` density.
- **Still not proved, because they are definitions:**
  - the sampling model (i.i.d. events from the normalised volume of the
    diamond between the tips);
  - the order map and the unlabelled order
    ([full_order_resolution_probe.md](full_order_resolution_probe.md));
  - that the model's order is `J^+` (Section 1);
  - the time orientation `Xi`, which the tips select (Remark 1.1);
  - that the coupling places the reference model in the wave's `y`
    chart, with the same tips (Remark 5.1).

  The last three were added after review, finding F13.
- **Standard background used** (list completed after review, F13):
  - an absolutely continuous function with non-negative derivative almost
    everywhere is non-decreasing;
  - the chain rule for a `C^1` function of an absolutely continuous one;
  - a diffeomorphism maps absolutely continuous curves to absolutely
    continuous curves;
  - Minkowski space's `J^+` is `Delta T >= |Delta z|`, for the same class
    of curves;
  - a positive conformal factor preserves causal character.

  Section 1's almost-sure statement uses only the explicit timelike curve
  and a measure-zero graph, not `I^+ = int J^+`.

## 8. Deterministic checks

Script: [independent_audit/wave_relation/checks.py](independent_audit/wave_relation/checks.py).
No sampling of orders. Seeded random draws are used only to test
inequalities and identities. The criteria are in the script's docstring
and were fixed before the full run.

**How the checks are built.** The curvature (W1) and geodesic (W3) checks
do not use Lemma 4's Christoffel symbols. They take fourth-order finite
differences of the 4x4 metric in all four coordinates and apply the
generic formulas.

**Run 1 failed** (script at `3bb1512`; output
[checks_output_run1.txt](independent_audit/wave_relation/checks_output_run1.txt)).
- **What failed.** 18 checks passed and 3 failed, all in W1 (Ricci). The
  worst case was q = 9, `lambda = 1`, `u = 0.45`, with `|R|/scale = 1.06`.
- **The cause** (corrected after review, finding F8). The failure is
  truncation error of the outer finite difference near the pole of `f`
  at `u = pi/(2 sqrt q)`. The step in `u` was `lambda/40` regardless of
  that pole.
  - At the worst case the outer stencil reaches `2h = 0.05`. The
    distance to the pole is 0.074.
  - The first version said the stencil "reached past the pole". That is
    wrong.
    - The metric `cos^2(sqrt(q) u) Sigma` is entire, so evaluating it
      past the pole is harmless.
    - `gamma^-1` is evaluated only at Christoffel centres, all at
      `u <= 0.50 < 0.524`.
    - The referee reproduced the failure (1.058 against 1.061) with the
      same outer stencil applied to the exact Christoffel symbols.
  - Commit `3d3f9ca`'s message and the script's Run history repeat the
    wrong mechanism. The script's docstring now carries a correction line.
- **The diagnosis.** The error converged as the step shrank. At the
  worst case, successive halvings cut it by 28.8, 17.9 and 16.4 times,
  approaching the 16 of fourth-order truncation. The values converged to
  the predicted ones: 0 for the vacuum and `2q` for the background. The
  defect was in the check, not in the formula.

**Run 2 passed** (step rule fixed at `3d3f9ca`; output
[checks_output.txt](independent_audit/wave_relation/checks_output.txt)).
- **The change.** The step in `u` became `min(lambda/80, d/200)` for the
  wave and `min(1e-3, d/200)` for the background, with `d` the distance
  from `u` to the pole. No tolerance was changed.
- **How the rule was chosen.** By trying two candidates on the W1 cases.
  The rejected `d/100` also passes. The referee re-ran the trials, whose
  outputs had not been kept, and added `d/400` (finding F10).

  | step rule | vacuum | control | background |
  |---|---|---|---|
  | `lambda/40` (run 1) | 1.06 | 0.967 | `2.8e-6` |
  | `min(lambda/80, d/100)` | `3.8e-7` | `3.5e-7` | `8.2e-7` |
  | `min(lambda/80, d/200)` (run 2) | `2.44e-8` | `2.20e-8` | `5.1e-8` |
  | `min(lambda/80, d/400)` (referee) | `1.3e-9` | `1.5e-9` | `3.3e-9` |

  The tolerances are `1e-5` for the vacuum and control and `1e-6` for the
  background.
- **All 21 checks pass** (numbers corrected after review, F11):
  - **W1:** the vacuum Ricci tensor is at most `2.44e-8` of the scale. The
    positive control (`s -> 1.1 s`) reproduces
    `R_uu = 2q - s^2/(2 lambda^2)` to `2.20e-8`. The background has
    `R_uu = 2q` to `5.1e-8`.
  - **W2:** the family algebra agrees to `1.5e-13` or better.
  - **W3:** 12 numerically integrated null geodesics end on
    `2 Delta v = Delta x^T M^-1 Delta x` to `4.9e-8`.
  - **W4:** 1600 random fixed-endpoint paths never beat the bound.
  - **W5:** the conformal chart holds to `2e-15`. Theorem 2 and the
    Minkowski relation classify 78,904 random pairs identically.
    Case-(b) pairs and near-null pairs are included.
  - **W6:** the x-chart and y-chart relations agree on 300 near-null
    pairs.

**What the checks do not test** (completed after review, finding F12).
- **Checks that cannot fail on the mathematics.** They can fail only
  through a code inconsistency.
  - **W4, "the optimal path attains it".** On a common quadrature it is
    the identity `p^T M p = Delta x^T M^-1 Delta x`.
  - **W4, "paths never beat the bound".** The ratio equals
    `1 + integral delta'^T gamma delta' / act` exactly: Cauchy-Schwarz
    with the remainder built in. The cross term is at most `1.3e-13`.
    The smallest ratio, 1.0106, reflects only the size of the
    perturbations.
  - **W2, "M = c m I - s J".** Both sides use the same quadrature nodes,
    so this is the inverse check summed with weights.
  - **W3, "u stays affine".** `Gamma^u` comes from finite differences of
    constant entries.
  - **W4, the Lemma 1 sign test.** It is a tautology, and it never tests
    the case `w^u = 0`.
  - **W6** repeats W2's identity on pairs near the boundary.
- **So no check in this script that could fail tests the "only if"
  direction for the wave.** That direction rests on the proof of
  Theorem 2. The referee supplied a test that could fail (review,
  Section 4, R2):
  - 600 random piecewise-linear causal curves, each segment certified
    causal;
  - 1007 of the segments pause in `u`, and 373 have zero velocity;
  - both a non-vacuum `beta` and the background are covered.

  None violates Theorem 2, and none decreases `u`.
- **The checks that could fail did test something.**
  - W1 recomputes the curvature from finite differences of the metric.
  - W3 recomputes the geodesics.
  - W5 compares Theorem 2 with the Minkowski relation.
- **Every check is at grid points or random draws.** The statements rest
  on the proofs of Sections 2 to 6.
## 9. Review and the corrections made after it

[independent_audit/wave_relation/review.md](independent_audit/wave_relation/review.md).
An AI referee in a fresh context reviewed this note as committed at
`1ba2e6e`.
- It also reviewed the script and both runs.
- It re-ran `checks.py`, and the output matched line for line.
- It derived the Ricci tensor exactly, by a different method.
- It tested Theorem 2 directly on piecewise-linear causal curves,
  including ones that pause in `u`.
- Its scripts and outputs are archived in
  [independent_audit/wave_relation/referee/](independent_audit/wave_relation/referee/).

**Verdicts.**
- **Mathematics: CORRECT WITH MINOR FIXES.**
  - Lemma 1, Theorem 2 (both directions), Remark 2.1 and Lemmas 3, 4 and 5
    were re-derived and found correct.
  - The class of curves is adequate.
  - The conformal argument of Lemma 5 is complete.
- **Script and runs: CORRECT WITH MINOR FIXES.**
  - Fixing the step rule after run 1 was legitimate. No tolerance
    changed, and the rejected candidate also passes.
  - The note's account of run 1's mechanism was wrong.

**Corrected after the review.** No statement of a theorem or lemma changed.

| finding | correction |
|---|---|
| F1 | Section 1: the almost-sure agreement of `J^+` and the chronological relation on sampled events is written out |
| F2 | Remark 1.1: the two time orientations, and why the tips select `Xi` |
| F3 | Remark 2.2: the slab condition for an extension |
| F4 | Remark 2.1: the boundary ray at `u = u_1` |
| F5 | Lemma 4's proof lists every zero factor, including all `Gamma^u` |
| F6 | notation: the vector field is `Xi`, the curve `C`, the unimodular matrix `Sigma`. `T` stays the matrix in Section 4 and the coordinate in Section 6; `Q` here is a matrix, not design.md's scalar |
| F7 | Remark 5.1: the coupling's `y` identification, and the density's bounds near the closed diamond |
| F8 | Section 8: run 1 failed through truncation error near the pole, not by crossing it; the script's docstring has a correction line |
| F9 | Section 8: the halving factors 28.8, 17.9, 16.4 |
| F10 | Section 8: the step-rule table, including the rejected `d/100`, which also passes |
| F11 | Section 8: `2.44e-8` and `1.5e-13`; draws test identities too |
| F12 | Section 8: every check that cannot fail on the mathematics is named; the referee's piecewise-linear test is credited for the "only if" direction |
| F13 | Section 7: the definitions (the order is `J^+`, the orientation, the `y` identification) and the standard facts used |

**The step the referee trusts least.** It was Section 1's sentence that the
model's order is `J^+` and agrees almost surely with the chronological
relation. That is now written out (F1). A close second was the `y`
identification, now Remark 5.1. Both are statements about the model's
definition, not about the geometry.

**Not checked by the referee:**
- whether any code in the programme samples orders with `I^+` rather
  than `J^+`;
- extensions beyond the Rosen chart;
- the rest of one_point_strip_bound.md and quantized_order_bound_r2.md;
- the earlier `check_a_exact_family.py`;
- an independent reproduction of W3;
- whether its piecewise-linear minimisation converges to 1 at high
  frequency. It relies on the explicit curve for that.