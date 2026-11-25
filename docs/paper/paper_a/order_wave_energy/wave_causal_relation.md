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
  - Section 6 lists deterministic checks.
- **Authorship and review.** AI-written. Awaiting a fresh-context AI review.
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

`T = d_u + d_v` has `g(T, T) = -2`, so it is timelike; it is declared
future-directed. A nonzero causal vector `w` is future-directed iff
`g(w, T) < 0`.

A **future-directed causal curve** is a locally absolutely continuous
curve `c : [0, 1] -> S` whose velocity is, at almost every parameter
value, zero or a future-directed causal vector. Piecewise-`C^1` curves are
included. Write `X <= Y` (`X` causally precedes `Y`) if `X = Y` or such a
curve runs from `X` to `Y`. This is the relation `J^+`. The order of the
model is this relation on distinct sampled events. Its null boundary has
measure zero, so for sampled events `J^+` and the chronological relation
agree almost surely.

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
`g(w, T) = -(w^u + w^v)`. Causal means `g(w, w) <= 0`, which is the
inequality. Since `gamma > 0`, it gives `w^u w^v >= 0`, so `w^u` and
`w^v` do not have opposite signs. With `w^u + w^v > 0`, both are
non-negative. If `w^u = 0`, the inequality gives
`w_x^T gamma w_x <= 0`, so `w_x = 0`, and then `w^v > 0`. []

## 3. The causal relation of a Rosen plane wave

**Theorem 2.** Let `X = (u_1, v_1, x_1)` and `Y = (u_2, v_2, x_2)` be
points of `S`. Write `Delta v = v_2 - v_1` and `Delta x = x_2 - x_1`.
Then `X <= Y` iff one of the following holds:

- (a) `u_1 < u_2` and `2 Delta v >= Delta x^T M(u_1,u_2)^-1 Delta x`;
- (b) `u_1 = u_2`, `Delta x = 0` and `Delta v >= 0`.

In particular, `X <= Y` never holds when `u_1 > u_2`.

*Proof of "only if".* Let `c(tau) = (u, v, x)(tau)` be a future-directed
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
  is null with `g(., T) = -Delta v < 0`, so the curve is future-directed.
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
holds since `(gamma^-1)' = -gamma^-1 gamma' gamma^-1`. So the boundary
`2 Delta v = Delta x^T M^-1 Delta x` of `J^+(X)` in `S` is swept by the
null geodesics from `X`.

**Remark 2.2 (which spacetime).** Theorem 2 concerns `S`, the whole
domain of the Rosen chart. A causal curve between two points with `u` in
a subinterval stays in that subinterval. So the relation between points of
the model's window `|u| <= 1/2` is decided inside the window. If the wave
is viewed as part of a larger plane-wave spacetime, such as its
Brinkmann form, the same holds as long as `u` extends to a function that
does not decrease along future causal curves. This note does not use or
prove any extension.

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

For `gamma = a^2 Gamma` with `det Gamma = 1`, this becomes
`R_uu = -2 a''/a - (1/4) tr Q^2`, with `Q = Gamma^-1 Gamma'`. For the
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
  - Every `Gamma^a_vb` is zero, and so is every `Gamma^a_uu`.
  - Nothing depends on `v` or `x`.
  - Every product in `R_ij`, `R_ui`, `R_uv` and `R_vv` contains one of
    these zero factors, so these components vanish.
  - `R_uu = -d_u Gamma^a_au - Gamma^i_uj Gamma^j_ui`, which is the
    formula above.
- **Reduction for `gamma = a^2 Gamma`.** Here
  `P = 2 (a'/a) I + Q`, and `tr Q = (log det Gamma)' = 0`. So
  `tr P' = 4 (a'/a)'` and `tr P^2 = 8 (a'/a)^2 + tr Q^2`. Then use
  `(a'/a)' + (a'/a)^2 = a''/a`.
- **The circular family.** Write `Gamma = R D R^T`, with
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
    ([full_order_resolution_probe.md](full_order_resolution_probe.md)).
- **Standard background used:** absolutely continuous functions and
  conformal invariance of causal structure.

## 8. Deterministic checks

Script: [independent_audit/wave_relation/checks.py](independent_audit/wave_relation/checks.py).
No sampling of orders. Seeded random draws are used only to test
inequalities. The criteria are in the script's docstring and were fixed
before the full run.

**Results.** To be filled in after the run.
