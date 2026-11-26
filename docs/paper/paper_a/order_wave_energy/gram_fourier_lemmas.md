# Argument 3's outline steps, written out

2026-11-26. **Proofs of the steps that
[gram_fourier_noise_r2.md](gram_fourier_noise_r2.md) (here "r2") left at
outline level.**

- **What it does.** r2 Section 8 lists the steps of argument 3 that were
  still at outline level, with the audit's referee reports supplying them.
  This note writes each one as a proof:
  - the index-counting proof of the noise bound (5), with the U-statistic
    moment and Hoeffding bounds it uses;
  - the canonical-map estimates (6) and (7), tips included;
  - the kernel coupling (8), uniform in the fixed endpoint, and its moment
    form (9);
  - the closeness of wave and background moments for every probe, which r2
    cited;
  - the population expansion behind the signal bound (11), including the
    two assertions r2 made there;
  - the anchor selection of r2 Section 6.

  It also re-checks the anchor-guard constants that (11) uses, in exact
  rational arithmetic (Section 8, G9).
- **Authorship and review.** AI-written. Reviewed by a fresh-context AI
  referee, and corrected after it (Section 9). No human has checked it.
- **What it is not.**
  - It adds no statistic and changes no statement of r2.
  - It is asymptotic, for one family: the exact circular vacuum plane wave,
    known `0 < q <= 1`, fixed-N i.i.d. physical-volume sampling.
  - The bracket `4/3 <= alpha_low <= alpha_high <= 2`
    ([AUDIT.md](independent_audit/AUDIT.md) Section 2) is unchanged.
  - Constants are explicit where that is cheap. They are far from any
    reachable size (AUDIT.md Section 5).

## 0. Map, imports and conventions

| r2 step (Section 8) | here |
|---|---|
| index count of (5) | Lemmas B, C; Proposition D |
| U-statistic moment and Hoeffding bounds | Lemma A |
| tip estimates (6), (7) | Lemmas E, F |
| (8), uniform over fixed endpoints in a compact interior set `K_(s_0)` | Lemmas G, H, I(a), (b), (d) |
| closeness of moments for every probe (cited in r2 (R3) and Section 5) | Lemma I(c), (e), (f); Remark I1 |
| (9), for `p = 2, 4` (the only cases r2 uses; r2 states every `p >= 2`) | Lemma J |
| population expansion behind (11), with "anchor-phase terms vanish" and "calibration errors are `o(lambda^2)`" | Lemma K, Corollary K2, Proposition L |
| anchor selection (r2 Section 6) | Lemma M |
| anchor-guard constants used by (11) | Section 8, G9 |

**Imported, proved elsewhere (AI-reviewed).**
- [WCR] = [wave_causal_relation.md](wave_causal_relation.md):
  - Theorem 2, the causal relation of a Rosen plane wave;
  - Lemma 3, the circular family in the `y` chart;
  - Lemma 5, the background's conformal chart and density.
- [QOB] = [quantized_order_bound_r2.md](quantized_order_bound_r2.md):
  - Lemma 2.1, `|beta| <= C_E lambda min(l, lambda)`;
  - Lemma 2.2 and Corollary 2.3, the sandwich;
  - Lemma 2.4, `T` is a common time function.
- From r2:
  - the raw moments (Section 1.1);
  - the calibration `Sigma` and its regular domain (Section 1.2);
  - the statistic (1a) and (1b);
  - the properties (R1) to (R4);
  - the identity T1 to T4 (Section 5);
  - the scan (Section 7).

**Imported, not proved here.**
- The definitions listed in [WCR] Section 7: the sampling model, the order
  map, the order being `J^+`, the time orientation and the `y`
  identification.
- That `Sigma` inverts the background population map on its regular
  domain ([null_coordinate_challenge.md](null_coordinate_challenge.md)
  Section 2, certified there). Remark K1 checks its pair-fraction formula
  independently.
- The probe weight and anchor tolerances of r2 Section 1.3 and
  [gram_fourier_challenge.md](gram_fourier_challenge.md) Section 5.

**Conventions.**
- `C` is a constant that may depend on `q`, on `K` where `K` appears, and
  on the fixed cutoffs, weight and tolerances. It never depends on
  `lambda`, the phase `phi`, `n`, `n_P`, `omega`, the references or the
  anchor triple. It may change from line to line.
- "Uniform" means uniform in all of these.

## 1. Setting

- **The window.** `0 < q <= 1`, `k = sqrt q`, `a(u) = cos(k u)`,
  `f = a^-2`, `f_max = sec^2(k/2) <= sec^2(1/2) < 1.2985`.
  - `m(s, t) = int_s^t f`, so `m(s, t) >= t - s`.
  - `m(u) = m(-1/2, u)` and `M = m(1/2)`.
- **The wave in the `y` chart** ([WCR] Lemma 3, with `B = -E`). For
  `s < t` put
  `beta(s, t) = (2 lambda k / c_lam) e^(i phi) int_s^t f(w) e^(iw/lambda) dw`,
  with `c_lam = sqrt(1 + 4 q lambda^2)`.
  - `B(s, t)` is the symmetric trace-free matrix with complex
    representative `beta`, that is `[[Re beta, Im beta], [Im beta, -Re beta]]`.
    Its eigenvalues are `+/- |beta|`, and `B^2 = |beta|^2 I`.
  - `M_lam(s, t) = m(s, t) I - B(s, t)`.
  - [QOB] Lemma 2.1: `|beta(s, t)| <= C_E lambda min(t - s, lambda)`, with
    `C_E = 2 k (4 f_max - 2)`.
  - The volume density is `c_lam a^2` in `(u, v, y)` for the wave and `a^2`
    for the background. The constant `c_lam` cancels on normalising.
- **Relations** ([WCR] Theorem 2). For `u_1 < u_2` put the pair slack
  `kappa_lam(X, Y) = Delta v - Delta y^T M_lam(u_1, u_2)^-1 Delta y / 2`.
  - `X <=_lam Y` iff `kappa_lam >= 0`, when `u_1 < u_2`.
  - When `u_1 = u_2`: iff `Delta y = 0` and `Delta v >= 0`.
  - Never when `u_1 > u_2`.

  The background is the case `lambda = 0`: `kappa_0 = Delta v - |Delta y|^2/(2 Delta m)`.
- **Tips and diamonds.** `p_- = (-1/2, 0, 0)`, `p_+ = (1/2, 1, 0)`.
  - `D_lam = {X : p_- <=_lam X <=_lam p_+}`.
  - For `|u| < 1/2` put
    `eta_lam(X) = kappa_lam(p_-, X) = v - y^T M_lam(-1/2, u)^-1 y / 2` and
    `zeta_lam(X) = kappa_lam(X, p_+) = 1 - v - y^T M_lam(u, 1/2)^-1 y / 2`.
  - By Theorem 2, a point with `|u| < 1/2` is in `D_lam` iff
    `eta_lam >= 0` and `zeta_lam >= 0`. The slices `u = +/- 1/2` have
    measure zero.
- **The Minkowski chart** ([WCR] Lemma 5, [QOB] Section 1).
  - `T = (m + v)/sqrt 2`, `Z = (m - v)/sqrt 2`, `z = (Z, y)`.
  - The background relation is `T_Y - T_X >= |z_Y - z_X|`.
  - Put `gt(X, Y) = |T_Y - T_X| - |z_Y - z_X|`. This is the null margin
    without orientation. Changing `X` or `Y` by a Euclidean displacement
    `d` in `(T, Z, y)` changes `gt` by at most `sqrt 2 d`.
- **Tip distance.** `s(X) = min(u + 1/2, 1/2 - u)` and
  `e(s) = C_E lambda min(1, lambda/s)`.
- **Standing smallness.** `lambda <= lambda_1 = 1/(4 C_E)`, so
  `C_E lambda <= 1/4`. At `q = 1`, `lambda_1 = 0.039`.

## 2. Moments and concentration of averages and U-statistics

References `R_1, ..., R_n` are i.i.d. A **one-reference coordinate** is
`eps = (1/n) sum_j kappa(R_j)` with `E kappa(R) = 0` and `|kappa| <= 1`.
A **two-reference coordinate** is the centred U-statistic
`eps = U - EU`, `U = binom(n,2)^-1 sum_(j<k) h(R_j, R_k)`, with `h`
symmetric and `0 <= h <= 1`. Its **Hoeffding decomposition** is

\[
\varepsilon=\frac2n\sum_j h_1(R_j)+\binom n2^{-1}\sum_{j<k}h_2(R_j,R_k),
\]

with `h_1(r) = E h(r, R) - Eh` and `h_2(r, r') = h(r, r') - h_1(r) - h_1(r') - Eh`.
Then `E h_1(R) = 0` and `E h_2(r, R) = 0` for every `r`, by direct
computation. Also `|h_1| <= 1` and `|h_2| <= 2`.

**Lemma A (moments and tails).**
- (i) **Lone index.** Let `Z` be a bounded function of the references that
  is a product of factors, each a function of one or two references. If
  some index `j` occurs in exactly one factor, and that factor has mean
  zero in `R_j` for every value of its other argument, then `E Z = 0`.
- (ii) **Moments.** For every integer `r >= 1` there is `c_r` with
  `E eps^(2r) <= c_r n^-r` for both kinds of coordinate. Hence, for a
  vector of 18 such coordinates,
  `E |eps|^(2r) <= 18^r c_r n^-r`.
- (iii) **Tails.** If `eps` is the centred mean of `n` independent
  `[0, 1]`-valued variables, then `Pr(|eps| >= x) <= 2 exp(-2 n x^2)`.
  For a two-reference coordinate,
  `Pr(|eps| >= x) <= 2 exp(-2 floor(n/2) x^2)`.

*Proof.*
- (i) Condition on every reference except `R_j`. The other factors are
  then constants, and the one factor containing `R_j` has conditional mean
  zero. Boundedness permits the conditioning.
- (ii) One-reference: `E eps^(2r) = n^(-2r) sum E prod_(l=1)^(2r) kappa(R_(j_l))`.
  By (i) a tuple contributes only if each index occurs at least twice, so
  it has at most `r` distinct indices. There are at most `B_(2r) n^r` such
  tuples (`B` a Bell number), each with expectation at most 1. So
  `E eps^(2r) <= B_(2r) n^-r`.
  - The same count works for any bounded centred kernel, with the bound's
    power in the constant.
  - Two-reference: `(x + y)^(2r) <= 2^(2r-1) (x^(2r) + y^(2r))`. The first
    part is a one-reference coordinate (kernel `2 h_1`). In the second,
    `2r` canonical factors occupy `4r` slots. By (i) a surviving tuple has
    at most `2r` distinct indices. So it contributes at most
    `binom(n,2)^(-2r) B_(4r) n^(2r) 2^(2r) <= C n^(-2r)`.
  - The vector bound is the power-mean inequality.
- (iii) The first is Hoeffding's inequality. For the second, Hoeffding's
  representation writes `U` as the average over all permutations `pi` of
  `W_pi = (1/m') sum_(l=1)^(m') h(R_(pi(2l-1)), R_(pi(2l)))`,
  `m' = floor(n/2)`. Each `W_pi` is a mean of `m'` independent
  `[0, 1]`-valued variables. By convexity of `exp`,
  `E e^(theta (U - EU)) <= max_pi E e^(theta (W_pi - EW)) <= e^(theta^2/(8 m'))`.
  The Chernoff bound gives the tail, and the other side is symmetric. []

## 3. The background noise bound (5)

Throughout this section the references and probes are background events:
i.i.d. `mu_0`, the normalised `a^2 du dv d^2y` on `D_0`. The anchor triple
is fixed, which is the conditioning of r2 Section 1.1. The kernels below
use the background relation.

**Lemma B (bounded variation along null generators).** Fix `(v, y)` and
let `X(u) = (u, v, y)`.
- (i) For `u <= u'`, `X(u) <=_0 X(u')`.
- (ii) For any fixed point `r`, `u -> 1{r <_0 X(u)}` is non-decreasing and
  `u -> 1{X(u) <_0 r}` is non-increasing.
- (iii) For fixed `r, r'` and anchor `A`, the ancestor-pair kernel
  `1{r < X} 1{r' < X} 1{r ~ r'}` and the interval kernel
  `1{A < r} 1{r < X}` are non-decreasing in `u`.
- (iv) Their centred versions and Hoeffding projections, as functions of
  `u` with the reference arguments fixed, have:
  - one-reference centred kernels: variation at most 2, sup at most 1;
  - `h_1`: variation at most 2, sup at most 1;
  - `h_2`: variation at most 4, sup at most 2.
- (v) Anchor kernels do not depend on the probe, so they are constant in
  `u`.

*Proof.*
- (i) [WCR] Theorem 2(a) with `gamma = a^2 I`: `Delta v = 0 >= 0 = |Delta y|^2/(2 Delta m)`.
- (ii) and (iii) follow from (i) by transitivity of `J^+`, and because a
  product of non-negative non-decreasing functions is non-decreasing.
- (iv) Write `g(u; r) = E_(R') h(X(u); r, R')`. This is non-decreasing
  with values in `[0, 1]`, and so is the population mean `a(u) = E h`.
  Then:
  - `h_1 = g - a` is a difference of two such functions;
  - `h_2 = h - g(r) - g(r') + a` is a signed sum of four.

  Centred one-reference kernels are handled like `h_1`. []

**Lemma C (one oscillatory probe integral).** Let `phi` be `C^1` and
supported in a compact subset `K_phi` of the interior of `D_0`. Fix
functions `kappa_1, ..., kappa_L` on `D_0` (any fixed references) whose
restrictions to every line `u -> (u, v, y)` have variation at most `V` and
sup at most `S`. Then for `omega >= 1`,

\[
\Big|\int_{D_0}\phi(x)\prod_{l=1}^L\kappa_l(x)\,e^{-i\omega u_x}\,dx\Big|
\le\frac{|\pi_{v,y}K_\phi|}{\omega}
\Big(S^L\sup_{v,y}\|\partial_u\phi\|_{L^1(du)}+L\,V\,S^{L-1}\sup|\phi|\Big),
\]

where `dx = du dv d^2y` and `|pi_(v,y) K_phi|` is the three-dimensional
volume of the projection of `K_phi` to `(v, y)`.

*Proof.* Use Fubini, with lines at fixed `(v, y)`. On each line,
`b(u) = phi prod kappa_l` vanishes outside a compact interval, and it has
bounded variation:
`Var b <= S^L int |d_u phi| + sup|phi| * L V S^(L-1)`, by the product rule
for variations. Stieltjes integration by parts gives
`int e^(-i omega u) b du = (i omega)^-1 int e^(-i omega u) db`, with no
boundary terms. So `|int| <= Var b / omega`. Integrate over `(v, y)`. []

**Proposition D (the noise bound (5)).** Let the anchor triple be
arbitrary. Let `n >= 2` estimation references and `n_P` probes be i.i.d.
`mu_0`, independent of each other and of the anchors. Let `F_omega` and
`S_n(omega)` be as in r2 (1b), with `A_reg` from r2 (1a), and `u~` smooth
on `R^3`. Then for every `K >= 1` and `omega >= 1`,

\[
\|S_n(\omega)\|_2\le C_K\Big[\sum_{k=1}^K\omega^{k-2}n^{-k/2}
+\frac{\omega^K}{n^{(K+1)/2}}+\frac1{\sqrt{n\,n_P}}\Big].
\tag{D}
\]

With `n_P >= n` the last term is at most `1/n`. That is r2 (5) after the
absorption described there (into the `k = 2` term for `K >= 2`, and into
the remainder `omega/n` for `K = 1`).

*Proof.* Write `mu_0(x)` for the population vector of a probe at `x`,
`eps(x) = mu_hat(x) - mu_0(x)`, and
`Sbar = int F(mu_0(x) + eps(x)) mu_0(dx)`. By r2 (R1), every point used
below lies in `[0, 1]^18`, which is inside `Q`.

*Step 1: the probe average.* Given references and anchors, the probes are
i.i.d. So `E |S_n - Sbar|^2 = E Var_X F(mu_hat(X)) / n_P <= E |F(mu_hat(X))|^2 / n_P`.
- By r2 (R2), `A_reg(mu_0(x)) = 0`.
- `A_reg` is Lipschitz on `Q`, so `|F(mu_hat(x))| <= C |eps(x)|`.
- Lemma A(ii) gives `E |eps|^2 <= C/n`.

So `||S_n - Sbar||_2 <= C (n n_P)^(-1/2)`.

*Step 2: Taylor expansion without paying for the top power of `omega`.*
Fix `x` and put `G(t) = J(t) e^(-i omega u~(mu_0 + t eps))`, where
`J(t) = int_0^1 DA_reg(mu_0 + s t eps)[eps] ds`. Then
`F(mu_0 + t eps) = t G(t)`, because `A_reg(mu_0) = 0`.
- The order-`k` Taylor term of `F(mu_0 + eps)` is
  `T_k = G^((k-1))(0)/(k-1)!`.
- The remainder is `R = G(1) - sum_(j<K) G^((j))(0)/j!`.
- `J^((i))(t) = int_0^1 s^i D^(i+1) A_reg(mu_0 + s t eps)[eps^(i+1)] ds`,
  so `|J^((i))| <= C |eps|^(i+1)`.
- The `j`-th derivative of the phase factor is the factor times a
  polynomial in `omega D^l u~ [eps^l]`, bounded by `C (1 + omega)^j |eps|^j`.
- Hence `|G^((K))| <= C_K omega^K |eps|^(K+1)` for `omega >= 1`, and
  `|R| <= C_K omega^K |eps|^(K+1)`.
- By Lemma A(ii) with `r = K + 1`,
  `|| int R dmu_0 ||_2 <= (int E|R|^2 dmu_0)^(1/2) <= C_K omega^K n^(-(K+1)/2)`.

*Step 3: the coefficients have degree at most `k - 1` in `omega`.*
- `T_k(x) = sum_(|alpha| = k) c_alpha(x) eps(x)^alpha`, with
  `c_alpha = d^alpha F(mu_0(x)) / alpha!`.
- By the Leibniz rule, `d^alpha (A e^(-i omega u~))` is `e^(-i omega u~)`
  times a polynomial in `omega` of degree `k`. Its degree-`k` coefficient
  is `A (-i)^k prod d u~`, which vanishes at `mu_0(x)` by (R2).
- Every derivative of `F` vanishes at a point outside `supp A_reg`, and
  also on its boundary, by continuity from the open complement. Where
  `mu_0(x)` is in `supp A_reg`, `u~(mu_0(x)) = Sigma_u(mu_0(x)) = u_x`,
  since `Sigma` inverts the background population map.

So

\[
c_\alpha(x)=e^{-i\omega u_x}\sum_{j=0}^{k-1}\omega^j\varphi_{\alpha,j}(x),
\]

with `phi_(alpha,j)(x) = P_(alpha,j)(mu_0(x))` for smooth `P` that vanish
outside `supp A_reg`.

*Step 4: the coefficients are `C^1` with uniform bounds.*
- `phi_(alpha,j)(x) != 0` forces `mu_0(x)` into `supp A_reg`.
- Then `w o Sigma != 0`, so `x` lies in the compact set
  `K_w = {x : (u, v, R^2)(x) in supp w}`. This set is inside the interior
  of `D_0` (probes have `eta >= 0.617`, `zeta >= 0.150`; r2 Section 1.3).
- Also every `b_0(A_i, x) >= b_flat/2 > 0`, so `x` is chronologically after
  each anchor.
- On a neighbourhood of that compact set, `x -> mu_0(x)` is smooth:
  - `p_0`, `f_0` and `a_0` are smooth on the interior (Corollary K2 with
    `lambda = 0`);
  - `b_0(A, x) = alpha(u_x - u_A) kappa_0(A, x)^2 / gamma` is smooth where
    `kappa_0 > 0` and `u_x > u_A`.
- So `phi_(alpha,j)` is `C^1` with compact support in the interior.
- Its sup and `d_u` bounds are uniform over anchor triples:
  - a non-zero coefficient needs `mu_0,anc` in `supp chi_anc`, that is,
    the anchors in the closed outer set, which is compact;
  - for triples outside it, every `T_k` vanishes identically;
  - on the support, `b_0(A_i, x) = alpha(u_x - u_A) kappa_0^2 / gamma >= b_flat/2`
    with `kappa_0 <= 1`, so `u_x - u_(A_i) >= alpha^-1(gamma b_flat/2) > 0`
    uniformly;
  - `mu_0` and `d_u mu_0` are jointly continuous in (anchors, `x`) on the
    compact set {anchors in the outer set, `x` in `K_w` with that gap}, so
    they are bounded there. (Added after review, finding F12.)

*Step 5: the index count.* Expand each factor `eps_(alpha_l)(x)` into its
Hoeffding parts: one slot for an average or `h_1` part, two slots for an
`h_2` part. Fix one choice of parts, with `L` slots, `k <= L <= 2k`.
- Its normalisation is at most `4^k n^-L`.
- The term is `N sum_t Z_t`, where `t` runs over index tuples and
  `Z_t = int c_alpha(x) prod_l kappa_l(x; R_(t_l)) mu_0(dx)`.
- Here `kappa_l` is the chosen part: a centred one-reference kernel, `h_1`
  or `h_2`, or an anchor kernel (constant in `x`).

Then:
- `E |sum_t Z_t|^2 = sum_(t, t') E[Z_t conj(Z_t')]`.
- If an index occurs once in `(t, t')`, the term is zero. This is Lemma
  A(i), applied inside the two `x` integrals by Fubini; each part has mean
  zero in each reference argument.
- Surviving pairs have at most `L` distinct indices, so there are at most
  `B_(2L) n^L` of them.
- For every realisation, `|Z_t| <= C omega^(k-2)`:
  - apply Lemma C to each `phi_(alpha,j) a^2 / Z_0` (the density of
    `mu_0` in `(u, v, y)` is `a^2/Z_0`);
  - the variation and sup constants come from Lemma B;
  - sum over `j <= k - 1`, with `omega^j / omega <= omega^(k-2)`.

So `||N sum_t Z_t||_2 <= 4^k n^-L (B_(2L) n^L)^(1/2) C omega^(k-2) <= C n^(-k/2) omega^(k-2)`,
since `L >= k`. There are finitely many `alpha` and choices of parts. This
is r2 (4).

*Step 6.* Add Steps 1, 2 and 5. []

**What the index count covers.** It covers:
- repeated indices within a product, such as the diagonal of `eps_p^2`;
- products mixing average, `h_1`, `h_2` and anchor parts, including cycles
  through `h_2`;
- shared references across all probes and anchors.

There is no bias term, because every coordinate is centred (the
U-statistic `a_hat` is unbiased; r2 Section 1.1).

## 4. The canonical map: (6), (7) and the shifts

For `|u| < 1/2` put:
- `M_L = M_lam(-1/2, u)` and `M_R = M_lam(u, 1/2)`;
- `beta_L = beta(-1/2, u)` and `beta_R = beta(u, 1/2)`, with matrices
  `B_L` and `B_R`;
- `m = m(u)`, `sigma = m/M`, and `W_0 = m (M - m)/M`;
- `E_L = B_L / m` and `E_R = B_R / (M - m)`, with norms `eps_L` and
  `eps_R`;
- `P = (I + E_L)/(1 - eps_L^2)` and `S = (I + E_R)/(1 - eps_R^2)`.

Since `(m I - B_L)(m I + B_L) = (m^2 - |beta_L|^2) I`,
`M_L^-1 = P/m` and `M_R^-1 = S/(M - m)`. Put

\[
Q=M_L^{-1}+M_R^{-1},\qquad W_0Q=(1-\sigma)P+\sigma S,\qquad
\mathcal T=(W_0Q)^{-1/2}.
\]

The **canonical map** is `Psi_lam : (u, v, y) -> (u, v~, z)`, with
`z = T^-1 y` and `v~ = eta_lam + |z|^2/(2m)`. It is the map of
null_coordinate_challenge.md Section 4.

**Lemma E (the canonical map).** Let `lambda <= lambda_1`.
- (i) **Bijection.** `Psi_lam` is a smooth bijection of `{|u| < 1/2}` onto
  itself that preserves `u`, `eta` and `zeta`: `eta_0 o Psi = eta_lam` and
  `zeta_0 o Psi = zeta_lam`. Hence it maps `D_lam ∩ {|u| < 1/2}` onto
  `D_0 ∩ {|u| < 1/2}`. Its Jacobian determinant is `1/det T`.
- (ii) **The tip control.** `eps_L, eps_R <= e(s) <= 1/4`. Moreover
  `W_0 Q = alpha_0 I + E`, where

  \[
  \alpha_0=\frac{1-\sigma}{1-\epsilon_L^2}+\frac{\sigma}{1-\epsilon_R^2},\qquad
  E=\frac{(1-\sigma)E_L}{1-\epsilon_L^2}+\frac{\sigma E_R}{1-\epsilon_R^2},
  \]

  and `E` is trace-free. So
  `det T = (alpha_0^2 - |E|^2)^(-1/2)`. This gives

  \[
  |\det\mathcal T-1|\le4\,e(s)^2,\qquad \|\mathcal T-I\|\le4\,e(s).
  \tag{6}
  \]

- (iii) **Shifts.** On `D_0`, `|z|^2 <= 2 W_0 <= 2 f_max s`. Write
  `delta y = y - z` and `delta v = v - v~` at corresponding points. Then

  \[
  |\delta y|\le4e(s)\sqrt{2f_{\max}s}\le C\lambda\min(\sqrt s,\lambda/\sqrt s),
  \qquad
  |\delta v|\le\frac{12\,C_Ef_{\max}}{M}\,\lambda^2 .
  \]

  The displacement of the point in the `(T, Z, y)` chart is
  `d(X) <= |delta v| + |delta y|`. So

  \[
  d(X)\le d(s):=C_d\big(\lambda^2+\lambda\min(\sqrt s,\lambda/\sqrt s)\big)
  \le 2C_d\lambda^{3/2}.
  \]

- (iv) **Interior first order.** On `s >= s_0 > 0`, `eps_L, eps_R <= C lambda^2`,
  and

  \[
  \mathcal T-I=-\tfrac12\big[(1-\sigma)E_L+\sigma E_R\big]+O(\lambda^4),\qquad
  \delta v=\frac1{2M}z^T(E_L-E_R)z+O(\lambda^4).
  \]

*Proof.*
- (i) The past slack of the image is `v~ - |z|^2/(2m) = eta_lam` by
  definition. For the future slack, note that `1/m + 1/(M-m) = 1/W_0` and
  `|z|^2 = y^T W_0 Q y`. So
  `1 - v~ - |z|^2/(2(M-m)) = 1 - eta_lam - y^T Q y / 2 = 1 - v - y^T M_R^-1 y / 2 = zeta_lam`.
  - At fixed `u`, the map `(v, y) -> (v~, z)` is
    `z = T^-1 y`, `v~ = v + (a function of y)`, with block-triangular
    Jacobian `[[1, *], [0, T^-1]]`.
  - Its inverse is `y = T z`, `v = v~ - |z|^2/(2m) + z^T T M_L^-1 T z / 2`.
  - The diamonds correspond by Section 1.
- (ii) **The bound on `eps`.** `eps_L = |beta_L|/m <= C_E lambda min(l_L, lambda)/l_L`,
  with `l_L = u + 1/2 <= m`. This is at most
  `C_E lambda min(1, lambda/s)`. The same holds on the right.
  - `W_0/m = 1 - sigma` and `W_0/(M - m) = sigma` give the form of
    `W_0 Q`. A sum of trace-free symmetric matrices is trace-free, and
    `alpha_0 I + E` has eigenvalues `alpha_0 +/- |E|`.
  - With `eps <= e <= 1/4`:
    `0 <= alpha_0 - 1 <= (16/15) e^2` and `|E| <= (16/15) e`.
  - So `x = alpha_0^2 - |E|^2` lies in
    `[1 - 1.14 e^2, 1 + 2.3 e^2]`. The mean value theorem for
    `x^(-1/2)` gives `|det T - 1| <= (1/2)(1 - 1.14/16)^(-3/2) 2.3 e^2 < 4 e^2`.
  - The eigenvalues of `T` are `(alpha_0 +/- |E|)^(-1/2)`. Both lie
    within `4e` of 1, because `alpha_0 +/- |E|` lies in
    `[1 - 1.07 e, 1 + 1.07 e + 1.07 e^2]` and `e <= 1/4`.
- (iii) **The transverse shift.** On `D_0`, `eta_0, zeta_0 >= 0` give
  `|z|^2/(2 W_0) <= 1 - eta_0 - zeta_0 <= 1`. Also
  `W_0 <= min(m, M - m) <= f_max s`. Then
  `|delta y| = |(T - I) z| <= 4 e(s) |z|`.
  - **The longitudinal shift.** Write `X = W_0 Q = P + sigma (S - P)`.
    Then
    `T M_L^-1 T - I/m = (1/m)[X^(-1/2) P X^(-1/2) - I] = -(sigma/m) X^(-1/2)(S - P) X^(-1/2)`.
  - `||S - P|| <= ||S - I|| + ||P - I|| <= 2(eps_L + eps_R)`, since
    `(eps + eps^2)/(1 - eps^2) = eps/(1 - eps) <= 2 eps`.
  - `||X^-1|| = 1/(alpha_0 - |E|) <= 3/2`.
  - So, with `|z|^2 <= 2 W_0`,
    `|delta v| = |z^T (T M_L^-1 T - I/m) z|/2 <= (1/2)(sigma/m)(3/2) 2 (eps_L + eps_R) 2 W_0 = 3 sigma (1 - sigma)(eps_L + eps_R)`.
  - Finally,
    `sigma eps_L <= (f_max l_L / M) C_E lambda min(1, lambda/l_L) <= f_max C_E lambda^2 / M`,
    and symmetrically `(1 - sigma) eps_R <= f_max C_E lambda^2 / M`. So
    `|delta v| <= 6 f_max C_E lambda^2 / M`; the stated 12 is generous.
  - The `(T, Z, y)` displacement is `(delta v/sqrt 2, -delta v/sqrt 2, delta y)`.
- (iv) For `s >= s_0`, `eps <= C_E lambda^2 / s_0`.
  - Expand `X = I + (1 - sigma) E_L + sigma E_R + O(lambda^4)`, then
    `X^(-1/2)`.
  - Also `T M_L^-1 T = (1/m)(I + tau)(I + E_L)(I + tau) + O(lambda^4) = (1/m)(I + E_L + 2 tau) + O(lambda^4)`,
    with `tau = T - I`. This gives
    `I + sigma (E_L - E_R)` over `m`, and `sigma/m = 1/M`. []

**Lemma F (the canonical one-point laws, (7)).** Let `nu_lam` be the law of
`Psi_lam(X)` for `X` uniform on `D_lam` (normalised wave volume). Then
`nu_lam` has density `a^2 det T / Z_lam` on `D_0` in `(u, v~, z)`, with
`Z_lam = int_(D_0) a^2 det T`. Also `Z_0 = pi gamma`, `gamma = alpha(1)`.
For `lambda <= min(lambda_1, 1/2)`,

\[
\mathrm{TV}(\nu_\lambda,\mu_0)\le\frac{4f_{\max}C_E^2}{\gamma}\,
\lambda^4\Big(1+2\log\frac1{2\lambda}\Big).
\tag{7}
\]

*Proof.*
- The density is Lemma E(i) and [WCR] Lemma 3(v). The constant `c_lam`
  cancels.
- The slice of `D_0` at fixed `u` is
  `{(v~, z) : |z|^2/(2m) <= v~ <= 1 - |z|^2/(2(M - m))}`. Its volume is
  `int (1 - |z|^2/(2 W_0))_+ d^2z = pi W_0`. Then
  `Z_0 = pi int a^2 W_0 du = pi Lambda_0(-1/2, 1/2) = pi alpha(1)`, by
  Remark K1.
- For unnormalised densities `p` and `r`,
  `TV(p/|p|, r/|r|) <= |p - r|_1 / |r|_1`. So, using `a <= 1`, (6), and
  `W_0 <= f_max s` on both halves,

  \[
  \mathrm{TV}\le\frac{\pi}{Z_0}\int a^2W_0|\det\mathcal T-1|\,du
  \le\frac{8\pi f_{\max}C_E^2\lambda^2}{Z_0}\int_0^{1/2}\min(1,\lambda^2/s^2)\,s\,ds .
  \]

- The last integral is `lambda^2/2 + lambda^2 log(1/(2 lambda))`.
- The logarithm is genuine: near a tip, `det T - 1` is about
  `-(1/2) eps_L^2`, which is of order
  `lambda^4 |e^(is/lambda) - 1|^2 / s^2` (audit report 3b, Section 3.1). []

## 5. The kernel coupling: (8), (9), and closeness for every probe

Couple a wave event `X^lam` with the canonical point `X = Psi_lam(X^lam)`.
Its physical position is `xi_X = Psi_lam^-1(X)`. A **mismatch** at a
canonical pair `(X, R)` means that the wave relation of `(xi_X, xi_R)`
differs from the background relation of `(X, R)` in at least one
orientation.

**Lemma G (strip).** For `lambda <= lambda_1`, every point `P` and every
`w >= 0`,

\[
\nu_\lambda\{R:|\tilde g(P,R)|\le w\}\le C_Sw,\qquad
\mu_0\{R:|\tilde g(P,R)|\le w\}\le C_Sw .
\]

One can take `C_S = 4 (2/Z_0) pi M (M + 1)/(2 sqrt 2)`.

*Proof.*
- In the `(T, Z, y)` chart built from canonical coordinates, the density
  of `nu_lam` is `a^4 det T / Z_lam`. By (6) with `e <= 1/4`,
  `det T` lies in `[3/4, 5/4]` and `Z_lam >= (3/4) Z_0`. So the density is
  at most `2/Z_0`.
- For fixed `z_R`, `|gt| <= w` confines `T_R` to two intervals of length
  `2w`.
- The `z`-projection of `D_0` has `Z` in an interval of length
  `(M + 1)/sqrt 2` and `|y|^2 <= M/2`. []

**Lemma H (relation comparison).** Let `X, R` be in `D_0`.
- (1) **Any pair.** A mismatch implies
  `|gt(X, R)| <= c_0 lambda^2 + sqrt 2 (d(X) + d(R))`, with
  `c_0 = 2 sqrt 2 C_E`.
- (2) **Pairs separated in `u`.** If `u_R <= u_X - l_0` with `l_0 > 0`, a
  mismatch implies
  `|kappa_0(R, X)| <= C_(l_0) (lambda^2 + |delta y_X| + |delta y_R|)`,
  for `lambda <= 2 l_0`. Here `kappa_0` is evaluated at the canonical
  points. The same holds with the roles of `X` and `R` exchanged.

*Proof.*
- (1) Let `gt'` be `gt` at the physical points `(xi_X, xi_R)`, in the common
  chart. Always `|gt - gt'| <= sqrt 2 (d(X) + d(R))`.
  - **No orientation flip.** `Psi_lam` preserves `u`, and both relations
    point from smaller to larger `u` when `u_X != u_R`. So the two pairs
    can be related only in the same orientation, and only the presence of
    the relation can differ. The slice `u_X = u_R` is a null set.
  - If the wave relation of `(xi_X, xi_R)` differs from the background
    relation of the same physical points, then `|gt'| <= c_0 lambda^2`, by
    [QOB] Corollary 2.3 (both orientations, via Lemma 2.4). Hence
    `|gt| <= |gt'| + sqrt 2 (d(X) + d(R))`.
  - If the background relations of `(xi_X, xi_R)` and of `(X, R)` differ,
    then `gt'` and `gt` lie on opposite sides of 0. So
    `|gt| <= |gt - gt'| <= sqrt 2 (d(X) + d(R))`.
  - A mismatch is one of these two cases. (The orientation sentence and
    the explicit step were added after review, finding F5.)
- (2) Only `R <= X` is possible, under both relations. Put
  `Delta m = m(u_R, u_X) >= l_0`.
  - **Matrix part.** Since `|beta| <= C_E lambda^2 <= lambda/4 <= Delta m / 2`,
    the eigenvalues of `M_lam^-1 - I/Delta m` are
    `+/- |beta| / (Delta m (Delta m -/+ |beta|))`, so
    `||M_lam^-1 - I/Delta m|| <= 2 |beta| / Delta m^2 <= 2 C_E lambda^2 / l_0^2`.
  - **Points part.** `|Delta y'|^2 - |Delta z|^2 = 2 Delta z . (delta y_X - delta y_R) + |delta y_X - delta y_R|^2`,
    and all of `|z|`, `|delta y|` are bounded on `D_0`.
  - So `|kappa_lam(xi_R, xi_X) - kappa_0(R, X)|` is at most:
    - `|delta v_X| + |delta v_R|`;
    - plus `C lambda^2 / l_0^2`;
    - plus `C (|delta y_X| + |delta y_R|) / l_0`.
  - By Lemma E(iii), `|delta v| <= C lambda^2`.
  - A mismatch means the two slacks have opposite signs. []

**Lemma I (mismatch probabilities and closeness of moments).** Fix
`s_0 in (0, 1/4]` and put `K_(s_0) = {X in D_0 : s(X) >= s_0}`. For
`lambda <= min(lambda_1, s_0)` (so that Lemma H(2) applies with
`l_0 = s_0/2`; added after review, finding F6), with constants depending on
`s_0`:
- (a) `sup_(X in K_(s_0)) nu_lam{R : mismatch(X, R)} <= C lambda^2`;
- (b) `(nu_lam x nu_lam){(R, R') : mismatch(R, R')} <= C lambda^2`;
- (c) `sup_(X in D_0) nu_lam{R : mismatch(X, R)} <= C lambda^(3/2)`.

Now couple each wave reference with a background reference by a maximal
coupling of `nu_lam` and `mu_0`, so that the canonical positions agree
except with probability `TV(nu_lam, mu_0)`. Let the fixed points be
probes or anchors at canonical positions. Then:
- (d) For every raw indicator kernel of r2 Section 1.1, with its fixed
  points in `K_(s_0)`, the kernels disagree with probability at most
  `C lambda^2`. This covers past, future, the ancestor pair, the interval
  `1{A < R < X}`, and the anchor kernels. This is r2 (8). Hence
  `|mu_lam(X) - mu_0(X)| <= C lambda^2` for `X` and the anchors in
  `K_(s_0)`.
- (e) For every `X in D_0`, the probe coordinates satisfy
  `|mu_lam,probe(X) - mu_0,probe(X)| <= C lambda^(3/2)`.
- (f) **Locating a point from its wave state.** Let `K` be a compact subset
  of `D_Sigma`. There are `lambda_K > 0` and `C_K` such that, for
  `lambda < lambda_K`: if `mu_lam,probe(X)` lies in `K`, then `X` is in the
  interior of `D_0`, and its canonical state `(u, m, v~, |z|^2)(X)` lies
  within `C_K lambda^(3/2)` of `Sigma(mu_lam,probe(X))`. Consequently:
  - **anchors:** if the wave anchor numbers `iota_lam` lie in the outer
    box, each canonical anchor `u_i` is within `0.002 + C lambda^(3/2)` of
    its frozen value. So the anchors are in `K_(0.14)` for small `lambda`,
    and (d) applies to them;
  - **probes:** if `(w o Sigma)(mu_lam,probe(X)) != 0`, the canonical
    state of `X` is within `C lambda^(3/2)` of `supp w`. So `X` is in the
    neighbourhood `N` of Remark I1 for small `lambda`.

  (Added after review, finding F1. As committed, Proposition L and (R3)
  applied (d) to the anchors before showing that they are in `K_(s_0)`.)

*Proof.*
- (a) **Interior references.** Take `s(R) >= s_0/2`. Then `d(X)` and
  `d(R)` are at most `C lambda^2` (Lemma E(iii) with `s` bounded below).
  By Lemma H(1), a mismatch needs `|gt| <= C lambda^2`, which has mass at
  most `C lambda^2` by Lemma G.
  - **Tip references.** Take `s(R) < s_0/2`. Then `|u_R - u_X| >= s_0/2`,
    and Lemma H(2) applies with
    `w(s) = C (lambda^2 + lambda min(sqrt s, lambda/sqrt s)) <= C lambda^2 (1 + s^(-1/2))`.
  - For fixed `(u_R, z_R)`, `kappa_0` is affine in `v~_R` with slope
    `+/- 1`. So a mismatch confines `v~_R` to an interval of length
    `2 w(s_R)`.
  - The density is at most `2/Z_0` and the `z`-disc has area at most
    `2 pi f_max s`. So the mass is at most
    `C int_0^(s_0/2) s lambda^2 (1 + s^(-1/2)) ds <= C lambda^2`.
  - The large shift `lambda^1.5` lives where there is almost no mass.
- (b) A mismatch needs `|gt(R, R')| <= c_0 lambda^2 + 2 sqrt 2 max(d(R), d(R'))`.
  Bound by the two cases, condition on the point with the larger shift,
  and use Lemma G:
  `<= 2 C_S (c_0 lambda^2 + 2 sqrt 2 int d dnu_lam)`. Also
  `int d dnu_lam <= C int_0^(1/2) s (lambda^2 + lambda min(sqrt s, lambda/sqrt s)) ds = O(lambda^2)`.
- (c) Use `d <= 2 C_d lambda^(3/2)` everywhere in Lemma H(1), then Lemma G.
- (d) Each kernel is a product of at most three relation indicators: two
  with a fixed endpoint, and at most one between two references. Apply a
  union bound with (a) and (b), and add `TV(nu_lam, mu_0) = O(lambda^4 log(1/lambda))`
  (Lemma F). A moment differs by at most the disagreement probability of
  its kernel.
- (e) The same, with (c) in place of (a).
- (f) Fix a closed neighbourhood `K'` of `K` inside `D_Sigma` on which
  `Sigma` is Lipschitz.
  - By (e), `mu_0,probe(X)` is within `C lambda^(3/2)` of
    `mu_lam,probe(X)`, so it is in `K'` for small `lambda`.
  - Points of `K'` have `p, f > 0`, so `X` is interior. `Sigma` inverts the
    background population map, so `Sigma(mu_0,probe(X))` is the canonical
    state of `X`.
  - The Lipschitz bound gives the distance.
  - **For the anchors,** take `K` to be the set of `(p, f, a)` whose
    `Sigma` state lies in the closed anchor box. That set is compact in
    `D_Sigma`: anchor states there have `eta >= 0.126` and `zeta >= 0.638`
    (argument-3 revision review, item 1).
  - **For the probes,** take `K = supp(w o Sigma)`. []

**Lemma J (moment form, (9)).** Let `X` and the anchors be in `K_(s_0)`.
Only `p = 2` and `p = 4` are proved; r2 states (9) for every `p >= 2` but
uses only these two (scope noted after review, finding F7).
Let `eps_lam(X)` and `eps_0(X)` be the centred errors built from the
coupled reference lists. If `n lambda^2 >= 1`, then for `p = 2, 4`,

\[
\|\varepsilon_\lambda(X)-\varepsilon_0(X)\|_p\le C\lambda/\sqrt n .
\tag{9}
\]

*Proof.*
- **One-reference coordinates.** `eps_lam - eps_0` is the centred mean of
  `D_j = kappa_lam(X; R_j^lam) - kappa_0(X; R_j^0)`. These are i.i.d. with
  `|D| <= 1` and `E D^2 <= pi_X <= C lambda^2` (Lemma I(d)).
  - Put `Dbar = D - ED`, so `|Dbar| <= 2`, `E Dbar^2 <= pi_X` and
    `E Dbar^4 <= 4 pi_X`.
  - Then `E (sum Dbar_j)^4 = n E Dbar^4 + 3 n (n - 1)(E Dbar^2)^2 <= 4 n pi_X + 3 n^2 pi_X^2`.
  - So the fourth norm of the mean is at most
    `(4 pi_X)^(1/4) n^(-3/4) + 3^(1/4) pi_X^(1/2) n^(-1/2)`.
  - The first term is `C (lambda/sqrt n)(n lambda^2)^(-1/4) <= C lambda/sqrt n`.
  - The second norm is `(E Dbar^2/n)^(1/2) <= C lambda/sqrt n`.
- **Two-reference coordinates.** The difference of the two U-statistics is
  the U-statistic of the kernel
  `D(rho, rho') = h_lam(R^lam, R'^lam) - h_0(R^0, R'^0)` on the i.i.d.
  coupled pairs `rho_j = (R_j^lam, R_j^0)`.
  - By Hoeffding's representation (Lemma A(iii)) and Minkowski's
    inequality, its centred `p`-norm is at most that of a mean of
    `floor(n/2)` i.i.d. copies of `D`.
  - `E D^2 <= C lambda^2` (Lemma I(d)), so the one-reference bound applies
    with `n` replaced by `floor(n/2)`. []

**How (10) now stands.** r2 Section 5 states the exact identity
`S_n - E F(mu_lam) = T1 + T2 + T3 + T4`, and bounds each term from (5),
(7), (8), (9), (R3), and the closeness of moments for probes outside the
active set. Each of those inputs is now proved:
- (5) is Proposition D;
- (7) is Lemma F;
- (8) is Lemma I(d);
- (9) is Lemma J;
- (R3)'s two bullets are Lemma I(d), (e) and (f), and Remark I1 below;
  (f) places the anchors and every probe with non-zero weight in
  `K_(0.14)`;
- T4's law change is Lemma F, by Minkowski's integral inequality:
  `|| int G(x)(nu_lam - mu_0)(dx) ||_2 <= sup_x ||G(x)||_2 * 2 TV`.

The term bounds themselves are r2's, and were reviewed there.

**Remark I1 (probes outside the active set; r2 Section 5).** Let `N` be a
compact neighbourhood, inside `K_(0.14)`, of the active set
`{x : (u, v, R^2)(x) in supp w}`.
- **Claim.** There is `eta_N > 0` such that, for `x in D_0 \ int N` and
  every `mu` in `Q` whose probe part lies within `eta_N` of
  `mu_0,probe(x)`, `A_reg(mu) = 0`.
- **Proof of the claim.**
  - `x -> mu_0,probe(x)` is continuous on the closed diamond, and
    `D_0 \ int N` is compact. So it suffices that `mu_0,probe(x)` never
    meets the compact set `supp(w o Sigma)` there.
  - If it did, then `p, f > 0`, so `x` is interior. Then
    `Sigma(mu_0,probe(x)) = (u, m, v, R^2)(x)` would be in `supp w`, so `x`
    would be in the active set, which is inside `int N`.
- **Consequence.** With Lemma I(e), `A_reg(mu_lam(X)) = 0` for `X` outside
  `N` once `C lambda^(3/2) < eta_N/2`.
  - T2 and T3 are then non-zero only if `|eps_lam| >= eta_N/2` or
    `|eps_0| >= eta_N/2`. By Lemma A(iii) and a union bound over the 18
    coordinates, that has probability at most `C e^(-c n)`.
  - Their contribution is `O((1 + omega)^2 e^(-c n/2))` in `L^2`.

## 6. The signal: the population expansion behind (11)

**Lemma K (exact interval volumes in a Rosen plane wave).** Let `gamma` be
any Rosen profile, with `M(s, t) = int_s^t gamma^-1`, and put

\[
\Lambda(s,t)=\int_s^t\sqrt{\det\gamma(w)}\,
\sqrt{\frac{\det M(s,w)\det M(w,t)}{\det M(s,t)}}\,dw ,\qquad
\Lambda_2(s,t)=\int_s^t\Lambda(w,t)\sqrt{\det\gamma(w)}\,
\sqrt{\frac{\det M(s,w)\det M(w,t)}{\det M(s,t)}}\,dw .
\]

For `P_1 <= P_2` with `u_1 < u_2`, write `kappa = kappa(P_1, P_2)` and
`I = {P : P_1 <= P <= P_2}`. Then:
- (i) `vol(I) = pi kappa^2 Lambda(u_1, u_2)`;
- (ii) the ordered pairs `R < R'` in `I` have measure
  `(pi^2 kappa^4 / 6) Lambda_2(u_1, u_2)`;
- (iii) the probability that two independent uniform points of `I` are
  related is `Lambda_2 / (3 Lambda^2)`. It depends only on `(u_1, u_2)`.

*Proof.*
- **Composition.** For `u_1 < u < u_2` and `P = (u, v, x)`,

  \[
  \kappa(P_1,P)+\kappa(P,P_2)=\kappa-\tfrac12(x-x^*)^TQ_u(x-x^*),\qquad
  Q_u=M(u_1,u)^{-1}+M(u,u_2)^{-1},
  \]

  for some `x*`. This holds because the minimum over `x` of
  `d_1^T M_1^-1 d_1 + d_2^T M_2^-1 d_2` with `d_1 + d_2 = Delta x` is
  `Delta x^T (M_1 + M_2)^-1 Delta x`, with Hessian `2 Q_u`, and
  `M_1 + M_2 = M(u_1, u_2)`.
- **Slice of `I`.** By Theorem 2, at fixed `(u, x)` the slice of `I` is
  the `v`-interval where both slacks are non-negative. Its length is
  `h(x) = (kappa - (x - x*)^T Q_u (x - x*)/2)_+`. The slices
  `u = u_1, u_2` are null.
- **The integrals.** `int_(R^2) h^j d^2x = 2 pi kappa^(j+1) / ((j + 1) sqrt(det Q_u))`.
  Also `det Q_u = det M(u_1, u_2) / (det M(u_1, u) det M(u, u_2))`. With
  the volume element `sqrt(det gamma) du dv d^2x`, the case `j = 1` gives
  (i).
- (ii) Integrate (i) for `I(R, P_2)` over `R` in `I`.
  - At fixed `(u_R, x_R)`, as `v_R` runs over its slice,
    `s_2 = kappa(R, P_2)` runs over `[0, h(x_R)]`, with `dv_R = -ds_2`.
    So `int s_2^2 = h^3/3`.
  - The case `j = 3` gives `pi kappa^4 / (6 sqrt(det Q_u))` per slice.
  - Integrate against `Lambda(u_R, u_2) sqrt(det gamma(u_R))`.
- (iii) Comparable means `R < R'` or `R' < R`. Equality has probability
  zero. So the probability is `2 * (ii) / (i)^2`. []

**Remark K1 (the background).** For `gamma = a^2 I`:
- `a(w)^2 m(s, w) m(w, t)/m(s, t) = sin(k(w - s)) sin(k(t - w)) / (k sin(k(t - s)))`.
- So `Lambda_0(s, t) = (1 - k l cot(k l)) / (2q) = alpha(l)`, with
  `l = t - s`.
- Hence `b_0(A, X) = (kappa_0)_+^2 alpha(u_X - u_A) / alpha(1)`, which is
  r2's identity (positive part added after review, F2), and `p_0 = eta^2 alpha(u + 1/2)/gamma`, which is (2) of
  null_coordinate_challenge.md.
- `Lambda_(2,0)/(3 Lambda_0^2)` equals that document's
  `p_0(q l^2) = C_1 C_3 / (3 C_2^2)` (check G1, to 1e-10). The flat limit
  is `1/10`. This independently confirms the pair-fraction formula that
  `Sigma` inverts.

**Corollary K2 (wave moments in canonical coordinates).** Let `X` be in
`K_(s_0)` (canonical), with physical point `xi_X`, and let `lambda <= lambda_1`.
- (i) **Exact forms.**
  - `p_lam(X) = eta^2 Lambda_lam(-1/2, u) / Lambda_lam(-1/2, 1/2)`.
  - `f_lam(X) = zeta^2 Lambda_lam(u, 1/2) / Lambda_lam(-1/2, 1/2)`.
  - `a_lam(X) = p_lam(X)^2 Lambda_(2,lam)(-1/2, u) / (3 Lambda_lam(-1/2, u)^2)`.
  - For `A` with `u_X - u_A >= l_0`:
    `b_lam(A, X) = (kappa_lam(xi_A, xi_X))_+^2 Lambda_lam(u_A, u_X) / Lambda_lam(-1/2, 1/2)`.
    The positive part matters: an unrelated pair has `b = 0`, since Lemma K
    needs `P_1 <= P_2`. (Corrected after review, finding F2. Every use below
    has pair slack at least 0.1.)

  The background has the same forms with `lambda = 0`, at the canonical
  points. Here `eta` and `zeta` are the same for both, by Lemma E(i).
- (ii) **Comparison.** For `t - s >= l_0`,

  \[
  \Big|\frac{\Lambda_\lambda(s,t)}{c_\lambda}-\Lambda_0(s,t)\Big|
  +\Big|\frac{\Lambda_{2,\lambda}(s,t)}{c_\lambda^2}-\Lambda_{2,0}(s,t)\Big|
  \le C\lambda^4\big(1+\log(1/\lambda)\big).
  \]

- (iii) **Consequence.** The probe coordinates `(p, f, a)` of
  `mu_lam(X)` and `mu_0(X)` differ by `O(lambda^4 log(1/lambda))`,
  uniformly on `K_(s_0)`. On compact subsets of the interior of `D_0` that
  lie in `K_(s_0)`, in particular on `N` and at the anchors, the calibrated
  states agree to the same order:
  `Sigma(mu_lam,probe(X)) = (u, m(u), v~, |z|^2) + O(lambda^4 log(1/lambda))`.
  Also
  `eta^(0)_AX = kappa_lam(xi_A, xi_X)(1 + O(lambda^4 log(1/lambda)))`.

*Proof.*
- (i) A reference precedes `xi_X` iff it lies in `I(p_-, xi_X)`, and
  `vol D_lam = vol I(p_-, p_+) = pi Lambda_lam(-1/2, 1/2)`, since
  `kappa(p_-, p_+) = 1`. Apply Lemma K, and the identity
  `a = p^2 P_2(I(p_-, X))`.
- (ii) In the `y` chart, `sqrt(det gamma) = c_lam a^2` and
  `det M_lam(s, t) = m(s, t)^2 (1 - rho(s, t))`, with
  `rho = |beta|^2/m^2 <= C_E^2 lambda^2 min(1, lambda^2/l^2) <= 1/16`.
  - The integrand of `Lambda_lam / c_lam` is `a^2 (m_1 m_2 / m)` times
    `r = sqrt((1 - rho_1)(1 - rho_2)/(1 - rho))`. Here `r` lies in
    `[(1 - rho_1)(1 - rho_2), 1 + rho]`, so
    `|r - 1| <= rho_1 + rho_2 + rho`.
  - With `m_1 m_2 / m <= f_max min(w - s, t - w)` and
    `Phi(l) = int_0^l x min(1, lambda^2/x^2) dx`, which is at most
    `lambda^2 (1/2 + log(l/lambda)^+)`, this gives, for any interval of
    length `l`,
    `|Lambda_lam / c_lam - Lambda_0| <= 3 f_max C_E^2 lambda^2 Phi(l)`.
    This is `O(lambda^4 (1 + log(1/lambda)))`.
  - For `Lambda_2`, split the integrand difference into two parts.
    - The first is `(Lambda_lam(w, t)/c_lam - Lambda_0(w, t)) r`. Since
      `int_0^l x Phi(x) dx = O(lambda^2 (1 + log(1/lambda)))`, it
      contributes `O(lambda^4 (1 + log(1/lambda)))`.
    - The second is `Lambda_0(w, t)(r - 1)`. With
      `Lambda_0(w, t) <= f_max (t - w)^2 / 4`, the three `rho` terms
      contribute `O(lambda^2 Phi)`, `O(lambda^4)` and `O(lambda^4)`.
  - `Lambda_0(s, t) = alpha(t - s) >= alpha(l_0) > 0`, so the relative
    errors are of the same order.
- (iii) `c_lam` cancels in every ratio. On a compact subset of the
  interior, the background `(p, f, a)` lie in a compact subset of
  `D_Sigma`, and `Sigma` is Lipschitz on a neighbourhood of it. (The
  restriction to the interior was added after review, finding F3: `Sigma`
  is undefined where `eta = 0` or `zeta = 0`.)
  - The background calibration returns `(u, m, v~, |z|^2)` at the canonical
    point. Indeed `v = 1 - sigma + sigma eta - (1 - sigma) zeta = v~` and
    `R^2 = 2 W_0 (1 - eta - zeta) = |z|^2` (r2 Section 1.2).
  - `eta^(0) = (b gamma / alpha(u_hat_X - u_hat_A))^(1/2)`, with the exact
    form (i) and `u_hat - u = O(lambda^4 log(1/lambda))`. []

**Proposition L (the signal bound (11)).** Let `0 < q <= 1`, with the
frozen probe weight `w`. Condition on an anchor triple whose wave
population anchor numbers `iota_lam` lie in the three-quarter set of r2
Section 1.3. There are `lambda_2 > 0` and `c_q > 0` such that, for every
`lambda < lambda_2`, every phase `phi` and every `|omega - 1/lambda| <= 1/8`,

\[
\big|\mathbb E_{X\sim\nu_\lambda}F_\omega(\mu_\lambda(X))\big|\ge c_q\lambda^2 .
\tag{11}
\]

Explicitly, `c_q = (sqrt q / 4) 0.0101 int W f du`, with
`W(u) = (a(u)^2/Z_0) int int w(u, v~, |z|^2) dv~ d^2z`.

*Proof.* Write `z_i`, `u_i`, `m_i = m(u_i)` for the canonical anchors.
- **The anchors and the probes are in `K_(0.14)`.**
  - `iota_lam` is in the three-quarter box. By Lemma I(f) the canonical
    anchors are in `K_(0.14)`, and every probe with non-zero weight is in
    `N`. Then, by Lemma I(d) and the smoothness of `iota`, the canonical
    anchor numbers
    `iota^0 = (u_i, v~_i, |z_i|^2, z_i . z_j)` differ from `iota_lam` by
    `O(lambda^2)`. So for small `lambda` they lie in the outer box, where
    the certificate bounds of Section 8 (G9) hold.
  - In particular `u_i` is within 0.002 of `(-0.35, -0.25, -0.15)`, and
    `u - u_i >= 0.248` for every active probe.
- **Step 1: reduce to the canonical frame.** On the active set, (R3) gives
  `A_reg(mu_lam(X)) = (w H)(mu_lam(X))`. Off a compact neighbourhood of
  it, `F_omega(mu_lam(X)) = 0` (Remark I1). Then make three replacements:
  - `nu_lam` by `mu_0`: the cost is `2 TV sup|wH| = O(lambda^4 log * lambda^2)`
    (Lemma F, and `|H| = O(lambda^2)` below);
  - the calibrated phase `u_hat` by `u`:
    `omega |u_hat - u| = O(lambda^3 log)`, so the cost is
    `O(lambda^5 log)` (Corollary K2(iii));
  - `w(Sigma(mu_lam))` by `w(u, v~, |z|^2)`: the cost is
    `O(lambda^4 log * lambda^2)`.
- **Step 2: expand `H`.** Write `ell^0` for the coefficients of the
  canonical Gram matrix. Then `sum_i ell^0_i z_i = 0`, with `ell^0_3 = 1`.
  - By Corollary K2(iii), every calibrated `g` equals
    `(|z_A|^2 + |z_X|^2)/2 - Delta m (Delta v~ - kappa_lam(xi_A, xi_X)) + O(lambda^4 log)`.
  - The canonical background slack is
    `kappa_can = Delta v~ - |Delta z|^2/(2 Delta m)`, and
    `(|z_A|^2 + |z_X|^2)/2 - Delta m (Delta v~ - kappa_can) = z_A . z_X`.
  - So `g_AX = z_A . z_X + Delta_A(X) + O(lambda^4 log)`, with
    `Delta_A(X) = Delta m (kappa_lam - kappa_can)`.
  - Expand `kappa_lam` with `y = (I + tau) z`, `v = v~ + delta v`, and
    `M_lam(u_A, u)^-1 = (I + B_AX/Delta m)/Delta m + O(lambda^4)`. By
    Lemma E(iv),

  \[
  \Delta_A(X)=\Delta m(\delta v_X-\delta v_A)-\Delta z\cdot(\tau_Xz_X-\tau_Az_A)
  -\frac{\Delta z^TB_{AX}\Delta z}{2\Delta m}+O(\lambda^4),
  \]

  where `Delta z = z_X - z_A` and `B_AX = B(u_A, u)` (check G6).
  - Similarly the calibrated anchor Gram entries are
    `z_i . z_j + O(lambda^2)`, so `ell = ell^0 + delta ell` with
    `delta ell = O(lambda^2)`, a constant in `X`.
  - Hence
    `H = (sum_i delta ell_i z_i) . z_X + sum_i ell^0_i Delta_i(X) + O(lambda^4 log)`.
- **Step 3: average over the transverse angle.** At fixed
  `(u, v~, |z_X|)`, the weight, the density `a^2/Z_0`, the phase and the
  slice of `D_0` are all invariant under rotating `z_X`. So average over
  the angle:
  - terms linear in `z_X` average to zero;
  - so do the trace-free quadratic forms `z_X^T E z_X`, since
    `<z^T E z> = |z|^2 tr E / 2`. This covers `delta v_X` and
    `z_X . tau_X z_X` (trace-free up to `O(lambda^4)`, by Lemma E(iv)), and
    `z_X^T B_AX z_X`.
  - Write `B_AX = B_L(u) - B_L(u_A)`. What survives is

  \[
  \langle H\rangle(u)=\mathrm{Re}\big[\beta_L(u)D_A(u)\big]+\tilde S(u)+O(\lambda^4\log(1/\lambda)),
  \qquad D_A(u)=-\frac12\sum_i\frac{\ell^0_i\,\overline{\zeta_i}^{\,2}}{m(u)-m_i},
  \]

  with `zeta_i = z_(i,1) + i z_(i,2)`. This uses
  `z^T B z = Re(beta conj(zeta)^2)`. The term
  `S~(u) = sum_i ell^0_i [z_i^T B_L(u_i) z_i / (2(m(u) - m_i)) - (m(u) - m_i) delta v_i - z_i . tau_i z_i]`
  collects the anchor-phase terms. Each is an `O(lambda^2)` constant times
  a function of `u` that is smooth with `lambda`-independent derivatives
  on `[0.1, 0.3]`.
- **Step 4: the local term.** Two integrations by parts give

  \[
  \beta_L(u)=\frac{2k}{c_\lambda}e^{i\phi}\big(-i\lambda^2f(u)e^{iu/\lambda}\big)+\kappa_{\lambda,\phi}+O(\lambda^3),
  \]

  uniformly on `[0.1, 0.3]` (check G6). Here `kappa_(lam,phi) = O(lambda^2)`
  is a constant, the endpoint term at `u = -1/2`. `Re[kappa D_A]` joins
  `S~`.
- **Step 5: the Fourier integral.** `E = int W(u) <H>(u) e^{-i omega u} du + O(lambda^5 log)`,
  and `W >= 0` is smooth and supported in `(0.1, 0.3)`.
  - The `S~` part is at most `||(W S~)''||_1 / omega^2 = O(lambda^4)`, by
    two integrations by parts. These are r2's "anchor-phase terms vanish".
  - In `Re[Y e^(iu/lambda)] = (Y e^(iu/lambda) + conj(Y) e^(-iu/lambda))/2`,
    the conjugate part oscillates at `omega + 1/lambda`. Its smooth
    amplitude is `O(lambda^2)`, so it contributes `O(lambda^4)`.
  - The `O(lambda^3)` and `O(lambda^4 log)` remainders contribute
    `O(lambda^3)`.
  - What remains is

  \[
  E=-\frac{ik\lambda^2}{c_\lambda}e^{i\phi}\int W(u)f(u)D_A(u)e^{-i\delta u}\,du+O(\lambda^3),
  \qquad\delta=\omega-1/\lambda .
  \]

  This is the boxed limit (4) of the challenge, re-derived. With a
  separable weight, `W f` is proportional to `w_u`.
- **Step 6: no cancellation.**
  - Rotating the canonical frame by `theta` multiplies `D_A` by
    `e^(-2 i theta)`. A reflection conjugates it, which amounts to
    replacing `delta` by `-delta`. So
    `|int W f D_A e^(-i delta u)|` can be computed in the certificate's
    frame (`a_1 = (sqrt G_11, 0)`, `a_2` with positive second coordinate).
  - There `Im D_A >= 0.0106978` and `|D_A| <= 0.0424` on `[0.1, 0.3]`
    (G9). Since `W f >= 0`,
    `|int W f D_A e^(-i delta (u - 0.2))| >= int W f (Im D_A - |D_A| |delta| 0.1) >= 0.0101678 int W f`.
  - Hence `|E| >= (k/c_lam) 0.0101678 lambda^2 int W f - C lambda^3`. For
    small `lambda`, `c_lam <= 2` and the remainder is at most half the main
    term. []

The two assertions of r2 Section 6 are now proved, in sharper form:
- **anchor-phase terms vanish:** they contribute `O(lambda^4)` (Step 5);
- **calibration errors are `o(lambda^2)`:** they are
  `O(lambda^4 log(1/lambda))` (Corollary K2(iii)).

## 7. Anchor selection

r2 Section 6 selects anchors as follows.
- In a block, `3 ceil(log m)` candidates form `ceil(log m)` disjoint
  triples.
- Pilot raw anchor moments are estimated from `n' >= c m` pilot references,
  disjoint from everything else.
- A triple, with a labelling, is accepted if its pilot intrinsic numbers
  `iota_hat` lie in the inner box (half tolerances). The first accepted
  triple is used.
- The outer box forces `u_1 < u_2 < u_3` with gaps of at least 0.096, so at
  most one labelling of a triple can be accepted.

**Lemma M.** Under the wave, as `m -> infinity` and `lambda -> 0`:
- (a) `Pr(no triple is accepted) -> 0`;
- (b) `Pr(the accepted triple has iota_lam outside the three-quarter box) -> 0`.

Under the background the same holds with `iota_0`. The null argument does
not need it (r2 (R2)).

*Proof.*
- **The regular set.** Let `Omega` be the preimage of the closed outer box
  under `iota`, in raw-moment space. It is compact inside the region where
  `iota` is smooth: the anchor states lie in `D_Sigma` and every
  `b_ij >= b_low` (r2 (R1); the argument-3 revision review, item 1). Let
  `L_iota` be a Lipschitz constant of `iota` on a closed
  `eta_1`-neighbourhood of `Omega`.
- **Pilot accuracy.** Put
  `eta = min(eta_1, tol_min / (8 L_iota))`. By Lemma A(iii), a triple's
  pilot raw moments are all within `eta` of their population values with
  probability at least `1 - 30 exp(-2 floor(n'/2) eta^2)`.
  - There are 15 of them: `p, f, a` for each anchor, and `b` for all six
    ordered pairs, so that whichever labelling is accepted is covered.
  - The factor was 24 before review, finding F4.
- (b) If a triple is accepted and its pilot moments are `eta`-accurate,
  then its population moments lie within `eta` of `Omega`. So `iota_lam`
  is within `L_iota eta <= tol/8` of `iota_hat`, which is in the half box.
  A union bound over `ceil(log m)` triples gives a failure probability of
  at most `30 ceil(log m) exp(-2 floor(n'/2) eta^2) -> 0`.
- (a) **Good positions have positive probability.** Let `U'` be the set of
  canonical triples with `iota_0` in the one-eighth box. It contains a
  neighbourhood of the frozen positions, where `iota_0 = iota*`, so
  `mu_0^(x3)(U') = 2 p_1 > 0`.
  - For small `lambda`, `nu_lam^(x3)(U') >= 2 p_1 - 3 TV >= 3 p_1 / 2`.
  - On `U'`, `iota_lam = iota_0 + O(lambda^2)` (Lemma I(d)), so `iota_lam`
    is in the quarter box.
  - With `eta`-accurate pilots, `iota_hat` is then within `tol/8` of
    `iota_lam`, so in the half box, and the triple is accepted.
  - The positions of different triples are independent, but the pilot
    references are shared. So bound the two causes separately:
    `Pr(none) <= (1 - 3 p_1/2)^(ceil(log m)) + 30 ceil(log m) exp(-2 floor(n'/2) eta^2)`,
    which tends to 0. []

As r2 says, this gives no practical sample size. The probability `p_1` is
tiny (audit report 3b, F4: about `7e-38` per candidate triple).

## 8. Deterministic checks

[independent_audit/argument3_lemmas/checks.py](independent_audit/argument3_lemmas/checks.py).
No order is sampled. Seeded random draws are used only to place test
points. Pass criteria are in the script's docstring, written before the
full run.

- **G1.** `Lambda_0 = alpha` by quadrature. The pair fraction
  `Lambda_(2,0)/(3 Lambda_0^2)` against `p_0 = C_1 C_3/(3 C_2^2)` (mpmath
  `0F1`). The flat value `1/10`.
- **G2.** Lemma K for the wave:
  - the composition identity at random points;
  - the slice length `h(x)` from the two slacks of Theorem 2 against the
    closed form;
  - the slice area by a two-dimensional grid against
    `pi kappa^2 / sqrt(det Q)`.
- **G3.** Corollary K2(ii): the stated bound
  `3 f_max C_E^2 lambda^2 Phi(l)` holds, and the ratio to
  `lambda^4 log(1/lambda)` stays bounded.
- **G4.** Lemma E:
  - `zeta` preserved;
  - the closed form of `det T`;
  - the bounds (6) and the shift bounds, on a `u` grid reaching
    `s = 1e-7`, for several `lambda`, `phi` and `q`.
- **G5.** Lemma F: the numerical
  `(pi/Z_0) int a^2 W_0 |det T - 1|` is below the bound (7).
- **G6.** The interior expansion of `Delta_A(X)`: residual over
  `lambda^4` bounded. The `beta_L` expansion: residual over `lambda^3`
  bounded.
- **G7.** The angular bookkeeping of Proposition L, Step 3:
  `<sum ell^0_i Delta_i> - Re[beta_L D_A] - S~` over `lambda^4` bounded,
  with exact `Delta_i`.
- **G8.** The Fourier coefficient. The exact angle-averaged
  `sum ell^0_i Delta_i` against `W e^(-i omega u)`, compared with the
  main term of Step 5: the relative difference falls as `lambda` falls.
- **G9.** Exact rational re-certification of the anchor-guard constants
  that (11) and r2 use, on the outer box:
  - both coordinates of `a_3`;
  - `|ell_1|` and `|ell_2|`;
  - `|a_2x a_2y|`;
  - `Im D_A`, `|D_A|` and the grid coefficient;
  - the anchor-pair slack and the active-probe slack.

  This uses a rational enclosure of `sec^2(1/2)`.
- **G10.** Lemma B on a grid: past indicators along background null
  generators never decrease. Lemma A and Lemma J: the fourth-moment
  identity, in exact arithmetic.

**Results.** The run at `14c309f`, from a clean worktree:
[checks_output.txt](independent_audit/argument3_lemmas/checks_output.txt).
**31 pass, 0 fail.** System Python 3.11 with NumPy and mpmath, CPU, under a
minute. Before it, G1, G9 and G10 had been run in full while debugging, and
had passed (the script's docstring says so).

| check | result |
|---|---|
| G1 | `Lambda_0 = alpha` to `3.5e-13`; pair fraction equals `p_0` to `1.2e-15`; flat value `0.099999997` |
| G2 | composition identity to `6.7e-16`; slice length to `1.7e-16`; slice area to relative `1.4e-8` |
| G3 | the `Lambda` bound holds; deviation over `lambda^4 (1 + log(1/lambda))` at most 6.4. For `Lambda_2`: 0.42, 0.50, 0.56 at `1/lambda` = 32, 64, 128 |
| G4 | `zeta` preserved to `6e-16`; `det T` closed form to `1e-15`; the four bounds hold with largest ratios 0.19, 0.11, 0.033, 0.073 |
| G5 | (7) holds. The quantity computed is the L1 upper bound `(pi/Z_0) int a^2 W_0 |det T - 1|` on TV, not TV itself (relabelled after review, finding F9). Over `lambda^4 log(1/lambda)` it is about 10.5 at `q = 0.25` and about 41 at `q = 1`. The referee's exact TV at `lambda = 1/512` gives 8.07 and 30.9 |
| G6 | `Delta` residual over `lambda^4` at most 8.1; `beta_L` residual over `lambda^3` at most 4.1 |
| G7 | angular-bookkeeping residual over `lambda^4` at most 2.3 |
| G8 | relative difference from the main term 0.194, 0.0151, 0.0047 at `1/lambda` = 64, 128, 256, for both `delta = 0` and `delta = 1/8` |
| G9 | all ten bounds hold (values below) |
| G10 | no decrease in 4830 switch-ons; the fourth-moment identity is exact |

The certified G9 values are:
- `a_3x >= 0.0796257` and `a_3y >= 0.0793745`;
- `|ell_1|, |ell_2| <= 1.0094340`;
- `|a_2x a_2y| <= 2.00626e-5`;
- `Im D_A >= 0.0107107`, using `sec^2(1/2)` in place of 1.30;
- `|D_A| <= 0.0423907`;
- grid coefficient `>= 0.0101808`;
- anchor-pair slack `>= 0.0289167`;
- active-probe slack `>= 0.1001972`.

**Comparison with the author's certificate** (rewritten after review,
finding F10):
- **Values that agree:** `|ell| <= 1.009434`, `|D_A| <= 0.0423907`, the
  anchor-pair slack and the active-probe slack. `a_3x` agrees with the
  value audit report 3b quotes (0.079626).
- **`a_3y`.** The certificate records 0.0793741; G9 gets 0.0793745. Both
  are valid lower bounds. G9's is computed as
  `(G_23 - |G_12| G_13 / G_11) / sqrt(G_22)` with outward bounds; the
  certificate used a slightly cruder enclosure.
- **`Im D_A` and the grid coefficient** differ by input, not by error:
  - G9 uses the rational bound on `sec^2(1/2)` (about 1.2985) where the
    certificate used 1.30;
  - G9 uses the corrected `|a_2x a_2y| <= 2.00626e-5`, where the
    certificate used `2e-5`.
- **`ell`.** r2's "1.0094" is the certificate's 1.009434 rounded down by
  `3.4e-5`. With 1.009434, r2's chain gives 0.01069785 (referee r5), so
  r2's 0.0106978 stands.
- **The bound 0.0106978 that Proposition L quotes** is r2's corrected value.
  G9's sharper 0.0107107 is not needed.

**What the checks do not test.**
- G1 to G8 are numerical, at finitely many `q`, `lambda`, phases and
  points. They test the identities and the size of the remainders, not the
  proofs' uniformity statements.
- G8 tests the core term of (11): the exact angle-averaged
  `sum ell^0 Delta` at one `(v~, |z|)`. It does not include the
  calibration, the law change or the `delta ell` term, which the proof
  bounds separately.
- No check exercises Proposition D, Lemma I or Lemma M on data, since
  nothing is sampled. Their content is the proofs.

## 9. Review

[independent_audit/argument3_lemmas/review.md](independent_audit/argument3_lemmas/review.md).

**The referee.**
- A fresh-context AI referee, given the note at `6257276`, its sources and
  its checks as claims to test.
- Only short deterministic checks. Its eight scripts are archived in
  [referee/](independent_audit/argument3_lemmas/referee/), with outputs.

**Verdict: CORRECT WITH MINOR FIXES.**
- Every numbered result is true, and every result the note relies on
  holds.
- One finding is a GAP (F1).
- One statement is false as stated (F2), with no effect downstream.
- The rest are minor: constants, scope, wording and the note's account of
  itself.

**What the referee's own checks found.**
- **The signal (Proposition L).** It used the exact map and slacks,
  including `delta ell`, and took the calibration as exact.
  - The relative error against the Step 5 main term is 0.19, 0.015,
    0.0047 and 0.0009 at `1/lambda` = 64, 128, 256 and 512. The same holds
    in a rotated and in a reflected frame.
  - It reproduces the challenge's limit (4), 0.017294 + 0.008246i.
  - The actual coefficient is 0.0184, against the certified floor 0.0102.
- **Lemma K** by direct quadrature, independent of the composition
  identity: (i) to `1.5e-10`, (ii) to `4e-14`.
- **Lemma E** down to `s = 1e-10` and `q = 0.01`.
  - The one tight bound is `eps <= e(s)`, which reaches 0.993 of `e(s)`; it
    comes from [QOB] Lemma 2.1.
  - The other constants are generous.
- **Lemma F:** the exact TV lies below (7) with ratio at most 0.026.
- **Lemma H:** 24,000 seeded near-cone pairs gave 2,981 mismatches, all
  within the bound.
- **The `Lambda_2` deviation** has rate `lambda^4 log(1/lambda)` (slopes
  near -1 per unit `log(1/lambda)`).

**Corrected after the review**, in the text above:

| finding | grade | correction |
|---|---|---|
| F1 | GAP | Lemma I(f) added: a point's canonical state is located from its wave state. With it, the canonical anchors are in `K_(0.14)` before Lemma I(d) is applied to them, and every probe with non-zero weight is in `N`. Proposition L and the (R3) line of Section 5 now cite it. |
| F2 | false as stated | `b = (kappa)_+^2 Lambda / Lambda(-1/2, 1/2)`, in Corollary K2(i) and Remark K1. |
| F3 | minor | Corollary K2(iii) restricted to compact subsets of the interior, which is all that is used. |
| F4 | minor | Lemma M: 15 pilot moments, factor 30 instead of 24. |
| F5 | minor | Lemma H(1): no orientation flip, and the explicit step `|gt| <= |gt'| + sqrt 2 (d_X + d_R)`. |
| F6 | minor | Lemma I stated for `lambda <= min(lambda_1, s_0)`. |
| F7 | minor (scope) | Section 0 map; Lemma J proves `p = 2, 4` only. |
| F8 | minor (overclaim) | The preamble said "reviewed" while Section 9 was pending. It now says reviewed and corrected. |
| F9 | minor | The G5 row now says that the quantity is the L1 bound on TV, and gives the referee's exact TV. The script's label is corrected in its docstring only; outputs are unaffected. |
| F10 | minor | Section 8 now says which certificate values agree and why the others differ. |
| F11 | minor (wording) | Section 0: (8) is uniform over endpoints in `K_(s_0)`. |
| F12 | minor | Proposition D, Step 4: joint continuity, and the uniform lower bound on `u_x - u_A`. |

No inequality of r2 changes, and no statement of this note's lemmas
changes except F2's positive part and F6's range of `lambda`.

**The referee's least-trusted step** was F1, now written as Lemma I(f).
Second was F12's uniformity, now written out.

**The review did not cover:**
- r2's own T1 to T4 bounds, (R1) to (R4) and the scan, beyond reading them;
- the imported [WCR] and [QOB] results;
- the certificate that `Sigma` inverts the background map;
- the fully calibrated population signal (its r4 takes the calibration as
  exact; Corollary K2 is checked separately);
- `q` outside {0.01, 0.25, 1}, and `lambda` below 1/512.

## 10. What is now written, and what remains

- **Written here:** every step that r2 Section 8 listed as outline level.
  (9) is written for `p = 2, 4`, the cases r2 uses.
- **Still imported:**
  - the definitions of [WCR] Section 7;
  - the background calibration's inverse and its certificate (in
    null_coordinate_challenge.md);
  - r2's own reviewed parts: (R1) to (R4), the bounds on T1 to T4, and the
    scan.
- **Re-checked, not re-derived:** the anchor-guard constants are
  re-certified in exact arithmetic by G9, an independent computation, and
  by the referee's r5. The author's certificate file is not edited.
- **Unchanged:**
  - the bracket `4/3 <= alpha_low <= alpha_high <= 2`;
  - every statement is asymptotic, with constants far from any reachable
    size;
  - nothing here shows that 2 is necessary;
  - argument 2's source document is still unrevised. (Pointer added later
    on 2026-11-26, after the review by the session that started this audit session: it has since been
    revised, in [local_null_challenge_r2.md](local_null_challenge_r2.md);
    see AUDIT.md Section 15.)
  - the review is AI-level. No human has checked it.
