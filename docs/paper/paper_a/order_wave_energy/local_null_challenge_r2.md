# Using shared-reference cancellation to improve the sufficient region, revision 2: the audit's items carried in and the steps written out

2026-11-26. **Revision of [local_null_challenge.md](local_null_challenge.md).**
- **The original is unchanged.** It was frozen at `7dec69a`, and no later
  commit touched it. SHA-256 of the original as it stands:
  `3ed347fa4eb357170e70688618d1908484297de02d97c0ecaffe6a515ef3dbc7`.
  Its companion [local_null_results.md](local_null_results.md) is unchanged
  too (`996f8dcf...`). Section 6 here says what that file should say about
  audit item 4.
- **What this revision does.**
  - It carries the independent audit's five items for argument 2 into the
    text: [AUDIT.md](independent_audit/AUDIT.md) Section 3, argument 2,
    items 1 to 5, and the referee report
    [proof2_shared_reference/report.md](independent_audit/proof2_shared_reference/report.md),
    findings 1 to 8.
  - It writes the original's sketched steps as lemmas with proofs:
    - the confidence bounds and their union bound (Section 2);
    - the witness geometry, with explicit constants (Section 3);
    - the uniform bias bound, the step the audit asked a human to check
      first (Section 4);
    - the passage from a population witness to a certified rejection,
      which the referee trusted least (Section 4).
  - The test itself is not changed. Section 0 maps every change.
- **Authorship and review.** AI-written. The original derivation is by one
  AI author. The audit's AI referee graded it CORRECT with `lambda -> 0`
  added. This revision follows and cites that referee's re-derivations, and
  imports lemmas from three companion notes that were each reviewed by a
  fresh-context AI referee (Section 8). A fresh-context AI referee reviewed
  this revision and graded it CORRECT WITH MINOR FIXES. It was corrected
  after the review (Section 9). No human has checked any of it.
- **What it is not.** It is not a human-reviewed or published theorem. It
  is asymptotic. It covers one family: the exact circular vacuum plane
  wave with fixed known `0 < q <= 1`, and fixed-N i.i.d. physical-volume
  sampling. It does not change the bracket
  `4/3 <= alpha_low <= alpha_high <= 2` (AUDIT.md Section 2).

For the same fixed known `0<q<=1` circular vacuum family and fixed-N i.i.d.
physical-volume protocol, the result is

\[
\boxed{\lambda\to0\ \text{ and }\ N\lambda^3/\log N\longrightarrow\infty
\quad\Longrightarrow\quad
\operatorname{TV}(P_{N,\lambda},P_{N,0})\longrightarrow1.}
\tag{1}
\]

`lambda -> 0` is new (audit item 1). A sequence with `lambda` fixed and
`N -> infinity` satisfies the original's hypothesis and is not covered:
every expansion below is a small-`lambda` expansion. Precisely
(Theorem 4.4):
- **(level)** under the background, for **every** `N >= 4`, the test
  rejects with probability at most `N^-2`. This part is not asymptotic. It
  holds for the exact test, and for any implementation whose enclosures are
  sound ((H-enc)(i), Section 4);
- **(power)** under the wave, `lambda -> 0` and `N lambda^3 / log N -> infinity`
  make the rejection probability tend to one, uniformly in the phase.

The previous sufficient condition used exponent 4. The difference here is
to bound the **difference** between two reconstructed coordinates with their
shared sampling noise retained, and to inspect all probe pairs. No
wavelength, phase, embedding, or generator-selected pair is supplied to the
observer.

## 0. What changed against the frozen original

| where | change | audit item / referee finding |
|---|---|---|
| preamble, (1) | `lambda -> 0` added; the result split into a non-asymptotic level and an asymptotic power (Theorem 4.4) | item 1; finding 1 |
| Section 1 | rewritten: the model and its imports; Lemma 1.1 (`a = p^2 c`); Lemma 1.2 (nested kernels); the calibration and its domain `D_F` stated | finding 8 (notation) |
| Section 2 | `B_N` is now any integer with `32 N^4 <= 2^B_N <= 64 N^4`. That is the least such `B_N`, or the implementation's `bit_length(32 N^4)`, which is one more when `N` is a power of two | item 3; finding 3 |
| Section 2 | the 2022 citation replaced by Hoeffding (1963) Section 5, written out as Lemma 2.2 | finding 6 |
| Section 2 | written as Lemmas 2.1 to 2.4 and Proposition 2.5; the level holds for every `N >= 4`, skipped pairs included | |
| Section 2 | Remark 2.6: the rectangular domain `R` forces `c = a/p^2` into a window 4.1% wide at `q = 1` and 1.0% at `q = 0.25`. The test is not changed | item 4; finding 4 |
| Section 3 | rewritten: Lemma 3.1 (the leading matrix (5)), with `u` in a compact subset of the open interval and `tau` in any compact subset of `(0, 2 pi)`; Lemma 3.2 (the slack of a short pair); Lemma 3.3 (the witness set, explicit constants); Lemma 3.4 (the masses (6)); Lemma 3.5 (the count (7)) | finding 5 |
| Section 3 | the matrix of the original's (4), `F_lambda`, renamed `Phi_lam`; the reference count `M` renamed `n_R`; the boxes renamed `box_b`, `box_d` | finding 8 |
| Section 4 | the imported bias expansion replaced by Lemma 4.1, which derives the uniform bound from the companion note's exact interval volumes, importing [GFL] Corollary K2(ii) | item 2; finding 2 |
| Section 4 | Lemmas 4.2, 4.3 and Theorem 4.4: increments `O(lambda)`, the domain `R` eventually inside `D_F`, the width (8), and the assembly. One power of `log N`, through the Bernstein term only | item 5; findings 2, 4, 7 |
| Section 4 | the implementation's enclosures stated as a hypothesis, (H-enc), which the audit did not review | referee Section 5, item 4 |
| Section 5 | the original's Section 5, unchanged apart from an italic note at its head | |
| Section 6 | new: which of the frozen witnesses the test, as specified, can ever evaluate. On the coverage event, ten of 12 are skipped at every reference count, including the `q = 0.25` rows at `r <= 512`, which the audit record's item 4 omits | item 4 |
| Sections 7 to 9 | new: deterministic checks, status and imports, review | |
| (after review) preamble, Sections 0 to 2, 4, 6, 7, 9 | the corrections of Section 9: (H-enc) in three clauses, Lemma 4.2(ii)'s constant, Lemma 4.1's range and constant, notation, and reporting | review |
| (later) Sections 4, 6, 8, 9 | corrections after the review by the session that started this audit session: (H-enc) separated from the imports of (1), the bias constant and its extrapolation caveat, and the account of the audit record (Section 9) | parent review |

Deterministic checks:
[independent_audit/argument2_revision/checks.py](independent_audit/argument2_revision/checks.py).

## 1. Setting, moments and nested kernels

### 1.1 The model

Notation and imports follow
[gram_fourier_lemmas.md](gram_fourier_lemmas.md) (here [GFL]) Section 1.
- `0 < q <= 1`, `k = sqrt q`, `a(u) = cos(k u)`, `f(u) = a(u)^-2`,
  `f_max = sec^2(k/2)`.
- `m(s, t) = int_s^t f`, `m(u) = m(-1/2, u)`, `M = m(1/2)`, `t(u) = m(u)/M`.
- The wave in the `y` chart: `M_lam(s, t) = m(s, t) I - B(s, t)`, where
  `B(s, t)` is the symmetric trace-free matrix with complex representative
  `beta(s, t) = (2 lambda k / c_lam) e^(i phi) int_s^t f(w) e^(iw/lambda) dw`
  and `c_lam = sqrt(1 + 4 q lambda^2)`. [QOB] Lemma 2.1 gives
  `|beta(s, t)| <= C_E lambda min(t - s, lambda)`, with
  `C_E = 2k(4 f_max - 2)`.
- The relation ([WCR] Theorem 2): for `u_1 < u_2`, `X <=_lam Y` iff the
  pair slack `kappa_lam(X, Y) = Delta v - Delta y^T M_lam(u_1, u_2)^-1 Delta y / 2`
  is non-negative. The background is `lambda = 0`. Both relations are
  `J^+` of a spacetime, hence transitive.
- The diamond `D_lam` between the tips `p_- = (-1/2, 0, 0)` and
  `p_+ = (1/2, 1, 0)`, with slacks `eta_lam(X) = kappa_lam(p_-, X)` and
  `zeta_lam(X) = kappa_lam(X, p_+)`.
- **Sampling.** `N` i.i.d. points from the normalised volume of `D_lam`;
  the observer receives the unlabelled order. These are the model's
  definitions ([WCR] Section 7).
- **Standing smallness.** `lambda <= lambda_1 = 1/(4 C_E)`. Then
  `c_lam <= sqrt(65/64) < 1.008`, because
  `4 q lambda_1^2 = 1/(16 (4 f_max - 2)^2) <= 1/64`.

**Notation** (completed after review, finding F9). Some symbols stay
overloaded, for continuity with the original and with [GFL]; the form of
each use separates them.
- `f(u) = a(u)^-2` is a function of `u`. The future fraction always
  carries a subscript or sits in an argument list: `f_X`, `f_0`,
  `F(p, f, a)`.
- `c_X` is a pair probability, `c_lam = sqrt(1 + 4 q lambda^2)`, and
  `c_L`, `c_W`, `c'` are constants.
- The wave matrix always has arguments or a subscript: `B(s, t)`, `B_L`,
  `B_tot`. The confidence integer of Section 2 is `B_N`.
- `s(X)` and `s_0` are tip distances. The past segment length in the
  calibration is `l_p`.
- The reference count is `n_R`; `M = m(1/2)`.
- The original's matrix `F_lambda` is `Phi_lam`; `F` is the calibration.

### 1.2 Moments

Split the order by an independent uniform random permutation into
`n_R = floor(N/2)` references and `P = N - n_R` probes. Given the split, all
`N` points are still i.i.d. from the normalised volume, so the references
are i.i.d. and independent of the probes. For a probe `X` put

\[
p_X=\Pr(R<X),\quad f_X=\Pr(X<R),\quad
a_X=\Pr(R<X,\ R'<X,\ R\sim R'),
\]

with `R, R'` independent references and `~` meaning comparable. The
estimates are `p_hat = (1/n_R) sum_j 1[R_j < X]`, the same for `f_hat`, and
the complete U-statistic
`a_hat = binom(n_R, 2)^-1 sum_(j<l) 1[R_j < X, R_l < X, R_j ~ R_l]`.
Its denominator is **not** the number of ancestor pairs.

**Lemma 1.1.** Under either hypothesis, `a_X = p_X^2 c_X`, where
`c_X = P_2(I(p_-, X))` is the probability that two independent uniform
points of the interval `I(p_-, X)` are comparable.

*Proof.* A point `R` of `D_lam` precedes `X` iff `p_- <= R <= X`, since
`X <= p_+` makes `R <= p_+` automatic. So given `R < X`, `R` is uniform on
`I(p_-, X)`. Two independent such points are both in it with probability
`p_X^2`, and then are comparable with probability `c_X`. []

[GFL] Lemma K(iii) gives `c_X = Lambda_2 / (3 Lambda^2)` over `(-1/2, u_X)`.
It depends only on `u_X`, under both hypotheses.

**Lemma 1.2 (nested kernels).** If `X <=_lam Y`, then for all `R, R'`:
- `1[R < X] <= 1[R < Y]`;
- `1[Y < R] <= 1[X < R]`;
- `1[R < X, R' < X, R ~ R'] <= 1[R < Y, R' < Y, R ~ R']`.

So the three increments `p_Y - p_X`, `f_X - f_Y` and `a_Y - a_X` have
kernels with values in `{0, 1}`, and each kernel's variance is at most its
mean.

*Proof.* Transitivity of `<=_lam`. A `{0,1}` variable `xi` has
`Var xi <= E xi^2 = E xi`. []

The empirical increments are exactly the differences of the estimates, for
example `(1/n_R) sum_j (1[R_j < Y] - 1[R_j < X]) = p_hat_Y - p_hat_X`. The
same references must be used at `X` and `Y`.

### 1.3 The calibration

Let `p0(z)` be the isotropic pair probability of
[null_coordinate_challenge.md](null_coordinate_challenge.md) Section 2, and
`alpha(l) = (1 - k l cot(k l))/(2q)`, `gamma = alpha(1)`. Put

\[
\mathcal F(p,f,a)=V(p,f,a/p^2;q),\qquad
V=1-t(l_p)+t(l_p)\sqrt{p\gamma/\alpha(l_p)}-(1-t(l_p))\sqrt{f\gamma/\alpha(1-l_p)},
\tag{2}
\]

where `l_p = (p0^-1(c)/q)^(1/2)` and `c = a/p^2`. Here `l_p` is the past
segment length `u + 1/2` (the original writes `s`), `t(l_p)` stands for
`t(u)` at `u = l_p - 1/2`, and `1 - t(l_p) = t(1 - l_p)` by the symmetry of
`f`.
- **Domain.** `D_F = {p > 0, f > 0, p0(q) < a/p^2 < 1/10}`. It is open.
  The certificate of null_coordinate_challenge.md Section 2 (also
  re-derived by the audit referee, report Section 3a) gives
  `d log p0 / dz <= -1/50` on `[0, 1]`, so `p0' < 0` there, and `q <= 1`.
  So `p0^-1` is smooth on `(p0(q), 1/10)`, and `F` is `C^infinity` on
  `D_F`.
- **Background identity** (imported, re-derived by the audit referee and in
  [GFL] Remark K1 and the proof of Corollary K2(iii)). For `X` in the
  interior of `D_0`, the background moments `theta_0(X) = (p_X, f_X, a_X)`
  lie in `D_F`, and `F(theta_0(X)) = v(X)`.
- **Background monotonicity.** If `X <_0 Y`, then `v_Y >= v_X`. This is
  [WCR] Theorem 2 with `lambda = 0`.

`q <= 1` is used only through the certificate: here, and in Lemma 4.2,
which needs `p0` decreasing on `[0, q]`. (Corrected after review, finding
F8: the text first said that Sections 3 and 4 do not use it.)

## 2. Simultaneous difference confidence bounds

The rule is the original's, with `B_N` fixed as follows.
- Choose an integer `B_N` with `32 N^4 <= 2^B_N <= 64 N^4`. The least `B_N` with
  `2^B_N >= 32 N^4` qualifies, and so does the implementation's
  `bit_length(32 N^4)`, which is one more than the least when `N` is a power
  of two (check N3). Either way `B_N <= 6 + 4 log_2 N`.
- For a mean with effective sample size `kk`, the Hoeffding radius is
  `r_H = sqrt(B_N/(2 kk))`. Here `kk = n_R` for `p` and `f`, and
  `kk = floor(n_R/2)` for `a`.
- For a non-negative quantity with empirical mean `d_hat`, put
  `d_U = min(1, d_hat + r_H)` and

\[
r_D=\min\Big(r_H,\ \sqrt{2d_U B_N/kk}+2B_N/(3kk)\Big).
\tag{3}
\]

- Apply (3) to each individual moment of each probe, and to each of the
  three non-negative increments of each observed probe pair `X < Y`.

**Lemma 2.1 (the radius (3)).** Let `xi_1, ..., xi_kk` be i.i.d. with
values in `[0, 1]`, mean `d` and variance `sigma^2 <= d`. Let `H` be the
Hoeffding event `|xi_bar - d| <= r_H` and `Bn` the Bernstein event
`|xi_bar - d| <= sqrt(2 sigma^2 B_N/kk) + 2B_N/(3kk)`. Then
`Pr(not H) <= 2 e^-B_N`, `Pr(not Bn) <= 2 e^-B_N`, and on `H ∩ Bn`,
`|xi_bar - d| <= r_D`.

*Proof.*
- Hoeffding: `Pr(|xi_bar - d| > r) <= 2 exp(-2 kk r^2)`, which is `2 e^-B_N`
  at `r = r_H`.
- Bernstein, with `|xi - d| <= 1`:
  `Pr(|xi_bar - d| >= r) <= 2 exp(-kk r^2 / (2 sigma^2 + 2r/3))`. The
  exponent equals `B_N` at
  `r* = B_N/(3kk) + sqrt(B_N^2/(9kk^2) + 2 sigma^2 B_N/kk)`, and
  `r* <= sqrt(2 sigma^2 B_N/kk) + 2B_N/(3kk)`.
- On `H`, `d <= d_hat + r_H` and `d <= 1`, so `sigma^2 <= d <= d_U`. On
  `Bn` the deviation is then at most `sqrt(2 d_U B_N/kk) + 2B_N/(3kk)`, and on
  `H` it is at most `r_H`. []

The Bernstein event is defined with the **true** variance, so it is a
fixed event. Only its radius is bounded through the data, on `H`.

**Lemma 2.2 (complete order-two U-statistics).** Let `R_1, ..., R_n` be
i.i.d., `h` a symmetric kernel with values in `[0, 1]`,
`U = binom(n, 2)^-1 sum_(j<l) h(R_j, R_l)`, and `kk = floor(n/2)`. Then
Lemma 2.1 holds for `xi_bar = U`, with this `kk`, `d = E h(R_1, R_2)` and
`sigma^2 = Var h(R_1, R_2)`.

*Proof* (Hoeffding 1963, Section 5).
- For a permutation `pi` of `{1, ..., n}` put
  `V_pi = (1/kk) sum_(i=1)^kk h(R_pi(2i-1), R_pi(2i))`. Each `V_pi` is a
  mean of `kk` i.i.d. copies of `h(R_1, R_2)`.
- `U = (1/n!) sum_pi V_pi`: every `pi` contributes `kk` pairs, and by
  symmetry every unordered pair `{j, l}` occurs equally often, `n! kk / binom(n, 2)`
  times (check N5).
- For every real `s`, convexity of `exp` gives
  `E e^(sU) <= (1/n!) sum_pi E e^(s V_pi) = E e^(s V_id)`.
- The Hoeffding and Bernstein tails of Lemma 2.1 are Chernoff bounds, that
  is bounds on `E e^(s(V_id - d))`. They therefore hold for `U`. []

This replaces the original's citation of a 2022 theorem, which neither the
original's author nor the audit consulted (finding 6). The variance in the
Bernstein bound is the full kernel variance, at most `d` for a `{0, 1}`
kernel.

**Lemma 2.3 (simultaneous coverage).** Under either hypothesis, with
probability at least `1 - (3/4) N^-2`, the conclusion of Lemma 2.1 holds
simultaneously for every individual moment of every probe and every
increment of every observed probe pair.

*Proof.*
- Condition on the split and on the probes. The references are then
  i.i.d. from the normalised volume and independent of the probes. The set
  of observed probe pairs is fixed.
- Each individual moment and each increment is a mean of i.i.d. `{0, 1}`
  variables or a complete order-two U-statistic with a `{0, 1}` kernel
  (Lemma 1.2). Lemmas 2.1 and 2.2 apply to each.
- There are at most `6P <= 6N` events for the probes and
  `6 binom(P, 2) <= 3N^2` for the pairs, since a strict order relates each
  unordered pair in at most one orientation. Each fails with probability at
  most `2 e^-B_N <= 2^(1-B_N) <= 1/(16 N^4)`.
- The union bound gives failure at most
  `(6N + 3N^2)/(16 N^4) <= (3/4) N^-2` for `N >= 1`. Dependence among the
  events, reuse of references, and later selection of a pair do not
  matter. Removing the conditioning does not change the bound. []

The bound is generous: the referee computed coverage at least
`1 - 0.09 N^-3.77` (finding 3). Nothing here uses that.

For an observed pair put `b = theta_X`, `d = theta_Y - theta_X`, with
`theta = (p, f, a)`. The `f` component of `d` is negative. Let `box_b` be
the individual-moment box for `b`, `box_d` the signed increment box from
(3), and `hull(0, box_d)` the smallest box containing `0` and `box_d`. Put

\[
R=\text{box}_b+\text{hull}(0,\text{box}_d).
\]

Skip the pair if `R` is not inside `D_F`. Otherwise let `L_i` bound
`|d_i F|` and `H_ij` bound `|d_i d_j F|` on `R`. Let `e_b_i` and `e_d_i`
be the radii of `box_b` and `box_d`, and `D_j` bound `|d_j|` on `box_d`.

**Lemma 2.4 (propagation, (4)).** If `b` is in `box_b` and `d` in `box_d`,
and `F` is `C^2` on a neighbourhood of `R`, then

\[
\big|[\mathcal F(\hat b+\hat d)-\mathcal F(\hat b)]
 -[\mathcal F(b+d)-\mathcal F(b)]\big|
\le\sum_iL_ie_{d,i}+\sum_{i,j}H_{ij}e_{b,i}D_j .
\tag{4}
\]

*Proof.*
- Split the left side as
  `[F(b_hat + d_hat) - F(b_hat + d)] + [G(b_hat) - G(b)]`, with
  `G(x) = F(x + d) - F(x)`.
- The segment from `b_hat + d` to `b_hat + d_hat` lies in
  `b_hat + box_d`, which is inside `box_b + box_d`, inside `R`. So the
  first bracket is at most `sum_i L_i e_d_i`.
- `G(b_hat) - G(b) = int_0^1 grad G(b_tau) . (b_hat - b) dtau`, with
  `b_tau = b + tau (b_hat - b)` in `box_b`, and
  `d_i G(b_tau) = int_0^1 sum_j d_i d_j F(b_tau + sigma d) d_j dsigma`.
- The points `b_tau + sigma d` lie in `box_b + hull(0, box_d) = R`,
  because `sigma d` is in the convex box `hull(0, box_d)`. So the second
  bracket is at most `sum_(i,j) H_ij e_b_i D_j`. []

The `hull(0, box_d)` part is needed: the referee's negative control
(report, s2 F3) violates (4) with bounds taken on `box_b` alone.

`F(b_hat + d_hat) = F(theta_hat_Y)`, since `d_hat = theta_hat_Y - theta_hat_X`.
Evaluate the central difference `F(theta_hat_Y) - F(theta_hat_X)` with
certified intervals and add the radius in (4). Reject if its upper endpoint
is negative for any observed `X < Y`.

**Proposition 2.5 (level).** Under the background, for every `N >= 4`, the
test rejects with probability at most `(3/4) N^-2 <= N^-2`. This needs the
enclosures to be sound ((H-enc)(i), Section 4). It does not need them to be
tight.

*Proof.* `N >= 4` makes `n_R >= 2` and `floor(n_R/2) >= 1`. On the event
of Lemma 2.3, take any observed pair `X <_0 Y` that is not skipped. Then
`b` is in `box_b`, `d` is in `box_d`, and `R` is inside `D_F`. Both probes
are almost surely interior, so by Section 1.3 the interval contains
`F(theta_Y) - F(theta_X) = v_Y - v_X >= 0`, and its upper endpoint is not
negative. Skipping only removes rejections. Data-dependent skipping is
harmless because the coverage event does not depend on it. []

### Executable derivative certification

*This subsection is the original's text, unchanged. Its `s` is `l_p` of
Section 1.3.*

The new source reuses the rational intervals, alternating-series calibration,
and inverse monotonicity certificate. For c=a/p^2, write
`z=p0^-1(c)` and `s=sqrt(z/q)`. Then

\[
\mathcal F=1-t(s)+A(s)\sqrt p-A(1-s)\sqrt f,
\qquad A(s)=t(s)\sqrt{\gamma/\alpha(s)}.
\]

The source differentiates this expression twice, including the implicit
inverse derivatives `z_c=1/p0'(z)` and `z_cc=-p0''(z)/p0'(z)^3`.
All values and derivatives are enclosed on R, not sampled at its centre.
Any box touching p=0, f=0, or an inverse-calibration endpoint is skipped.
On compact interior population sets these derivatives are uniformly bounded
for fixed q, and the empirical boxes eventually remain in such a set.
The constants may be large; no finite-N 90%-power price is claimed.

**Remark 2.6 (the calibration window; audit item 4).** `R` is a box, so on
`R` the ratio `c = a/p^2` takes every value in
`[a_min / p_max^2, a_max / p_min^2]`. All of it must lie in
`(p0(q), 1/10)`.
- That window is narrow: 4.1% of `1/10` wide at `q = 1`
  (`p0(1) = 0.095919`) and 1.0% at `q = 0.25` (`p0(0.25) = 0.099031`)
  (check N1).
- With unlimited references the boxes shrink to points, and `R` tends to
  the box with corners `theta_X` and `theta_Y`. Its `c` range is then
  `[a_X/p_Y^2, a_Y/p_X^2]`, which is wider than the window allows unless
  `p_Y/p_X` is close to one. On the coverage event of Lemma 2.3, where
  `R` contains `theta_X` and `theta_Y` and hence the limit box, a pair whose
  limit range leaves the window is skipped at **every** reference count.
  Section 6 applies this to the frozen witnesses.
- Asymptotically it is harmless: for the witness pairs of Section 3, `R`
  shrinks like `lambda` around a point at fixed distance from the boundary
  of `D_F` (Lemma 4.3).
- The constants are not uniform as `q -> 0`, where the window closes.
- The referee suggested two repairs: bound the derivatives on the path set
  `{b_tau + sigma d}`, or calibrate in `(p, f, c)`. Either would change the
  frozen test, so neither is adopted here.

## 3. Short null-neighbour pairs have a larger signal-to-noise ratio

### 3.1 The canonical map and the leading matrix

Use the canonical map `Psi_lam : (u, v, y) -> (u, v~, z)` of [GFL]
Lemma E, which is the slice map of
[null_coordinate_challenge.md](null_coordinate_challenge.md) Section 4.
- It is a bijection of `D_lam ∩ {|u| < 1/2}` onto `D_0 ∩ {|u| < 1/2}` that
  preserves `u`, `eta` and `zeta`.
- Its inverse is `y = T(u) z` and `v = v~ + dv(u, z)`, with
  `dv(u, z) = z^T Gamma(u) z / 2` and `Gamma = T M_L^-1 T - I/m`. Here
  `M_L = M_lam(-1/2, u)` and `T` is as in [GFL] Lemma E. Note that `dv`
  does not depend on `v~`.
- **The matrix of the original's (4).** Put `B_L(u) = B(-1/2, u)` and
  `B_tot = B(-1/2, 1/2)`, and

\[
\Phi_\lambda(u)=-\frac{B_L(u)-t(u)B_{\rm tot}}{m(u)\,[M-m(u)]}.
\]

  This is the original's `F_lambda`, with `E_lambda = -B` ([WCR] Lemma 3).
  On `s(u) >= s_0`, [GFL] Lemma E(iv) gives `Gamma = (E_L - E_R)/M + O(lambda^4)`,
  and `(E_L - E_R)/M = -Phi_lam` because `B(u, 1/2) = B_tot - B_L(u)`. So

\[
v=\tilde v-\tfrac12z^T\Phi_\lambda(u)z+O(\lambda^4),\qquad
y=z+O(\lambda^2),
\tag{5'}
\]

  uniformly on `s >= s_0`. This is the original's challenge (5),
  null_coordinate_challenge.md.

For a complex number `w` write `[w]` for the symmetric trace-free matrix
`[[Re w, Im w], [Im w, -Re w]]`. It has eigenvalues `+/- |w|`, and the map
`w -> [w]` is linear and preserves the operator norm.

**Lemma 3.1 (the leading matrix, (5)).** Fix a compact `U` in `(-1/2, 1/2)`
and a compact `T_set` in `(0, 2 pi)`. There are `C` and `lambda_U > 0` such
that, for `lambda < lambda_U`, every phase `phi`, every `u` in `U` and every
`tau` in `T_set`,

\[
\Phi_\lambda(u+\tau\lambda)-\Phi_\lambda(u)=\lambda^2[\ell(u,\tau)]+\mathcal R,
\qquad \|\mathcal R\|\le C\lambda^3,
\tag{5}
\]

\[
\ell(u,\tau)=\frac{2ik f(u)}{c_\lambda\,m(u)[M-m(u)]}\,
e^{i(u/\lambda+\phi)}\,(e^{i\tau}-1).
\]

Moreover `|l(u, tau)| = 4k f(u) |sin(tau/2)| / (c_lam m (M - m))`, which is
at least `c_L = 16 k s_T / (1.008 M^2)` with `s_T = min over T_set of |sin(tau/2)|`.

*Proof.* Put `u' = u + tau lambda`. Choose `lambda_U` so that `u'` stays in
a compact `U'` inside `(-1/2, 1/2)`. Write
`Phi(u') - Phi(u)` as the sum of three terms.
- **The main term**, `-[B_L(u') - B_L(u)] / (m(u)(M - m(u)))`. Integrating
  by parts,

  \[
  \int_u^{u'}fe^{iw/\lambda}dw=-i\lambda f(u)e^{iu/\lambda}(e^{i\tau}-1)
  -i\lambda\big(f(u')-f(u)\big)e^{iu'/\lambda}+i\lambda\int_u^{u'}f'e^{iw/\lambda}dw .
  \]

  The last two terms are at most `2 tau lambda^2 max|f'|`. Multiplying by
  `2 lambda k e^(i phi) / c_lam` gives the representative of
  `B_L(u') - B_L(u)` as `-lambda^2 (2ik/c_lam) f(u) e^(i(u/lambda + phi))(e^(i tau) - 1) + O(lambda^3)`.
  Dividing by `-m(M - m)` gives `lambda^2 [l]`.
- **`-B_L(u') [1/(m'(M - m')) - 1/(m(M - m))]`.** `|beta(-1/2, u')| <= C_E lambda^2`
  by [QOB] Lemma 2.1. The bracket is at most `C tau lambda`, since
  `1/(m(M - m))` has bounded derivative on `U'`. So this is `O(lambda^3)`.
- **`B_tot [t(u')/(m'(M - m')) - t(u)/(m(M - m))]`.** The same argument
  gives `O(lambda^3)`.
- **The lower bound.** `|e^(i tau) - 1| = 2 |sin(tau/2)|`, `f >= 1`,
  `m (M - m) <= M^2 / 4`, and `c_lam < 1.008`. []

So `(pi/2, 3 pi/2)` is not needed: any compact subset of `(0, 2 pi)` gives
a positive `c_L` (finding 5). The remainder constant is large; the referee
measured an error of about `400 lambda` relative to the leading matrix at
`q = 1` (report, s3 B1). The construction below uses
`T_set = [pi/2, 3 pi/2]`, the original's choice, so `s_T = 1/sqrt 2`.

### 3.2 The slack of a short pair

**Lemma 3.2.** Let `U`, `T_set` and `lambda < lambda_U` be as in
Lemma 3.1, and take canonical points `X = (u_1, v~_1, z_1)` and
`Y = (u_2, v~_2, z_2)` in `D_0`, with `u_1` in `U`, `u_2 = u_1 + tau lambda`,
`tau` in `T_set`, and `w = z_2 - z_1` with `|w| <= lambda`. Let
`xi_X = Psi_lam^-1(X)` and `xi_Y` be the physical points. Then

\[
\kappa_\lambda(\xi_X,\xi_Y)=(\tilde v_2-\tilde v_1)
-\frac{\lambda^2}2z_1^T[\ell(u_1,\tau)]z_1-\mathcal A+\epsilon,
\qquad |\epsilon|\le C(\lambda^3+\lambda^2|w|),
\]

where the action `A = Delta y^T M_lam(u_1, u_2)^-1 Delta y / 2` satisfies

\[
0\le\mathcal A\le\frac{2}{3\tau\lambda}|\Delta y|^2,\qquad
|\Delta y|\le(1+4e)|w|+C_T\lambda^2|z_1|,
\]

with `e = C_E lambda^2 / s_0` and `s_0` the least tip distance on `U'`.
The constants depend on `q`, `U` and `T_set` only.

*Proof.*
- **The longitudinal part.** `Delta v = Delta v~ + dv(u_2, z_2) - dv(u_1, z_1)`.
  Write
  `dv(u_2, z_2) - dv(u_1, z_1) = z_1^T (Gamma_2 - Gamma_1) z_1 / 2 + w^T Gamma_2 (z_1 + z_2)/2`.
  The second term is at most `C lambda^2 |w|`, because
  `||Gamma|| = O(lambda^2)` and `|z| <= 1` on `D_0`. By (5'),
  `Gamma_2 - Gamma_1 = -(Phi(u_2) - Phi(u_1)) + O(lambda^4)`, and Lemma 3.1
  replaces this by `-lambda^2 [l] + O(lambda^3)`.
- **The transverse part.** `Delta y = T(u_2) w + (T(u_2) - T(u_1)) z_1`.
  [GFL] Lemma E(ii) gives `||T - I|| <= 4e`. By Lemma E(iv),
  `T - I = -[(1 - sigma) E_L + sigma E_R]/2 + O(lambda^4)`, and the bracket
  has `u`-derivative `O(lambda)`, because `|d beta_L / du| = (2 lambda k / c_lam) f`.
  So `||T(u_2) - T(u_1)|| <= C_T lambda^2`.
- **The action.** `m(u_1, u_2) >= u_2 - u_1 = tau lambda`, since `f >= 1`.
  `|beta(u_1, u_2)| <= C_E lambda min(tau lambda, lambda) <= (tau lambda)/4`,
  since `C_E lambda <= 1/4`. So `M_lam(u_1, u_2) >= (3/4) tau lambda I`,
  and `A <= |Delta y|^2 / ((3/2) tau lambda)`. []

### 3.3 The witness set

Fix `U = [-1/4, 1/4]` and `T_set = [pi/2, 3 pi/2]`, with
`rho_1 = 1/10`, `rho_2 = 1/5`. From Lemma 3.1, `c_L = 8 sqrt 2 k / (1.008 M^2)`.
Put

\[
g_1=\frac{c_L\rho_1^2}{16},\qquad g_2=\frac{c_L\rho_1^2}{8},\qquad
c'=\Big(\frac{3\pi c_L\rho_1^2}{128}\Big)^{1/2}.
\]

At `q = 1`: `c_L = 9.40`, `g_1 = 0.00588`, `g_2 = 0.0118`, `c' = 0.0832`.
At `q = 0.25`: `c_L = 5.38`, `g_1 = 0.00336`, `g_2 = 0.00672`,
`c' = 0.0629` (check N8 prints them).

The **witness set** `W_lam` is the set of canonical pairs `(X, Y)` with:
- `u_1` in `U` and `tau = (u_2 - u_1)/lambda` in `T_set`;
- `eta_0(X)` in `[3/10, 2/5]`;
- `|z_1|` in `[rho_1, rho_2]`, with the direction of `z_1` within `pi/6`
  of the line spanned by an eigenvector of `[l(u_1, tau)]` for the eigenvalue
  `-|l|`;
- `|z_2 - z_1| <= c' lambda^(3/2)`;
- `v~_2 - v~_1` in `[-g_2 lambda^2, -g_1 lambda^2]`.

Let `K = {X in D_0 : |u| <= 0.3, eta_0 >= 1/4, zeta_0 >= 2/5}`. It is a
compact subset of the interior of `D_0`.

**Lemma 3.3 (witness pairs are related).** There is `lambda_3 > 0`,
depending on `q` only, such that for `lambda < lambda_3`, every phase and
every `(X, Y)` in `W_lam`:
- (i) `X` and `Y` are in `K`;
- (ii) `xi_X <_lam xi_Y`, with `kappa_lam(xi_X, xi_Y) >= (3/64) c_L rho_1^2 lambda^2`;
- (iii) `v~_2 - v~_1 <= -g_1 lambda^2`.

*Proof.*
- (i) For `X`: `|u_1| <= 1/4` and `eta_0 >= 3/10`. Also
  `W_0 = m (M - m)/M >= min(m, M - m)/2 >= 1/8` on `U`, so
  `zeta_0 = 1 - eta_0 - |z_1|^2/(2 W_0) >= 1 - 2/5 - 4/25 = 0.44`. For
  `Y`: `u`, `v~` and `z` move by `O(lambda)`, so the slacks move by
  `O(lambda)`.
- (ii) Write `z_1 = |z_1| (cos psi, sin psi)` and let `psi_-` be the
  eigen-direction for `-|l|`. Then
  `z_1^T [l] z_1 = -|l| |z_1|^2 cos(2(psi - psi_-))`. With
  `|psi - psi_-| <= pi/6` this gives
  `-(lambda^2/2) z_1^T [l] z_1 >= (lambda^2/4) c_L rho_1^2`.
- By Lemma 3.2, with `|w| <= c' lambda^(3/2)` and `tau >= pi/2`,
  `A <= (2/(3 tau)) c'^2 lambda^2 (1 + o(1)) <= (c_L rho_1^2/32) lambda^2 (1 + o(1))`,
  and `|epsilon| = O(lambda^3)`.
- So `kappa_lam >= (-1/8 + 1/4 - 1/32) c_L rho_1^2 lambda^2 - o(lambda^2) = (3/32) c_L rho_1^2 lambda^2 - o(lambda^2)`.
  For small `lambda` this is at least half that. Since `u_1 < u_2`,
  [WCR] Theorem 2 gives `xi_X <_lam xi_Y`.
- (iii) is the definition. []

**Lemma 3.4 (the masses, (6)).** Let `X` and `Y` be independent uniform
points of `D_lam`, and `X~ = Psi_lam(X)`, `Y~ = Psi_lam(Y)`. Put
`pi_lam = Pr((X~, Y~) in W_lam)`. For `lambda < lambda_3` there are
`0 < c_W <= C_W`, depending on `q` only, with

\[
c_W\lambda^6\le\pi_\lambda\le C_W\lambda^6,
\qquad
\sup_x\Pr((x,\tilde Y)\in W_\lambda)+\sup_y\Pr((\tilde X,y)\in W_\lambda)\le C_W\lambda^6 .
\tag{6}
\]

*Proof.*
- **The density.** By [GFL] Lemma F, `X~` has density `a^2 det T / Z_lam`
  on `D_0` in `(u, v~, z)`. With `|det T - 1| <= 4 e^2 <= 1/4`, also
  `Z_lam` lies in `[(3/4) Z_0, (5/4) Z_0]`, with `Z_0 = pi gamma`. So on
  `K` the
  density lies between `d_- = (3/5) cos^2(0.3 k) / Z_0` and
  `d_+ = (5/3)/Z_0`.
- **The volume of `W_lam`.** Use `(u_1, eta_0, z_1)` for `X`, since
  `dv~ = d eta_0` at fixed `(u, z)`, and `(tau, w, v~_2 - v~_1)` for `Y`,
  with `du_2 = lambda dtau`. For each `(u_1, tau)`, the admissible
  directions of `z_1` form two arcs of total angle `2 pi/3`. By Fubini,

  \[
  \mathrm{vol}(W_\lambda)=\tfrac12\cdot\tfrac1{10}\cdot\pi\cdot\lambda\cdot
  \tfrac13\pi(\rho_2^2-\rho_1^2)\cdot\pi c'^2\lambda^3\cdot(g_2-g_1)\lambda^2 ,
  \]

  where `1/2 = |U|`, `1/10` is the length of the `eta` range and `pi = |T_set|`.
- So `d_-^2 vol <= pi_lam <= d_+^2 vol`, which is order `lambda^6`.
- **Conditional masses.** For fixed `x`, the `Y`-section has volume at
  most `pi lambda * pi c'^2 lambda^3 * (g_2 - g_1) lambda^2`. For fixed
  `y`, the `X`-section requires `u_1` in `u_y - lambda T_set`, `z_1` within
  `c' lambda^(3/2)` of `z_y`, and `v~_1` in `v~_y + [g_1, g_2] lambda^2`,
  with the same volume bound. Multiply by `d_+`. []

**Lemma 3.5 (the count, (7)).** Let `W` be the number of ordered pairs
`(i, j)` of distinct probes with `(X~_i, X~_j)` in `W_lam`. For
`lambda < lambda_3` and `P >= 2`,

\[
\mathbb EW=P(P-1)\pi_\lambda,\qquad
\operatorname{Var}W\le P^2\pi_\lambda+4P^3\pi_\lambda C_W\lambda^6,
\qquad
\Pr(W=0)\le\frac{4}{c_WP^2\lambda^6}+\frac{16C_W}{c_WP}.
\tag{7}
\]

So `N lambda^3 -> infinity` makes `Pr(W = 0) -> 0`, since `P >= N/2`.

*Proof.* Given the split, the probes are i.i.d. uniform on `D_lam`. Write
`W = sum I_ij`.
- `Var I_ij <= pi_lam`, over `P(P - 1)` terms.
- `I_ij I_ji = 0`, because `W_lam` requires `u_j > u_i`. So those
  covariances are at most zero.
- Pairs of pairs sharing exactly one index number at most `4 P^3`. For
  them `E[I_ij I_kl] <= pi_lam C_W lambda^6`, by conditioning on the shared
  point and (6).
- Disjoint pairs are independent.
- Chebyshev, with `P - 1 >= P/2`:
  `Pr(W = 0) <= Var W / (E W)^2`. []

The witness pairs are pairs of probes related in the wave order, so the
observer sees them. The observer checks every observed probe pair and does
not receive the witness selection. Using only `N` disjoint pairs would miss
the `N^2` factor.

**Remark 3.6 (why `v~` is not monotone; the referee's view).** The
`u`-lines are null generators: `(u, v, y) -> (u + tau lambda, v, y)` is a
causal move with slack exactly zero. Along it, `v~ = v + y^T Phi_lam(u) y / 2 + ...`
changes by `(lambda^2/2) z^T [l] z`, which has either sign. Lemma 3.3 makes
this quantitative.

## 4. Remove the remaining population coordinates

**Lemma 4.1 (the uniform bias bound; audit item 2).** There are `C_b` and
`lambda_b <= lambda_1`, depending on `q` only, such that for
`lambda < lambda_b`, every phase and
every canonical `X~` in `K`, with physical point `xi = Psi_lam^-1(X~)`,

\[
|\theta_\lambda(\xi)-\theta_0(\tilde X)|_\infty
+\big|\mathcal F(\theta_\lambda(\xi))-\tilde v(\tilde X)\big|
\le C_b\lambda^4\big(1+\log(1/\lambda)\big).
\]

*Proof.* This is [GFL] Corollary K2, applied on `K`, which lies in
`K_(0.2)`.
- K2(i) writes the wave moments exactly through Lemma K:
  `p = eta^2 Lambda_lam(-1/2, u)/Lambda_lam(-1/2, 1/2)`, the same for `f`,
  and `a = p^2 Lambda_2,lam / (3 Lambda_lam^2)`. The background has the same
  forms at the canonical point, with the same `eta`, `zeta`, by Lemma E(i).
- K2(ii) bounds `|Lambda_lam/c_lam - Lambda_0|` and the same for
  `Lambda_2` by `C lambda^4 (1 + log(1/lambda))` on intervals of length at
  least `l_0 = 0.2`.
- `theta_0(K)` is a compact subset of `D_F` (Lemma 4.2), `F` is Lipschitz
  on a neighbourhood of it, and `F(theta_0(X~)) = v~(X~)` by the background
  identity, since `v = v~` in the background. The moment bound puts
  `theta_lam(xi)` in that neighbourhood once `lambda < lambda_b`. (The range
  `lambda < lambda_b` replaces `lambda <= lambda_1` after review, finding
  F4. Only small `lambda` is used.) []

**The form the audit asked for.** The audit's item 2 states the bound on
the scalar kernel of the chain hierarchy:
`0 <= K_0 - Kbar_lam <= C lambda^2 min(l, lambda)^2 / l`, with
`Kbar_lam = K_0 sqrt(1 - rho)`. Corollary K2's integrand estimate is the
same estimate: `rho(s, t) = |beta(s, t)|^2 / m(s, t)^2 <= C_E^2 lambda^2 min(1, lambda^2/l^2)`,
and `|r - 1| <= rho_1 + rho_2 + rho` for the factor `r` it multiplies.
Only this upper bound is used. The sharp logarithmic coefficient of
[full_order_chain_asymptotics.md](full_order_chain_asymptotics.md) is not
needed.

**The constant is large.** The audit referee measured
`sup |bias| / (lambda^4 log(1/lambda))` up to about `2.6e4` at interior
`u`, and about `4e4` at `u = -0.2` (finding 2). On all of `K` the constant
`C_b` is about `2.5e5`, ten times the first figure.
- This revision's referee fitted `bias / lambda^4 = A + B log(1/lambda)`
  at `(u, eta, zeta) = (-0.3, 0.6, 0.4)` (review finding F3, script r4).
  It found `B` about `2.4e5` (`q = 0.25`) and `2.5e5` (`q = 1`), with `A`
  negative.
- So the ratio `|F bias| / (lambda^4 log(1/lambda))` is still rising
  towards `B`. It is `1.37e5` and `1.42e5` at `lambda = 1/2560`; those
  values are not the constant.
- The two-point bias term of Theorem 4.4 falls below the witness signal
  `g_1 lambda^2` at about `lambda = 3e-5` for `q = 0.25` and `4e-5` for
  `q = 1`. Depending on the form taken for the constant, the range is
  `2.5e-5` to `3.4e-5`, and `3.4e-5` to `4.6e-5`.
- **These are extrapolations.** The single-log law was checked only for
  `lambda >= 1/2560`, at one `u`.
- **They mark where the proof's guarantee begins, not how the test
  behaves.** On the frozen rows, the actual difference of the two points'
  biases is 0.9 to 3.1 times `lambda^4 log(1/lambda)` at `r = 2048`
  (check N10).

(Corrected after the review by the session that started this audit session. The text first quoted the
value `1.4e5` at `lambda = 1/2560` as the constant.)

**Lemma 4.2 (increments and the domain).** There are `C_theta`, `r_0 > 0`
and a compact `K_theta` inside `D_F`, depending on `q` only, such that for
`lambda < min(lambda_3, lambda_b)` and every `(X~, Y~)` in `W_lam`:
- (i) `|theta_lam(xi_Y) - theta_lam(xi_X)|_inf <= C_theta lambda`;
- (ii) every set within sup-distance `r_0` of `theta_lam(xi_X)` lies in
  `K_theta`;
- (iii) `F` is `C^infinity` on a neighbourhood of `K_theta`. Put
  `Lbar = max |grad F|` and `Hbar = max |Hess F|` there.

*Proof.*
- **`theta_0(K)` is a compact subset of `D_F`.** On `K`, `l_p = u + 1/2`
  is in `[0.2, 0.8]`, so `c_0 = p0(q l_p^2)` lies in
  `[p0(0.64 q), p0(0.04 q)]`. That is strictly inside `(p0(q), 1/10)`,
  because `p0` is strictly decreasing. Also `p_0 >= (1/4)^2 alpha(0.2)/gamma`
  and `f_0 >= (2/5)^2 alpha(0.2)/gamma`. `theta_0` is continuous, so
  `theta_0(K)` is compact, at some distance `r_* > 0` from the complement of
  `D_F`. Take `r_0 = r_*/4` and `K_theta` the closed `2 r_0`-neighbourhood
  of `theta_0(K)`.
- (i) `theta_0 = (eta^2 alpha(l_p)/gamma, zeta^2 alpha(1 - l_p)/gamma, p_0^2 p0(q l_p^2))`
  is smooth on a neighbourhood of `K`. The canonical coordinates of `X~`
  and `Y~` differ by `|Delta u| <= (3 pi/2) lambda`,
  `|Delta v~| <= g_2 lambda^2` and `|Delta z| <= c' lambda^(3/2)`, so
  `|Delta eta|, |Delta zeta| <= C lambda`. Hence
  `|theta_0(Y~) - theta_0(X~)| <= C lambda`. Add Lemma 4.1 at both points.
- (ii) By Lemma 4.1, `theta_lam(xi_X)` is within `r_0` of `theta_0(K)` for
  small `lambda`, so a set within `r_0` of it lies within `2 r_0` of
  `theta_0(K)`. (Corrected after review, finding F2. The first statement,
  for boxes of half-width `r_0` containing a point within `r_0` of
  `theta_lam(xi_X)`, reached `4 r_0 = r_*` and could touch the boundary of
  `D_F`. Lemma 4.3 uses (ii) only for `R`, which lies within `r_0` of
  `theta_lam(xi_X)`.)
- (iii) Section 1.3. []

**Hypothesis (H-enc) on the implementation.** The mathematical test
evaluates `F(theta_hat_Y) - F(theta_hat_X)` and bounds `L_i`, `H_ij` on `R`.
An implementation must satisfy three clauses. (The first two were added
after review, finding F1.)
- (i) **Soundness.** Every enclosure contains the true value, and the
  derivative bounds dominate the true suprema of `|d_i F|` and
  `|d_i d_j F|` on `R`. The level (Proposition 2.5) needs only this.
- (ii) **No spurious skip.** A pair whose `R` lies in `K_theta` is
  evaluated. This is a real condition: the implementation also returns no
  interval when an enclosure of `p0(q)` or of `s` fails, or when the
  calibration of an endpoint fails (`derivative_boxes`,
  `difference_interval`, `calibrated`).
- (iii) **Bounded and tight.** The derivative bounds used on any `R` inside
  `K_theta` are at most fixed constants `L'`, `H'`, and the enclosure of the
  central difference has width `o(lambda^2)` on such pairs.

The power needs all three. Exact evaluation with the exact suprema
satisfies all three, with `L' = Lbar` and `H' = Hbar`.

**(1) does not depend on (H-enc).** (1) is a statement about the total
variation between two laws, so one test suffices, and the exact test
satisfies all three clauses. (H-enc) is the condition for carrying
Theorem 4.4 over to the code's output. It matters only when the frozen
outputs of Section 5 are quoted as certified. (This paragraph was added
after the review by the session that started this audit session.) The original's
rational-interval implementation
(`experiments/positive_control/order_wave_local_null.py`, precision
`40 + bit_length(N)` bits, derivative enclosures over `R`) is designed to
satisfy them, but neither the audit nor this revision's referee reviewed
its rounding, enclosure or skip logic (referee report Section 5, item 4;
review Section 7). Its enclosures are very
loose at finite size: `local_null_results.md` records a Hessian bound near
`1.5e12` before the eight-band refinement.

**Lemma 4.3 (box sizes and the width (8)).** Let `lambda -> 0` with
`rho_N = N lambda^3 / log N -> infinity`. Take the exact test, or an
implementation satisfying (H-enc). On the event
of Lemma 2.3, for all large `N`, every observed pair of probes `(X, Y)`
with `(X~, Y~)` in `W_lam` satisfies:
- `R` lies in `K_theta`, so the pair is evaluated, by (H-enc)(ii);
- the radius in (4) is at most

\[
C\left(\sqrt{\frac{\lambda\log N}{N}}+\frac{\log N}{N}
+\lambda\sqrt{\frac{\log N}{N}}\right)=o(\lambda^2).
\tag{8}
\]

*Proof.*
- `kk >= floor(floor(N/2)/2) >= N/5` for `N >= 20`, and
  `B_N <= 6 + 4 log_2 N`. So every `r_H <= C epsilon_N`, with
  `epsilon_N = sqrt(log N / N) = lambda^(3/2) rho_N^(-1/2) = o(lambda^(3/2))`.
- **Individual boxes.** Their radii are at most `r_H`, so
  `e_b_i <= C epsilon_N`.
- **Increments.** Each true increment is at most `C_theta lambda`
  (Lemma 4.2(i)). On the coverage event `d_hat <= d + r_H`, so
  `d_U <= C_theta lambda + 2 r_H <= 2 C_theta lambda` for large `N`. Then
  `e_d_i <= r_D <= sqrt(4 C_theta lambda B_N/kk) + 2B_N/(3kk) <= C (sqrt(lambda log N/N) + log N/N)`,
  and `D_j <= |d_j| + 2 e_d_j <= 3 C_theta lambda`.
- **The domain.** `box_b` lies within `2 r_H` of `theta_X`, and
  `hull(0, box_d)` within `|d_i| + 2 r_D` of `0`, in each coordinate. So
  `R` lies within `C_theta lambda + 4 r_H` of `theta_lam(xi_X)`, which is
  less than `r_0` for small `lambda`. By Lemma 4.2(ii), `R` is inside `K_theta`, hence inside `D_F`.
- **The radius.** By (H-enc)(iii),
  `sum L_i e_d_i + sum H_ij e_b_i D_j <= 3 L' max_i e_d_i + 9 H' max_i e_b_i max_j D_j`.
  By the three bullets above this is at most (8).
- **Each term is `o(lambda^2)`.** Divided by `lambda^2`, the three terms
  are `rho_N^(-1/2)`, `lambda / rho_N` and `(lambda/rho_N)^(1/2)`. []

**One power of `log N` (audit item 5).** `B_N ~ 4 log_2 N` enters (8) through
the Bernstein term `sqrt(2 d_U B_N/kk)` with `d_U ~ lambda`. Requiring it to
be `o(lambda^2)` is exactly `N lambda^3 / log N -> infinity`. The witness
count needs only `N lambda^3 -> infinity`, and the bias needs nothing.

**Which term dominates (finding 7).** At every accessible `lambda` the
third term of (8), `H D e_b`, dominates, because `H/L` is about `2e4`. At
`r = 1/lambda = 2048` the referee found `L e_d = 3.3e-9` against
`H D e_b = 1.6e-8` at `N = 7.6e23`. Certifying that witness needs about
`1e24` events with ideal constants, and a reference count between `1e25`
and `1e26` with the implemented enclosures, against `lambda^-3 = 8.6e9`
(AUDIT.md Section 5).

**Theorem 4.4.** Fix known `0 < q <= 1`.
- (a) **Level.** Under the background, for every `N >= 4`, the exact test,
  and any implementation satisfying (H-enc)(i), rejects with probability at
  most `N^-2`.
- (b) **Power.** Under the wave, if `lambda -> 0` and
  `N lambda^3 / log N -> infinity`, the exact test, and any implementation
  satisfying (H-enc)(i) to (iii), rejects with probability tending to one,
  uniformly in the phase.

Hence `TV(P_(N,lam), P_(N,0)) >= Pr_lam(reject) - Pr_0(reject) -> 1`, which
is (1). The test is randomised only through the split, which is independent
of the order, so it is a test on the order's law.

*Proof.*
- (a) is Proposition 2.5.
- (b) The coverage event of Lemma 2.3 has probability at least
  `1 - N^-2` under the wave as well: Lemma 2.3 uses only that the references
  are i.i.d. On it, and on `{W >= 1}`, take a witness pair `(X, Y)`.
  - It is observed (Lemma 3.3(ii)) and evaluated (Lemma 4.3).
  - By (4), the upper endpoint is at most
    `F(theta_Y) - F(theta_X) + 2 x (8) + o(lambda^2)`, the last term from
    (H-enc)(iii).
  - By Lemma 4.1 at both points and Lemma 3.3(iii),
    `F(theta_Y) - F(theta_X) <= (v~_2 - v~_1) + 2 C_b lambda^4 (1 + log(1/lambda)) <= -g_1 lambda^2 + o(lambda^2)`.
  - So for large `N` the upper endpoint is negative, and the test rejects.
  - `Pr(reject) >= 1 - N^-2 - Pr(W = 0) -> 1`, by Lemma 3.5.
- **Uniform in phase.** No constant above depends on `phi`. In fact the
  order's law does not depend on the phase, which is a rotation of the
  transverse plane. []

**What is used where.**
- Known `q` is essential: `t`, `alpha`, `gamma` and `p0(q l_p^2)` all
  depend on it.
- `q <= 1` is used only for the certified inverse (Section 1.3 and
  Lemma 4.2).
- Both endpoints of a witness pair must be probes, because the increment
  bounds need references independent of both. The split costs constant
  factors only.
- The constants are not uniform as `q -> 0` (Remark 2.6).
- "Eventually" is far away. For the original's witness geometry the referee
  estimated `lambda` below about `1/1000` at `q = 1` and `1/4000` at
  `q = 0.25`, and `N` around `1e24` or more (report Section 4). On the
  set `W_lam` above, the window needs `lambda` near `2^-17` (check N9,
  Section 6), and the bias bound needs `lambda` below about `3e-5`
  (`q = 0.25`) or `4e-5` (`q = 1`), by extrapolation (after Lemma 4.1).

The test needs no frequency-dependent threshold or phase estimate.

## 5. Frozen deterministic challenge and limits

*This section is the original's text, unchanged. Section 6 adds what audit
item 4 implies for the frozen witnesses.*

Implementation: `experiments/positive_control/order_wave_local_null.py`.
Reuse the exact WaveFixture, its scalar hierarchy, and the old independent
five-state Rosen ODE cross-check; do not fork a second geometry generator.

- Freeze q in `{0.25,1}`, r in `{64,128,256,512,1024,2048}`, phase 0.37.
  At each rung set `u1=1/5`, `u2=1/5+pi/r`, canonical v=0.5 and transverse
  radius 0.2. Use the positive eigenvector of `F(u1)-F(u2)`. This is a
  fixed construction, with no phase/u search or reselection.
- Construct the canonical pair, record its exact causal margin, and lower
  the second v_tilde by half that margin when it is positive. Retain failures.
  Pass only the population (p,f,c) values to the inverse calibration, and
  record the recovered decrease divided by lambda^2. Record the three
  unconditional moment increments divided by lambda as a variance control.
- At q=1,r=2048 retain the exact same pair, double quadrature order and
  tighten the primitive. Then check that pair again using the existing
  independent five-state Rosen kernel and direct action quadrature. Freeze
  the independent-formulation criterion as positive recovered decrease and
  less than 1% relative change, as in the previous witness check.
- Rational controls test first and second derivative enclosures against
  independently evaluated finite differences of the calibrated value,
  using whole perturbation boxes. Check the signed difference interval on
  nearby background population moments, and observer input/count invariants.
- At q=1,r=2048 also feed rounded population moments, with an explicitly
  synthetic reference count `10^30`, into the certified difference routine.
  Record its interval without treating this as an observed order or a power
  result. This checks the algebraic decision path at nearly exact moments;
  it supplies no evidence that such a sample was collected or is affordable.

These are constructed population witnesses and numerical consistency checks.
They are not random-sprinkling classification accuracy, a finite-N powered
campaign, or an independent researcher proof audit. The implementation has
cubic reference-chain counting and potentially many rational pair bounds;
no efficient large-N runtime is established.

The new exponent 3, if this argument survives audit, still does not match
the [necessary powers below 4/3](quantized_order_bound.md). No claim is made
at either critical endpoint, for unknown q, or for arbitrary wave families.

## 6. Which frozen witnesses the test can ever evaluate (audit item 4)

[local_null_results.md](local_null_results.md) reports 12 constructed
population witnesses, 10 with a positive population decrease. Audit item 4
says that the test as specified skips some of them at every reference
count, and that the results file should say so. Check N2 computes, from
the archived population moments
([local_null_evidence/rows.jsonl](local_null_evidence/rows.jsonl)), the
limit range `[a_X/p_Y^2, a_Y/p_X^2]` of Remark 2.6 for all 12 rows.

| `r = 1/lambda` | `q = 0.25`, window `(0.099031, 0.1)` | `q = 1`, window `(0.095919, 0.1)` |
|---:|---|---|
| 64 | `[0.074560, 0.132767]` skipped | `[0.072984, 0.131380]` skipped |
| 128 | `[0.085922, 0.115252]` skipped | `[0.084392, 0.113799]` skipped |
| 256 | `[0.092414, 0.107175]` skipped | `[0.090916, 0.105712]` skipped |
| 512 | `[0.095890, 0.103299]` skipped | `[0.094410, 0.101836]` skipped |
| 1024 | `[0.097688, 0.101401]` skipped | `[0.096219, 0.099940]` evaluated |
| 2048 | `[0.098604, 0.100462]` skipped | `[0.097139, 0.099002]` evaluated |

- **On the coverage event, ten of the 12 rows are skipped at every
  reference count.** They are all
  six `q = 0.25` rows and the `q = 1` rows with `r <= 512`. Only the `q = 1`
  rows at `r = 1024` and `2048` can ever be evaluated.
- So of the ten witnesses with a positive population decrease, the test as
  specified could certify at most two, even with unlimited references. The
  "10 of 12" count of local_null_results.md is a count of population
  witnesses, not of pairs this test could certify.
- **A correction to the audit record.** AUDIT.md Section 3, argument 2,
  item 4 lists the skipped `q = 0.25` witnesses as those at
  `r = 1024, 2048`.
  - That item follows the wording of finding 4 in the audit referee's
    report (its findings table).
  - The same report's Section 6.2 has the full list, `r <= 2048`, and the
    computation agrees with it.
  - The audit referee's script s6 X4 computed only `(0.25, 1024)`,
    `(0.25, 2048)` and `(1, 2048)`.
  - AUDIT.md's earlier wording is not changed. Its Section 15, written with
    the corrections after review, records the correction. Pointers to
    Section 15 were added at the end of item 4 and at the end of
    Section 14's line on argument 2.
  - (At `cfe4b8e` this sentence named a section not yet written: review
    finding F5. After the review by the session that started this audit session, it was corrected
    again: it had said that the earlier sections were not edited, and it
    did not say where item 4's wording came from.)
- The three rows the referee computed (s6 X4) are reproduced to `4.3e-7`.
- This does not affect (1). The witness set of Section 3 enters the window
  as `lambda -> 0` (Lemma 4.3). Check N9 shows how far away that is: on
  `W_lam`, all sampled limit boxes are inside the window only at
  `lambda = 2^-17` (`q = 1`) and `2^-20` (`q = 0.25`), and at
  `lambda = 2^-11` only 41.5% (`q = 1`) and 0.1% (`q = 0.25`) are.

**What local_null_results.md should say.** Add, when it is next revised:
"The test as specified can evaluate only the `q = 1` witnesses at
`r = 1024` and `2048`, given enough references. On its coverage event, its rectangular domain puts the other
ten outside the calibration window at every reference count
(local_null_challenge_r2.md
Section 6). The asymptotic claim is unaffected."

## 7. Deterministic checks

Script: [independent_audit/argument2_revision/checks.py](independent_audit/argument2_revision/checks.py).
No order is sampled. Seeded draws only place test points. The criteria are
in the script's docstring and were committed before the run.

| check | what |
|---|---|
| N1 | `p0(1)`, `p0(0.25)` and the window widths of Remark 2.6 |
| N2 | the limit `c` range of the 12 frozen witnesses (Section 6); the nestedness signs of their increments |
| N3 | `B_N`: `32 N^4 <= 2^B_N <= 64 N^4`, `B_N <= 6 + 4 log_2 N`, and the union bound of Lemma 2.3, for the least `B_N` and the implementation's `B_N` |
| N4 | Lemma 2.1 by exact binomial enumeration |
| N5 | Lemma 2.2: the averaging identity, and `E e^(sU) <= E e^(sV)`, by exact enumeration |
| N6 | Lemma 2.4 on seeded boxes, with an explicit function whose derivative bounds are exact |
| N7 | Lemma 3.1: the residual of (5) over `lambda^3`, and the sign of (5') |
| N8 | Lemma 3.3: every sampled witness pair related, with the stated floor |
| N9 | Lemma 4.2: increments over `lambda`, `c_0(K)` inside the window, and the scale of `lambda` at which the limit boxes enter the window |
| N10 | information only: the archived rows' `|decrease - margin| / (lambda^4 log(1/lambda))` |

**Results.** One full run, after the criteria were committed at `5ffcdd9`:
**26 of 26 pass**
([checks_output.txt](independent_audit/argument2_revision/checks_output.txt)).
The run took about 5 seconds on CPU.
- **N1.** `p0(1) = 0.0959192` and `p0(0.25) = 0.0990315`: window widths
  4.081% and 0.969%.
- **N2.** All 12 rows have the nested signs. Section 6 gives the ranges.
- **N3.** 3094 values of `N`, both choices of `B_N`.
- **N4.** The worst exact failure probability is 0.121 of `4 e^-B_N`.
- **N5.** The identity is exact for `n = 2` to 7. The largest
  `E e^(sU) / E e^(sV)` is 0.961.
- **N6.** The largest ratio of the two sides of (4) is 0.739.
- **N7.** The residual of (5) over `lambda^3` is at most 148 (`q = 0.25`)
  and 281 (`q = 1`), and decreases from `1/256` to `1/4096`. On
  `[pi/2, 3 pi/2]`, `min |l| / c_L = 1.008`, the factor `1/c_lam` the bound
  gives away. The sign of (5') holds: `|dv + z^T Phi z/2| / lambda^4` is
  at most 0.32.
- **N8.** All 632 sampled pairs are related at every `lambda`, including
  `1/512`, with both points in `K`. Over all three `lambda` the least
  `kappa/lambda^2` is 0.0098 (`q = 0.25`, at `1/512`) and 0.0171 (`q = 1`),
  against floors 0.0025 and 0.0044. (Corrected after review, finding F6:
  the text first gave 0.0101 for `q = 0.25`.) N8's samples miss the binding
  corner. The referee's adversarial pairs (script r3: `u_1` near 0,
  `tau = pi/2`, `|z_1| = rho_1`, offset `pi/6`, worst direction of `w`)
  reach 0.00557 and 0.00974. That is still above the floors, and about 1.1
  times the leading-order floor `(3/32) c_L rho_1^2` (0.00504 and
  0.00881).
- **N9.** `|Delta theta_0|_inf / lambda` is between 2.8 and 3.2 at every
  `lambda`. `c_0(K)` lies inside the window with margins `3.5e-4` and
  `3.8e-5` at `q = 0.25`, and `1.5e-3` and `1.5e-4` at `q = 1`. The
  fractions of limit boxes inside the window are in Section 6.
- **N10** (information). `|decrease - margin| / (lambda^4 log(1/lambda))`
  falls from about 360 at `r = 64` to 3.1 (`q = 0.25`) and 0.9 (`q = 1`) at
  `r = 2048`. That is consistent with Lemma 4.1. Its bound's constant is
  far larger (Lemma 4.1).

**What the checks do not cover** (added after review, finding F10).
- Section 6's central claim is printed by N2 but not gated. The referee's
  r1 reproduces it independently, with `p0` computed by quadrature.
- No check touches Lemmas 3.4, 3.5, 4.2(ii), 4.3 or Theorem 4.4. They were
  checked by hand, by the author and the referee.
- N10 measures the difference of the two points' biases on the frozen rows,
  not Lemma 4.1's bound. The referee's r4 measures the bias on `K`.
- N9 uses background moments only.
- Some checks can catch only code errors: N1's series agreement, N3
  (arithmetic), N5(b) (Jensen's inequality on one kernel) and N9(ii)
  (monotonicity of `p0`). N6 tests Lemma 2.4 on a synthetic function, with
  no negative control; the audit referee's s2 F3 is one.

## 8. Status, imports and what remains

**Imported, proved elsewhere (each reviewed by a fresh-context AI referee).**
- [WCR] = [wave_causal_relation.md](wave_causal_relation.md): Theorem 2
  (the causal relation), Lemma 3 (the `y` chart).
- [QOB] = [quantized_order_bound_r2.md](quantized_order_bound_r2.md):
  Lemma 2.1 (`|beta| <= C_E lambda min(l, lambda)`).
- [GFL] = [gram_fourier_lemmas.md](gram_fourier_lemmas.md): Lemma E (the
  canonical map), Lemma F (the canonical one-point law), Lemma K and
  Corollary K2 (exact interval volumes and the bias bound).

**Imported, not proved here.**
- The model's definitions ([WCR] Section 7).
- The calibration's background identity `F(theta_0) = v` and the
  monotonicity certificate of `p0` on `[0, 1]`
  ([null_coordinate_challenge.md](null_coordinate_challenge.md) Section 2;
  re-derived by the audit referee, report Section 3a; the pair-fraction
  formula is re-derived in [GFL] Remark K1).

These two, with the reviewed imports above, are all that (1) uses.

**Not an import of (1): the hypothesis (H-enc) on the implementation.** (1)
holds for the exact test. (H-enc) is needed only to carry Theorem 4.4 over
to the code's output, for example when the frozen outputs of Section 5 are
quoted as certified. No referee has checked the code against it. (Separated
from the imports after the review by the session that started this audit session.)

**What this revision does not change.**
- The bracket `4/3 <= alpha_low <= alpha_high <= 2`. Exponent 3 is not the
  bracket's upper end; argument 3 gives 2.
- Every statement is asymptotic, with constants far from any reachable
  size (Section 4, AUDIT.md Section 5).
- The frozen test, its design and its results.
- The review is AI-level. No human has checked it.

## 9. Review

[independent_audit/argument2_revision/review.md](independent_audit/argument2_revision/review.md).

**The referee.**
- A fresh-context AI referee. It was given this revision at `cfe4b8e`, its
  imports, the frozen original, the audit's items for argument 2 and the
  audit referee's report, as claims to test. It was told not to read the
  other reviews.
- Only short deterministic checks. Its five scripts are archived in
  [referee/](independent_audit/argument2_revision/referee/), with outputs.

**Verdict: CORRECT WITH MINOR FIXES.**
- Every numbered result is correct as stated, or after a minor fix.
- Each import is used within its stated hypotheses. In particular, `K`
  lies in `K_(0.2)`, every interval used has length at least 0.2, and `p`,
  `f`, `a` are exactly the forms of [GFL] Corollary K2(i).
- No step needs a new idea. All ten findings are MINOR.
- All five audit items are addressed. Section 6 is correct, and the audit
  record's item 4 does omit the `q = 0.25` rows at `r <= 512`.

**What the referee's own checks found.**
- **r1:** `p0` computed from [GFL] Remark K1's interval kernel by
  quadrature, not from the series: `p0(1) = 0.09591918` and
  `p0(0.25) = 0.09903149`. The Section 6 table is reproduced to `4.8e-7`.
- **r2:** the exact laws of complete U-statistics in 384 cases. The `r_D`
  failure probability is at most 0.080 of `4 e^-B_N`. The negative control
  with `kk = C(n, 2)` violates its bound by a factor of 29.7.
- **r3:** an independent canonical map, and adversarial witness pairs down
  to `lambda = 1/32768`. All are related and in `K`. The least
  `kappa/lambda^2` is 0.00557 (`q = 0.25`) and 0.00974 (`q = 1`), and
  `|kappa - leading| / lambda^3` stays at most 5.5.
- **r4:** the bias on `K`, wave against background, from exact interval
  volumes. It behaves as `lambda^4 (A + B log(1/lambda))` with stable
  slopes; the constants are in finding F3. The criterion of its first
  version failed for `F` at both `q`, by a factor of 2.04, before the
  asymptotic regime. That output is kept, and the script's docstring
  explains the change.
- **r5:** exact counts behind N9's fractions: 2000 of 2000 at `2^-17` for
  `q = 1`, and 1994 of 2000 for `q = 0.25`.

**Corrected after the review**, in the text above:

| finding | grade | correction |
|---|---|---|
| F1 | minor | (H-enc) now has three clauses: (i) soundness, (ii) no spurious skip, (iii) bounded and tight. The level (preamble, Proposition 2.5, Theorem 4.4(a)) is stated for the exact test or a sound implementation. |
| F2 | minor (constant slip) | Lemma 4.2(ii) is restated for sets within `r_0` of `theta_lam(xi_X)`, with `K_theta` the closed `2 r_0`-neighbourhood. Lemma 4.3 and Theorem 4.4 are unaffected. |
| F3 | minor (practicality) | After Lemma 4.1: the bias constant on all of `K`, about `2.5e5` (the limit of the log slope; `1.4e5` at `lambda = 1/2560` and still rising), and its consequence, `lambda` below about `3e-5` (`q = 0.25`) or `4e-5` (`q = 1`), by extrapolation. "What is used where" says so too. (This row first said `1.4e5`; corrected after the review by the session that started this audit session.) |
| F4 | minor (scope) | Lemma 4.1 holds for `lambda < lambda_b`, and Lemma 4.2 for `lambda < min(lambda_3, lambda_b)`. |
| F5 | minor (false status sentence) | Section 6's sentence about AUDIT.md Section 15 is now true: that section was written with these corrections. A note records that it was not true at `cfe4b8e`. |
| F6 | minor | Section 7, N8: the least slack is 0.0098 at `q = 0.25`, and the referee's adversarial minima are added. N8 itself is not changed after its run. |
| F7 | minor (wording) | "On the coverage event" is added in Remark 2.6, Section 6 and Section 0. |
| F8 | minor (wording) | Section 1.3: `p0' < 0` comes from the certificate, and `q <= 1` is also used in Lemma 4.2. |
| F9 | minor | A notation paragraph in Section 1.1. The confidence integer is now `B_N` and the past segment length `l_p`. The `Z_lam` range and the Section 0 wording on Lemma 4.1 are fixed. |
| F10 | minor | Section 7: what the checks do not cover. |

No inequality of the proof changes, except Lemma 4.2(ii)'s constant (F2) and
Lemma 4.1's range of `lambda` (F4). The test, (1), and the level and power
statements are unchanged, apart from the conditions on the implementation
(F1).

**The referee's least-trusted step** is Lemma 4.1 feeding Theorem 4.4.
- The theorem subtracts the bias bound at the two points separately.
- With the constant measured on `K`, about `2.5e5`, that term exceeds the
  witness signal `g_1 lambda^2` until `lambda` is about `3e-5`
  (`q = 0.25`) or `4e-5` (`q = 1`).
- The referee states two caveats: it took [GFL] Corollary K2(ii) as given,
  and it checked the single-log law only down to `lambda = 1/2560`, at one
  `u`. So these values are extrapolations. They mark where the proof's
  guarantee begins, not how the test behaves. (These caveats were added
  after the review by the session that started this audit session.)
- Check N10 suggests that the difference of the two biases is far smaller
  than either bias. A bound on that difference would bring "eventually"
  closer, but no such bound is proved here.
- The referee trusts the asymptotic logic.

**Corrected after the review by the session that started this audit session.** Later the same day,
the session that had started the audit session reviewed this revision at
`faf6262`, and the user approved its corrections. It judged all three of
its points to leave (1) and the bracket unchanged. The corrections are
wording and records only; no inequality and no part of the test changes.
- **The audit record (Section 6).** The earlier wording of AUDIT.md is not
  changed, but pointers were added in two places. Item 4 followed the
  wording of finding 4 of the audit referee's report; the same report's
  Section 6.2 had the full list. The sentence for local_null_results.md
  now says "can evaluate only ..., given enough references".
- **(H-enc) (Sections 4 and 8).** (1) holds for the exact test, so (H-enc)
  is not an import of (1). Theorem 4.4(b) is stated for the exact test and
  for any implementation satisfying (H-enc). (H-enc) is the condition for
  carrying the theorem over to the code's output.
- **The bias constant (after Lemma 4.1, F3 row, least-trusted step).** It
  is about `2.5e5`, the limit of the log slope. `1.4e5` is its value at
  `lambda = 1/2560`, and it is still rising. The crossover is about `3e-5`
  (`q = 0.25`) and `4e-5` (`q = 1`), by extrapolation. It marks where the
  proof's guarantee begins, not how the test behaves.

**The review did not cover:**
- the proofs of the imports ([GFL] E, F, K, K2; [WCR] Theorem 2 and
  Lemma 3; [QOB] Lemma 2.1), only their statements;
- the `p0` certificate itself, beyond monotonicity on a grid;
- the implementation's soundness, rounding, precision and skips;
- Lemmas 3.4 and 3.5 numerically, which were checked by hand only;
- phase uniformity numerically (phase 0.37 only, plus the analytic
  argument);
- geometry beyond `lambda = 1/32768`, and the bias beyond `1/2560`;
- local_null_results.md and the provenance of Section 5.
