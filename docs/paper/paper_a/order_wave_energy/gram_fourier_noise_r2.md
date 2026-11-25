# Shared-reference Fourier noise, revision 2: the regularised statistic written out

2026-11-26. **Revision of [gram_fourier_noise.md](gram_fourier_noise.md).**
- **The original is unchanged.** It is frozen: its SHA-256 is recorded in
  the anchor-guard certificate at `f1072f0`.
- **What this revision does.** It writes the regularised statistic down as
  one object (Section 1). It also makes the corrections the independent
  audit asked for ([AUDIT.md](independent_audit/AUDIT.md) Section 3,
  argument 3, items 1 to 7) and defines `B` in Section 4. Everything else
  is the original text, which Section 0 maps.
- **Authorship and review.** AI-written. The original derivation is by one
  AI author and was audited by AI referees. This revision is reviewed in
  [independent_audit/argument3_revision/](independent_audit/argument3_revision/).
- **What it is not.** It is not a human-reviewed or published theorem.
  It is asymptotic only, and it covers one family: the exact circular
  plane wave, known `0 < q <= 1`, fixed-N i.i.d. physical-volume sampling.
- **Status.** The original's preamble said the programme record remained
  "insufficient below 4/3, sufficient at 3". It also said not to promote
  the upper endpoint on the strength of this draft. Since then the AI-only
  audit has stated the bracket `4/3 <= alpha_low <= alpha_high <= 2`, with
  its conditions (AUDIT.md Section 2). This revision writes argument 3's
  statistic and its wave-side chain in full. Section 8 lists what remains
  at outline level.

The proposed result is: for each fixed integer K>=1 there is a randomized,
order-only, wavelength-free test, with known fixed 0<q<=1, separating the exact
circular vacuum family from its effective background whenever

\[
N\asymp\lambda^{-\alpha},\qquad \alpha>2+2/K.                 \tag{1}
\]

Consequently the family would cover each fixed power above 2. It does not
establish a critical-window law at 2, show that 2 is necessary, or give a single
test uniform in a growing K. The impossibility exponent is still below 4/3.
The scan cap is the only binding requirement for the `2/K` (Section 7).

## 0. What changed against the frozen original

| where | change | audit item |
|---|---|---|
| Section 1 | rewritten: the 18 raw moments listed; the regularised statistic `A_reg` defined with explicit cutoffs; properties (R1) to (R4) proved; the anchor Gram block renamed `G_A` | items 1, 2 |
| Section 2, last sentence | the reason for (3) corrected | item 6 |
| Section 3 | (5) states `n_P >= c n` and the conditioning on the anchor triple, and lists `1/n` once; the elementary check holds for general `omega` with the factor `1 - sinc^2(omega/2)` | item 6 |
| Section 4 | `B` defined (trace-free, not "positive") | notation; requested with the audit corrections |
| Section 5 | rewritten as the identity `S_n - E F(mu_lam) = T1 + T2 + T3 + T4`, with the out-of-active-set remark | item 4 |
| Section 6 | the false identity `|a_2x a_2y| = |G12|` corrected; certificate value 0.0106978; "anchor pair slacks" defined as `b_ij > 0` | items 3, 5 |
| Section 7 | null vote bound `(K C_K / log n)^2`; no anchor-acceptance event needed under the null; scan-cap remark | items 2, 6, 7 |
| Section 8 | status updated | |
| preamble | replaced; the original's status sentences summarised under "Status"; a sentence on the scan cap after (1) | |
| Section 2 | a note that `B(u)` in (3) is generic | |
| Section 3 | the simplified-bound sentence reworded around the absorbed `1/n` | item 6 |
| Section 6 | "outer good set" replaced by "three-quarter set, where `chi_anc = 1`"; the left side of (11) and "together with (R3)" added | items 1, 2 |
| Section 7 | "and regularized calibration being realized" dropped, since (1a) is now written; the guard paragraph rewritten, with the "retain a failure" rule kept | |
| (after review) Sections 1, 3, 4, 5, 6, 7, 8 | the corrections of Section 9 | review |

Deterministic checks of the new numbers:
[independent_audit/argument3_revision/checks.py](independent_audit/argument3_revision/checks.py).

## 1. Raw moments, the regularised statistic, and the pointwise zero

### 1.1 Raw moments

Condition on an anchor triple `A1, A2, A3`, selected from disjoint pilot
events as in Section 6. References `R_j` and probes `X_i` are independent of
that selection. Use the bounded order kernels

\[
p_X=E1_{R<X},\quad f_X=E1_{X<R},\quad
a_X=E[1_{R<X}1_{R'<X}1_{R\sim R'}],\quad
b_{AX}=E1_{A<R<X}.
\]

Here ~ means comparable, not a link or a null relation. Estimate a by the
unbiased order-two U-statistic over **all reference pairs**, before forming
`c_hat=a_hat/p_hat^2`. In particular, a conditional ancestor-pair fraction
has a different finite-n denominator and must not be silently substituted
in this expansion. Anchor p,f,a and anchor-to-anchor b use the same references.

For a probe X the raw-moment vector `mu(X)` has **18** coordinates:

- probe: `p_X, f_X, a_X` (3);
- anchors: `p_{A_i}, f_{A_i}, a_{A_i}` for `i = 1, 2, 3` (9);
- anchor to anchor: `b_12, b_13, b_23`, with `b_ij = b_{A_i A_j}` (3);
- anchor to probe: `b_1X, b_2X, b_3X`, with `b_iX = b_{A_i X}` (3).

The 12 anchor coordinates, the second and third groups, are written
`mu_anc`. `mu_0(X)` and `mu_lam(X)` are the population vectors under the
background and the wave, and `epsilon(X) = mu_hat(X) - mu(X)` is the
centred empirical error. Every coordinate of `epsilon` is an average or an
order-two U-statistic.

### 1.2 Calibration and the Gram functional

Let `Sigma : (p, f, a) -> (u, m, v, R^2)` be the background state
calibration of `gram_fourier_challenge.md` Section 1:
- `c = a / p^2`, `u` from `c` by inverting the background chain
  calibration, `m = m(u)`, `t = m / M`;
- `eta` and `zeta` from `p` and `f`;
- `v = 1 - t + t eta - (1 - t) zeta`;
- `R^2 = 2 m (M - m)(1 - eta - zeta) / M`.

Its regular domain `D_Sigma` is the open set of `(p, f, a)` with `p, f > 0`,
`c` strictly inside the range of the background chain calibration,
`eta, zeta > 0` and `eta + zeta < 1`. On `D_Sigma` the map `Sigma` is
smooth with nonsingular derivative, for fixed `q > 0`.

For an ordered pair with `u_j > u_i` and `b_ij >= 0`, put

\[
\eta^{(0)}_{ij}=\sqrt{b_{ij}\gamma/\alpha(u_j-u_i)},\qquad
g_{ij}=\frac{R_i^2+R_j^2}{2}-(m_j-m_i)(v_j-v_i-\eta^{(0)}_{ij}),
\]

with `alpha` the background interval coefficient and `gamma = alpha(1)`.

**Notation.**
- `alpha(.)` with an argument is the interval coefficient
  (`alpha(L) = (1 - x cot x)/(2q)`, `x = sqrt(q) L`). Without an argument,
  `alpha` is the sampling exponent in (1).
- `c = a/p^2`. The wave normalisation `sqrt(1 + 4 q lambda^2)` is written
  `c_lam` (Section 4).
- `K` is the Taylor order. `K_lambda` and `K_0` in (8) are order kernels.

- Let `G_A` be the 2 by 2 anchor Gram block, with diagonal `R_1^2, R_2^2`
  and off-diagonal `g_12`. The challenge called this block `B`. It is
  renamed here so that it cannot be confused with the trace-free matrix
  `B` of Section 4.
- Where `G_A` is positive definite, put
  `ell = (-(g_13, g_23) G_A^-1, 1)` and

\[
H(\mu)=\ell_1 g_{1X}+\ell_2 g_{2X}+g_{3X}.
\]

### 1.3 The regularised statistic

**Anchor cutoff.** Write `iota(mu_anc)` for the 12 intrinsic anchor
numbers `(u_i, v_i, R_i^2, g_ij)`, `i < j`, and `iota*` for their values at
the frozen triple of the challenge.
- Use the outer tolerances `tau_k` of Section 6: 0.002 for the `u_i` and
  `v_i`, and 0.00002 for the six Gram entries `R_i^2` and `g_ij`.
- Use the standard smooth step
  `s(x) = e^{-1/x} / (e^{-1/x} + e^{-1/(1-x)})` on `(0, 1)`, with `s = 0`
  for `x <= 0` and `s = 1` for `x >= 1`. Put `psi(x) = s(4 (1 - |x|))`.
  This is even and smooth, with `psi = 1` on `[-3/4, 3/4]` and `psi = 0`
  outside `(-1, 1)`.
- Put

\[
\chi_{\rm anc}(\mu_{\rm anc})=\prod_{k=1}^{12}\psi\big((\iota_k-\iota^*_k)/\tau_k\big)
\]

where `iota` is defined, and `chi_anc = 0` elsewhere.

**Interval cutoff.** Let `b_*(q) > 0` be the infimum of the background
population `b` over:
- the three anchor pairs, with anchors in the outer set (support of
  `chi_anc`);
- the three anchor-probe pairs, with probes in the active set (support of
  `w`).

An explicit computable lower bound is
`b_low(q) = s_A^2 alpha(0.096) / gamma`, with `s_A = 0.028916...`, because:
- `b = (slack)^2 alpha(Delta u) / gamma`;
- Section 6 certifies anchor-pair slacks of at least `s_A` and active-probe
  slacks of at least 0.1002;
- on these sets every u-gap is at least 0.096, and `alpha` increases with
  the gap.

At `q` = 1, 0.25 and 0.01, `b_low` = 7.18e-6, 7.58e-6 and 7.70e-6. The
corresponding bound for anchor-probe pairs is about 80 times larger.

Put `b_flat = b_low(q) / 2 <= b_*(q) / 2` and
`chi_b(b) = s(2 b / b_flat - 1)`. This is smooth, with `chi_b = 0` on
`(-infinity, b_flat / 2]` and `chi_b = 1` on `[b_flat, infinity)`.

**Probe weight.** `w(u, v, R^2)` is the challenge's separable bump. Its
support lies in the interior of `Sigma(D_Sigma)`: probes are strictly
inside the diamond, with `eta >= 0.617` and `zeta >= 0.150`. So
`w o Sigma` has compact support inside `D_Sigma` and extends by zero to a
smooth function of `(p_X, f_X, a_X)`.

**The statistic.**

\[
A_{\rm reg}(\mu)=\chi_{\rm anc}(\mu_{\rm anc})
\prod_{1\le i<j\le3}\chi_b(b_{ij})\prod_{i=1}^{3}\chi_b(b_{iX})\,
(w\circ\Sigma)(p_X,f_X,a_X)\,H(\mu),                          \tag{1a}
\]

with `A_reg = 0` wherever one of the cutoff factors vanishes. The cutoff
acts on all 18 raw moments:
- on the 12 anchor coordinates through `chi_anc`, and also through
  `chi_b(b_ij)`;
- on the three anchor-probe intervals through `chi_b(b_iX)`;
- on the probe's `(p, f, a)` through `w o Sigma`.

**Phase.** Let `u~` be any smooth bounded function of `(p_X, f_X, a_X)`
equal to `Sigma_u` on an open neighbourhood of the probe projection of
`supp A_reg`.

**The score.** For `omega >= 1`,

\[
F_\omega(\mu)=A_{\rm reg}(\mu)\,e^{-i\omega\tilde u(\mu)},\qquad
S_n(\omega)=\frac1{n_P}\sum_{X}F_\omega(\hat\mu(X)).              \tag{1b}
\]

This single `F_omega` is the object of (5), of (10) and of (11) below.

### 1.4 Properties

**(R1) Smooth, and bounded on the box that matters.** Let
`Q = [-1, 2]^18`.
- `A_reg` is smooth on `R^18`.
- `supp A_reg ∩ Q` is compact.
- On `Q`, every derivative of `A_reg` is bounded, and every derivative of
  `F_omega` of order `j` is bounded by `C_j (1 + omega)^j`.

Every point the argument evaluates lies in `[0, 1]^18 ⊂ Q`, because raw
moments and their empirical versions are probabilities:
- population vectors;
- empirical vectors;
- segments between them.

The bounds are **not** global. Large `b_iX` is not cut off, and `H` grows
like `sqrt(b_iX)` there: at the frozen anchors, `H` is about 0.93 at
`b_3X = 1` and 10 at `b_3X = 100` (review E1). Nothing in the proof leaves
`Q`.

*Proof.* On the support of the product of cutoffs:
- `chi_anc != 0`, so the anchors are in the outer set. This is a compact
  set where the anchor states lie in `D_Sigma`, `u_1 < u_2 < u_3` with
  gaps of at least 0.096, and `G_A` is positive definite (Section 6).
- `w o Sigma != 0`, so the probe state is in `supp w`. In particular
  `u_X >= 0.1`, so `u_X - u_3 >= 0.248`.
- `chi_b != 0`, so all six `b >= b_flat / 2 > 0`. Hence every square root
  in `eta^(0)` is smooth there, and every `alpha(u_j - u_i)` is bounded
  below.

So `H` is smooth on a neighbourhood of the support, and the product is
smooth. Within `Q` the support is closed and bounded, so all derivatives
are bounded on `Q`. `u~` is smooth there, and `F` vanishes off the
support. Also, `alpha(L) = (1 - x cot x)/(2q)` with `x = sqrt(q) L <= 1`
is increasing in `L`. This is the monotonicity used for `b_low`.

**(R2) The background zero, for every anchor triple.** For every
background probe X and **every** anchor triple,

\[
A_{\rm reg}(\mu_0(X))=0.                                       \tag{2}
\]

*Proof.* If any cutoff factor vanishes, there is nothing to prove.
Otherwise:
- all six population `b` are positive. A population interval probability
  `b_ij` is positive exactly when `i` is chronologically before `j`; a
  null-related pair has `b = 0`. So `A1 < A2 < A3 < X`, chronologically
  and in that order;
- X is interior;
- `G_A` is positive definite.

The background identity then holds exactly (audit report 3a, item 1a):
`b_AX = alpha(u_X - u_A) slack^2 / gamma`, so `g_AX = y_A . y_X`. Three
plane vectors with `G_A` nonsingular satisfy `sum ell_i y_i = 0`, hence
`H = 0`.

No "good anchor" condition is used, and the probe need not be after anchor
3 for the bound to hold: if it is not, `b_3X = 0` and the cutoff vanishes.
The challenge's frozen statistic, which has only `w`, does not have this
property: a probe that is not after anchor 3 can have `H = 0.0029`
(audit report 3a, F1).

**(R3) The population signal is unchanged.** Under the alternative,
condition on an anchor triple whose population anchor numbers lie in the
three-quarter set. By Section 6, an accepted triple does so with
probability tending to one. There is `lambda_0 > 0` such that, for
`lambda < lambda_0`:
- for every probe whose wave population state lies in `supp w`,
  `A_reg(mu_lam(X)) = (w H)(mu_lam(X))`;
- for every probe, `|A_reg(mu_lam(X))| <= C lambda^2`.

*Proof.*
- `chi_anc = 1` on the three-quarter set.
- Wave population `b` differ from background ones by `O(lambda^2)`
  uniformly on the compact active and anchor sets, so they are at least
  `b_* - C lambda^2 >= b_flat`, where `chi_b = 1`.
- For the second bullet, `A_reg` is Lipschitz on `Q` and vanishes at
  `mu_0(X)`:
  - inside a compact neighbourhood of the active set,
    `|mu_lam - mu_0| = O(lambda^2)` (Section 4 and (8));
  - outside it, `A_reg` vanishes on a fixed neighbourhood of `mu_0(X)`.
    That neighbourhood contains `mu_lam(X)` once
    `sup_X |mu_lam(X) - mu_0(X)|` is small.

**Cited, not derived here.** The uniform closeness of wave and background
population moments for **every** probe, tips included, with anchors in a
compact interior set, is taken from audit report 3b (Sections 3.2 and 3.4,
`O(lambda^1.5)` near the tips). (8) itself is stated in Section 4 only for
active probes.

**(R4) Fixed and computable.** The cutoffs are fixed functions of the raw
moments with fixed numerical constants: `tau`, `psi`, `b_flat(q)` and `w`.
Known `q` enters `b_flat` and `Sigma`. Nothing is chosen after the data.

Every derivative up to a fixed order is bounded, with constants depending
on q, K and these compact sets. There is no assertion of uniformity as
q->0.

## 2. A bounded-variation lemma specific to this background

Fix (v,y) and vary u. In the background this curve is a future null generator.
Each past indicator increases at most once; each future indicator decreases
at most once. An anchor-interval indicator is a past indicator multiplied by
the constant 1_{A<R}. The raw ancestor-pair kernel increases at most once.
All therefore have uniformly bounded variation in u.

For the order-two kernel h, write its Hoeffding decomposition

\[
U_nh-Eh=\frac2n\sum_jh_1(R_j)+
\binom n2^{-1}\sum_{j<k}h_2(R_j,R_k).
\]

Taking an expectation of a bounded-variation function does not increase its
variation bound. Thus h_1 and the canonical h_2 retain uniform variation,
and h_2 is centered in each reference argument. Anchor kernels are constant
in the probe u. Products of any fixed number of these kernels also have
bounded variation. Multiplying by a smooth coefficient and the background
density preserves the bound. Integration by parts as a Stieltjes integral
gives, for every fixed reference-index configuration,

\[
\left|\int e^{-i\omega u}B(u)\,du\right|
\le\frac{|B(\hbox{left})|+|B(\hbox{right})|+\operatorname{Var}B}{\omega}.
                                                               \tag{3}
\]

(Here `B(u)` is a generic bounded-variation function. It is not the matrix
of Section 4.)

All integrations in v,y are over a fixed compact set. The probe weight is
strictly inside the diamond, so moving domain boundaries introduce no
uncontrolled surface term.

**Why (3) is needed.** The original said that independent noise at every
probe "would lose (3)". That is misleading: independent references for
each probe would give a smaller variance. The bound is needed because
shared references add a common-mode term across probes, and that is the
term (3) controls (audit report 3a, F3).

## 3. Taylor terms, repeated indices, and the remainder

In the order-k Taylor term of F around mu_0, the coefficient has degree at
most k-1 in omega. The degree-k term would multiply A_reg(mu_0), which is
zero by (2). For any monomial, expand each empirical error into its
order-one projection and its canonical order-two part. Suppose this
produces L reference slots in the product, k<=L<=2k, with normalization
comparable to n^(-L).

To bound its squared expectation, use two copies, with L+L' slots. If a
reference index appears only once across those slots, conditioning on all
other indices kills the term, even before the u integrals. Consequently a
surviving equality pattern has at most (L+L')/2 distinct indices. There are
only finitely many patterns at fixed k. Their normalization and counts give
at most C_k n^(-(L+L')/2), hence at most C_k n^(-k).

For each surviving configuration, (3) applies separately in the two probe
integrals. Thus the L2 norm of the integrated degree-k term is bounded by

\[
C_k\omega^{k-2}n^{-k/2},\quad\omega\ge1.                    \tag{4}
\]

This argument includes repeated indices within a product, bias terms and
cycles involving canonical order-two kernels. It does not assume independent
edges or discard the repeated-index terms. For k=1 it is 1/(omega sqrt(n));
for k=2 it is 1/n.

The remainder requires care: directly bounding the (K+1)-st derivative of
F would appear to cost omega^(K+1). Instead, on the segment mu_0+t epsilon,

\[
A_{\rm reg}(\mu_0+t\epsilon)
=t\int_0^1 DA_{\rm reg}(\mu_0+st\epsilon)[\epsilon]ds=tJ(t).
\]

Expand J(t) exp(-i omega u~(mu_0+t epsilon)) through degree K-1, then multiply
by t. Its order-K remainder is bounded by
`C_K omega^K |epsilon|^(K+1)`. Bounded-kernel U-statistic moment bounds give
`E|epsilon|^(2K+2)<=C_K n^(-(K+1))`; this can also be obtained by the same
index-centering count without invoking an asymptotic normal approximation.

**The probe average.** The anchor triple is fixed by conditioning. It is
independent of the references and the probes, and (2) holds for every
triple. Use `n_P >= c n` independent probes. Then the conditional variance
of the probe average is at most `E|F|^2 / n_P <= C E|epsilon|^2 / n_P <=
C / n^2`, because by (2) `|A_reg(mu_0 + epsilon)| <= C |epsilon|`. Its L2
contribution is `1/n`. For `K >= 2` that is the order of the `k = 2` term,
so it is absorbed there. For `K = 1` it is absorbed into the remainder
`omega / n`, since `omega >= 1`. Combining terms yields, conditionally on the anchor
triple and uniformly over triples:

\[
\|S_n(\omega)\|_2\le C_K\left[
\sum_{k=1}^K\omega^{k-2}n^{-k/2}
+\frac{\omega^K}{n^{(K+1)/2}}\right].                          \tag{5}
\]

For 1<=omega<=n^(K/(2K+2)), the sum is bounded by a constant times
`1/(omega sqrt(n)) + 1/n`, and the remainder is no larger than the first
term. The three-term candidate bound in the challenge is therefore valid on
the frequency range the test uses. No claim of that simplified bound for
arbitrary omega is needed.

**Elementary check of the mechanism.** Let X, R be uniform on [0, 1], with
`H_hat(X) = F_n(X) - X`. With n references and m independent probes the
exact variance is

\[
\Big(1-\frac1m\Big)\frac{1-\operatorname{sinc}^2(\omega/2)}{n\,\omega^2}
+\frac1{6nm},\qquad \operatorname{sinc}x=\frac{\sin x}{x}.
\]

At `omega = 2 pi j` the sinc term vanishes and this is the original's
formula, which is exact only there (audit report 3a, F5; checked by
quadrature in `checks.py`). It follows by integrating the empirical-CDF
error once by parts. It explains both leading terms of (5), but is not a
validation of the nonlinear causal-set instrument.

## 4. Canonical density and kernel coupling, including the tips

**Definition of B.**
- `B(s, t)` is the trace-free symmetric 2 by 2 matrix whose complex
  representative is
  `beta(s,t) = (2 lambda k / c_lam) e^{i phi} int_s^t f(u') e^{i u'/lambda} du'`,
  with `k = sqrt(q)` and `c_lam = sqrt(1 + 4 q lambda^2)`. It is not
  positive: being trace-free, it has eigenvalues `±|beta|`.
- The exact wave action is `Mt(s, t) = Delta m(s, t) I - B(s, t)`, and
  `|beta| <= C lambda min(Delta u, lambda)`.
- At a point with null coordinate u, put `B_left = B(-1/2, u)`,
  `B_right = B(u, 1/2)`, `M_left = m I - B_left` and
  `M_right = (M - m) I - B_right`.

Let s be distance in u to the nearer tip. The oscillatory integral bound
above implies a relative matrix error at most

\[
e(s)=C\lambda\min(1,\lambda/s).
\]

The canonical slice map y=T_lambda(u)z has
`T=(W_0(M_left^-1+M_right^-1))^(-1/2)`. Its first-order perturbation is
trace-free. Matrix Taylor expansion, uniformly where e<=C lambda, gives

\[
|\det T-1|\le C e(s)^2.                                    \tag{6}
\]

The background slice volume is proportional to W_0~s. For s<lambda the
integrated error is O(lambda^2 int_0^lambda s ds)=O(lambda^4); for s>lambda
it is O(lambda^4 int_lambda^1 ds/s). Normalization changes by the same order.
Thus the canonical one-point laws satisfy

\[
\|\nu_\lambda-\mu_0\|_{TV}\le C\lambda^4\log(1/\lambda).
                                                               \tag{7}
\]

This is stronger than the common-physical-chart O(lambda^2) estimate because
the canonical map preserves both tip slacks. It does not follow from an
interior-only expansion.

For later order-kernel estimates, the transverse shift obeys
`|delta y|<=C lambda min(sqrt(s),lambda/sqrt(s))`. Its L1 average is
O(lambda^2). The longitudinal shift is uniformly O(lambda^2). For the latter,
use either
`delta v=(z^T T M_left^-1 T z-|z|^2/m)/2`
or, from T(M_left^-1+M_right^-1)T=W_0^-1 I, the corresponding right expression.
Choose the expression with the long denominator; |z|^2<=Cs removes the tip
singularity. These are uniform bounds, not merely formal first variations.

Two independently sampled reference points have background-cone disagreements
under this coordinate change on O(lambda^2) mass. Indeed the Minkowski null
margin is Lipschitz in both points, and its conditional strip probability is
uniformly at most C times the strip width. Split the crossing event into the
strip attributable to either endpoint, condition on that endpoint, and
integrate its displacement. The wave-versus-background cone difference in
the common physical chart adds O(lambda^2), by the existing uniform sandwich.

The same estimate holds with one endpoint fixed in the compact anchor/probe
region. Split the other endpoint into an interior region, where all shifts
are O(lambda^2), and the tips. At a tip its u-separation from the fixed point
is bounded below; the action is then uniformly Lipschitz in transverse
coordinates and longitudinal slack. Its strip width is bounded by
`C lambda^2/sqrt(s)`, which is integrable. Hence, uniformly for active X,

\[
\Pr\{K_\lambda(X,R,R')\ne K_0(X,R,R')\}\le C\lambda^2       \tag{8}
\]

for each raw indicator kernel. For ancestor pairs, apply a union bound to
the two ancestor indicators and the reference-reference comparison. Couple
the reference locations using (7); its additional mismatch probability is
smaller than (8).

## 5. The wave-side comparison, written as an identity

Let `epsilon_lam` and `epsilon_0` be centred errors from the coupled
reference lists of Section 4. For every fixed `p >= 2`, (8), boundedness
and the independent-pair representation of an order-two U-statistic imply,
when `n lambda^2 >= 1`,

\[
\|\epsilon_\lambda(X)-\epsilon_0(X)\|_p\le C_p\lambda/\sqrt n.
                                                               \tag{9}
\]

For completeness, a bounded centered difference kernel with variance
O(lambda^2), averaged over independent pairs, has p-th norm bounded by
`C_p[lambda/sqrt(n)+lambda^(2/p)/n^(1-1/p)]`. The second term is no larger
when n lambda^2>=1. Averaging over pairings and Jensen transfers this bound
to the U-statistic. The order-one case is the same calculation. Only
`p = 2, 4` are used. In every block of the test `n lambda^2 -> infinity`,
because `alpha > 2`.

Let X range over wave probes, with law `nu_lam` in canonical coordinates.
With `F = F_omega` from (1b), so that `F(mu_0) = 0`, the following is an
exact identity:

\[
\begin{aligned}
S_n-E_XF(\mu_\lambda(X))
&=\underbrace{\big[\mathrm{avg}_XF(\mu_\lambda)-E_XF(\mu_\lambda)\big]}_{T_1}
+\underbrace{\mathrm{avg}_X\big[F(\mu_\lambda+\epsilon_\lambda)-F(\mu_\lambda+\epsilon_0)\big]}_{T_2}\\
&\quad+\underbrace{\mathrm{avg}_X\big[F(\mu_\lambda+\epsilon_0)-F(\mu_\lambda)-F(\mu_0+\epsilon_0)\big]}_{T_3}
+\underbrace{\mathrm{avg}_XF(\mu_0(X)+\epsilon_0(X))}_{T_4}.
\end{aligned}
\]

The derivatives satisfy, on `Q` (by (R1); every point used lies in
`[0, 1]^18`), `|DF| <= C(1 + omega|A_reg|)` and
`|D^2F| <= C(1 + omega + omega^2 |A_reg|)`. Throughout,
`omega = 1/lambda + O(1)`. Condition on an anchor triple as in (R3).

- **T1.** `|F(mu_lam)| <= C lambda^2` by (R3), with independent probes:
  `O(lambda^2 / sqrt(n))`.
- **T2.**
  - On the segment, `|A_reg| <= C(lambda^2 + |epsilon_lam| + |epsilon_0|)`.
  - So the term is at most
    `C(1 + lambda + omega(|epsilon_lam| + |epsilon_0|)) |epsilon_lam - epsilon_0|`.
  - By (9) with `p = 2, 4`, and `omega ~ 1/lambda`, this is
    `C(lambda/sqrt(n) + 1/n)`.
- **T3.**
  - It equals
    `int int D^2F(mu_0 + t delta + s epsilon_0)[delta, epsilon_0]`, with
    `delta = mu_lam - mu_0 = O(lambda^2)` and
    `|A_reg| <= C(|delta| + |epsilon_0|)` there.
  - So it is at most
    `C(omega + omega^2 lambda^2 + omega^2 |epsilon_0|) lambda^2 |epsilon_0|`.
  - At `omega = 1/lambda + O(1)` that is
    `C(lambda |epsilon_0| + |epsilon_0|^2) = O(lambda/sqrt(n) + 1/n)`.
- **T4.**
  - The conditional probe fluctuation is `O(1/n)`.
  - Changing the probe law from `nu_lam` to `mu_0` costs at most
    `C lambda^4 log(1/lambda) / sqrt(n)`, by (7).
  - The remainder `E_{X ~ mu_0} F(mu_0(X) + epsilon_0(X))` is bounded by
    (5).

Only T4's remainder uses the oscillatory gain, and it is entitled to it:
- background kernels of i.i.d. references;
- probes integrated against `mu_0`;
- the background phase.

Wave kernels are never asked for the gain. The wave-specific terms T2 and
T3 sit at `lambda/sqrt(n)`, the order of the leading background term, so
(8) has no slack against the exponent 2.

**Probes outside the active set.** The score sums over all probes. For X
outside a compact neighbourhood of the active set, `A_reg` vanishes on a
fixed neighbourhood of `mu_0(X)`. `mu_lam(X)` lies in it for small
`lambda`; that uniform closeness is cited from audit report 3b, as in
(R3). T2 and T3 are then nonzero only if `|epsilon| >= const`. By
Hoeffding that has probability `exp(-c n)`, uniformly in X.

Together with (5), at the target frequency the total score error is

\[
O_{L^2}\left(\lambda/\sqrt n+1/n+
\omega^K n^{-(K+1)/2}\right).                              \tag{10}
\]

There is **no** assumption that all N coupled sample coordinates agree, and
hence no unnecessary restriction `N lambda^4 log -> 0`.

## 6. The signal survives intrinsic anchor selection

For general anchor vectors a_i and ell_i with sum ell_i a_i=0, define

\[
D_A(u)=-\frac12\sum_{i=1}^3
\frac{\ell_i(a_{i1}-ia_{i2})^2}{m(u)-m_i}.
\]

This is only a proof coordinate system. Rotate it away using the anchor Gram
matrix; the test never receives it. For the frozen triple, it reduces to D
in the challenge, whose imaginary part is strictly positive. Its complex
integral is therefore nonzero. The same holds uniformly in a sufficiently
small open neighbourhood, including its rotations and reflections. Reflection
reverses the complex helicity representation; it changes which conjugate
coefficient a positive frequency selects, while preserving its magnitude.

Such a neighbourhood is specified intrinsically by u_i,v_i,R_i^2 and g_ij.
For a concrete compact neighbourhood, use outer tolerances 0.002 for the
three u_i and v_i, and 0.00002 for all Gram entries, around the frozen
triple. This is the outer set of `chi_anc`.

**The bound on `a_2x a_2y`, corrected.**
- In the orientation `a_1 = (sqrt(G11), 0)`, `a_2` has positive second
  coordinate, `a_2x = G12 / sqrt(G11)` and
  `a_2y = sqrt(G22 - G12^2 / G11)`. So

\[
a_{2x}a_{2y}=G_{12}\sqrt{G_{22}/G_{11}-G_{12}^2/G_{11}^2},
\]

  which equals `G12` only if `G22 = G11`.
- On the box `G11, G22` in `[0.00638, 0.00642]`, `|G12| <= 0.00002`, its
  supremum is `2.0063e-5`, not `2e-5`.
- The elementary Gram formulas still bound `a_3`'s two coordinates below
  0.079, and `|ell_2| <= 1.0094` (the certificate's bound; the original
  quoted the rounder `1.02`).
- For `q <= 1`, `m(u) - m_3 < 1.30 * 0.452` on `0.1 <= u <= 0.3`, and
  `m(u) - m_2 > 0.348`.
- So `Im D_A(u) >= 0.0106978 > 0.01` throughout that interval. With
  `1.02` in place of `1.0094` the same chain gives `0.0106972`, also above
  0.01.

The frozen anchor-guard certificate (`f1072f0`) records
`0.01069797`, overstated by `1.8e-7` because it used `2e-5`. It is not
edited; the corrected value is the one above. Every target survives:
- `Im D_A > 0.01`, with margin `7e-4`;
- the grid coefficient `0.0101679 > 0.009`;
- `|D_A| < 0.05`, which uses only the diagonal Gram entries.

**Anchor and probe slacks, defined.** By the definition of `g_ij`, the
order-derived slack is

`Delta v - (R_i^2 + R_j^2 - 2 g_ij) / (2 Delta m) = eta^(0)_ij = sqrt(b_ij gamma / alpha)`.

- The certificate bounds this below by `0.0289` on the outer set. So
  `b_ij > 0`, and every anchor pair is timelike. This is what "anchor pair
  slacks stay positive" means.
- Active probes have slacks of at least `0.1002`, so they are after all
  three anchors.
- The anchor Gram block `G_A` stays positive definite.

These bounds can be checked with rational interval arithmetic in the full
instrument; they are not numerical results of this run.

Take disjoint candidate triples from O(log n) candidate events. Pilot moment
estimates use disjoint O(n) references. Accept the first triple in an inner
neighbourhood with half these tolerances.
- A still smaller neighbourhood has positive physical probability.
- Bounded U-statistic concentration and a union bound over O(log n)
  candidates imply that some triple is accepted with probability tending
  to one.
- Every accepted triple has its population anchor numbers in the
  three-quarter set, where `chi_anc = 1`, with probability tending to one.
- Wave population moments converge uniformly to the canonical background
  moments there, so the same conclusion holds under the alternative.

This argument does not supply a practical sample size: the good-set
probability and inverse-calibration constants may be very small. Under the
null the acceptance event is not needed at all, by (R2).

The population expansion in the challenge, including angular integration,
together with (R3), now yields uniformly over these anchors and wave phase

\[
|E_\lambda[F_\omega(\mu_\lambda(X))]|=|E_\lambda[wH e^{-i\omega u}]|\ge c_q\lambda^2,
\quad |\omega-1/\lambda|\le1/8,                             \tag{11}
\]

for sufficiently small lambda. The O(1) frequency offset cannot cancel the
integral: |D_A|<0.05 on the same outer set, and centering the Fourier phase at
u=0.2 bounds its change by (0.1/8)*0.05<0.001, below the imaginary lower bound
0.0106978. Smooth terms involving anchor phases vanish on integration by
parts; one-point calibration errors are o(lambda^2). Taking the magnitude
respects the previous even-likelihood symmetry and does not resurrect a
forbidden signed, orientation-fixed score.

## 7. A scan with no wavelength input or unknown threshold constant

Here a simpler asymptotic construction replaces the optional Monte Carlo
calibration discussed in the challenge. It is a new derivation, not a
re-scoring of that population experiment.

Randomly relabel events and split them into L=ceil(16 log N) equal blocks,
discarding the remainder. In each block reserve 3 ceil(log m) candidates,
then divide the other m events between pilot references, estimation references
and probes. Put n equal to the smallest of these three list sizes, so that
`n_P >= n`. If the block is too small or has no accepted anchor triple, its
vote is zero.

For fixed K use the same deterministic grid in every block:

\[
\Omega_n=\{j/4:4\le j\le\lfloor4n^{K/(2K+2)}\rfloor\}.
\]

At each frequency let a block vote if its score magnitude exceeds

\[
t_n(\omega)=\log n\left[
1/(\omega\sqrt n)+1/n+\omega^K/n^{(K+1)/2}\right].            \tag{12}
\]

Reject the background if a strict majority votes at any common grid point.
All quantities in this definition are order counts, known q, deterministic
functions or independent relabeling randomness. There is no true wavelength,
phase, embedding, labelled generator anchor, or coordinate oracle input.

**Under the null.**
- By (5), which holds conditionally on any anchor triple, and by Markov's
  inequality, each fixed-frequency block vote has probability at most
  `(K C_K / log n)^2`.
- No acceptance event is needed: a block with no accepted triple votes
  zero, which only lowers the vote probability.
- So, uniformly on the grid, the vote probability is eventually below 1/4.
- Blocks are independent, so Hoeffding gives a majority probability of at
  most `exp(-L/8) <= N^-2`.
- The grid has fewer than `4 sqrt(N)` points. A union bound therefore makes
  the **whole scan**'s size tend to zero.

The logarithm in (12) eventually dominates the unknown fixed `C_K`; this is
an asymptotic statement, not a finite-N significance guarantee.

**Under the alternative, given (1).** `1/lambda` is eventually strictly
below the scan cap, and a grid point is within 1/8 of it. Equations (10) to
(12) give a signal-to-noise ratio tending to infinity, including the
logarithmic threshold and block splitting. The block vote tends to one, and
the majority rejects with probability tending to one. This proves the
proposed implication (1) for the statistic (1a) and (1b), conditional on
the lemmas as stated.

**Remark (scan cap).** The cap `n^(K/(2K+2))` is the only binding
requirement on the `2/K`. A grid up to `sqrt(n)` keeps the null control by
(5). Given (10), it would yield `alpha > 2 + 2/(K+1)` (audit report 3a,
F7). The statement (1) is left as it is.

Finite precision must also be part of the instrument: compute the regularized
score with a certified absolute error, for example n^-4, and vote only when its
lower magnitude bound exceeds (12). This error is negligible compared with
the threshold throughout the scan range. Fixed double precision and a fixed
root tolerance do not implement this asymptotic requirement.

The concrete geometric guards were frozen with the original document and
verified by the anchor-guard certificate at `f1072f0`. It used Fraction
interval arithmetic and the certified cosine series at q=1, whose derivative
bound controls every 0<q<=1.
- **The `Im D_A` correction was not re-certified.** Section 6 computes it
  by hand from the certificate's recorded terms; it has not been redone in
  interval arithmetic.
- **The rule for any re-certification stands,** as the original put it:
  it must retain a failure without enlarging the neighbourhood or weakening
  its target inequalities.

This algebraic control is not a noise-rate or finite-data-power test.

## 8. Status, remaining work and prior art

- **What this revision completes.** Argument 3's statistic is now written
  in full, as one object for which (5), (10) and (11) hold. So are its
  properties (R1) to (R4), the wave-side chain (the identity T1 to T4),
  and the corrected constants. The review in
  [independent_audit/argument3_revision/](independent_audit/argument3_revision/)
  graded this CORRECT WITH MINOR FIXES, which are applied (Section 9).
- **What is still at outline level**, as audited, with the referee reports
  filling in the steps:
  - the index-counting proof of (5);
  - the tip estimates (6) to (9);
  - the uniform closeness of wave and background moments for every probe,
    which is cited;
  - the population expansion behind (11), which rests on the challenge's
    predicted limit (derived by referee 3b) and on two assertions (anchor
    phase terms vanish; calibration errors are `o(lambda^2)`);
  - the anchor-guard certificate, which is the author's code. One of its
    steps was wrong and was corrected here by hand, not re-certified;
  - the standard U-statistic moment and Hoeffding bounds, cited.
- **On the bracket.** `4/3 <= alpha_low <= alpha_high <= 2` is the audit's
  statement (AUDIT.md Section 2), not something new here. Its lower end is
  not newly written out: `quantized_order_bound.md` has not been revised
  since the audit, so the audit's items for argument 1 are still missing
  from its text. Among them is linking the constants, without which the
  claim fails at `L = 1`.
- **What it does not change.**
  - All statements are asymptotic. The constants put any finite-N test out
    of reach (AUDIT.md Section 5).
  - Nothing here shows that 2 is necessary.
- **The next required deliverable for an instrument**, as before: the
  complete order-only count implementation with fixed regularization,
  anchor guards, certified increasing precision and the block scan. Its
  nontrivial positive controls must exercise the same data path; a zero
  score on a small sample that never found anchors is not evidence of
  detection power. Any numerical accuracy or power campaign must be sized
  and frozen before sampling.
- **The sharp lower bound remains a separate obstacle.** This proof does
  not exclude more informative geometry-fitting instruments between powers
  4/3 and 2.

Order-derived embeddings, distances and Gram algebra are established tools,
not our novelty claim: see [Henson (2006)](https://arxiv.org/abs/gr-qc/0601069)
and [Johnston (2025)](https://arxiv.org/abs/2502.09701). The possible contribution
would be the oscillatory statistical bound in this exact high-frequency
experiment and eventually a matching resolution law. A bounded search on
2026-09-24 did not retrieve this complete result, which is not a novelty proof.
For general U-statistic background see
[Gine, Latala and Zinn (2000)](https://arxiv.org/abs/math/0003228); no theorem
there has been substituted for the model-specific variation and coupling
arguments above.

## 9. Review and the corrections made after it

[independent_audit/argument3_revision/review.md](independent_audit/argument3_revision/review.md).
An AI referee in a fresh context, given this revision, the frozen original,
the challenge, AUDIT.md, both referee reports on the original and the
certificate. It ran only short deterministic checks.

**Verdict: CORRECT WITH MINOR FIXES.** `A_reg` (1a) is one well-defined
object, and (5), (10) and (11) hold for it. The identity T1 to T4 is
exact. AUDIT.md's items 1 to 7 for argument 3 are all addressed. Text that
Section 0 does not mark as changed is word-identical to the original, and
the original still matches the certificate's hash.

Corrected after the review, in the text above:

- **E1, an error as stated.** (R1) claimed compact support and globally
  bounded derivatives on `R^18`. That is false, because large `b_iX` is
  not cut off and `H` grows like `sqrt(b_iX)`. It now states smoothness on
  `R^18` and compact support and bounded derivatives on
  `Q = [-1, 2]^18`, which contains every point the proof uses. Section 5
  reads "on `Q`" instead of "globally". Nothing downstream changes.
- **M1.** Section 0 lists every change, including the preamble's status
  sentences. The re-certification rule of the guard paragraph is restored.
- **M2.** The cutoffs are named: the standard `exp(-1/x)` smooth step.
- **M3.** `|ell_2| <= 1.0094`, which gives 0.0106978. The original's 1.02
  gives 0.0106972. Both exceed 0.01.
- **M4.** For `K = 1` the probe term is absorbed into the remainder
  `omega / n`.
- **M5.** "Chronologically related, in that order", in (R2).
- **M6.** (R3) and Section 5 now state that they condition on a triple in
  the three-quarter set. They also say that the uniform closeness of wave
  and background moments for every probe, tips included, is cited from
  audit report 3b, not derived here.
- **M7.** Notation: `c_lam`, `alpha(.)` against the exponent `alpha`, and
  `K` against the kernels `K_lambda`.
- **The Section 8 headline**, which overclaimed, now says what is written
  in full and what is still outline. It also notes that the bracket's lower
  end (argument 1) has not been revised in its own text.

The referee's least-trusted remaining step is the coupling uniform in X.
- T2 and T3 sit exactly at `lambda/sqrt(n)` by (8), with no slack against
  the exponent 2.
- The "every probe" extension rests on a citation.
- (8) at the tips is still the original's sketch. Its proof is in audit
  report 3b, Section 3.2.