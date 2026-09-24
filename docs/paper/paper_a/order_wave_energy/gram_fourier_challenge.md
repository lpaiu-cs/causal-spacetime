# Anchor Gram identities and an oscillatory order statistic

2026-09-24. A new candidate after the linear cone-residual obstruction.
The established bracket remains exponents below 4/3 versus sufficient 3.
The proposed route could cover every fixed sampling power above 2, but that
improvement is **not yet claimed**: the stochastic and observation-removal
steps below must be completed and audited. The sharp boundary remains the
full objective, even if this candidate supplies another improved upper bound.

## 1. A pointwise identity expressed entirely in order probabilities

Keep fixed known 0<q<=1. Use the previous p,f,a moments, with c=a/p^2,
and recover u, m, eta, zeta as in the null-coordinate construction. In addition
to v recover the transverse squared radius

\[
R_X^2=\frac{2m_X(M-m_X)}M(1-\eta_X-\zeta_X).
\]

For a timelike pair A<X, the order probability
`b_AX=Pr(A<Z<X)` is its interval volume divided by the whole diamond volume.
With alpha and gamma=alpha(1) from the background, set

\[
\eta_{AX}^{(0)}=\sqrt{b_{AX}\gamma/\alpha(u_X-u_A)},\quad
g_{AX}=\frac{R_A^2+R_X^2}{2}
 -(m_X-m_A)(v_X-v_A-\eta_{AX}^{(0)}).
\]

Under the background, g_AX is exactly the transverse inner product y_A.y_X.
Choose three mutually timelike anchors A1<A2<A3 and subsequent probes X.
Build the anchor Gram matrix G with diagonal R_A^2 and off-diagonal g_AB.
If its leading 2 by 2 block B is positive definite, set

\[
\ell=(-G_{3,\{1,2\}}B^{-1},1),\qquad
H(X)=\ell_1g_{A_1X}+\ell_2g_{A_2X}+g_{A_3X}.
\tag{1}
\]

In two transverse dimensions, H(X)=0 for **every** background probe.
Only order probabilities enter (1); no spatial orientation is supplied.
This identity is distinct from inserting a coordinate-oracle forcing into
a conditional expectation. Its empirical error comes from shared reference
counts, and the detector will use the magnitude of an oscillatory average.

## 2. The candidate statistic

Split events into anchor candidates, pilot references, estimation references,
and probes. Select a well-conditioned anchor triple using only pilot moment
estimates. Estimate p,f,a and all required b_AX with the same estimation
references. For smooth compactly supported w define

\[
S_n(\omega)=\frac1{n_P}\sum_X
 w(\widehat u_X,\widehat v_X,\widehat R_X^2)
 \widehat H(X)e^{-i\omega\widehat u_X}.
\tag{2}
\]

The proposed test uses |S_n| (or its square), not a signed component in an
external transverse frame. Wavelength and phase must eventually be removed
by a deterministic frequency grid and a whole-test calibration. A score with
fixed supplied anchor indices is only a component of that complete observer.

The crucial candidate noise bound, for each fixed integer K, is

\[
\|S_n(\omega)\|_{L^2(P_0)}\le C_K\left[
\frac1{\omega\sqrt n}+\frac1n+
\frac{\omega^K}{n^{(K+1)/2}}\right].
\tag{3}
\]

Reason to investigate (3): at fixed (v,y), varying u follows a background
null generator. Every reference past/future/interval indicator has at most
one jump in u. Products of finitely many centered kernels have uniformly
bounded variation. Integration against exp(-i omega u) therefore supplies
1/omega. In a Taylor expansion of (2), the order-k term has at most
omega^(k-1), because the population residual H is identically zero. Expanding
its squared moment in reference indices, centering kills a term with an
index occurring only once. This suggests order-k noise bounded by
`C_K omega^(k-2) n^(-k/2)`; the independent probe average contributes 1/n.

This is a proof programme, not an invocation of a generic U-statistic bound
that automatically proves (3). The order-two ancestor kernel, repeated-index
terms, smooth cutoff, anchor-selection conditioning, and Taylor remainder
must all be accounted for. Moment inequalities for canonical U-statistics
are relevant background, e.g.
[Gine, Latala and Zinn (2000)](https://arxiv.org/abs/math/0003228).

## 3. A nonvanishing population Fourier coefficient

Use the canonical slice coordinates (u,v_tilde,z) from the earlier proof.
To distinguish signs from that document, write
`B_lambda(u)=(2 lambda sqrt(q)/c) J(-1/2,u)`, so
`Mtilde(s,t)=Delta m I-[B_lambda(t)-B_lambda(s)]`.
On a compact interior interval its complex trace-free representation is

\[
B_\lambda(u)=-2i\sqrt q\lambda^2 f(u)
 e^{i(u/\lambda+\phi)}
 +\text{a smooth endpoint term}+O(\lambda^3),\quad f=a^{-2}.
\]

Take canonical anchor transverse vectors z1=(d,0), z2=(0,d), z3=(d,d),
at three distinct past u values. Then ell=(-1,-1,1) at leading order.
For probes averaged over transverse angle, the local high-frequency part
of g_AX is

\[
-\frac{z_A^TB_\lambda(u)z_A}{2(m(u)-m_A)}.
\]

Terms from the canonical coordinate change either have zero angular average
or depend on anchor phases with a smooth coefficient in probe u. The change
of ell also multiplies a background inner product whose angular mean is zero.
Consequently, putting

\[
D(u)=\frac{d^2}{2}\left[\frac1{m(u)-m_1}-\frac1{m(u)-m_2}\right]
       +i\frac{d^2}{m(u)-m_3},
\]

the local part of the angle-averaged H is Re[D(u) B_lambda(u)]. For separable
radial weight w_u(u)w_v(v)w_t(|z|^2), divide the population S by its weight
mass. The predicted limit is

\[
\boxed{\frac{S(1/\lambda)}{\lambda^2\,\mathbb E w}
\longrightarrow -i\sqrt q e^{i\phi}
\frac{\int w_u(u)D(u)\,du}{\int w_u(u)a(u)^2\,du}.}
\tag{4}
\]

The imaginary part of D has one strict sign on the probe interval, so the
coefficient cannot vanish. Rotating the canonical anchors by theta changes
the displayed phase to phi-2 theta without changing its magnitude. Thus an
open set of anchor triples has a signal at least c lambda^2 uniformly in
phase. Such triples can in principle be selected through their order-derived
Gram entries; their true coordinates must not be passed to the observer.

## 4. What remains before claiming exponent 2

If (3) and its alternative counterpart are valid, fixed-frequency noise is
O(1/(omega sqrt(n))) on successively larger ranges
`omega<=n^(K/(2K+2))`. Signal lambda squared then dominates at omega=1/lambda
for every fixed `n~lambda^(-alpha)` with alpha>2, choosing K large enough.

The following steps are required, not optional qualifications:

1. Complete the centered-kernel expansion proving (3), including every
   repeated-index term and a valid global smooth extension of the calibration.
2. Bound the wave-side error. A possible route is the canonical-slice coupling:
   its one-point density differs by O(lambda^4 log(1/lambda)) in L1, while
   individual order kernels can differ on O(lambda^2) mass. Prove these
   statements with endpoint weights, rather than assuming uniform interior
   estimates extend to the tips. For 2<alpha<3 the full-sample density cost
   would vanish and centered kernel differences would have size lambda/sqrt(n).
3. Specify an order-only anchor-selection rule with a positive-probability
   compact good set, pilot certification, and success probability tending to 1.
4. Remove frequency input and control the whole test. One possible construction
   uses O(log N) independent blocks, constant-size **independent null Monte
   Carlo rank tests** at each frequency within each block, then majority voting
   over a polynomial grid. Exchangeability supplies finite-sample per-block
   size without unknown derivative constants; the block probability bound
   must precede the grid union bound. Do not claim that a single-frequency
   variance estimate controls a scan.
5. Implement that complete observer with increasing numerical precision and
   perform appropriate finite-data validation. Population quadrature alone
   cannot satisfy this step or establish the sharp boundary.

The contrast between an even likelihood and a complex signal is not a
contradiction: the anchor-relative signal has an unknown phase; its modulus
uses even joint information. The actual finite-data distribution, not a
nonzero oracle residual, must provide the discrimination claim.

## 5. Frozen population challenge

Freeze this document and the corresponding source before new geometry runs.
Use CPU, the existing WaveFixture and scalar hierarchy, and no new dependency.

- q=1, r in `{64,128,256,512}`, metric phase 0.37. Canonical anchors have
  u=(-0.35,-0.25,-0.15), v=(0.15,0.25,0.35), transverse vectors above with
  d=0.08. Rotate all transverse anchor vectors by theta in `{0,pi/8}`:
  eight rows, retained without reselection.
- The weight is the product of `exp(-1/(1-x^2))` on |x|<1, zero otherwise,
  with centers/half-widths u=(0.2,0.1), v=(0.7,0.05), radius squared
  (0.02125,0.01875). Apply the weight and Fourier phase to **reconstructed**
  values. Integrate over the encompassing canonical ranges u=[0.05,0.35],
  v=[0.60,0.80], radius squared=[0.0001,0.06].
- Gauss orders are max(96,ceil(r/4)+48) for u and 16 for v and radius squared.
  Four transverse angles integrate the quadratic Gram residual exactly at
  each u,v,radius; the reconstructed weight is angularly independent.
- Record the complex normalized S/lambda^2, prediction (4), weight mass,
  anchor conditioning, and minimum causal slack over nonzero-weight probes.
  Record failure of a calibration or anchor regularity without changing the
  points or filtering out the case.
- At r=512,theta=0 tighten the primitive and double all quadrature orders.
  Require less than 1% relative change for this numerical consistency check.
  There is no preregistered finite-data accuracy or power claim.
- Algebraic self-checks reconstruct background states and Gram inner products
  from exact population moments, verify the zero residual, and compare angular
  averaging of a quadratic expression with a dense angular grid.

The immediate question is whether (4), including its complex phase and
normalization, survives the exact geometry calculation. A pass does not yet
upgrade the established sufficient region; the proof and observer work in
Section 4 remain necessary.
