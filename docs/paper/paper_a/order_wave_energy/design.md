# Can unresolved vacuum waves leave an intrinsic energy signal?

2026-09-23. **Pre-execution derivation and deterministic challenge.** This is a
new discovery-oriented question, not an extension of the frozen P14 accuracy
verdict. Written and implemented by the same AI; no independent proof audit.
Success, priority, and publication readiness are undecided.

The approved question is whether finite causal observations retain a universal
effective energy and a predictive resolution law after losing microscopic wave
shape. Merely reproducing Isaacson--Burnett backreaction is not a discovery.
If the stronger question fails, salvaging a small safe publication is not the
default objective. The first attack below tests the proposed mechanism; it does
not certify a finite-sample instrument or a gravitational stress tensor.

## 1. Exact vacuum family, with energy held finite

On a positive Rosen patch in four dimensions, take

\[
g_\lambda=-2du\,dv+a_\lambda^2 S_{\lambda,AB}dx^A dx^B,
\quad S=S^T>0,\quad\det S=1.
\]

For \(\gamma=a^2S\), direct contraction of the connection gives

\[
R_{uu}=-\tfrac12\operatorname{tr}(\gamma^{-1}\gamma'')
 +\tfrac14\operatorname{tr}[(\gamma^{-1}\gamma')^2]
=-2a''/a-\tfrac14\operatorname{tr}[(S^{-1}S')^2].
\]

Indeed, \(\gamma^{-1}\gamma'=2(a'/a)I+T\), with
\(T=S^{-1}S'\), \(\operatorname{tr}T=0\), and
\(R_{uu}=-\tfrac12(\operatorname{tr}\gamma^{-1}\gamma')'
-\tfrac14\operatorname{tr}[(\gamma^{-1}\gamma')^2]\).
All other Ricci components vanish. Thus choose

\[
a''=-Q_\lambda a,\qquad Q_\lambda=\tfrac18\operatorname{tr}(T^2).
\tag{1}
\]

For linear polarization, \(S=\mathrm{diag}(e^{2\beta},e^{-2\beta})\),
\(Q=\beta'^2\). Set \(\beta_\lambda=\lambda b f(u/\lambda+\phi)\),
where \(f\) is smooth periodic and \(q=b^2\langle f'^2\rangle\).
The bounded oscillatory forcing in (1) averages to \(q\). With fixed common
initial data on a compact positive patch, integral-equation averaging gives
\(a_\lambda\to a_0\), \(a_0''=-qa_0\). The limit metric has

\[
g_0=-2du\,dv+a_0^2(dx^2+dy^2),\qquad R_{uu}[g_0]=2q.
\tag{2}
\]

This is the known high-frequency vacuum/nonvacuum-limit mechanism, not a new
Einstein solution principle. The metric oscillation is \(O(\lambda)\); its
curvature is generally \(O(\lambda^{-1})\). The existing **fixed curvature
amplitude** high-frequency obstruction is a different scaling and does not
forbid (2). In the Einstein interpretation with \(G_{ab}=8\pi G T_{ab}\),
\(T_{uu}^{\rm eff}=q/(4\pi G)\); its attribution to waves uses knowledge of
the underlying vacuum family, not causal order alone.

For a polarization check, let
\(S=R(\theta)\mathrm{diag}(e^{2\beta},e^{-2\beta})R(\theta)^T\).
Then
\[
Q=\beta'^2+\sinh^2(2\beta)\theta'^2.
\]
The circular family
\(\theta=(u/\lambda+\phi)/2\),
\(\beta=\tfrac12\operatorname{asinh}(2\lambda\sqrt q)\)
has \(Q=q\) exactly, hence \(a_\lambda=a_0\) at every frequency.

## 2. An observable using only causal relations

For \(s<t\), define
\[
M(s,t)=\int_s^t\gamma^{-1}(u)du,\qquad
K(s,t)=a(s)a(t)\sqrt{\det M(s,t)}.
\tag{3}
\]
The null action is
\(\tfrac12(x_t-x_s)^T M(s,t)^{-1}(x_t-x_s)\).
Equation (3) is the square root of the Brinkmann Jacobi determinant: endpoint
zweibeins contribute \(a(s)a(t)\). The positive Rosen transverse action also
ensures there is no conjugate point inside the chart; an arbitrary large
diamond outside that patch is not covered.

Reuse the already derived [exact chain hierarchy](../window_response/plane_wave_chain_hierarchy.md):
\(C_1=K\), \(C_{m+1}(s,t)=\int_s^t C_m(s,u)K(u,t)du\). In four dimensions,
\[
V=\pi h^2 C_2/K,\qquad
P_2=\frac{KC_3}{3C_2^2},\qquad
P_3=\frac{K^2C_4}{15C_2^3}.
\tag{4}
\]
Here \(P_2\) is the probability that two uniform interior events are related;
the finite-order estimator is simply the relation count divided by
\(\binom N2\). Coordinates, profile, and energy are absent from that statistic.
They remain in the forward generator used to validate it. No blind inverse
campaign has yet been run.

In the limiting background, \(K_0(s,t)=\sin(\sqrt q(t-s))/\sqrt q\).
For \(L=t-s\),
\[
C_m=\frac{L^{2m-1}}{(2m-1)!}\,
{}_0F_1(m+\tfrac12;-qL^2/4),
\quad
P_2=\frac1{10}\left[1-\frac4{105}qL^2+O(q^2L^4)\right].
\tag{5}
\]
The same expression follows by taking the inverse Laplace transform of
\((p^2+q)^{-m}\); \(C_m=(-\partial_q)^{m-1}K_0/(m-1)!\).
Also \(P_3=(1/350)[1-4qL^2/45+O(q^2L^4)]\).

Define the dimensionless order statistic and its clock-normalized version:
\[
\mathcal Q=\frac{105}{2}(1-10P_2),\qquad
\mathcal E_{\lambda,L}=\mathcal Q/L^2.
\tag{6}
\]
For axial tips \(\Delta u=L,\Delta v=L/2\), the central proper time is
\(L\), the central observer has \(du/d\tau=1\), and
\(\mathcal Q\to L^2R(U,U)=2qL^2+O(q^2L^4)\).
Order alone supplies \(\mathcal Q\), not an absolute length or energy unit.
Other observers see the appropriate null-energy projection; this experiment
does not reconstruct the complete tensor or wave direction.

In each fixed smooth vacuum metric, the Ricci term in the small-diamond
relation fraction vanishes. In the averaged background it does not. Therefore
the *population* readout has noncommuting resolution limits:
\[
\boxed{
\lim_{L\to0}\lim_{\lambda\to0}\mathcal E_{\lambda,L}=2q,
\qquad
\lim_{\lambda\to0}\lim_{L\to0}\mathcal E_{\lambda,L}=0.}
\tag{7}
\]
The limits do not say an arbitrarily small finite sprinkling is accurate.
At fixed density, the number of events shrinks with the window; sampling noise
and an actual clock/density convention remain separate constraints.

There is also a finite-data obstruction to identifying the *source* of the
energy. For any fixed N on a compact positive patch, the complete finite-poset
law of this vacuum family converges in total variation to that of the null-matter
background (2). To see this, uniform positive-definite convergence of gamma
gives uniform convergence of M, the causal actions, and volume densities away
from coincident u values. Axial diamonds lie in one bounded coordinate domain.
Diamond boundaries and pairwise null boundaries have product-volume measure
zero. Dominated convergence therefore applies to the integral defining the
probability of each labelled N-event order. There are finitely many such
orders; summing their probability differences proves total-variation convergence.
The N=0 case is trivial. Passing to unlabelled orders cannot increase the distance.
This statement is for fixed N; it supplies no rate when N grows with frequency.
It makes the wave-versus-matter ambiguity stronger than merely matching P2.

## 3. Predictive finite-resolution law (weak waves, phase averaged)

This is the higher-risk extension beyond reproducing the known limit. Expand
the exact kernel to quadratic order in the wave amplitude. Writing
\(a=1+A+O(\beta^4)\), with \(A''=-\beta'^2\), affine integration
constants cancel from \(K\). For \(\ell=t-s\),
\[
\delta K(s,t)=2\ell\operatorname{Var}_{[s,t]}(\beta)
 -\int_s^t(u-s)(t-u)\beta'(u)^2du.
\tag{8}
\]
This also follows by expanding each diagonal entry of \(M\) in (3).
If \(\beta=\sum_n b_nP_n(2(u-s)/\ell-1)\), Legendre orthogonality gives
\(\delta K=-\ell\sum_{n\ge2}[n(n+1)-2]b_n^2/(2n+1)\).
Constant and linear strain modes cancel; the remaining correction is negative.

For \(\beta=\epsilon\sin(\omega u+\phi)\), phase averaging and
\(q=\epsilon^2\omega^2/2\) give
\[
D(\ell):=\langle\delta K\rangle_\phi
=q\left[-\frac{\ell^3}{6}+\frac{2\ell}{\omega^2}
-\frac{4(1-\cos(\omega\ell))}{\omega^4\ell}\right].
\tag{9}
\]
Varying (4) around the flat kernel (not fitting any measured coefficient),
\[
\delta\log P_2=\frac{D(L)}L
+\frac{60}{L^5}\int_0^L(L-s)^3D(s)ds
-\frac{24}{L^3}\int_0^L(L-s)D(s)ds.
\tag{10}
\]
Let \(r=\omega L\). Then
\(\langle\mathcal E\rangle=2qF(r)+O(q^2L^2)\), where
\[
F(r)=1+\frac{105}{r^4}A(r),
\quad
A=1-\cos r+\int_0^1
\frac{36-156x+180x^2-60x^3}{x}(1-\cos(rx))dx.
\tag{11}
\]
For numerical stability at small \(r\), the runner evaluates (9)--(10).
In closed form,
\[
A=36[\gamma_E+\log r-\operatorname{Ci}(r)]-85-\cos r
+36\frac{\sin r}{r}+\frac{180-60\cos r}{r^2}
-120\frac{\sin r}{r^3}.
\]
Thus
\[
F(r)=\frac{23}{288}r^2+O(r^4),\qquad
F(r)=1+O(\log r/r^4)\quad(r\to\infty).
\tag{12}
\]
The small-window coefficient agrees with the repository's existing constant
vacuum calculation after averaging the local squared curvature. It is not an
independent new discovery. The formula also predicts an **overshoot above 1**
at intermediate resolution (test at \(r=2\pi,8\)): a monotone interpolation
between vacuum and effective energy is not promised.

For several Fourier components, phase averaging removes cross terms at this
order. The predicted readout is \(2\sum_n q_nF(nr)\). In particular,
the mixed profile below predicts \(0.8F(r)+0.2F(2r)\), with no new fit.
Equal energy therefore implies the same unresolved limit, **not** the same
finite-resolution curve. A claim of waveform-independent crossover would be
too strong. Circular polarization is included in the exact limit check, not
in the scalar linear-polarization transfer formula.

## 4. Frozen first challenge (diagnostic, not a confidence-interval verdict)

Code: `experiments/positive_control/order_wave_energy.py`. Commit this design
and source before execution, run clean, and keep all observations. No frozen
prior campaign is changed. SciPy/NumPy versions, full commit, source/design
hashes, raw rows, and numerical refinements are recorded. No GPU.

Use \(u\in[-1/2,1/2]\), \(a(-1/2)=\cos(\sqrt q/2)\),
\(a'(-1/2)=\sqrt q\sin(\sqrt q/2)\). The three families are sine,
mixed \(f=(\sin z+0.25\sin2z)/\sqrt{1.25}\), and circular. Both linear
families have \(\langle f'^2\rangle=1/2\); set \(b=\sqrt{2q}\).

1. Analytic engineering controls: flat \(P_2=1/10\); numerical convolution
   against (5); circular \(Q=q\) by independent differentiation; (12); a sine
   kernel at q=0.7, lambda=0.13, phase=0.37 against an independent trace-free
   Brinkmann Jacobi ODE; integral/closed-form agreement for (11).
2. Exact limit challenge: \(q=0.25,1,4\), all three polarizations/profiles,
   \(\lambda=1/8,1/16,1/32,1/64,1/128\), fixed phase 0.37. Compare with the
   entire nonlinear background (5), not the small-q approximation (6).
3. Transfer challenge: both linear profiles, \(q=0.0025,0.00125\),
   \(r=0.25,0.5,1,2,4,2\pi,8,16,64\), 12 uniform phases. Predict the
   phase mean from (11) without fitting. Retain per-phase variation. Halving q
   should reduce the leading truncation error, subject to numerical precision.
4. Fixed independent numerical refinements: mixed q=4,r=128 at phase 0.37;
   sine q=0.00125,r=0.25 and mixed q=0.00125,r=8 with 24 phases. Tighten
   DOP853 tolerances by 100, halve its maximum step, and use 384 quadrature
   nodes. Main quadratures are 96/192, with 384 automatically if their
   difference exceeds 2e-11. Report differences, not certified error bounds.
5. Evaluate analytic F on 160 log-spaced r in [0.01,128] to expose an overshoot
   or nonmonotonicity. This is a prediction plot, not 160 independent data.

The null directions are substantive: a wrong sign/coefficient or persisting
nonconvergence attacks the proposed mechanism; a transfer-law mismatch attacks
its predictive content. A difference between equal-energy finite-r waveforms
refutes a single universal crossover curve but was already predicted above.
Numerical controls failing makes the corresponding test unadjudicated. These
are deterministic calculations, so no binomial p-values or statistical 95%
accuracy claim will be attached to their quadrature discrepancies.

## 5. What could still defeat the scientific claim?

- Finite-order laws can converge to the same effective geometry for vacuum
  waves and actual null matter. The limit cannot identify the *origin* of the
  energy; attributing it to gravitational waves needs the vacuum-model premise.
- A fixed-order finite-N sample retains stochastic uncertainty. Convergence of
  population chain probabilities is not a recovered energy from one poset.
- A scalar observer projection and a positive plane-wave patch do not establish
  general polarization, localized packets, arbitrary curvature, singularities,
  or complicated global causal structure.
- Existing Ricci-from-chain work plus known backreaction might already entail
  much of (7). A literature non-retrieval is not proof of novelty. The candidate
  additional content is the finite-frequency transfer law, its spectral test,
  and eventually a genuinely intrinsic finite-data resolution limit.
- Uniform random thinning changes the sample density in the same geometry;
  it does not by itself create the effective Ricci tensor. Window scale and
  frequency must be distinguished from the event count.

Primary background: [Huneau--Luk, high-frequency Einstein solutions (2024)](https://arxiv.org/abs/2404.07659),
[Harte, nonlinear plane-wave optics (2015)](https://arxiv.org/abs/1502.03658),
[Roy--Sinha--Surya, small-diamond chain curvature (2013)](https://arxiv.org/abs/1212.0631).
These support the surrounding theory; they do not establish priority of the
specific transfer formula. Its derivation and numerical implementation need
independent scrutiny before any discovery or submission claim.
