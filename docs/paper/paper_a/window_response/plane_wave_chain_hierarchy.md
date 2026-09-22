# Exact chain hierarchy in plane-wave diamonds

2026-09-22. Post-study derivation draft, written after the local and global
volume investigations at `3ecca2b` and `0159730`. The frozen campaign is unchanged.
Proofs and numerical checks here were produced by the same AI; this is not an
independent expert review or a claim of priority.

## Research decision

The earlier volume checks alone do not establish a strong submission case.
The stronger candidate is an exact, nonperturbative reduction of **every fixed
chain abundance** for an arbitrary matrix-valued plane-wave profile to powers
of one scalar Volterra kernel. It gives a caustic scaling law and distinguishes
two otherwise confounded population limits: strong focusing in four dimensions
and increasing dimension in flat spacetime.

This is a general exact formula within the plane-wave class, not a formula for
all oscillating metrics. It requires solving a transverse linear ODE and scalar
integral recurrences; arbitrary profiles need not have elementary closed forms.
Its computational content is the elimination of all longitudinal and transverse
event coordinates, leaving recurrences whose quadrature dimension does not
increase with chain length. It is not merely an unevaluated spacetime-volume
definition. Publication readiness remains conditional on independent proof review
and a closer priority comparison, especially with curved-spacetime propagator
and chain-abundance literature.

## 1. Setting and scope

Let `d=n+2`, with `n>=1`, and consider the full, unquotiented plane wave

\[
 ds^2=-2\,du\,dw+x^T H(u)x\,du^2+dx^Tdx,\qquad H(u)=H(u)^T.
 \tag{1}
\]

The profile is smooth on a neighborhood of `[a,b]`. The transverse coordinates
range over all of `R^n`; there are no timelike walls or identifications. Vacuum
means `tr H=0`. Neither constant polarization nor commutativity of `H(u)` at
different times is assumed. General symmetric profiles are also covered.

Assume the Dirichlet quadratic form

\[
 Q[y]=\int_a^b(|y'|^2+y^THy)\,du
 \tag{2}
\]

is strictly positive for nonzero `y` vanishing at both endpoints. This is the
interval before the first conjugate point from `a`. It is stronger than simply
checking a nonzero determinant at the final endpoint: after multiple conjugate
points that determinant can again be nonzero. Positivity on the full interval
implies positivity on every subinterval by zero extension.

For initially axial tips take `p=(a,0,0)`, `q=(b,h,0)`, `h>0`, and orient time
by increasing `u`. The measure is `du dw d^n x`. Null generators with constant
`u` have measure zero and do not change the chain integrals.

Define the Jacobi matrix and the positive kernel

\[
 \partial_t^2 B(t,s)=H(t)B(t,s),\quad B(s,s)=0,\quad
 \partial_tB(s,s)=I,\qquad K(s,t)=\sqrt{\det B(t,s)}\quad(s<t).
 \tag{3}
\]

There are no determinant zeros before conjugacy, so the square root has the
positive branch continued from `B(t,s)~(t-s)I`. Define

\[
 C_1(s,t)=K(s,t),\qquad
 C_{r+1}(s,t)=\int_s^t C_r(s,u)K(u,t)\,du.
 \tag{4}
\]

Here `C_r` denotes a continuum convolution, **not** a random chain count.
We use `N_k` for the latter. All kernels are nonnegative; associativity follows
from Tonelli's theorem. This is a two-argument Volterra kernel, not an ordinary
translation-invariant convolution except for constant profiles.

## 2. Exact theorem for all chain lengths

Let `I_k` be the physical-volume integral over ordered interior events
`p<X_1<...<X_k<q`. The two tips are not counted among the `k` events. Then

\[
 \boxed{I_k(a,b;h)=
 \frac{(2\pi)^{kn/2}h^{kd/2}}{\Gamma(kd/2+1)}
 \frac{C_{k+1}(a,b)}{K(a,b)}}\qquad(k\ge1).
 \tag{5}
\]

In particular `V=I_1`. For `k` independent uniform points in the diamond, let
`P_k` be the probability that they form a total order, allowing all `k!` temporal
labelings. For fixed `N` points, or for a Poisson sprinkling of density `rho`,

\[
 P_k=\frac{k!I_k}{V^k},\qquad
 \mathbb E_N N_k=\binom Nk P_k,\qquad
 \mathbb E_\rho N_k=\rho^k I_k.
 \tag{6}
\]

These are chains of causal relations, not chains of links, maximal chains, or
inclusive intervals. The complete joint distribution of arbitrary order motifs
is not determined here.

### Proof: action, slack simplex, and determinant composition

For a causal path parametrized by `u`, (1) gives

\[
 w'\ge\tfrac12(|x'|^2+x^THx).
 \tag{7}
\]

Strict positivity of (2) makes the stationary transverse action its unique
minimum. Write the symplectic transfer matrix on one segment as
`F=[[A,B],[C,D]]`. Its classical action is

\[
 S(x_0,x_1)=\tfrac12 x_0^TB^{-1}Ax_0-x_0^TB^{-1}x_1
             +\tfrac12 x_1^TDB^{-1}x_1.
 \tag{8}
\]

Fix ordered `u` coordinates `a=u_0<u_1<...<u_k<u_{k+1}=b`. Sum (8) over
segments with `x_0=x_{k+1}=0`. This gives `S_total=X^TMX/2`, where `X` has
`nk` coordinates and `M` is a positive block tridiagonal Hessian. The `k`
intermediate `w` integrals are a slack simplex of volume

\[
 \frac{[h-X^TMX/2]_+^k}{k!}.
 \tag{9}
\]

Changing variables to the unit ball and performing its beta integral gives

\[
 \int_{\mathbb R^{nk}}\frac{[h-X^TMX/2]_+^k}{k!}\,dX
 =\frac{(2\pi)^{nk/2}h^{k+nk/2}}
        {\Gamma(k+nk/2+1)\sqrt{\det M}}.
 \tag{10}
\]

The crucial identity is

\[
 \det M=\frac{\det B(b,a)}
 {\prod_{i=0}^{k}\det B(u_{i+1},u_i)}.
 \tag{11}
\]

For two segments the intermediate Hessian is
`M=D_1 B_1^{-1}+B_2^{-1} A_2`. Direct multiplication, in the displayed order,
gives

\[
 B_2 M B_1=B_2D_1+A_2B_1=B_{21}.
 \tag{12}
\]

Taking determinants proves (11) for one intermediate event. Eliminating further
intermediate coordinates by Schur complements composes neighboring actions;
the same identity then proves (11) inductively. No matrices were commuted.
Substituting (11) into (10) and integrating the ordered `u` simplex yields (5).

As a separate normalization check, Laplace transforming (9) in `h` gives
`exp(-lambda X^TMX/2)/lambda^(k+1)`. Its Gaussian integral is
`(2 pi)^(kn/2) lambda^(-kd/2-1)/sqrt(det M)`. Inverse transformation yields
the same gamma factor and power of `h` in (5).

For nonaxial tips replace `h` by

\[
 h_{\rm eff}=w_q-w_p-S_{a,b}(x_p,x_q)>0.
 \tag{13}
\]

Translate the intermediate transverse positions by their classical path.
Stationarity removes cross terms, the Hessian is unchanged, and the Jacobian
is one. Equations (5)-(6) follow with `h_eff`. In particular normalized chain
probabilities depend on the profile and `[a,b]`, but not on this positive slack.

## 3. Four-dimensional consequences and flat checks

For `d=4`, the theorem becomes

\[
 I_k=\frac{(2\pi)^k h^{2k}}{(2k)!}\frac{C_{k+1}}K,
 \quad V=\pi h^2\frac{C_2}{K},\quad
 \boxed{P_2=\frac{KC_3}{3C_2^2},\quad
 P_3=\frac{K^2C_4}{15C_2^3}}.
 \tag{14}
\]

Consequently the population ratio

\[
 \boxed{\mathcal R=\frac{P_3}{P_2^2}
       =\frac35\frac{C_2C_4}{C_3^2}}
 \tag{15}
\]

does not explicitly contain the vanishing full-interval Jacobi determinant.

In flat spacetime, putting `nu=d/2`,

\[
 K(s,t)=(t-s)^{\nu-1},\quad
 C_r(a,b)=\frac{\Gamma(\nu)^r}{\Gamma(r\nu)}(b-a)^{r\nu-1}.
 \tag{16}
\]

Equations (5)-(6) reduce to the known flat result

\[
 P_k^{\rm flat}(d)=(k-1)!
 \left[\frac{\Gamma(d+1)}2\right]^{k-1}
 \frac{\Gamma(d/2)\Gamma(d)}
      {\Gamma(kd/2)\Gamma((k+1)d/2)}.
 \tag{17}
\]

In four dimensions `P_k=k!12^k/[(2k)!(2k+1)!]`, so
`P_2=1/10`, `P_3=1/350`, `P_4=1/29400`, and `R=2/7`.
This agrees with Roy-Sinha-Surya Eq. (9), after their `chi_k` is multiplied
by `k!`. Recovering their formula is a normalization control, not a new result.

## 4. Conditional degree moments and an exact finite-N variance

For a uniformly chosen point `X` at coordinate `u`, let `V_p(X)` and `V_f(X)`
be its past and future interval volumes within the diamond. Set

\[
 \alpha(u)=\frac{C_2(a,u)}{K(a,u)},\quad
 \beta(u)=\frac{C_2(u,b)}{K(u,b)},\quad
 \gamma=\frac{C_2(a,b)}{K(a,b)},\quad
 W(u)=\frac{K(a,u)K(u,b)}{K(a,b)}.
 \tag{18}
\]

For integers `r,s>=0`, an exact mixed moment is

\[
 \boxed{\mathbb E\left[(V_p/V)^r(V_f/V)^s\right]
 =\frac{2(2r)!(2s)!}{(2r+2s+2)!}
 \frac{\int_a^b W(u)\alpha(u)^r\beta(u)^s\,du}
      {\gamma^{r+s+1}}.}
 \tag{19}
\]

To prove it, at fixed `u` diagonalize the boundary Hessian and put
`z=x^TMx/2`. The transverse radial measure is `2 pi W(u) dz`.
The longitudinal slacks `eta,zeta>=0` obey `eta+zeta=h-z`, and the subinterval
volumes from (14) are `pi alpha eta^2`, `pi beta zeta^2`. Their integral is the
Dirichlet simplex beta integral
`h^(2r+2s+2)(2r)!(2s)!/(2r+2s+2)!`. Normalization by `V^(r+s+1)` yields (19).

Let `q(X)=(V_p+V_f)/V` and `zeta_1=Var(q(X))`. Then

\[
 \mathbb E q=P_2,\quad
 \mathbb E q^2=
 \frac{\int_a^b W[(\alpha^2+\beta^2)/15+\alpha\beta/45]du}{\gamma^3},
 \tag{20}
\]

and for the observed unordered relation fraction `f_N`,

\[
 \operatorname{Var}(f_N)=\frac{2P_2(1-P_2)}{N(N-1)}
 +\frac{4(N-2)\zeta_1}{N(N-1)}.
 \tag{21}
\]

The variance identity is standard U-statistic combinatorics. The new candidate
calculation is its exact evaluation for arbitrary plane-wave profiles. This
does not imply an exact Gaussian sampling law or certify a confidence interval.

## 5. A closed weak-field coefficient for every k

For the constant vacuum profile `H=diag(omega^2,-omega^2)`, write
`L=b-a`, `z=omega L`. Its kernel is

\[
 K(s,t)=\frac{\sqrt{\sin(\omega(t-s))\sinh(\omega(t-s))}}{\omega}
 =(t-s)[1-\omega^4(t-s)^4/180+O(\omega^8(t-s)^8)].
 \tag{22}
\]

The flat simplex of `r` gaps has normalized Dirichlet parameters `(2,...,2)`.
Its moment `E sum gap_i^4/L^4 = r(2)_4/(2r)_4` gives

\[
 \frac{C_r}{C_r^{\rm flat}}=
 1-\frac{z^4}{3(2r+1)(2r+2)(2r+3)}+O(z^8).
 \tag{23}
\]

Thus, for each fixed `k`,

\[
 \boxed{\frac{P_k}{P_k^{\rm flat}(4)}=1+b_kz^4+O(z^8),\quad
 b_k=\frac1{180}-\frac1{3(2k+3)(2k+4)(2k+5)}-\frac{k}{252}.}
 \tag{24}
\]

This recovers `V/V_flat=1+z^4/252+O(z^8)` and the previous
`P_2=1/10-23z^4/75600+O(z^8)`. It also gives
`b_3=-139/20790` and `R=(2/7)[1-5z^4/8316+O(z^8)]`.
The general coefficient is derived from the exact hierarchy, not a fit to the
original campaign. The remainder here is for fixed `k`, not uniform as `k->infty`.

## 6. Caustic scaling theorem

Keep `a` fixed and approach the first conjugate endpoint `b_c` from below.
Assume a regular zero

\[
 \det B(b,a)\sim c\delta^m,\quad c>0,\quad\delta=b_c-b\downarrow0,
 \tag{25}
\]

and a positive limiting effective slack. Here `m` is the determinant zero's
order; for the usual regular crossing it equals the conjugate multiplicity.
The same argument applies to a smooth family of profiles on a fixed interval
with this degeneration. Every proper subinterval remains strictly positive.
The kernels are continuous and bounded on the compact limiting triangle, so
dominated convergence gives

\[
 C_r(a,b)\longrightarrow c_r\in(0,\infty)\quad(r\ge2).
 \tag{26}
\]

Positivity follows by integrating over any interior simplex of proper gaps.
Combining (5), (6), and (25)-(26), for every fixed `k`,

\[
 \boxed{I_k\asymp\delta^{-m/2},\qquad
 P_k\asymp\delta^{m(k-1)/2}.}
 \tag{27}
\]

In particular the volume diverges, while normalized chain probabilities vanish
for `k>=2`. At and beyond the first caustic the axial diamond in these full
plane-wave models has infinite volume; a uniform fixed-N volume sample is then
undefined. Replacing determinants by their absolute values after conjugacy
does not repair the probability model. Truncating the transverse space would
define a new, window-dependent problem.

### Flat dimension becomes arbitrarily misleading

The flat relation fraction is

\[
 F(d)=\frac{\Gamma(d+1)\Gamma(d/2)}{2\Gamma(3d/2)}
 \sim\sqrt{3\pi d/2}\left(\frac2{3\sqrt3}\right)^d.
 \tag{28}
\]

Putting `a_0=log(3 sqrt(3)/2)`, the population flat calibration
`F(d_pop)=P_2` therefore satisfies

\[
 \boxed{d_{\rm pop}=
 \frac{m}{2a_0}\log(1/\delta)
 +\frac1{2a_0}\log\log(1/\delta)+O(1).}
 \tag{29}
\]

The physical spacetime dimension has not changed. Curved-space bias of the
Myrheim-Meyer estimator is already known; (29) is a specific nonperturbative
caustic asymptotic derived here, not the first observation of that bias.

### Three-chain ratio separates the two population limits

For a four-dimensional wave, (15) gives

\[
 \mathcal R\longrightarrow\frac35\frac{c_2c_4}{c_3^2}>0.
 \tag{30}
\]

For a flat diamond of dimension `d`, (17) instead gives

\[
 \mathcal R_{\rm flat}(d)=
 \frac{2\Gamma(d)\Gamma(3d/2)}{\Gamma(d/2)\Gamma(2d)}
 \sim\sqrt{8/3}\left(\frac34\right)^{3d/2}\longrightarrow0.
 \tag{31}
\]

Thus any sequence of flat dimensions matching the vanishing `P_2` near the
wave caustic fails to match `P_3/P_2^2`. This is a separation between these two
population limits, not a universal identification of curvature, topology, or
dimension from two numbers.

For the constant vacuum profile with `[a,b]=[0,1]`, `delta=pi-z`, deterministic
quadrature gives `R_limit~0.262779744688`. Representative values are:

| delta | P2 | flat-calibrated d | wave R | flat R at the same P2 |
|---:|---:|---:|---:|---:|
| 0.1 | 0.04352855 | 4.97565 | 0.26661147 | 0.18815262 |
| 0.001 | 0.00464733 | 7.52142 | 0.26282338 | 0.06301134 |
| 0.000001 | 0.000147066 | 11.34430 | 0.26277979 | 0.01214279 |

These decimals are reproducible numerical evaluations, not certified enclosures.
The limiting positivity and contrast with (31) follow analytically without
relying on their decimal values.

An independent closed-form control is the isotropic **nonvacuum** profile
`H=-omega^2 I_2`. Here `K(s,t)=sin(omega(t-s))/omega`, and Laplace transformation
gives, with `j_l` the spherical Bessel function,

\[
 C_r(a,b)=\frac{L^r j_{r-1}(\omega L)}{(2\omega)^{r-1}(r-1)!},
 \quad L=b-a.
 \tag{32}
\]

For integer `r` these reduce to elementary sine/cosine expressions. At the
first caustic the determinant has a double zero, so `P_k` vanishes as
`(pi-omega L)^(k-1)`, rather than the half-power sequence of the constant vacuum
example. Direct evaluation of (32) gives the exact limit
`R -> 2(15-pi^2)/45`. This additional case checks the multiplicity dependence;
it is not evidence that the isotropic profile satisfies the vacuum equations.

### The finite-sample cost cannot be hidden

At fixed `N`, a union bound gives
`Pr(at least one k-chain) <= binom(N,k) P_k`. Consequently a fixed-size sample
converges to an antichain in the caustic limit. The observed ratio is undefined
if there are no relations, and a ratio of counts is not an unbiased version of
the population ratio.

For seeing any three-chain with nonvanishing probability it is necessary that
`N^3 delta^m` not tend to zero. This is only a necessary condition, not a power
calculation or a sufficient condition for accurate estimation. At `delta=1e-6`
and the earlier campaign's `N=1200`, the expected three-chain count is only
about `1.63` (about `0.075` in the matching flat model). A new claim of practical
classification would require a frozen estimator, valid uncertainty, and power
sizing. No such campaign has been run here.

## 7. Numerical and provenance audit

Run from the repository root in the existing NumPy/SciPy environment:

```powershell
.venv/Scripts/python.exe docs/paper/paper_a/window_response/check_plane_wave_chains.py
```

The script writes `plane_wave_chain_audit.json`, binding its LF-normalized source
hash and software versions. It reads no sampled poset artifacts. Its checks are:

- Exact flat normalizations for `k=1,...,6`; general-dimensional comparison with
  (17) at `d=3,4,5,6,8`; rational beta-moment reduction of (24).
- Independently integrated segment transfer matrices versus full-interval
  determinants for `k=1,...,5`, including a noncommuting rotating polarization.
  The largest relative determinant discrepancy was about `1.6e-13`.
- An oscillating diagonal profile and a rotating profile with changing amplitude.
  Both obey `sup ||H||<=5<pi^2`, proving positivity on the unit interval by
  Poincare's inequality. Symplectic residuals and time reversal are checked.
- Grid refinements through 2048 panels; the last relative change across
  `P_1,...,P_6` was below `9.4e-9` in both nonconstant profiles. The initial
  512-to-1024 change for `P_6` exceeded the chosen `1e-7` numerical check, so the
  grid was refined rather than loosening that check. This is a deterministic
  development diagnostic, not a preregistered scientific accuracy gate.
- Cross-check of volume, `P_2`, and conditional-degree variance against the old
  independent Gauss quadrature at amplitudes `0,1/64,1,2`; all absolute
  differences were below `2.7e-12`.
- Near-caustic `k=1,...,4` scaling at five distances; the last refinement change
  was below `4.0e-8`. Separate nested adaptive quadrature verifies the limiting
  ratio. Neither refinement differences nor adaptive error estimates are
  asserted to be rigorous mathematical error bounds.
- The isotropic nonvacuum profile checks `C_1,...,C_7` against the independent
  closed Bessel expressions (32) and checks the double-zero caustic limit.

Quadrature thresholds are implementation checks. They do not replace the proof,
certify all profiles numerically, or imply an independent scientific replication.

## 8. What is and is not a publication contribution

Primary-source comparison, bounded to the papers and searches inspected on
2026-09-22:

| Primary source | Established result to credit | Candidate distinction here |
|---|---|---|
| [Roy, Sinha, Surya, arXiv:1212.0631](https://arxiv.org/html/1212.0631), Eqs. (7)-(9), (16), sections 3-6 | Nested chain integrals, flat chain means, small-diamond Ricci corrections, curvature-sensitive dimension estimators | Exact all-k matrix-profile reduction without a curvature expansion; caustic scaling |
| [Kambor, Nomaan, arXiv:2007.03835](https://arxiv.org/html/2007.03835), section 2 | Chain-based geometric inference, tested in dS2 and FLRW3 | Exact plane-wave hierarchy and its near-caustic population limit |
| [Harte, Drivas, arXiv:1202.0540](https://arxiv.org/html/1202.0540), section III | Jacobi matrices, conjugate structure, plane-wave causal geometry | Unweighted causal-chain integrals and the ratio in (30)-(31) |
| [Nomaan, Dowker, Surya, arXiv:1701.07212](https://arxiv.org/html/1701.07212), Eqs. (15), (18), (25)-(34) | Causal-matrix powers and chain expectations as spacetime convolutions; four-dimensional link-based propagators | Reduction of unweighted chains to a one-dimensional Jacobi-determinant kernel, not just a formal spacetime convolution |

The action/Jacobi machinery, Gaussian composition, nested-chain method, and
U-statistic variance identity are established tools. Their mere reuse is not
the novelty claim. In particular, the last source's Eq. (18) already states
the general convolution identity for chain expectations. Equation (5) evaluates
its transverse and longitudinal variables in a plane wave; writing a chain
expectation as a convolution on spacetime by itself would add nothing. Its
four-dimensional link propagator also contains an interval-emptiness factor,
`exp(-rho V)`, absent from our unweighted chain integral.

Targeted searches for plane-wave chain abundances, causal-set
Volterra formulas, and plane-wave propagator connections did not retrieve a
direct match, but also returned irrelevant Fourier-plane-wave material. This
limited retrieval is not proof that no matching result exists. In particular,
weighted Green-function compositions must be compared carefully with the
unweighted causal-order chains here before making a firstness claim.

The research now has a specific theoretical contribution candidate: (5), (24),
and the linked consequences (27)-(31). It is stronger than a list of corrected
volume coefficients, because it calculates an infinite family of observables
and supplies a testable failure mechanism for flat dimensional inference.
It still does not justify promising acceptance or submitting immediately.

The decisive next checks are (i) independent scrutiny of the noncommuting
determinant proof and the caustic hypotheses, and (ii) an equation-level priority
comparison with propagator/chain formulas. If those survive, a theorem-focused
paper can stand on the exact hierarchy and asymptotic contrast; a practical
classifier claim additionally needs new sampling work. Earlier weak-amplitude
holdout results cannot be repurposed as validation of that stronger claim.
