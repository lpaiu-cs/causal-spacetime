# Why the coordinate-coupling lower bound stops at N lambda

2026-09-24. Analytic probe and frozen deterministic check design. This is a
**barrier for a proof method**, not a solution of the order-only problem.
The argument below gives a sharp law for a strictly stronger experiment
which reveals coordinates. The main research goal remains active.
All new derivations are by the same author and await independent audit.

## 1. Two experiments that must not be conflated

Keep fixed `0<q<pi^2`, axial window `[-1/2,1/2]`, slack h=1, and the exact
circular vacuum family in [the original probe](full_order_resolution_probe.md).
Let `a(u)=cos(sqrt(q)u)`, `f=a^-2`, `m(s,t)=integral_s^t f`, and
`M=m(-1/2,1/2)`. Use the transverse chart `x=sqrt(c)y`, with
`c=sqrt(1+4q lambda^2)`, used in that coupling proof. In this chart

\[
\widetilde M_\lambda(s,t)=m(s,t)I+E_\lambda(s,t),\qquad
E_\lambda=-\frac{2\lambda\sqrt q}{c}\int_s^t f(u)T(u/\lambda+\phi)du,
\]

where `T(theta)=[[cos(theta),sin(theta)],[sin(theta),-cos(theta)]]`.
The trace-free E has eigenvalues `+e,-e`, with

\[
e_\lambda(s,t)=\frac{2\lambda\sqrt q}{c}
\left|\int_s^t f(u)e^{iu/\lambda}du\right|,
\qquad e/m\le 2\lambda\sqrt q/c<1.
\tag{1}
\]

An event has normalized density `a(u)^2/Vbar_lambda` on its wave diamond in
these coordinates; the constant c in physical volume cancels when N is fixed.
Denote this one-point law by mu_lambda, and its effective-background version
by mu_0. The earlier strip bound gives `TV(mu_lambda,mu_0)<=C lambda^2`.

The actual experiment supplies only the sampled finite order. Define an
auxiliary experiment Q_(N,lambda) which additionally supplies all the real
coordinates `(u,v,y)` and their association with the order's events. The
chart is supplied by a model oracle, including the above rescaling. This is
extra information; it cannot be obtained from the original data by assertion.

With q known, let Z count unordered event pairs whose supplied relation
differs from the relation predicted by the effective background at those
coordinates. Its null action is `|y_j-y_i|^2/(2m(u_i,u_j))`, with u ordered.
Rejecting when Z>0 has exactly zero false positives in exact arithmetic and
requires no wavelength or phase input **once the coordinate oracle is given**.
It is not an order-only test, or a claimed floating-point certification.

## 2. Sharp auxiliary law

There is an explicit `kappa(q)>0`, given below, such that the disagreement
probability of one independent unordered pair is

\[
p_\lambda=\kappa(q)\lambda^2+o(\lambda^2).
\tag{2}
\]

For every fixed first event, the conditional disagreement probability is
bounded by `C lambda^2`, uniformly on the bounded sampling domains. Indeed,
at longitudinal separation ell a crossing needs `|Delta y|^2<=C ell`, and
its v-strip width is bounded by
`C |Delta y|^2 lambda min(ell,lambda)/ell^2`. Integrating over the two
transverse displacement coordinates gives `C lambda min(ell,lambda)`;
integrating ell gives `C lambda^2`. The other v/position factors and the
normalized densities are bounded uniformly for this fixed positive patch.

Indicators for disjoint event pairs are independent. For pairs sharing one
event, conditioning on that event bounds their joint expectation by
`C lambda^4`. Consequently

\[
\mathbb E Z={N\choose2}p_\lambda,\qquad
\operatorname{Var} Z\le C(N^2\lambda^2+N^3\lambda^4).
\tag{3}
\]

Thus `N lambda -> infinity` makes `Pr_lambda(Z=0)->0` by Chebyshev's
inequality, while the existing coupling proves indistinguishability when
`N lambda -> 0` even for this augmented experiment.

There is also a critical-window statement. If `N lambda -> tau` for fixed
`0<tau<infinity`, the dependency graph has O(N) neighbours per pair and
`b1+b2=O(N^3 lambda^4)=o(1)`, `b3=0`. The Poisson approximation of
[Arratia, Goldstein and Gordon, Theorem 1](https://dornsife.usc.edu/larry-goldstein/wp-content/uploads/sites/221/2023/06/AGG-1.pdf)
therefore yields `Z -> Poisson(kappa(q) tau^2/2)` in total variation.
Their norm has a factor-two convention; only convergence is used here.

To relate this to the entire augmented law, write A for the coordinate
configurations where the two deterministic order maps agree. The overlap
of their joint densities is `integral_A min(dmu_lambda^N,dmu_0^N)`. Hence

\[
0\le \operatorname{TV}(Q_{N,\lambda},Q_{N,0})-\Pr_\lambda(Z>0)
\le\operatorname{TV}(\mu_\lambda^N,\mu_0^N)\le CN\lambda^2.
\]

Combining these statements gives, with TV in [0,1],

\[
\boxed{\operatorname{TV}(Q_{N,\lambda},Q_{N,0})
\longrightarrow 1-\exp[-\kappa(q)\tau^2/2]
\quad\text{if }N\lambda\longrightarrow\tau\in[0,\infty].}
\tag{4}
\]

The value at tau=infinity is 1 by (3), without requiring `N lambda^2 -> 0`.
No such equality has been established for P, the order-only law.

## 3. Pair boundary integral and an explicit coefficient

For fixed `s<t` set `A=m(-1/2,s)`, `b=m(s,t)`, `B=m(t,1/2)`, so
`A+b+B=M`. Integrating the two v-slacks and four transverse coordinates in
the background chain `tip < X < Y < tip` gives the unweighted slice volume

\[
C(A,b,B)=\frac{\pi^2h^4}{6}\frac{AbB}{M}.
\tag{5}
\]

Keep the outer segments A,B fixed when differentiating with respect to b.
Then `C_b'=pi^2 h^4 A B(A+B)/(6 M^2)`. This is a boundary integral: the
derivative of the middle indicator inserts `|Delta y|^2/(2b^2)` on its null
boundary. Additional membership constraints in the fixed background diamond
are inactive there except on a set of boundary measure zero, by transitivity.

To first order in E, the action shift is
`-Delta y^T E Delta y/(2b^2)`. The average absolute quadratic form of a
trace-free 2x2 matrix with eigenvalues +/-e on a circle is `2e|Delta y|^2/pi`.
Rotational invariance therefore makes the disagreement slice volume

\[
\frac{2e}{\pi}C_b'+o(e)
=\frac{\pi h^4e}{3M^2}AB(A+B)+o(e).
\tag{6}
\]

For s,t separated from each other and the tips, integration by parts in
(1) gives, uniformly on compact sets,

\[
e/\lambda^2=2\sqrt q\,|f(t)e^{it/\lambda}-f(s)e^{is/\lambda}|+O(\lambda).
\]

Rapid-phase averaging replaces the absolute value by

\[
\mathcal B(F,G)=\frac1{2\pi}\int_0^{2\pi}
\sqrt{F^2+G^2-2FG\cos\psi}\,d\psi
=\frac{2(F+G)}\pi\operatorname{E}\!\left(\frac{4FG}{(F+G)^2}\right),
\tag{7}
\]

where E uses the elliptic-integral **parameter**, not its square root.
This averaging is in the separation variable; no independence of phases at
finite lambda is assumed. The uniform slice strip bound from Section 2 lets
one discard tip strips and `t-s<eta` at cost `C eta lambda^2`, first take
lambda to zero, and then eta to zero. It controls the otherwise nonuniform
short-segment remainder in (6). Replacing mu_lambda by mu_0 in a disagreement
probability changes it by `O(lambda^4)`: integrate the signed one-point
measure difference against the uniform conditional strip bound, for each
coordinate in turn. Thus this substitution does not change (2).

The reference one-point volume in y coordinates is
`V0=pi h^2 C2(q)/K(q)`, where `K=sin(sqrt(q))/sqrt(q)` and C2 is its
two-fold segment convolution. Accounting for both temporal orientations,

\[
\boxed{\kappa(q)=\frac{4\pi\sqrt q\,h^4}{3V_0^2M^2}
\int_{-1/2<s<t<1/2}
a(s)^2a(t)^2 AB(A+B)\mathcal B(f(s),f(t))\,ds\,dt.}
\tag{8}
\]

This is positive and finite for every fixed `0<q<pi^2`; h cancels.
It also gives the analytic small-q control

\[
\lim_{q\downarrow0}\kappa(q)/\sqrt q=\frac{32}{5\pi^2}.
\tag{9}
\]

Here `mathcal B(1,1)=4/pi` and
`integral_(A,B>=0,A+B<=1) AB(A+B)dA dB=1/30`.
At q=0 itself the models coincide; (9) is a limiting coefficient, not a
nonzero detection result at q=0.

## 4. An exact finite-lambda check, with its sampling law stated

To avoid checking (8) only against another first-order calculation, compute
the **exact lost-relation mass for iid reference coordinates**: the pair is
causal under H0, but not under the wave middle action. This reference-coordinate
experiment is a deterministic integration check, not actual wave sprinkling.
Its leading coefficient is kappa/2; twice this loss has the same leading
term as full disagreement, but is not its exact finite-lambda value.

Put `chi=e/b`, `T=(A+B)/M`, `theta0=acos(chi)/2`, and

\[
d(\theta)=\frac{T\chi(\chi-\cos 2\theta)}{1-\chi^2},\qquad
L(\chi,T)=\frac2\pi\int_{\theta_0}^{\pi/2}
\frac{d(\theta)}{1+d(\theta)}d\theta.
\tag{10}
\]

The exact lost slice volume is `C(A,b,B)L(chi,T)`. A direct derivation shifts
the transverse centre, leaving total null slack
`h-alpha|Delta y|^2-beta|w|^2`, with
`alpha=M/[2b(A+B)]` and `beta=(A+B)/(2AB)`. Integrating w and the two
nonnegative v-slacks turns loss into
`pi h^4/(48 beta) integral[1/alpha-1/(alpha+d_action_positive)]dtheta`.
Dividing by (5) gives (10). In particular `L=T chi/pi+O(chi^2)`.

The normalized lost probability is

\[
\ell_\lambda^{\rm ref}=\frac2{V_0^2}\int a(s)^2a(t)^2
C(A,b,B)L(e/b,(A+B)/M)\,ds\,dt,
\quad 2\ell_\lambda^{\rm ref}/\lambda^2\longrightarrow\kappa(q).
\tag{11}
\]

### Frozen check design

Implementation: `experiments/positive_control/order_wave_oracle_barrier.py`.
Before execution, freeze this file and the implementation in Git. Use q in
`{0.25,1,4}`, r=1/lambda in `{32,64,128,256,512,1024}`. Evaluate (11) with
the existing complex-primitive ODE; use 32 Gauss nodes per separation panel
whose endpoints are multiples of `2 pi/r`, 64 midpoint nodes, and 24 angular
nodes on the exact positive sector. Compare its scaled value to (8), computed
with 128x128 quadrature. At q=1,r=1024 double every quadrature order and use
the tighter ODE. Freeze a relative refinement tolerance `1e-5`. Retain all
rows, including any numerical inconsistency or finite-frequency failure to
approach the coefficient. This is not a finite-N accuracy/power campaign.

Runnable controls cover (9), the isotropic derivative of (5), the positive
angular integral against its independent arctangent antiderivative, and
the primitive against direct segment quadrature. No Poisson simulation or
order-only classifier performance is inferred from this integration check.

## 5. Consequence for the unresolved research goal

The old `N lambda -> 0` proof compares a stronger, coordinate-augmented
experiment and then forgets coordinates. Equation (4) shows that this
stronger experiment really becomes distinguishable at `N~1/lambda`.
Consequently, improving the upper bound on **that augmented TV** cannot
prove an order-only impossibility region extending beyond this scale.

A correlated coupling of orders, a bound on latent-coordinate posteriors,
or another argument after marginalizing coordinates could still improve
the actual lower bound. This theorem does not prohibit those approaches,
prove that the order-only boundary exceeds `1/lambda`, or remove the gap
between necessity exponent 1 and sufficient exponent 4. It identifies the
precise point where the earlier proof throws away the question of interest.
Do not mark the discovery-sized goal complete on the basis of this result.
