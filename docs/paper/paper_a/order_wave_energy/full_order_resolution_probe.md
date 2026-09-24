# Next high-risk question: the full-order resolution law

Latest bound (2026-09-24): the [quantized-coordinate proof](quantized_order_bound.md)
improves the impossibility region to every fixed power `N=O(lambda^(-alpha))`
with `alpha<4/3`. The earlier sections below retain their derivation history.
The [shared-reference construction](local_null_challenge.md), with
[retained checks and failures](local_null_results.md), improves the known-q<=1
sufficient condition to `N lambda^3/log N -> infinity`. The sharp boundary
and an observer attaining it remain unresolved; these proofs await independent
researcher audit.

2026-09-24. Analytic probe after the pointwise pair-fraction inverse was stopped.
This is not a preregistered experiment or a claimed discovery. The potentially
substantial target is a **sharp event-count versus wavelength law for the full
order**, with a statistic that attains it. The estimate below is only one side
of that problem, and its proof has not received independent review.

## Concrete success and failure

Let P_(N,lambda) be the law of an N-event causal order in an exact vacuum wave,
and P_(N,0) the law in the effective null-matter geometry. A complete answer
would determine the scale N*(lambda) separating total-variation convergence
to zero from nontrivial or consistent distinguishability, in a stated wave
class. It would also identify a test requiring only the order, event count,
and declared clock protocol, with unseen phases and frequencies.

This would answer whether pair-count compression discards recoverable wave
information. A few successful classifiers or another fitted continuum response
are not that conclusion. If only the elementary coupling bound below is
obtained, the discovery-sized target remains unachieved; it is not automatically
turned into a small publication.

## Exact circular vacuum family

Fix q>0 and u in [-1/2,1/2] inside a positive Rosen patch,
a(u)=cos(sqrt(q)*u)>0. Let

\[
\gamma_\lambda=a^2R(\theta)
\operatorname{diag}(e^{2\beta},e^{-2\beta})R(\theta)^T,
\quad\theta=(u/\lambda+\phi)/2,
\quad\sinh(2\beta)=2\lambda\sqrt q.
\]

The metric is exactly vacuum by the earlier Q=q identity. Put
c=cosh(2 beta)=sqrt(1+4 q lambda^2), s=sinh(2 beta). For a segment [v,w],
\[
M_\lambda=c m I-sJ,\quad
m=\int_v^w a^{-2}du,\quad
J=\int_v^w a^{-2}
\begin{pmatrix}\cos(u/\lambda+\phi)&\sin(u/\lambda+\phi)\\
\sin(u/\lambda+\phi)&-\cos(u/\lambda+\phi)\end{pmatrix}du.
\]

Rescale transverse coordinates by x=sqrt(c)*y. This preserves the order and
the axial tips. The effective action matrix in y coordinates is
\(\widetilde M=M_\lambda/c=mI-(s/c)J\), while the physical volume density
is c*a^2. The constant c cancels on conditioning on N uniform events.
This coordinate coupling is used to prove a bound; coordinates are not inputs
to any proposed observing instrument.

## A joint-limit obstruction beyond fixed N

Write ell=w-v. Boundedness and one integration by parts give, uniformly in
phase and in subintervals of the fixed positive patch,
\[
\|J(v,w)\|\le C\min(\ell,\lambda),\qquad
\|\widetilde M-mI\|\le C\lambda\min(\ell,\lambda).
\tag{1}
\]
The constant may depend on q, the fixed axial slack and the positive-patch
margin, but not on lambda, N, phase or segment endpoints. Both inverse
matrices are uniformly elliptic relative to ell for sufficiently small lambda.
The null action difference is bounded by
\[
|\delta A|\le C|\Delta y|^2
\frac{\lambda\min(\ell,\lambda)}{\ell^2}.
\tag{2}
\]

For a causal-boundary crossing inside the bounded axial diamond, the
longitudinal slack is bounded and |Delta y|^2<=C*ell. At fixed endpoint u
values and one transverse location, integrate the longitudinal threshold
strip over the two-dimensional transverse displacement:
\[
\int_{|\Delta y|^2\le C\ell}|\delta A|d^2\Delta y
\le C\lambda\min(\ell,\lambda).
\tag{3}
\]
The remaining bounded longitudinal/position factors do not change this rate.
Integrating ell over (0,1) bounds the measure of pair-relation disagreements
by C*lambda^2. Applying the same estimate with one point at either axial tip
bounds the symmetric difference of the sampling diamonds by C*lambda^2.
All densities and the nonzero reference diamond volume have uniform bounds.
The normalized one-point sampling laws in common y coordinates therefore
differ in total variation by at most C*lambda^2.

Maximally couple the N event locations, and use a union bound over all pairs
when the locations agree. Taking a measurable causal-order map and then
forgetting labels cannot increase total variation. Consequently
\[
\boxed{\operatorname{TV}(P_{N,\lambda},P_{N,0})
\le\min\{1,\ C[N+N(N-1)]\lambda^2\}.}
\tag{4}
\]

In particular, **N*lambda -> 0 implies full-order indistinguishability**, even
when N tends to infinity. In the normalized window the frequency parameter is
r=1/lambda; N=o(r) is insufficient. This strengthens the previous fixed-N
observation for this specific exact circular family. It does not show that
N proportional to r suffices, and it does not determine the sharp exponent.
The same phase-uniform upper bound holds after mixing phases. Unconditional
known-density Poisson counts require a separate count coupling and are not
silently substituted for the fixed-N protocol here.

## Why the gap is worth testing, without claiming it is solved

In the separately ordered weak-amplitude, phase-averaged limit, the pair
response F(r)-1 falls as log(r)/r^4. A root-N raw pair-fraction comparison then
suggests an r^8/log(r)^2 event price at fixed weak signal coefficient. That
heuristic is **not** a uniform fixed-q theorem as r grows; nonlinear remainders
must be controlled before using it as an upper-rate claim.

Even without asserting that rate, the elementary N=o(r) impossibility bound
and the costly pair-based inverse leave a large unresolved interval. The next
attack should first sharpen the order-law bound or exhibit a full-order
statistic that preserves the missing information. Any candidate must survive
removal of matched flat coordinates, waveform labels and phase labels. Exact
wave and effective-matter sprinklings need separate physical-volume sampling.
Finite-event confidence and power requirements must be frozen before a campaign.

The estimate uses standard coupling and oscillatory-integral tools; applying
them is not by itself a new physical law. A matching upper/lower resolution
law, or a counterexample defeating any proposed universal law, is the intended
payoff. This probe gives a specific question and a first bound, not that payoff.

## Fixed-q sufficient side added on 2026-09-24

The subsequent [fixed-energy derivation](full_order_chain_asymptotics.md)
removes the weak-q-first loophole for the exact circular family. It obtains
`P2_lambda-P2_0=-A(q) lambda^4 log(1/lambda)+O_q(lambda^4)` with an explicit
strictly positive A(q), and constructs a wavelength-free pair-count test at
fixed known q. Thus `N lambda^8 log^2(1/lambda) -> infinity` is sufficient
for separation. The [frozen checks and limits](full_order_chain_results.md)
also explain why a fixed list of ordinary chain fractions does not improve
this rate for regular Gaussian tests. The sharp full-order boundary and an
observable attaining it remain unresolved; the research goal is still active.

## Conditional-order improvement added on 2026-09-24

The [null-coordinate construction](null_coordinate_challenge.md) and its
[results/audit](null_coordinate_results.md) subsequently improve the sufficient
region to `N lambda^4/log N -> infinity` for fixed known `0<q<=1`. The test
uses conditional order counts to reconstruct a background null time and
checks its monotonicity along causal pairs. It is wavelength-independent.
The frozen numerical evidence consists of constructed population witnesses,
not a finite-N power campaign. The necessity side is still only
`N lambda -> 0`; exponent 4 is not established as the optimal boundary.

## Subsequent impossibility improvement on 2026-09-24

The [coordinate-oracle calculation](coordinate_oracle_barrier.md) first showed
that the existing `N lambda` comparison is sharp when exact coordinates are
revealed, including a Poisson critical-window law. The subsequent
[quantized-coordinate argument and frozen checks](quantized_order_bound.md)
retain uncertainty inside small cells instead. For any fixed integer K>=2,
`N^(3K+1) lambda^(4K) -> 0` now suffices for full-order indistinguishability.
Thus every fixed power below 4/3 is insufficient, improving the former
exponent 1. This does not establish the logarithmic/critical boundary at
4/3 itself. In the range `N~lambda^(-alpha)`, `1<alpha<4/3`, the arguments
give a strict separation: coordinates plus order distinguish consistently,
while order alone does not. These proofs await independent audit and do not
complete the sharp full-order resolution goal.
