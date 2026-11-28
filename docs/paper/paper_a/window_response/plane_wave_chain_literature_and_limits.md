# Plane-wave chains: literature comparison and statistical limits

2026-09-22. Follow-up to `plane_wave_chain_hierarchy.md` at `0fb08fa`.
This is a bounded primary-source audit and a deterministic theory extension.
It changes no frozen experiment, estimator, gate, or previous result.
The derivations and checks are by the same AI, not independent expert review.

## 1. Decision after the reference audit

Keep the direction: **exact unweighted chain integrals in finite plane-wave
diamonds, and the limits of flat-dimensional interpretation near focusing**.
Narrow the contribution. Neither higher-chain geometric inference, nor the
idea that three-point probabilities resolve a two-point degeneracy, is new.
Nor are Jacobi propagators, Van Vleck determinants, Gaussian composition,
or multiplicity-dependent caustic singularities new.

No equation-identical treatment of the arbitrary-profile unweighted hierarchy
was retrieved in the sources below. This is a candidate distinction, not a
finding of priority. In particular, the hierarchy may be regarded as a useful
evaluation of known convolution machinery; that editorial risk remains even
if no identical formula is found. A first-Weyl-detection claim is excluded.

### Source ledger and equation-level overlap

Equation numbers refer to the linked HTML versions unless PDF pages are given.
Sources were inspected on 2026-09-22; abstract-only material and targeted
source inspection are explicitly distinguished from a complete proof audit.

| Primary source / inspected location | Established content and effect on our claim |
|---|---|
| [Roy, Sinha, Surya, *The Discrete Geometry of a Small Causal Diamond*, arXiv:1212.0631v1](https://arxiv.org/html/1212.0631), Eqs. (7)-(9), (16), (30)-(31), sections 5-6 | Flat chain means, nested interval integrals, all-length leading Ricci corrections and curvature-sensitive dimension inference are prior. Our candidate is an exact plane-wave evaluation, including vacuum Weyl effects beyond that expansion. |
| [Kambor, Nomaan X, *Manifold Properties from Causal Sets using Chains*, arXiv:2007.03835v3](https://arxiv.org/html/2007.03835), Eqs. (8)-(16), (59) | Consecutive-chain geometric estimators and density-normalized chain moments are prior. The ratio used here is a combination of such known normalized moments, not a new general statistical method. |
| [Aghili, Bombelli, Pilgrim, *Discrete Spacetime: a Web of Chains*, arXiv:1807.08701v2](https://arxiv.org/html/1807.08701), section 2, Eqs. (6)-(7) | The fixed-N falling-factorial correction and flat-dimensional chain abundance agree exactly after the index conversion below. They are not contributions of this project. |
| [Bombelli, *Statistical Lorentzian geometry and the closeness of Lorentzian manifolds*, arXiv:gr-qc/0002053v2](https://arxiv.org/html/gr-qc/0002053), sections II-IV, displayed Eqs. (26)-(33) | Already gives poset probability laws, discusses total variation, and constructs a rectangle/cylinder pair with equal two-point and unequal three-chain probabilities. The general two-to-three-point idea is a direct priority overlap. Our proposed distinction is the exact focusing/high-dimension asymptotic and its sampling limit. |
| [Nomaan X, Dowker, Surya, *Scalar Field Green Functions on Causal Sets*, arXiv:1701.07212v5](https://arxiv.org/html/1701.07212), Eqs. (15), (18), (25)-(34) | Powers of the causal matrix and their expected chain counts are already spacetime convolutions. Four-dimensional link propagators additionally involve empty intervals. Our unweighted chain reduction must be compared with the former; it is not a new convolution principle or a calculation of link counts. |
| [Harte, Drivas, *Caustics and wave propagation in curved spacetimes*, arXiv:1202.0540](https://arxiv.org/html/1202.0540), Eqs. (64), (67), (95), (98), (110), (116) | World functions, Jacobi determinants, simple/double conjugacy singularities, and the constant-profile sine/sinh example are prior. The candidate addition is their consequence for normalized chain abundances and dimension calibration, not a new classification of caustics. |
| [Díaz, Dorn, *Propagators and WKB-exactness in the plane wave limit of AdS x S*, arXiv:hep-th/0406240](https://arxiv.org/html/hep-th/0406240), section 3.1, Eqs. (13)-(19) | Exact scalar kernels with the Van Vleck square-root factor and flat functional form are prior. Their displayed plane wave is isotropic. Appendix D concerns unequal-radius product backgrounds; it is not an arbitrary-profile unweighted-chain theorem. Section 3 below makes the weight distinction explicit. |
| [*Gravitons in a gravitational plane wave*, EPJC (2024)](https://link.springer.com/article/10.1140/epjc/s10052-024-12986-1), sections 2 and 4, Eqs. (61)-(65) | Related exact-kernel/propagator treatment. It is not used as a normalization authority for our probability integrals; its propagator and density conventions have not been fully reconciled here. No unweighted-chain priority clearance is inferred from this inspection. |
| [Alfyorov-Shnyukov, *Weyl curvature from the Hasse diagram: a parameter-free bridge formula for causal sets*, Paper 7](https://inspirehep.net/files/32388edfda552e0b4cfb6697bda66d7b), section 2.3 (PDF p. 6), Theorem 5.6 (p. 19), section 6 | The coupled flat/curved Hasse statistic and conditional bridge to an electric Weyl invariant are related work. Its observable and assumptions differ from intrinsic unweighted chain probabilities. This prevents broad first-Weyl language; it does not by itself establish equality of the two estimands. |
| [Ahmed-Bombelli, *Indications of Weyl curvature in a causal set*, APS 2026 abstract](https://meetings-archive.aps.org/smt/2026/apr-r93/1/) | Abstract inspected: plane gravitational waves and causal-set geodesic deviation; it refers to a planned paper. No full paper was retrieved by the exact-title search. Detailed overlap remains unresolved, rather than cleared. |
| [Ghosh, *Sample-Complexity Bounds and an Exact Bulk-Sector Channel Obstruction for Curvature Estimation on Finite Causal Sets*, v4, 2026-08-14](https://zenodo.org/records/21931074), source sections 2, 6.1, 6.3-6.4, 7-8 | The conditional weak-curvature counting-cost argument, position-space Fisher bound and fixed-orientation bulk-sector chain obstruction are prior related results. Section 6.4 explicitly separates the light-cone/ordering contribution. Our full-domain plane-wave integral includes that contribution; (L9) is a two-law testing bound near conjugacy, not the same curvature-estimation bound or a refutation of it. |

The Zenodo browser fetch failed, but direct public API access succeeded.
The distributed `Complexity-of-Causet-Curvature_V4.zip` contains the inspected
`Sample-Complexity-of-Causet-Curvature.tex`, SHA-256
`8f948873e975e796014a155b98c20b615faa1e1907f35494347f1afe29b772be`.
Its section 6.1 labels the order/position Fisher-scaling argument heuristic;
our proof of (L9) does not assume it or H1-H3. Targeted source inspection does
not certify that paper's numerical or algebraic proofs.

The fixed-N comparison is literal, not just qualitative. This note's k counts
interior sampled events, whereas the chain between fixed tips in ABP Eq. (6)
has j=k+1 edges. Substitution gives

\[
 \langle c_{k+1}\rangle_{\rm ABP}
 =(N)_k\left[\frac{\Gamma(d+1)}2\right]^{k-1}
 \frac{\Gamma(d/2+1)\Gamma(d)}
 {\Gamma(kd/2+1)\Gamma((k+1)d/2)}
 =\binom Nk P_k^{\rm flat}(d).
 \tag{L1}
\]

The identity follows from Gamma recursion. It was also checked numerically
at d=2,3,4,6,10 and k=1,...,6, with N=37. The largest log difference was
`1.07e-14`. Changing chain-length conventions would otherwise manufacture a
spurious disagreement with the literature.

### Search boundary

Queries included `plane wave causal chains spacetime volume`, `plane waves
causal sets chain`, `plane wave chain abundance`, `plane wave causal chain
Volterra`, `causal set caustic dimension`, `Myrheim caustics`, `Myrheim conjugate
dimension`, `plane wave Riesz spacetime Green functions`, and `plane wave heat
kernel exact Van Vleck`, plus exact-title searches for the competing work.
Backward references led to the propagator and statistical-geometry papers.
Irrelevant Fourier-wave results were excluded. This is not an exhaustive
citation-network, dissertation, or unpublished-manuscript search. The APS
full-paper comparison and exhaustive propagator priority remain open. No
author has been contacted.

## 2. Scope and notation retained

Use exactly the finite-diamond hypotheses of the hierarchy note: a smooth
symmetric transverse profile H(u), full transverse space, no quotient or walls,
and a strictly positive Dirichlet quadratic form on the entire endpoint
interval. Nonzero final determinant alone is insufficient. Set d=n+2,
nu=d/2, L=b-a, and let h>0 be the effective endpoint slack after subtracting
the classical transverse action. Vacuum additionally requires tr H=0.

Write K(s,t)=sqrt(det B(t,s)), C_r=K to the r-th Volterra power, and

\[
 I_k=\frac{(2\pi)^{kn/2}h^{k\nu}}{\Gamma(k\nu+1)}
       \frac{C_{k+1}(a,b)}{K(a,b)},\qquad
 P_k=\frac{k!I_k}{I_1^k}.
 \tag{L2}
\]

These are ordered integrals and unordered-subset chain probabilities,
respectively. They count all chains, without a link or maximality condition.
The quantities C_r are continuum kernel powers, not sampled counts.

## 3. Exact comparison with Van Vleck weights

The known plane-wave Van Vleck determinant is

\[
 \Delta(s,t)=\frac{(t-s)^n}{\det B(t,s)}.
 \tag{L3}
\]

Define J_k by taking the same ordered event integral as I_k and multiplying
each of its k+1 consecutive edges, including both endpoint edges, by
sqrt(Delta). The extra weights depend only on the u endpoints, so the
transverse/slack integration used for (L2) still applies. For each edge,

\[
 K(s,t)\sqrt{\Delta(s,t)}=(t-s)^{\nu-1}.
\]

The remaining simplex integral is exactly the flat Dirichlet integral:

\[
 \boxed{J_k=
 \frac{(2\pi)^{kn/2}h^{k\nu}}{\Gamma(k\nu+1)K(a,b)}
 \frac{\Gamma(\nu)^{k+1}L^{(k+1)\nu-1}}
      {\Gamma((k+1)\nu)}
 =\sqrt{\Delta(a,b)}\,I_k^{\rm flat}(L,h).}
 \tag{L4}
\]

This identifies the cancellation explicitly, including for noncommuting H.
It is compatible with the familiar flat-form-times-Van-Vleck structure of
propagators. J_k is a weighted causal-indicator integral, not automatically
the heat kernel or Green function in the cited papers: their differential
equation, support/distribution, mass and boundary conditions still matter.
In particular, J_1 is not the diamond volume.

For the unweighted observable the residual kernel is K, not the flat power
of t-s. Thus the profile-dependent C_{k+1} in (L2) does not disappear.
This is an estimand distinction, not proof that (L2) has never been evaluated
before. Equation (L4) is treated as an explanatory corollary of established
Jacobi/Gaussian machinery, not as a new propagator result.

## 4. All higher chain ratios separate the population limits

In d=4 define R_k=P_k/P_2^(k-1) for each fixed integer k>=3. Direct
cancellation in (L2) yields

\[
 \boxed{R_k=\frac{k!2^k3^{k-1}}{(2k)!}
      \frac{C_{k+1}C_2^{k-2}}{C_3^{k-1}}.}
 \tag{L5}
\]

Under the caustic hypotheses of the hierarchy note, C_r tends to a finite
positive c_r for r>=2. Therefore every R_k has a finite positive limit.
This extends the previously worked three-chain case without choosing a
favorable chain length after seeing simulations.

For the known flat-dimensional Gamma formula, Stirling expansion at fixed k
gives

\[
 P_k^{\rm flat}(d)\sim A_k d^{(k-1)/2}e^{-\psi(k)d},\quad
 A_k=\Gamma(k)(\pi/2)^{(k-1)/2}\sqrt{k(k+1)/2},
\]
\[
 \psi(k)=\tfrac12[k\log k+(k+1)\log(k+1)]-k\log2.
 \tag{L6}
\]

Consequently

\[
 \boxed{R_k^{\rm flat}(d)\sim
 \frac{\Gamma(k)\sqrt{k(k+1)/2}}{3^{(k-1)/2}}
 e^{-\lambda_k d},\qquad
 \lambda_k=\psi(k)-(k-1)\psi(2)>0\quad(k\ge3).}
 \tag{L7}
\]

Proof of the strict sign: psi(1)=0 and
psi''(x)=(1/x+1/(x+1))/2>0, so its secant slope from 1 to k>2 is greater
than that from 1 to 2. For k=3 this recovers
lambda_3=(3/2)log(4/3) and prefactor sqrt(8/3).

For the constant vacuum profile H=diag(z^2,-z^2) on [0,1], the numerical
values at delta=pi-z=1e-6 are:

| Ratio | Value near first conjugacy | Flat large-d decay rate lambda_k |
|---|---:|---:|
| R3 | 0.262779789 | 0.431523109 |
| R4 | 0.0281844832 | 1.159281024 |
| R5 | 0.00152758963 | 2.114052276 |
| R6 | 0.0000476118610 | 3.253224584 |

The decimals are quadrature diagnostics, not certified limiting constants.
The positive-versus-zero limit follows from the proof for every fixed k,
with no uniform-in-k claim.

The real-valued inverse F^-1(P_2) is a calibration parameter. Noninteger d
does not define an ordinary d-dimensional Minkowski sampling law. To state
the contrast using actual geometries, take integer d tending to infinity and
choose a sequence of wave parameters with P_2=F(d). Such a sequence exists
by continuity and P_2 tending to zero with a positive leading coefficient.
Then (L5) and (L7) differ along this sequence of genuine geometries.

## 5. Exact full three-point law, not just its chain component

Bombelli's finite-poset-law framework already supplies the general question.
Here the plane-wave conditional-volume moments evaluate all five answers.
For a uniformly sampled X in the diamond, put A=V_p(X)/V and B=V_f(X)/V.
Let p=P_2, t=P_3, and a2=E[A^2], b2=E[B^2]. In the notation alpha, beta,
gamma, W of hierarchy Eq. (18), its moment formula gives

\[
 a2=\frac{\int_a^b W(u)\alpha(u)^2du}{15\gamma^3},\quad
 b2=\frac{\int_a^b W(u)\beta(u)^2du}{15\gamma^3},\quad
 E[AB]=t/6. \tag{L8}
\]

The unlabelled three-point order, retaining time orientation, has law

| Order | Probability |
|---|---|
| Antichain | 1-3p+3(a2+b2) |
| Exactly one related pair | 3p+t-6(a2+b2) |
| One minimum, two incomparable future points | 3b2-t |
| One maximum, two incomparable past points | 3a2-t |
| Three-chain | t |

To derive the forks, choosing the minimum in three ways gives 3E[B^2],
which includes each chain once; subtract t. Reverse time for the other fork.
The expected number of related pairs is 3p; normalization then gives the
remaining rows. This is known combinatorial bookkeeping applied to the exact
plane-wave moments, not a new general identity for posets.

For a flat four-dimensional diamond the vector is exactly
`(265,66,9,9,1)/350` in table order. In a varying profile the two fork
probabilities need not agree. The rotating-profile check gives approximately
`0.0253775332` and `0.0253410378`, while preserving normalization and the
pair marginal. No inference of a physical time arrow is made from that
profile-dependent asymmetry.

## 6. A finite-sample indistinguishability bound

Let mu_N and nu_N be two actual laws of the entire induced order on N sampled
interior points. Suppose the respective probabilities for a uniformly chosen
pair to be related are p and q. The unique N-point antichain A_N obeys

\[
 1-\mu_N(A_N)\le\binom N2p,\qquad
 1-\nu_N(A_N)\le\binom N2q.
\]

Since total variation is one minus the shared probability mass, its shared
antichain mass alone gives the sharper bound

\[
 \boxed{\operatorname{TV}(\mu_N,\nu_N)
 \le 1-\min\{\mu_N(A_N),\nu_N(A_N)\}
 \le\min\{1,\binom N2\max(p,q)\}.}
 \tag{L9}
\]

For equal prior probabilities, every classifier based on the entire order has
success probability at most (1+TV)/2. For matched p=q asymptotic to a positive
constant times delta^(m/2), all such classifiers approach chance whenever
N^2 delta^(m/2) tends to zero. This is an elementary consequence of the
antichain mass, stronger in scope than a warning about one noisy ratio.
It does not exclude classification when that scaling condition fails.

Thus N=o(delta^(-m/4)) cannot retain information in this matched limit.
Seeing even one k-chain additionally requires that
N^k delta^(m(k-1)/2) not tend to zero, by the first-moment bound.
Neither condition is sufficient for estimation or a 90%-power design.
The bound can be vacuous at practical N: for N=1200 and delta=1e-6 it exceeds
one before capping. No useful finite-N power claim is inferred there.

There is no contradiction between sections 4 and 6. Population ratios retain
a contrast while their constituent rare-event probabilities vanish. The
limits of infinite sampling and approaching the caustic cannot be interchanged
without a sampling-rate argument. At or beyond conjugacy the full-model
diamond has infinite volume and the uniform fixed-N law itself is undefined.

## 7. Reproducible checks and remaining work

Run `python docs/paper/paper_a/window_response/check_chain_limits.py` in the
existing experiment environment. It writes `chain_limits_audit.json` with
LF source hashes for itself and the reused `check_plane_wave_chains.py`.
The latter script and its previous JSON were not changed. The new checks are:

- ABP finite-N normalization at 30 dimension/length combinations.
- Stirling prefactors and rates for k=3,...,6 at integer d=100 and 1000.
- Exhaustive verification of the five-state reconstruction on all 19 labelled
  strict posets on three points, independent of the geometric quadrature.
- Flat exact three-point probabilities, a noncommuting rotating vacuum profile,
  and constant-vacuum caustic distances 1e-3 and 1e-6. Grid refinement from
  1024 to 2048 panels changes P1,...,P6 by less than `6.7e-8` relatively.
- Agreement between the C4 chain calculation and the separate star integrals
  through E[(A+B)^2]=a2+b2+t/3. The full N=3 laws satisfy (L9) against an
  actual flat four-dimensional law; this is not a fractional-dimension test.

These checks passed. Floating-point agreement is not a certified enclosure,
independent scientific replication, or new confirmatory experiment.

The theorem-focused contribution candidate is now: the arbitrary-profile
unweighted reduction, its complete fixed-length focusing asymptotics, and
the explicit distinction between population identifiability and finite-sample
information. The full three-point law is an additional computable corollary.
The general higher-chain idea and weighted-propagator machinery must be cited
as prior work, not advertised as discoveries.

Remaining work before a submission claim is an independent proof audit
(especially the noncommuting composition and first-conjugacy hypotheses),
resolution or explicit qualification of the open priority items, and a
coherent theorem-focused manuscript. A practical classifier claim would need
a separately frozen, power-sized sampling study. The existing manuscript PDF
and archived bundle predate this continuation.
