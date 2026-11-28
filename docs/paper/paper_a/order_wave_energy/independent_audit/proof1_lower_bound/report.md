# Independent referee report: quantized_order_bound.md (order-only lower bound)

Referee: independent AI referee, 2026-11-24. Scope: Sections 1-5 of
`docs/paper/paper_a/order_wave_energy/quantized_order_bound.md`, with the
definitions it imports from `full_order_resolution_probe.md` and
`coordinate_oracle_barrier.md`. Sibling referee folders were not read. Status
sentences in the sources were ignored; every formula was re-derived by hand
and, where possible, tested by a script in this folder.

(The referee returned this report as text; the audit orchestrator saved it
unchanged apart from this sentence.)

## 0. Overall verdict

**CORRECT.**

No false step and no step whose conclusion fails. The text is a compressed
proof: several constants and quantifiers are implicit and the Section 4 lemma
is a sketch. Each such place is a MINOR finding below, with the argument that
completes it. None changes a statement, an exponent, or the final claim.

1. Proof as written: correct; complete at research-note level. A reader must
   supply (a) the link between the edge constant L, the cell h and c0;
   (b) why an edge has one orientation for both kernels; (c) the routine
   parts of the smoothing lemma; (d) the Hellinger convention and
   "#components with an edge <= #uncertain pairs".
2. Theorem: believed true. For fixed K >= 1, 0 < q < pi^2 there is C_K(q) with
   TV(P_{N,lam}, P_{N,0}) <= C_K [ N lam^2 + (N lam^4)^(1/(2K+1))
   + (N^(3K+1) lam^(4K))^(1/(2K+1)) ] for all N >= 1, 0 < lam <= 1, uniformly
   in phase. No repair was needed.

CORRECT does not mean: finite-lam content (F7), anything at alpha = 4/3, or a
numerical confirmation of the conclusion itself (Section 5).

## 1. Findings table

| # | Location | Grade | What is wrong / missing | How checked |
|---|---|---|---|---|
| F1 | Sec. 2, (3)-(4), lines 60-94 | MINOR | (3) lacks its hypothesis (crossing set); constants not given. True with e <= C_E lam min(ell,lam), C_E = 2 sqrt(q)(4 sec^2(sqrt(q)/2) - 2), c0 = 2 sqrt2 C_E; on the wave null cone the margin is at most 2 sqrt2 e dv/(b+dv). Only e <= C delta and e < b are used. | Sec. 3.2; check_b_sandwich.py B1-B4: 0 violations in 9 x 1.1e6 reference pairs and 24 x 1.2e5 null-targeted pairs, separations to 1e-7 |
| F2 | Sec. 3, lines 111-115 | MINOR | "fixed large L" and "c0 delta << h" not linked; need (L - 1 - sqrt3) h >= c0 delta. | D1: 0 violations in 1.47e7 pairs at the tightest cell; D1n (L = 1): 154 violations |
| F3 | Sec. 3, lines 121-123 | MINOR | Orientation for both kernels asserted, not argued; follows from (4) on the backward margin or from T being a common time function. | D2, B4 |
| F4 | Sec. 4, lines 158-176 | MINOR | Lemma is a sketch, no constant. True with L = [(m+1) rbar + max_j rbar'_j]/h^2, rbar = prod sup(h rho_i), rbar'_j = sup(h^2 rho_j') prod_{l != j} sup(h rho_l). | Sec. 3.5; check_e_* scripts |
| F5 | Sec. 5, lines 218-227 | MINOR | (8) derived only for c1 delta <= h <= h0; true for all N, lam after enlarging C_K. Chosen h is not the minimiser of (7) when N h dominates (factor up to 1840 in my grid): rate only; the vanishing condition of (7) optimised over h is exactly (1). | X1-X4, R1-R9 |
| F6 | Sec. 5, lines 202-213 | MINOR | Hellinger convention unstated (need H^2 = sum (sqrt p - sqrt q)^2); "#components with an edge <= #uncertain pairs" unstated. | H1-H5 |
| F7 | Lines 8-22, 229-232 | MINOR (caveat) | Purely asymptotic: C_K contains (K+1)^(K-1) and 2^(K(K-1)/4); margin (4K - (3K+1) alpha)/(2K+1). alpha = 1.3 needs K >= 14, margin <= 0.0035. c0 = 6.4, 18.1, 132, 13532 for q = 0.25, 1, 4, 9. | R7 |
| F8 | Lines 24-29 | MINOR (scope) | "Strict separation" also needs kappa(q) > 0 from coordinate_oracle_barrier.md, not audited. MC fractions (0.32, 0.63, 1.21) lam^2 at lam = 1/32 agree with the quoted 0.322, 0.632, 1.138. | B4 |
| F9 | Line 9 | MINOR (remark) | K = 1 valid: returns N lam -> 0. | R4, R6 |

No ERROR, no GAP.

## 2. Verdicts on the items of the brief

(a) TRUE. (b) TRUE, explicit constants supplied. (c) TRUE; measured TV/lam^2 =
1.61, 3.28, 8.05 for q = 0.25, 1, 4. (d) TRUE under (L - 1 - sqrt3) h >=
c0 delta. (e) TRUE. (f) TRUE. (g) follows; K = 1 gives alpha < 1; alpha = 4/3
not covered for any K.

## 3. Re-derivations

Notation: a = cos(sqrt(q) u), f = a^-2, m(s,t) = int_s^t f, s = sinh(2 beta)
= 2 lam sqrt(q), c = cosh(2 beta), t = s/c < 1, delta = lam^2,
fmax = sec^2(sqrt(q)/2), h_ax = 1.

### 3.1 (a) Exact family
* For any Rosen metric the only non-zero Ricci component is
  R_uu = -tr(g^-1 g'')/2 + tr[(g^-1 g')^2]/4 (SymPy, general gamma). For the
  circular family with constant beta, R_uu = 2q - sinh^2(2 beta)/(2 lam^2),
  zero exactly at sinh(2 beta) = 2 lam sqrt(q) (40 digits). Same a(u) for
  every lam.
* gamma^-1 = a^-2 [c I - s T(u/lam + phi)], so M = c m I - s J, det gamma =
  a^4, and in the y chart Mtilde = m I + E, E = -(s/c) J symmetric trace-free
  with eigenvalues +-e, e = t |int f e^{iu/lam}|.
* A causal curve has u non-decreasing; parametrised by u, v' >= x'^T gamma
  x'/2, and Cauchy-Schwarz gives int x'^T gamma x' >= dx^T M^-1 dx with
  equality on x' = gamma^-1 p, a null geodesic. Curves between points of the
  slab stay in it, and the Rosen chart covers it because a > 0 (this is
  0 < q < pi^2). So X <_lam Y iff du > 0 and 2 dv >= dy^T Mtilde^-1 dy
  (or du = 0, dy = 0, dv >= 0), exactly.
* Background: with m' = a^-2, T = (m+v)/sqrt2, Z = (m-v)/sqrt2, g0 =
  a^2(-dT^2 + dZ^2 + dy^2), Ricci = 2q du du, relation dT >= |dz|, volume
  a^2 du dv d^2y = a^4 dT dZ d^2y. The reference domain is the Minkowski
  diamond between the tips (D0).

### 3.2 (b) Operator bound and sandwich
Integration by parts: |int_s^t f e^{iu/lam}| <= lam (f(s) + f(t) + TV f) <=
lam (4 fmax - 2), and trivially <= m(s,t) <= fmax (t - s). Hence
e <= C_E lam min(ell, lam), and e <= t b < b for every lam.

Let u1 < u2, b = m(u1,u2), Q = dy^T (bI+E)^-1 dy in
[|dy|^2/(b+e), |dy|^2/(b-e)], g = dT - |dz|, g (dT + |dz|) = 2 b dv - |dy|^2,
dT = (b + dv)/sqrt2.
* g >= c0 delta > 0: then dv > 0 and 2 b dv - |dy|^2 >= c0 delta (b+dv)/sqrt2
  >= 2 e dv once c0 delta >= 2 sqrt2 e, so |dy|^2 <= 2 dv (b - e), Q <= 2 dv.
* 2 dv >= Q: then dv >= 0, |dy|^2 <= 2 dv (b + e), so
  g >= -2 sqrt2 e dv/(b + dv) >= -2 sqrt2 e.
* u1 = u2: relation forces dy = 0, g = 0. u1 > u2: no forward relation, and
  g >= c0 delta is impossible.
So (4) holds for all ordered pairs with c0 = 2 sqrt2 C_E; the reverse
orientation by exchanging points. g is monotone in dv, so the smallest valid
c0 is sup |g|/delta on the wave null cone. Measured: 2.96, 7.16, 36.6 for
q = 0.25, 1, 4 (saturating, below c0 = 6.39, 18.1, 132.4); q = 9 still
rising (2141 at lam = 1/4096, c0 = 13532). The supremum is smaller for short
separations (proportional to ell/lam for ell << lam).

### 3.3 (c) One-point laws
Both densities are proportional to a^2 on different domains, so
TV(mu_lam, mu_0) <= max(A(D_lam \ D_0), A(D_0 \ D_lam))/V_0. For fixed (u,y)
both sections are v-intervals whose ends move by at most
|y|^2 e_-/(2 m_-(m_- - e_-)) and the upper analogue, on
|y|^2 <= 2 h_ax (m_- + e_-). Hence
A(symdiff) <= pi h_ax^2 (1+t)^2/(1-t) int a^2 (e_- + e_+) du
           <= 2 pi h_ax^2 (1+t)^2/(1-t) C_E lam^2.
Product measures give C N delta. Measured TV/lam^2 = 1.61, 3.28, 8.05,
stable over lam = 1/8 ... 1/128: the first term of (7) is of the right order.

### 3.4 (d) Marks, uncertain pairs, bad marks
* Revealing (cube index, labelled order) refines the unlabelled order; the
  mark law is common; given marks the locations are independent.
* Boundary cubes: <= C N h (measured 6.3 h per event at q = 1). In an
  interior cube the conditional T density lies between rho_min/(rho_max h)
  and rho_max/(rho_min h) with derivative <= sup|d_T rho|/(rho_min h).
* Non-edges: moving inside two cubes changes a margin by at most
  (1 + sqrt3) h; a non-edge has both centre margins beyond L h, so both true
  margins exceed (L - 1 - sqrt3) h >= c0 delta with fixed sign and (4) fixes
  the relation for both kernels. D1: 0 violations; all 11 297 genuine
  wave/background disagreements lie on edges.
* Uniform bound: an uncertain second cube forces
  | |T_Y - T_cX| - |z_Y - z_cX| | <= (L + (1+sqrt3)/2) h, two T-intervals per
  z_Y, so probability <= C_unc h, C_unc = 4 (L + (1+sqrt3)/2)
  vol_3(proj_z D)/V_0, uniformly in the first event and the grid (worst
  sampled first event 0.075 against 0.259 at h = 0.004). E U <= C N^2 h.
* Orientation: disjoint cube T-intervals fix sign(dT) and exclude the reverse
  relation for both kernels; otherwise |dT_c| < h and |dz_c| <= (L+1) h,
  probability O(h^4) per pair.
* Trees: a component with more than K vertices contains K+1 vertices spanned
  by a tree of uncertain pairs. Peel a leaf: its location is independent of
  the rest and the edge event depends only on leaf and parent, costing
  C_unc h conditionally. Induction gives (C_unc h)^K with no independence
  between edges; (K+1)^(K-1) trees and C(N,K+1) sets give (5). Measured
  scalings per halving of h: 3-vertex trees 3.9 and 3.8, 4-vertex trees 7.7
  and 7.3.

### 3.5 (e) Smoothing lemma
Statement. rho = prod rho_i on [0,1]^k, 0 <= rho_i <= R_i, rho_i Lipschitz
with constant R'_i; rows +-(e_j - e_i) over distinct pairs;
p(tau) = int rho prod 1{a_e.t >= tau_e}. Then p is C^1 and grad p is
Lipschitz with constant (m+1) prod R_i + max_j R'_j prod_{l != j} R_l.

Proof. Integrating t_j first, p(tau + eta e_e) - p(tau) =
-int_0^eta G_e(tau + s e_e) ds with G_e the slice integral of
F = rho x (other indicators) x (cube) on t_j = t_i +- tau_e.
(i) Other threshold: the change of G_e is the integral of rho over a slab of
width eta inside a (k-1)-dimensional section, transverse with unit Jacobian
whether or not the rows share a vertex: <= eta prod R.
(ii) Own threshold c: G_e(c) = int F(t_i, t_i + c, rest). For fixed "rest"
the admissible t_j form [lo, hi], from the cube faces of t_j and the other
rows through j; no other row contains both i and j, so lo, hi do not depend
on t_i. The indicator of t_i + c in [lo, hi] changes only on two t_i-intervals
of length eta: <= (2 prod R + R'_j prod_{l != j} R_l) eta.
(iii) Gershgorin. Dependence among three or more rows never enters.
Scaling t = h s gives L = C_K/h^2. The author's sketch is the same
computation and is right.

|grad p|^2 <= 2 L p: descent inequality at tau - grad p/L, using p >= 0 and
the bound on all of R^m. Then sqrt p is sqrt(L/2)-Lipschitz where p > 0; on a
segment with zeros, split at the first and last zero. Global, zero set
included.

Sandwich step: conditionally on spatial coordinates, the wave indicator of an
edge lies pointwise between 1{t_j - t_i >= theta + c0 delta} and
1{t_j - t_i >= theta - c0 delta}. For a pattern the product is squeezed, so
p(tau + c0 delta 1) <= p_lam(S) <= p(tau - c0 delta 1), with p_0(S) = p(tau)
inside; hence |sqrt p_lam - sqrt p_0| <= sqrt(L/2) 2 c0 delta sqrt(m) =
C_K delta/h. Sum over <= 2^m patterns; joint convexity of H^2 removes the
conditioning.

Evidence. Exact rationals: Hessian quotients never exceed 1.0000 (repeated),
0.9998 (mixed) at steps 1/8 ... 1/4096; |grad p|^2/(2p) reaches 2.9987 at
the cyclic-triangle corner (true L >= 3, my bound 4); sqrt p passes across
the zero set everywhere. Non-uniform densities: second differences below
0.24 L_rho. Exact wave family, one uncertain pair of cubes: |sqrt p_lam -
sqrt p_0| about 1.5e-2 while p_0 runs 0, 1.2e-3, 2e-2, 0.5; three rows have
p_0 = 0 and p_lam > 0, all under the bound. The author's control values
(9/2 eps^2 - 4 eps^3; H^2/eps^2 = 13.44677) are reproduced exactly.

### 3.6 (f) Assembly
With H^2 = sum (sqrt p - sqrt q)^2: TV <= H; H^2 subadditive over
independent factors; TV of (marks, order) is the mark-average of conditional
TV. On good marks TV(order | marks) <= sqrt(U) C_K delta/h and
E sqrt U <= sqrt(E U) <= sqrt(C N^2 h): this is (7).
Algebra: with h = (delta/N^K)^(2/(2K+1)), N^(K+1) h^K = N delta h^(-1/2) =
(N^(3K+1) delta^(2K))^(1/(2K+1)); N h = (N delta^2)^(1/(2K+1));
N^2 h^4 <= (N h)^3. Under the hypothesis with N >= 1: lam -> 0, h -> 0,
N h -> 0, N lam^2 -> 0, and h/delta -> infinity because
N lam^(2 - 1/K) <= N lam^(4K/(3K+1)), the exponent gap being
(2K+1)(K-1)/(K(3K+1)) >= 0.

### 3.7 Explicit form
With M = K(K-1)/2 and L_K from 3.5, for c0 delta/(L - 1 - sqrt3) <= h <= h0:
TV <= C_1 N delta + C_bd N h + C_amb N^2 h^4
      + (K+1)^(K-1) C(N,K+1) (C_unc h)^K
      + 2 c0 sqrt(2^M M L_K h^2/2) (delta/h) sqrt(C_unc N^2 h/2).

### 3.8 (g) Final statement
Main exponent (4K - (3K+1) alpha)/(2K+1), positive iff alpha < 4K/(3K+1);
increasing strictly to 4/3. Every fixed alpha < 4/3 is covered by
K = floor(alpha/(4 - 3 alpha)) + 1; not uniform as alpha -> 4/3.
* K = 1: valid, gives N lam -> 0; no Hellinger term needed.
* alpha = 4/3: exponent -4/(3(2K+1)) < 0 for every K; nothing proved.
  Sequences lam^(-4/3) g(lam) with g -> 0 sub-polynomially are not covered.
* The vanishing condition of (7) optimised over h is exactly (1).

## 4. Step trusted least
The Section 4 lemma for general K, specifically the uniform bound on the
repeated second derivative. The improvement from exponent 1 to 4/3 rests on
the per-component Hellinger cost O(delta^2/h^2) uniformly in the thresholds,
including zero reference probability, and this is where the text is a sketch.
Confidence high: elementary slicing proof with explicit constant, exact
tests to k = 5, flat |Delta sqrt p| in the exact family. Numerical coverage
stops at k = 5 (uniform) and k = 3 (non-uniform).

## 5. Not verified
1. The conclusion itself (TV between order laws) was never computed.
2. The lemma numerically for K >= 6; non-uniform densities for k >= 4.
3. The size of C_K(q); the limit q -> pi^2.
4. kappa(q) > 0 and the oracle-side statements.
5. Whether 4/3 is the true boundary.
6. A Poisson (random N) version.
7. Section 7 provenance. The recorded control numbers are reproduced.

## 6. Scripts and observed output
Deterministic, ASCII, PASS/FAIL lines, nonzero exit on FAIL; system python
(numpy, sympy, mpmath). 2026-11-24: 301 PASS, 0 FAIL, all exit codes 0.

| Script | Tests | PASS | FAIL | Seconds |
|---|---|---:|---:|---:|
| check_a_exact_family.py | A1-A5 | 19 | 0 | 12 |
| check_b_sandwich.py | B0-B4 | 82 | 0 | 22 |
| check_c_one_point_law.py | C0-C2 | 45 | 0 | 13 |
| check_d_bad_marks.py | D0-D6 | 23 | 0 | 5 |
| check_e_smoothing_lemma.py | E1-E7 | 48 | 0 | 67 |
| check_e_lemma_density.py | F1-F4 | 10 | 0 | 2 |
| check_e_wave_pair_component.py | G1-G4 | 56 | 0 | 32 |
| check_f_assembly_exponents.py | H1-H5, X1-X4, R1-R9 | 18 | 0 | 2 |

audit_common.py holds shared helpers; nothing is imported from the author's
code. Two non-errors worth recording: the stated h is not the optimiser of
(7) for small alpha (my first draft of test R8 assumed it was and failed; the
proof did not), and the edge constant L must exceed 1 + sqrt 3 (control D1n).
