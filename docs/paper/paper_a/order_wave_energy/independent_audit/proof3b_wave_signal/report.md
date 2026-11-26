# Referee report 3b: wave-side coupling, noise comparison and signal

Subject: gram_fourier_noise.md Sections 4-6 (eqs. (6)-(11)) and
gram_fourier_challenge.md Section 3 (boxed (4)). Reviewed 2026-11-24.
Independent referee; sibling audit folders not opened.

(The referee returned this report as text; the audit orchestrator saved it
unchanged apart from this sentence.)

## 1. Verdict

CORRECT, within scope, with (10) conditional on the background bound (5),
which I treated as a hypothesis and did not audit.

(i) Is the written argument complete and correct as written? Correct: every
step I checked is true. Complete: as a sketch yes; as a proof it needs two
short arguments written out (F2, F3; supplied in 3.4). If "compressed" counts
as GAP in your rubric, those are the only two places, both repaired. One
sentence is literally false (F1); it moves a certified number in the seventh
digit and no inequality.

(ii) Do I believe (10) and (11)? (11): yes, unconditionally. (10): yes,
conditional on (5). No repair needed.

## 2. Findings (all MINOR)

| # | Location | What | How checked |
|---|----------|------|-------------|
| F1 | noise Sec. 6 "|a_2x a_2y|=|G12|<=0.00002"; guard_check term ell_bound*off.hi/F(348,1000); results "Im D_A >= 0.01069797" | Identity false; supremum 2.0063e-5; certificate overstated by 1.8e-7; all targets survive (Sec. 4) | algebra; check_05 |
| F2 | noise Sec. 5 | Chain is prose only, identity T1-T4 never written; bounds stated "on the active compact set" while the score sums over all probes | re-derivation |
| F3 | noise Sec. 6 "anchor pair slacks stay positive" | Certified expression is the order-derived slack, so positivity means b_ij > 0; unstated, reads as circular | algebra; exact rational check in check_05 |
| F4 | noise Sec. 6 selection | Limit true; quarter-tolerance neighbourhood has probability ~7e-38 per candidate triple, so needs log m >> 1e37 | check_05, order of magnitude |
| F5 | "calibration errors are o(lam^2)" | True, constants extreme: max|u_hat-u| ~ 5.2e3 lam^4 log(1/lam) against |D| ~ 0.02; calibration shift of the coefficient is 16.3%, 0.57% of the limit at r = 256, 512 | check_04 |
| F6 | results table, r = 64 | Stored row converged only to 0.5%: mine -0.18061-0.03881i; author's routine at factor 4 gives -0.180576-0.038817i | own code; one read-only re-run of population(...) |

## 3. Re-derivations

### 3.0 Setting

a = cos(k u), k = sqrt(q), f = a^-2, m(u) = int_{-1/2}^u f, M = m(1/2).
Wave relation in the rescaled chart: v2 - v1 >= (1/2) dy^T Mt(u1,u2)^-1 dy,
Mt(s,t) = dm I - B(s,t), B trace-free symmetric. Store [[b1,b2],[b2,-b1]] as
b1 + i b2; then y^T B y = Re(beta conj(zeta)^2), zeta = y1 + i y2, and
beta(s,t) = (2 lam k/c) e^{i phi} int_s^t f e^{iu/lam} du, c = sqrt(1+4q lam^2).
|beta| <= C lam min(l, lam), dm >= l. Sign and helicity confirmed by
integrating gamma^-1 entrywise (check_01).

Plane-wave null rotations preserve relation and volume, giving exactly:
interval volume pi eta^2 (K*K)/K with K = a(s)a(t) sqrt(det Mt), and pair
probability K (K*K*K)/(3 (K*K)^2) (1/10 in flat space). Canonical coordinates:
z = (W0 Q)^(1/2) y, Q = M_L^-1 + M_R^-1, W0 = m(M-m)/M,
v_tilde = eta + |z|^2/(2m). The canonical image of the wave diamond is exactly
the background diamond, and the reconstructed state (u_hat, v_hat, R_hat^2)
is an exact function of (u, v_tilde, |z|^2), hence rotation invariant in z.
K_lam/K_0 - 1 ~ -(1/2)|beta|^2/dm^2 <= C lam^2 min(1, lam^2/l^2); integrated
against K_0 K_0 this gives one-point state errors O(lam^4 log(1/lam)).

### 3.1 Eqs. (6)-(7) and the shifts

E_L = B_L/m, E_R = B_R/(M-m), norms eps_L, eps_R, t = m/M. A trace-free
symmetric 2x2 matrix squares to a multiple of I, so
M_L^-1 = (I + E_L)/(m (1 - eps_L^2)) and

    W0 Q = alpha_0 I + E,
    alpha_0 = (1-t)/(1-eps_L^2) + t/(1-eps_R^2),
    E = (1-t) E_L/(1-eps_L^2) + t E_R/(1-eps_R^2)   (trace-free),
    det T = (alpha_0^2 - |E|^2)^(-1/2).

alpha_0 - 1 = O(eps_L^2 + eps_R^2) and |E| <= C(eps_L + eps_R), so
|det T - 1| <= C(eps_L^2 + eps_R^2) with no truncated expansion.
eps_L <= C lam min(1, lam/(u+1/2)), eps_R <= C lam min(1, lam/(1/2-u)); the
far side is <= C lam^2 <= e(s). T exists for every u once 2 lam k f_max < 1.
(6) holds, tips included.

(7): canonical density is a^2 det T/Z_lam on the background diamond (Jacobian
det T; the v-shift is a translation); slice volume pi W0 ~ s.

    s < lam:  int_0^lam lam^2 s ds = lam^4/2
    s > lam:  int_lam^{1/2} (lam^4/s^2) s ds = lam^4 log(1/(2 lam))

The log is genuine: near a tip det T - 1 = -(1/2) eps_L^2 and
eps_L^2 ~ lam^4 |e^{is/lam} - 1|^2/s^2.

Transverse: |T - I| <= C(eps_L+eps_R), |z| <= C sqrt(s), so
|dy| <= C lam min(sqrt s, lam/sqrt s); mean against slice volume
lam [int_0^lam s^(3/2) ds + lam int_lam s^(1/2) ds] = O(lam^2).
Longitudinal: P = m M_L^-1, S = (M-m) M_R^-1, W0 Q = P + t(S-P),

    T M_L^-1 T - I/m = (1/m)[(P+t(S-P))^(-1/2) P (P+t(S-P))^(-1/2) - I],

norm <= C t (eps_L+eps_R)/m. With |z|^2 <= 2 W0:
|dv| <= C t(1-t)(eps_L+eps_R) <= C lam^2 everywhere. "Long denominator" is
only a device; both expressions agree by T(M_L^-1+M_R^-1)T = I/W0 (verified
to 1e-15).

### 3.2 Eq. (8), anchors fixed

phi(P,P') = Minkowski null margin in the chart (m, v, y), 2-Lipschitz in each
point.
(a) Sandwich: wave and background relations of two points agree unless
|phi| <= c0 lam^2; from Mt^-1 = (dm I + B)/(dm^2 - |beta|^2) and
|beta| <= C lam min(l, lam); uniform as l -> 0.
(b) The canonical-to-physical map displaces a point by
d(xi) <= C(lam^2 + lam min(sqrt s, lam/sqrt s)).
So a mismatch forces |phi(xi,xi')| <= 2 d(xi) + 2 d(xi') + c0 lam^2.

- Both random: split by the larger displacement, condition on that point,
  Pr(|phi(xi,.)| <= w) <= C w uniformly, integrate: C(E d + lam^2).
- One endpoint fixed in the compact anchor/probe region: its own d is
  O(lam^2). Uniform in the fixed point.
- 1{A<R<X}, A and X fixed: R is interior, all shifts O(lam^2); union bound.
  R close to A or X in u is covered by (a).
- Ancestor pair: two fixed-endpoint indicators plus one both-random.

Why the tip stays O(lam^2) although the shift reaches lam^1.5: for the random
point at tip distance s, the u-gap to the fixed endpoint is bounded below, so
the margin is Lipschitz in transverse coordinates with a bounded constant. A
mismatch then requires v_tilde in a window of width C(lam^2 + |dy|)
<= C lam^2 (1 + s^(-1/2)) for each (u, z). The slice carries volume ~ s:

    int s * lam^2 (1 + s^(-1/2)) ds < infinity,

and the worst zone s < lam contributes int_0^lam s * lam sqrt(s) ds
= (2/5) lam^3.5. The large shift lives where there is almost no mass.
Position coupling by (7) adds lam^4 log(1/lam).

### 3.3 Eq. (9)

Coupled pairs rho_j = (R_j^lam, R_j^0) are iid; a mismatch kernel D has
|D| <= 1, E D^2 <= C lam^2. Rosenthal for the centred mean:
C_p [lam/sqrt(n) + lam^(2/p) n^(-(1-1/p))]; second term <= first iff
(n lam^2)^(1/2-1/p) >= 1. Order two: Hoeffding's average over pairings into
floor(n/2) independent pairs, then convexity of the L^p norm.

Exact, p = 4, Bernoulli(pi = lam^2): mu4 = n pi(1-pi)[1 + 3(n-2)pi(1-pi)], so

    || centred mean ||_4 = (lam/sqrt(n)) (3 + 1/(n lam^2))^(1/4) (1 + O(lam^2)),

<= sqrt(2) lam/sqrt(n) iff n lam^2 >= 1, unbounded ratio below. n lam^2 >= 1
enters twice: here, and in 1/n <= lam/sqrt(n) (T2, T3). Per block
n ~ lam^(-alpha)/log, alpha > 2, so n lam^2 -> infinity however close alpha
is to 2. Only p = 2, 4 are used.

### 3.4 Section 5: the chain and (10)

X = canonical coordinates of a wave probe (law nu_lam); mu_lam(X), mu_0(X)
population vectors there; eps_lam, eps_0 centred errors from coupled lists.
F = A e^{-i omega u}, F(mu_0) = 0. Exact identity:

    S_n - E_X F(mu_lam(X))
      = [avg_X F(mu_lam) - E_X F(mu_lam)]                         T1
      + avg_X [F(mu_lam+eps_lam) - F(mu_lam+eps_0)]               T2
      + avg_X [F(mu_lam+eps_0) - F(mu_lam) - F(mu_0+eps_0)]       T3
      + avg_X F(mu_0(X)+eps_0(X))                                 T4

|DF| <= C(1+omega|A|), |D^2F| <= C(1+omega+omega^2|A|), globally.

T1: |F(mu_lam)| <= C lam^2, independent probes: lam^2/sqrt(n). Pointwise.
T2: on the segment |A| <= C(lam^2+|eps_lam|+|eps_0|); term
  <= C(1 + lam + omega(|eps_lam|+|eps_0|)) |eps_lam-eps_0|; by (9), p = 2, 4:
  C lam/sqrt(n) + C omega n^(-1/2) lam n^(-1/2) = C(lam/sqrt(n) + 1/n).
  Pointwise; triangle inequality over probes.
T3: = int int D^2F(mu_0 + t delta + s eps_0)[delta, eps_0],
  delta = mu_lam - mu_0 = O(lam^2), |A| <= C(|delta|+|eps_0|) there. Term
  <= C(omega + omega^2 lam^2 + omega^2|eps_0|) lam^2 |eps_0|
  = C(lam|eps_0| + |eps_0|^2): C(lam/sqrt(n) + 1/n). Pointwise. The term
  omega A(mu_lam) eps = lam|eps| sits here, handled correctly.
T4: conditional probe fluctuation <= C/(n n_P)^(1/2) ~ 1/n; law change
  int ||F(mu_0+eps_0)(x)||_2 |nu_lam - mu_0|(dx) <= C lam^4 log(1/lam)/sqrt(n);
  remainder E_{X~mu_0} F(mu_0(X)+eps_0(X)), bounded by (5) + C/n.

Only T4's remainder uses the oscillatory gain, and it is entitled to it:
background kernels of iid mu_0 references, probes integrated against mu_0,
background phase. Wave kernels are never asked for the gain (they are not
monotone along canonical u-lines). The wave-specific terms are pointwise at
lam/sqrt(n), the order of the leading background term: (8) has no slack
against exponent 2. Sum: (10), with no whole-sample agreement needed.

Argument supplied for F2: for X outside an enlarged active set, A vanishes on
a fixed neighbourhood of mu_0(X), and mu_lam(X) is within O(lam^1.5) of it.
T2, T3 are then nonzero only if |eps| >= const, probability exp(-c n) by
Hoeffding, uniformly in X.

Argument supplied for F3 (Section 6): by definition of g_ij,
dv - (R_i^2 + R_j^2 - 2 g_ij)/(2 dm) = eta0_ij = sqrt(b_ij gamma/alpha).
The certificate bounds this from below by 0.0289 on the windows, hence
b_ij > 0 and the pair is timelike. Independent check: a non-timelike pair has
g_ij = (R_i^2+R_j^2)/2 - dm dv <= -0.0028, -0.0288, 0.0004 for pairs
(1,2), (1,3), (2,3), outside the windows around 0, 0.0064, 0.0064, any q >= 0.

### 3.5 Signal: challenge Section 3 and (11)

First order in lam^2, canonical coordinates, dz = z - z_A, tau = T - I
(trace-free at first order), dv = -(1/2) z^T F z, y = z + tau z,
Mt^-1 = (1/dm)(I + [B(u)-B(u_A)]/dm):

    g_AX = z_A.z + dm (dv_X - dv_A)
           - dz^T [B(u) - B(u_A)] dz/(2 dm)
           - dz.(tau_X z - tau_A z_A) + O(lam^4 log(1/lam)).

Angular averaging over the direction of z at fixed (u, v_tilde, |z|^2) is
exact for the observer-side weight and phase (3.0). It kills every term
linear in z and every z^T(trace-free)z. Survivors: the local term
-z_A^T B(u) z_A/(2(m - m_A)); z_A^T B(u_A) z_A/(2 dm) (anchor phase, smooth
in u); constants. The change of ell multiplies z_i.z, mean zero. So

    <H> = Re[beta(u) D_A(u)] + smooth + O(lam^4 log),
    D_A = -(1/2) sum_i ell_i conj(zeta_i)^2/(m(u) - m_i),

the author's D for the frozen triple.

B expansion (three integrations by parts):

    beta(-1/2,u) = (2k/c) e^{i phi} [e^{ix/lam}(-i lam^2 f + lam^3 f'
                   + i lam^4 f'')]_{x=-1/2}^{x=u},

leading local term -2 i sqrt(q) lam^2 f e^{i(u/lam+phi)}, endpoint term
smooth, remainder O(lam^3). Confirmed.

Boxed limit: at omega = 1/lam + delta, (1/2) beta D e^{-i omega u} has
non-oscillatory part -i sqrt(q) lam^2 f e^{i phi} D e^{-i delta u}; the
conjugate sits at 2/lam and smooth terms at 1/lam, both O(lam^infinity)
against a smooth compact weight. Probe density a^2 det T, a^2 f = 1:

    S(omega)/(lam^2 E w) -> -i sqrt(q) e^{i phi}
         int w_u D_A e^{-i delta u} du / int w_u a^2 du.

Sign and normalisation as boxed.

Non-cancellation. Im D_A = sum ell_i a_ix a_iy/(m - m_i) >= 0.0107 on the
neighbourhood (Sec. 4), |D_A| <= 0.0424.
- Weight: any W >= 0, separable or not: Im int W D_A >= 0.0107 int W.
- Offset |delta| <= 1/8: |e^{-i delta(u-0.2)} - 1| <= 0.0125, modulus
  >= (0.0107 - 0.0125*0.0424) int W = 0.0102 int W.
- Rotation by theta: factor e^{-2 i theta}.
- Reflection: D_A -> conj(D_A), modulus unchanged.
Reconstructed coordinates: omega(u_hat - u) = O(lam^3 log), 2.4e-4 rad at
r = 512; weight shift multiplies a first-order quantity. Both second order.

Selection: raw moments concentrate uniformly (Hoeffding, union over O(log n)
candidates); inner set compact in the regular calibration domain; disjoint
triples independent with p0 > 0, so acceptance probability
1 - (1-p0)^(ceil(log m)) -> 1; wave moments within O(lam^2). Correct as a
limit under both hypotheses (F3, F4).

## 4. F1 in full

Text: "|a_2x a_2y| = |G12| <= 0.00002". With a_1 = (sqrt(G11), 0):
a_2x = G12/sqrt(G11), a_2y = sqrt(G22 - G12^2/G11), so

    a_2x a_2y = G12 sqrt(G22/G11 - G12^2/G11^2).

Equal to G12 only if G22 = G11. On the box (G11, G22 in [0.00638, 0.00642],
|G12| <= 2e-5) the supremum is 2e-5 * 1.0031 = 2.0063e-5 (scan: 2.0043e-5).

Certificate: Im D_A >= 0.079626*0.079374/(1.30*0.452) - 1.0094*2e-5/0.348
= 0.010756 - 0.000058 = 0.01069797. Corrected second term uses 2.0063e-5:
bound 0.0106978, lower by 1.8e-7.

Surviving targets:
- Im D_A > 0.01: 0.0106978, margin 7.0e-4.
- Grid coefficient 0.0106978 - 0.0424/80 = 0.0101679 > 0.009 (code target),
  margin 1.2e-3; text comparison 0.001 < 0.01 untouched.
- |D_A| < 0.05: unaffected (uses only G_ii), 0.0424.
True infimum on scan: Im D_A = 0.01378, max |D_A| = 0.0260.

## 5. Numerical tables (q = 1 unless stated)

(6)-(7), tips to s = 1e-7:

| r | sup|detT-1|/e^2 | TV/lam^4 | sup|dv|/lam^2 | sup|dy|/bound |
|---|---|---|---|---|
| 32 | 16.80 | 83.06 | 4.47 | 2.66 |
| 128 | 9.08 | 136.34 | 4.58 | 3.07 |
| 512 | 8.00 | 192.66 | 4.70 | 3.18 |
| 1024 | 17.66 | 221.80 | 4.75 | 3.20 |

A priori constant for (6): 81.6. TV/lam^4 rises 26-29 per doubling:
TV ~ 40 lam^4 log(1/lam).

(8), Pr/lam^2 envelopes over probes and relative phase:

| Kernel | r~32 | r~128 | r~512 |
|---|---|---|---|
| past 1{R<X} | 2.02 | 1.79 | 1.93 |
| future 1{X<R} | 0.176 | 0.170 | 0.171 |
| interval 1{A<R<X}, anchors fixed | 1.48 | 1.54 | 1.80 |
| pair 1{R~R'} | 0.82 | 0.82 | 0.82 |
| ancestor pair, fixed X | 0.10 | 0.10 | 0.11 |
| tip part (s < 4 lam) | 0.093 | 0.011 | 0.0013 |

Exploratory interval run to r ~ 2048: 1.52, 1.87, 1.74, 1.86.

B expansion: max |B - lead|/lam^3 = 4.42 over 24 (r,u), r to 2048; local
coefficient -0.0032-1.9988i, 0.0017-1.9997i, 0.0004-2.0000i at r = 128, 512,
2048 (claimed -2i).

Signal, S/(lam^2 E w), phase 0.37, theta = 0, limit 0.01729433+0.00824625i:

| r | referee value | vs author's row | from limit | calibration off |
|---|---|---|---|---|
| 128 | -0.02276391-0.06884441i | 1.9e-5 | 453% | 1.5% |
| 256 | 0.01524720+0.00587606i | 3.5e-5 | 16.3% | 0.47% |
| 512 | 0.01739303+0.00821975i | 6.0e-7 | 0.53% | 0.089% |
| 1024 | 0.01729266+0.00825177i | not run | 0.030% | 0.043% |

At r = 512: theta = pi/8 0.56%; reflected 0.60%; offsets +-1/8 0.54%;
q = 1/4 1.2% (0.04% calibration off); phase 2.1 0.54%. Modulus 0.0191-0.0192
against certified floor 0.0106.

## 6. Least trusted, not verified, disclosure

Least trusted: the hand-off to (5) for the T4 remainder. (10) is only as good
as (5) applied to the same regularised F, uniformly over the anchor good set.
In scope, the estimate whose failure would do most damage is (8) being
exactly O(lam^2); verified analytically and numerically.

Not verified: (5) and Sections 1-3, 7; (10) end to end (analytic only, no
simulation of the finite-order score); existence and nesting of the concrete
smooth cutoffs; uniformity beyond r <= 1024 (signal), r <= 804 (kernels),
q in {1, 1/4}; selection rates (limit only); explicit constants (mine are
observations, not bounds).

Disclosure: thresholds test boundedness, not predicted constants, and were
fixed after exploratory runs. The first version of check_03 compared single
r values and reported one FAIL (past kernel 1.86 at r = 512 against 1.23 at
r = 32); cause is phase dependence (r = 45.76 already gives 2.02), so I
switched to envelopes over clusters of r. Exploratory runs not in scripts:
interval kernel to r ~ 2048; kernels and check_02 at q = 1/4; r = 64 with
refinement. Author's code: read for definitions; run once, read-only, for F6.
No existing file modified; no state-changing git command.

## 7. Scripts

numpy only, deterministic, ASCII; system python 3.11. 53 PASS, 0 FAIL.

| Script | PASS | FAIL | Time |
|---|---|---|---|
| check_01_primitive_and_B.py | 6 | 0 | 1 s |
| check_02_canonical_density.py | 8 | 0 | 1 s |
| check_03_kernel_mismatch.py | 8 | 0 | 68 s |
| check_04_population_signal.py | 19 | 0 | 67 s |
| check_04b_population_r1024.py | 3 | 0 | 95 s |
| check_05_anchor_neighbourhood.py | 9 | 0 | 1 s |

Support modules: wavegeom.py, popsignal.py. check_01, 02, 05 also pass under
the venv python.
