# Independent referee report: proof 2, shared-reference differences

Document under review: docs/paper/paper_a/order_wave_energy/local_null_challenge.md,
Sections 1-4 (claim: N lambda^3 / log N -> infinity implies
TV(P_{N,lambda}, P_{N,0}) -> 1, witnessed by an order-only test).

Referee date 2026-11-24. All derivations were redone from the Rosen metric; all
numerics use code written for this audit (audit_geometry.py), except script s6,
which deliberately feeds the referee's numbers into the author's routine.
Status sentences in the reviewed documents were ignored.

(The referee returned this report as text; the audit orchestrator saved it
unchanged apart from this sentence.)

## 1. Overall verdict

CORRECT -- in the asymptotic regime lambda -> 0, which the boxed statement
leaves implicit (finding 1).

1. Is the written proof complete and correct as written? Every displayed formula
   and inequality of Sections 1-4 was re-derived and confirmed; no false step
   found. The text is a detailed sketch. In two places its justification lives
   outside the document (the standing assumption lambda -> 0; the uniformity of
   the imported bias expansion). Both are closed below with short arguments whose
   ingredients already appear in the companion documents, so they are graded
   MINOR, not GAP.
2. Is the theorem true? I believe so: for fixed known 0 < q <= 1, lambda -> 0
   and N lambda^3 / log N -> infinity, the stated test has false-positive
   probability at most N^-2 and power tending to one, hence TV -> 1. The only
   repair is wording. For lambda bounded away from zero the document proves nothing.

Nothing in Sections 1-4 was left undecided.

## 2. Findings

| # | Location | Grade | What | How checked |
|---|---|---|---|---|
| 1 | Boxed (1); Sections 3-4 | MINOR (scope) | lambda -> 0 is not stated. Every expansion used is a small-lambda expansion. A sequence with lambda fixed and N -> infinity satisfies the hypothesis and is not covered. State "lambda -> 0 and". | Reading; finite-lambda failures in s4 and in the author's own r=64 rows. |
| 2 | Section 4, first paragraph (import of F(theta_X) = v_tilde(X) + O(lambda^4 log(1/lambda))) | MINOR | Uniformity in the endpoint position is asserted, not shown. It is true: 0 <= K_0 - Kbar_lambda <= C lambda^2 min(l,lambda)^2 / l uniformly, so sup_s int |dK| dt <= C lambda^4 log(1/lambda) and all three moments move by that much on compact interior sets (Section 6.1). The constant is large: sup |bias| / (lambda^4 log r) is 2.6e4 at u in {-0.15, 0.05, 0.25} and about 4e4 at u = -0.2 (the latter from an exploratory run, not in a script). | Own derivation; s4 e0 (constant 9.3 for q=1, 2.1 for q=0.25, flat in lambda), s4 e (16 bias(2r)/bias(r) in [1.17, 1.36]). |
| 3 | Section 2, "Choose integer B with 2^B >= 32 N^4" | MINOR | Must be the least such B (B = O(log N)), otherwise (8) fails. Also exp(-B) <= 2^-B is lossy: the union bound gives coverage >= 1 - 0.09 N^-3.77; B ~ ln(32 N^4) would do. | s2 E1, E2. |
| 4 | Section 2, rectangular domain R and "skip if the calibration is not defined throughout R" | MINOR (design; finite lambda only) | R pairs a of one endpoint with p of the other, so a/p^2 sweeps [a_X/p_Y^2, a_Y/p_X^2], which must lie inside (p0(q), 1/10). That window is 4.1% wide at q=1 and 1.0% wide at q=0.25. Consequence: the test as specified SKIPS the frozen q=0.25 witnesses at r=1024 and r=2048, and the q=1 witnesses at r <= 512, for every reference count. Asymptotically harmless (R shrinks like lambda), but not mentioned, and it means the "10 of 12" population witnesses are not all certifiable even with unlimited data. Constants are not uniform as q -> 0. Fix: bound derivatives on the path set {b_tau + sigma d}, or work in (p,f,c). | s6 X4 (author's routine returns no interval at synthetic count 10^40) plus corner computation; q=1 at r=512, 256 confirmed by a one-off run of the same routine (not in a script). |
| 5 | Section 3, eq. (5) | MINOR (confirmed) | Value, sign and factor 2 i sqrt(q) f(u) / (m (M-m)) are right. The O(lambda) remainder has a large constant (error/lambda about 350-450 for q=1): relative error 54% at r=64, 11% at r=256, 0.6% at r=4096. u must be in a compact subset of the OPEN interval. (pi/2, 3pi/2) is not needed; any compact subset of (0, 2 pi) gives a lower bound. | s3 A, B1, B2, C; the r=64 figure is from an exploratory run. |
| 6 | Section 2, citation of "Theorem 1.6 of the 2022 concentration paper" | MINOR | Not consulted and not needed. Hoeffding (1963, Section 5): averaging over permutations of disjoint pairs plus convexity of exp transfers any moment-generating-function bound, hence Bernstein with k = floor(M/2) and the full kernel variance. | Re-derived; s2 D3-D6 (including a negative control with the naive k = choose(M,2)). |
| 7 | Section 4, eq. (8) | MINOR (note) | Exponents are right. At every accessible lambda the THIRD term (H D e_b) dominates, because H/L ~ 2e4: at r=2048, L e_d = 3.3e-9 against H D e_b = 1.6e-8 at N = 7.6e23 (one-off computation). Certified detection of the r=2048 witness needs N ~ 1e24 with ideal constants and a reference count between 1e25 and 1e26 with the implemented enclosures, against lambda^-3 = 8.6e9. The author claims no finite-N power; this is consistent. | s4 INFO lines; s6 INFO line. |
| 8 | Notation | MINOR | f (future fraction and cos^-2), c (cosh 2 beta, chain probability, generic constant), M (reference count and total m), F (calibration and F_lambda), B (integer and boxes). | Reading. |

No ERROR. No GAP that needed a new idea.

## 3. Verdict on each requested step

### 3a. Section 1

- a_X = p_X^2 c_X: correct, under both hypotheses. Given Z_1, Z_2 < X
  (probability p_X^2), the two points are iid uniform in I(p_tip, X), whose
  comparability probability is c_X by definition. The empirical a with
  denominator choose(M,2) is the complete U-statistic of the kernel
  1[Z_1 < X, Z_2 < X, Z_1 ~ Z_2], which is what makes its increments
  U-statistics with nested kernels. (s1 M2: Monte Carlo ratio 1.04 with the
  right denominator; the ancestor-pair denominator is off by a factor of about 100.)
- "Under the background, F equals the exact null coordinate v": correct.
  I re-derived the three calibration formulas of null_coordinate_challenge.md
  eq. (2) from the metric: for an interval with past slack eta,
  Vol = pi eta^2 (K*K)/K (complete the square in the transverse variable),
  and P(comparable) = K C_3 / (3 C_2^2) (the ordered two-chain measure is
  (pi^2 eta^4 / 6) C_3 / K). With K_0 = sin(k l)/k this gives
  alpha = (1 - k l cot k l)/(2q) and p0 = C_1 C_3/(3 C_2^2), p0(0) = 1/10.
  Then 1 - t + t eta - (1-t) zeta = v identically (s1 S2-S5, N2). p0 is
  strictly decreasing on [0,1]: d log p0/dz <= -0.0381 on a 2001-point grid
  (the author certifies -0.0355). Monte Carlo confirms p, f at q=1 and
  p0(6) = 0.0625 (s1 M1, M3; an extra 8 x 20000-event run gave 0.06251 +- 0.00019).
- Nestedness for X < Y, both hypotheses: correct. It uses only transitivity.
  Z < X implies Z < Y; Y < Z implies X < Z; the pair kernel of X is dominated
  by that of Y. The exact wave relation dv >= dy^T Mt^-1 dy / 2 is the causal
  relation of the Rosen patch (the action is a convex free-particle functional,
  no conjugate points) and is transitive; s5 W1b, W2b, W4 find 0 violations.
- "Variance bounded by mean": correct; the kernels take values in {0,1}.

### 3b. Section 2

- Radius (3): correct. Bernstein with true variance sigma^2 has exact root
  t* = B/(3k) + sqrt(B^2/(9k^2) + 2 sigma^2 B/k) <= sqrt(2 sigma^2 B/k) + 2B/(3k).
  The Bernstein event is defined with the TRUE variance (a fixed event of
  probability >= 1 - 2 exp(-B)); on the Hoeffding event sigma^2 <= d <= d_U,
  so r_D is valid on the intersection. Exact binomial enumeration over 280
  cases: failure at most 0.15 x 4 exp(-B) (s2 A-C).
- Transfer to order-two U-statistics, k = floor(M/2): correct (finding 6).
- Union bound: correct. (6N + 6N^2) 2 exp(-B) <= N^-2 reduces to 12N <= 20N^2.
  The count is generous (only about N/2 probes and N^2/8 pairs). Conditioning
  on the probes makes the set of related probe pairs fixed, and the references
  are iid given the probes because the split is an independent random permutation.
- Propagation bound (4): correct, and R is the right domain. The segment from
  b_hat + d to b_hat + d_hat lies in B_b + B_d; the points b_tau + sigma d lie
  in B_b + hull(0, B_d). s2 F1-F2 confirm it on 20000 random boxes; F3 shows
  the bound FAILS if derivative bounds are taken on B_b alone, so the
  hull(0,B_d) part is necessary.
- False positives under the background, every N, including skipped pairs:
  at most N^-2 for every N >= 4, non-asymptotically. Skipping can only remove
  rejections. Data-dependent skipping is harmless: on the coverage event every
  non-skipped pair has theta_X, theta_Y in R, F is the genuine inverse there,
  and the interval contains v_Y - v_X >= 0. s6 X3 confirms this on the
  author's routine with exact background moments.

### 3c. Section 3

- Leading matrix (5): correct, including sign (finding 5). Derivation:
  E_lambda(u') - E_lambda(u) = -(2 lambda sqrt q / c) e^{i phi} int_u^{u'} f e^{i w/lambda} dw
  exactly; the smooth coefficients 1/D, t/D change by O(lambda) and multiply
  O(lambda^2) quantities.
- Lower bound: |e^{i tau} - 1| = 2 |sin(tau/2)| >= sqrt 2 on [pi/2, 3pi/2];
  numerically the exact norm is at least 0.995 of 2 sqrt(2q) f/(m(M-m)) at
  r=4096 (s3 C).
- Positive margin of order lambda^2: correct. The slice map preserves both
  slacks exactly (s3 D1, error 1.7e-16), v_wave = v_tilde - z^T F z/2 + O(lambda^4),
  y = z + O(lambda^2), Jacobian 1 + O(lambda^4). For equal z and equal v_tilde
  the exact margin is -lambda^2 z^T L z / 2 + O(lambda^3)
  (|margin - prediction| / lambda^3 <= 11 for |z| = 0.2), and 13 of 36
  directions (>= 1/3) have margin >= 0.2 |L| |z|^2 lambda^2.
  A more transparent way to see it (not in the text): the u-lines are null
  generators, so (u, v, y) -> (u + tau lambda, v, y) is a causal move with
  margin exactly 0, while v_tilde = v + y^T F_lambda(u) y / 2 + ... changes by
  lambda^2 z^T L z / 2. So v_tilde is not monotone along a causal curve (s3 F).
- O(lambda^3) action correction: correct. |dy| / lambda^2 <= 2.0 and
  action / lambda^3 <= 0.57 over r = 64..4096 (s3 E2).
- Neighbourhood and pair mass (6): correct. With explicit constants
  (q=1: kappa = 0.0189, c' = 0.179; q=0.25: kappa = 0.0108, c' = 0.136)
  every one of 1200 sampled pairs per rung, with second transverse point in the
  ball of radius c' lambda^(3/2) and v_tilde gap in
  (-kappa lambda^2, -kappa lambda^2 / 2), is causally related, with margin at
  least 2.2 kappa lambda^2 (derived floor 1.5 kappa lambda^2), for
  r = 64..16384 (s3 G). The mass is lambda x lambda^3 x lambda^2.
- Conditional upper mass and (7): correct. W counts pairs in the constructed
  set, so the conditional mass is at most the box volume. Pairs of pairs
  sharing two indices contribute <= C N^2 lambda^6, sharing one index
  <= C N^3 lambda^12, disjoint ones nothing. A scaled-box Monte Carlo of the
  full local violating set gives 1.115, 1.124, 1.124 (times lambda^6) at
  r = 256, 1024, 4096 (s3 H).
- Pr(W = 0) -> 0: correct by Chebyshev,
  Var W / (E W)^2 <= C (1/(N lambda^3)^2 + 1/N).

### 3d. Section 4

- Uniform expansion at these endpoints: holds (finding 2, Section 6.1).
  Uniformity can fail only when a sub-segment length s or 1-s tends to 0, when
  eta or zeta tends to 0, or as q -> 0 (ds/dc = 1/(2 q s p0')). None occurs on
  the witness set. The exact wave moment formulas themselves were validated by
  brute-force sampling at large amplitude (s5): volumes agree at |z| < 1 where
  the background formulas are off by 13 to 64 standard errors; the wave pair
  probability is 0.0682 +- 0.0005 against formula 0.0674 and background 0.0625.
- |theta_Y - theta_X| = O(lambda): correct; simplest route is
  theta_wave(X) = theta_bg(u, eta, zeta) + O(lambda^4 log(1/lambda)) with
  (u, eta, zeta) moving by O(lambda). Numerically |d|_1 / lambda is 2.03 (q=1)
  and 2.07 (q=0.25), stable from r=64 to r=1024, with the nestedness signs (s4 f).
- Width (8) and o(lambda^2): correct. With rho = N lambda^3 / log N the three
  terms divided by lambda^2 are rho^-1/2, lambda/rho, sqrt(lambda/rho) (s2 G1).
- Final logic: sound. Witness pairs are probe pairs related in the wave order,
  hence observed. On the coverage event (valid under the wave too, because the
  concentration bounds are distribution-free) every witness pair is eventually
  not skipped and its upper endpoint is at most
  -kappa lambda^2 / 2 + O(lambda^4 log(1/lambda)) + 2 x (8) < 0. All constants
  are uniform in phase (a transverse rotation; the order law does not depend on
  it), in position over the compact witness set, and in lambda -> 0. They are
  NOT uniform in q -> 0.
- Population check: with my own quadrature of the exact moments I reproduce
  the author's 12-row table to 4.6e-7 and the r=2048 witness
  (margin 3.86450424167e-8, decrease 3.864466e-8). On an independent grid
  (other u, tau, direction, phase) all 90 exact margins are positive and F
  decreases for all 18 pairs at r=512 and r=1024; at r=64 none of 18 decrease,
  at r=128 13 fail, at r=256 5 fail (s4 c1, c2, d1-d4).

### 3e. Bookkeeping

- Power of log N: exactly one. It enters only through B ~ 4 log_2 N in the
  Bernstein term sqrt(2 d_U B / k) with d_U ~ lambda; requiring this to be
  o(lambda^2) is literally N lambda^3 / log N -> infinity. The witness count
  needs only N lambda^3 -> infinity; the bias needs nothing.
- Known q <= 1: known q is essential (t, alpha, gamma, p0(q s^2) all depend on
  it). q <= 1 is used only for the certified inverse (series enclosures and the
  monotonicity certificate on z in [0,1]); Sections 3-4 do not use it.
  Numerically p0 stays strictly decreasing up to z = 9.59.
- Probes and the split: both endpoints of a witness pair must be probes,
  because the increment bounds need a reference sample independent of both.
  There are N - floor(N/2) probes, so E W ~ c N^2 lambda^6 / 4; the reference
  count floor(N/2) costs a factor 2 in variance. Constants only.

## 4. The step I trust least

Section 4, the passage from a population witness to a certified rejection. It
rests on (i) a bias bound imported through two earlier documents, and (ii) the
unstated requirement that the rectangular domain R stays inside the calibrated
chain window. I verified both, but this is where every finite-lambda failure
sits and where the constants are worst: bias constant 2.6e4 to 4e4, |grad F| up
to 1.6e4, |Hessian F| up to 3.4e8, and a calibrated window only 1% to 4% wide.
"Eventually" therefore means lambda below about 1/1000 (q=1) or 1/4000 (q=0.25)
for this witness geometry, and N around 1e24 or more. The asymptotic statement
itself I trust.

## 5. What I could not verify

1. The literal statement for lambda not tending to zero. Not proved, not refuted.
2. Any rejection on a sampled order. No simulable N lets the test reject (or
   even evaluate) a pair, so the power side has only population-level evidence,
   the author's and mine.
3. The cited 2022 theorem (not consulted; the classical route suffices).
4. The author's rational-interval implementation as a whole. I checked its
   derivative enclosures at one point, its decision path on one witness and four
   background pairs (s6), not its rounding logic line by line, nor whether
   bits = 40 + N.bit_length() is enough precision.
5. Uniformity beyond sampled points: analytically it follows from Section 6.1;
   numerically I covered u in [-0.3, 0.3], q in {0.25, 1}, r <= 2048 for
   moments and r <= 16384 for geometry.
6. Section 5 of the document (frozen design, hashes, provenance) and the sibling
   lower-bound documents were not reviewed.

## 6. Suggested text repairs

### 6.1 Uniform bias (replaces the import in Section 4)

Kbar_lambda = K_0 sqrt(1 - eps |I|^2 / m^2), eps = 4 q lambda^2/(1 + 4 q lambda^2).
Because |I| <= m, |I| <= C min(l, lambda) and m >= c l uniformly in the
endpoints, 0 <= K_0 - Kbar_lambda <= K_0 eps |I|^2/m^2 <= C lambda^2 min(l,lambda)^2 / l.
Hence sup_s int |dK(s,t)| dt <= C lambda^4 (1/2 + log(1/lambda)), and
|dK| <= C lambda^4 / l_0 on segments of length at least l_0. Expanding
C_2 = K*K and C_3 = K*K*K multilinearly with K_0 <= 1 gives
|dC_2|, |dC_3| <= C lambda^4 log(1/lambda) uniformly over sub-segments of
length at least l_0. So alpha_lambda / gamma_lambda and c_lambda differ from
their background values by O(lambda^4 log(1/lambda)) uniformly for
s, 1-s >= l_0, and V is C^1 on the corresponding compact subset of its domain
because |p0'| >= 0.0018. Only this upper bound is needed; the sharp logarithmic
coefficient of full_order_chain_asymptotics.md is not.

### 6.2 Other edits

- Boxed (1): add lambda -> 0.
- Section 2: "the least integer B with 2^B >= 32 N^4".
- Section 2: say that R rectangular forces [a_X/p_Y^2, a_Y/p_X^2] into
  (p0(q), 1/10), or replace R by the path set.
- Section 3: "compact subset of (-1/2, 1/2)"; optionally add the
  null-generator remark of Section 3c.
- local_null_results.md: note that the q=0.25 population witnesses (r <= 2048)
  and the q=1 ones with r <= 512 are skipped by the test as specified.

## 7. Scripts and observed output

All in this folder. Deterministic, fixed seeds, exit code nonzero on any FAIL.
audit_geometry.py is the shared helper module (no checks of its own).

| Script | Interpreter | Time | Result |
|---|---|---|---|
| s1_background_calibration.py | system python | 21 s | 12 PASS, 0 FAIL |
| s2_confidence_bounds.py | system python | 2 s | 22 PASS, 0 FAIL |
| s3_leading_matrix_and_margin.py | system python | 11 s | 27 PASS, 0 FAIL |
| s4_population_witness.py | system python | 58 s | 18 PASS, 0 FAIL |
| s5_strong_wave_montecarlo.py | system python | 83 s | 6 PASS, 0 FAIL |
| s6_author_routine_crosscheck.py | E:/lab/causal-spacetime/.venv/Scripts/python.exe | 67 s | 4 PASS, 0 FAIL |

Total 89 PASS, 0 FAIL. During development three of my own pass criteria were
wrong and were corrected (not the author's claims): a negative control in s2
that used too small a B; ratio tests in s3 that ignored the oscillating
O(lambda) coefficient; and an arbitrary threshold on the Jacobian constant.

Selected output lines:

    s1 PASS S4 1 - t + t eta - (1-t) zeta = v in the background (any transverse radius)
    s1 PASS N1b d log p0/dz <= -1/50 on [0,1] ... max over grid -0.03809524 (author certifies -0.03555)
    s1 PASS M1 background p ... p formula 0.094944 MC 0.095385 (z=2.13); f formula 0.019977 MC 0.019961 (z=-0.17)
    s1 PASS M3 background pair probability at q=6 matches p0(6), not 1/10 | p0(6) = 0.06249, MC 0.06491
    s2 PASS C  exact binomial: data-driven r_D of eq.(3) failure <= 4 exp(-B) | worst ratio to 4exp(-B) 1.546e-01
    s2 PASS D5 negative control: naive k=choose(M,2) violates its nominal budget | B=10: observed 0.00220 vs nominal 0.00018
    s2 PASS E1 union bound ... worst N^2*failure: paper count 3.193e-03, actual probe/pair count 5.506e-04
    s2 PASS F3 negative control: bounds on B_b alone (not R) are violated | 4051/20000 violations
    s3 PASS B1 q=1.0: eq.(5) leading matrix ... sup error r=256,1024,4096: 1.7199, 0.3457, 0.0944; error/lambda: 440, 354, 387
    s3 PASS B2 q=1.0: sign and factor i of (5) ... opposite sign 1.996, real prefactor 1.409
    s3 PASS E1 q=1.0: ... |margin-pred|/lambda^3 at r=64,256,1024,4096: 11.03, 8.80, 4.70, 7.19
    s3 PASS G1 q=1.0: every pair of the neighbourhood is causally related for r >= 1024 | min margin/(kappa lambda^2): r=64: 3.049, r=256: 2.912, r=1024: 2.479, r=4096: 2.417, r=16384: 2.234
    s3 PASS H  conditional mass of violating partners / lambda^6 is stable ... 1.11477, 1.12400, 1.12366
    s4 PASS c1 independent replication of the 12-row frozen table | max |referee - published| = 4.63e-07
    s4 PASS c2 ... margin 3.86450424167e-08 (published 3.86450424167e-08); F decrease 3.86446618039e-08 (published 3.86446596945e-08)
    s4 PASS e0 q=1.0: |Kbar_lambda - K_0| <= C lambda^2 min(l,lambda)^2/l ... 9.270, 9.084, 9.152, 8.971
    s4 INFO q=1.0: pairs whose calibrated coordinate does NOT decrease, out of 9: r=64: 9, r=128: 6, r=256: 2, r=512: 0, r=1024: 0
    s4 INFO q=0.25: ... r=64: 9, r=128: 7, r=256: 3, r=512: 0, r=1024: 0
    s4 PASS e  q=1.0: ... 16*bias(2r)/bias(r) in [1.172, 1.357]; sup |bias|/(lambda^4 log r): r=64: 14319, ..., r=1024: 25803
    s4 INFO event count at which the ideal-constant version of radius (4) drops below half the r=2048 signal: N ~ 7.61e+23
    s5 PASS W1 q=1.0, lambda=0.5 (eps=0.707): ... p 0.085229 vs MC 0.085245 (z=0.11); f 0.025893 vs MC 0.025840 (z=-0.66); background formula off by z=-64.3, 1.0
    s5 PASS W3 q=6, lambda=0.3: ... wave formula 0.06744, background p0 0.06249, MC 0.06822 (empirical se 0.00054)
    s6 PASS X2 author's routine certifies the witness from the referee's moments (count 10^30) | [-3.879062863e-08, -3.849869999e-08]
    s6 INFO synthetic reference count needed ...: 10^24: upper +9.22e-08; 10^25: upper +3.62e-09; 10^26: upper -2.50e-08
    s6 PASS X4 ... q=0.25, r=1024: [0.097688, 0.101401] vs (0.099031, 0.1) -> SKIPPED; q=0.25, r=2048: [0.098604, 0.100462] -> SKIPPED; q=1.0, r=2048: [0.097139, 0.099002] vs (0.095919, 0.1) -> evaluated

## 8. Independence and limits of this audit

- The geometry code derives the causal relation, the slice map, F_lambda, the
  kernel hierarchy and the calibration from the metric; its primitive is a
  composite Gauss-Legendre rule, not the author's ODE. It reproduces the author's
  archived moments to 7e-15.
- The Monte Carlo checks use only the pairwise causal criterion and the density
  a(u)^2; they do not use the volume or chain formulas they test. Their
  precision is about 0.1% to 0.3% for volumes and about 1% for pair
  probabilities, so they catch structural errors, not lambda^4 effects.
- q = 6 appears in two Monte Carlo checks only because the pair formula is far
  from 1/10 there; it is outside the theorem's range.
- The sibling folders under independent_audit/ were not read. No existing file
  was modified and no state-changing git command was run.
