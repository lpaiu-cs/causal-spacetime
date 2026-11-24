# Independent referee report: gram_fourier_noise.md Sections 1-3 and 7

Referee: independent AI referee, 2026-11-24. Scope: the background noise
bound (5) and the frequency scan of
`docs/paper/paper_a/order_wave_energy/gram_fourier_noise.md`. Sections 4-6
(wave-side coupling, signal) were refereed separately; their displayed
conclusions (10) and (11) are treated here as hypotheses. Sibling referee
folders were not read.

(The referee returned this report as text; the audit orchestrator saved it
unchanged apart from this sentence.)

# Part A: summary

## Overall verdict: CORRECT

This covers Sections 1-3 and 7 only. Every step was re-derived; no false step
was found and nothing needed more than a routine argument.

- **Which statistic:** (5) is proved for the regularised statistic Section 1
  postulates. That statistic is not written down; the one frozen in
  `gram_fourier_challenge.md` eq. (2) has only the weight `w` and is not
  covered.
- **Conditional:** theorem (1) also needs (10), (11) and the Section 6 anchor
  selection, which this referee did not review.
- **(i) Written proof:** correct, and complete as an outline. The omitted
  pieces are routine and are supplied below.
- **(ii) Belief:** yes. (5) holds for all `omega >= 1` with `n_P >= n`; the
  scan delivers exactly `alpha > 2 + 2/K`, strictly.

## Findings (none ERROR or GAP)

| # | Location | Grade | What is wrong | How checked |
|---|---|---|---|---|
| F1 | Sec. 1 extension | MINOR | The regularised statistic is postulated, not constructed. The cutoff must act on all 18 raw moments, not through `w` alone. | Derivation; `s1` I6 (`H = 0.0029` for a probe not after anchor 3) |
| F2 | Sec. 1, Sec. 7 | MINOR (strengthening) | With the cutoff in `b`, (2) and (5) hold for every anchor triple. "Good anchor" and the pilot-failure term are unnecessary under the null. The one-line argument is missing. | Derivation |
| F3 | Sec. 2, last sentence | MINOR (wording) | "Independent noise ... loses (3)" is misleading: independent references give a smaller variance. Shared references add the common-mode term that (3) controls. | `s2` D: 1.188e-3 vs 2.276e-3 |
| F4 | Eq. (5) | MINOR | `1/n` is listed twice; `n_P >= c n` and the anchor conditioning are only implied. | Derivation |
| F5 | Elementary check | MINOR | Exact only for `omega = 2 pi j`; general `omega` carries `1 - sinc^2(omega/2)`. | `s2` sympy, 5 Monte Carlo cases, largest z 1.29 |
| F6 | Sec. 7 null vote | MINOR | The bound is `(K C_K/log n)^2`, not `C_K/(log n)^2`. Immaterial. | `s6` L4 |
| F7 | Scan cap, eq. (1) | MINOR | The cap is the only binding requirement. A grid up to `sqrt(n)` keeps the null control and, given (10), yields `alpha > 2 + 2/(K+1)`. | `s6` L4b, L5, L6 |
| F8 | Sec. 7 "eventually" | MINOR (practical) | With `C_K = c_q = 1`, `1/lambda` enters the scan at `lambda ~ 1e-40..1e-43` for `alpha = 2+2/K+0.1`. At entry `lambda^2 < t_n` in 4 of 12 cases. | `s6` L6 |
| F9 | Constants | MINOR (practical) | At q=1, `sd(u_hat) ~ 74.5/sqrt(n)`: about 2e6 references before `u_hat` is inside the window. | `s1` Monte Carlo |
| F10 | Remainder | observation | In the analogue the remainder term of (5) is never attained; noise stays near `c/(omega sqrt n)` up to `omega = n/4`. | `s5` B |

## Fitted exponents (1+1 analogue, L2 norm, bootstrap SE)

| Quantity | Claim | Fitted |
|---|---|---|
| Order 1 | `omega^(<=-1) n^(-1/2)` | `omega^(-1.140 +- 0.010) n^(-0.480 +- 0.009)` |
| Order 2 | `omega^(<=0) n^(-1)` | `omega^(-0.051 +- 0.014) n^(-1.009 +- 0.022)` |
| Order 3 | `omega^(<=1) n^(-3/2)` | `omega^(0.733 +- 0.018) n^(-1.476 +- 0.034)` |
| Remainder after K=3 | `omega^(<=3) n^(-2)` | `omega^(1.107 +- 0.038) n^(-1.865 +- 0.048)` |
| Full statistic, first-order regime | `omega^(<=-1) n^(-1/2)` | `omega^(-1.127 +- 0.033) n^(-0.502 +- 0.014)` |
| Probe term, `n_P = n` | `n^(-1)` | `n^(-1.016 +- 0.059)` |

Two 95% intervals narrowly miss the claimed value: order 1 by 0.002 and the
remainder by 0.04. Free-intercept fits contain both (`-0.486 +- 0.009`,
`-1.912 +- 0.050`). The referee reads this as misfit of the joint power law,
because order 1 is also checked without a fit: simulated variance over the
deterministic influence-function variance is 0.861 to 1.124 across 15 points,
worst z 2.52.

## Least-trusted step

The interface around the regularised statistic (F1). The null bound needs
only smoothness and compact support. But the same cutoff must be the one for
which (10) and (11) are proved, and it must equal one on the wave population
vectors. Within the text, "a zero extension then preserves (2) for every
background probe" is true only because of the cutoff in `b`.

## Could not verify

- (5) for the real 3+1 instrument by simulation (F9); only the mechanism in
  1+1 and the exact inputs in 3+1.
- The constants `C_K`, `c_q`, and when the asymptotics begin.
- (10), (11) and Section 6, which were out of scope.
- The `s`-dependence of `c = p0(q s^2)` at q=1 numerically: Monte Carlo
  agrees but cannot separate it from 0.1. The formula is confirmed at q=6.

## Scripts and final-run output

All are in this folder. Each is deterministic, ASCII, and exits 0.

| Script | PASS | FAIL | Time |
|---|---|---|---|
| `s1_background_identity.py` | 15 | 0 | 37 s |
| `s2_elementary_ecdf.py` | 12 | 0 | 8 s |
| `s3_toy_validation.py` | 10 | 0 | 29 s |
| `s4_toy_taylor_orders.py` | 11 | 0 | 65 s |
| `s5_toy_full_statistic.py` | 7 | 0 | 67 s |
| `s6_scan_logic.py` | 9 | 0 | 2 s |
| `toy_common.py` (module) | - | - | - |

After that run the referee normalised line endings to LF; `s2` and `s6` were
re-run in full and the rest only re-compiled.

# Part B: report body

## B1. Step-by-step verdicts

**1a. Section 1.**

- **Identity is exact.** The metric is `a^2(-2 dm dv + dy^2)` with
  `dm = du/a^2`. The interval section at fixed `u` has volume
  `pi slack^2 m1 m2/(m1+m2)`. With `a = cos(ku)`,
  `a^2 m1 m2/(m1+m2) = sin(k(u-u_A)) sin(k(u_X-u))/(k sin(k(u_X-u_A)))`,
  which depends on differences only. So
  `b_AX = alpha(u_X-u_A) slack^2/gamma` and `g_AX = y_A.y_X`. Three plane
  vectors give `sum ell_i y_i = 0`, hence `H = 0`.
- **Conditions used:** `X` interior, `A1<A2<A3<X` strictly timelike, `B`
  nonsingular.
- **Numerical check (`s1`):** interval law to 3.8e-12; nine order
  probabilities against 6e5 sprinkled events, largest z 1.97;
  `max|g - y.y| = 1.1e-12` and `max|H| = 5.2e-11` over 400 random chains;
  `max|H| <= 2.6e-13` for the frozen anchors at q = 1, 0.25, 0.01.
- **Extension is legitimate.** Take `A_ext = chi * w * H` with `chi` a smooth
  cutoff in the full raw-moment vector. Inside its support all `b > 0`, so
  the configuration is a timelike chain and (2) holds. This is true for every
  anchor triple (F2).
- **U-statistic versus ratio is consistent.** `a_hat` is unbiased, so
  `epsilon` is centred. The conditional ancestor-pair fraction equals
  `c_hat (1-1/n)/(1-1/(n p_hat))`; it would break exact centring by a smooth
  `O(1/n)` factor but not the bound.

**1b. Section 2.**

- **Monotonicity holds for all four kernels.** `u -> (u,v,y)` is null and
  future-directed, so pasts are nested. `s1` I7a: zero violations over 60
  lines x 1500 references, 19443 jumps seen. A transverse-line control shows
  59 references switching on and off.
- **Projections:** `h_1` has variation at most 2, canonical `h_2` at most 4,
  uniformly.
- **(3) is the right Stieltjes bound** and applies per reference
  configuration with constants independent of the references.
- **Coefficient, density, weight and boundary cooperate.** The coefficient
  carries `w` at the true `(u,v,R^2)`. Its `u`-support lies in `[0.1,0.3]`,
  strictly inside the diamond (`eta >= 0.617`, `zeta >= 0.150`), so endpoint
  terms vanish.

**1c. Section 3.**

- **Degree `k-1`:** only the `j = k` Leibniz term has degree `k`, and it
  multiplies `A(mu_0) = 0`. The phase is the true background `u`.
- **Index counting:** every slot is centred given the rest, so a lone index
  kills the term. Survivors number at most `C n^((L+L')/2)`, each bounded by
  `omega^(d+d') (C/omega)^2`. Diagonal terms and bias are inside the second
  moment.
- **Remainder:** `F(mu_0 + t eps) = t G(t)` and
  `|G^(K)| <= C omega^K |eps|^(K+1)`.
- **Moment bound and probe term:** standard;
  `E|S_n - Sbar|^2 = E Var_X(F)/n_P`. `s3` V6 gives ratios 1.02 to 1.13
  against that formula.
- **Domination:** for `omega <= sqrt(n)` each sum term is at most the first;
  on the scan range the remainder is too.

**1d. Elementary check.** The diagonal is `1/(6n)`. Off-diagonal, integration
by parts gives `(1-|phi|^2)/(n omega^2)` with `phi(2 pi j) = 0`. Confirmed.

**1e. Section 7.**

- Blocks are independent.
- Null vote probability is at most `(K C_K/log n)^2`, uniformly on the grid.
- Majority bound `exp(-L/8) <= N^-2` is right; exact binomial tails checked
  for `L <= 400`.
- Union over fewer than `4 sqrt(N)` points gives size at most `4 N^(-3/2)`.
- Critical powers under the alternative: 2, 2, `(2K+4)/(K+1)`, and
  `(2K+2)/K` for the scan range. The last exceeds the third by
  `2/(K(K+1))`. No term of (12) fails to be beaten by `lambda^2`, logs and
  block splitting included.

**1f. Silent assumptions.** None fails. Anchor selection is independent of
estimation references and probes; anchor moments sit inside `epsilon`. Finite
probes cost `1/sqrt(n n_P)` and no more. The 3+1 calibration is
single-valued.

## B2. The analogue

- **Model:** iid events on the unit square in null coordinates, with shared
  references.
- **Moments:** `p`, `f`, `a` (order-two U-statistic over all pairs), `b`
  (anchor interval), and anchor `p`, `f`.
- **Calibration:** square-root map to `(u_hat, v_hat)`, plus
  `0.15 (2a/p^2 - 1)` added to `u_hat` so the phase depends on the
  U-statistic.
- **Functional:** `H = (u_hat-uA_hat)(v_hat-vA_hat) - b + 0.3 (2a/p^2 - 1)`,
  zero on the population.
- **Order separation:** counts are exact; `G(t) = int F(mu_0 + t eps) dX`
  separates Taylor orders by Richardson combination.
- **First-order law:** deterministic quadrature gives
  `omega sigma_1 = 0.01556, 0.01551, 0.01551, 0.01553` at
  `omega = 256..2048`, a pure `1/omega`.
- **Scan range:** the full statistic equals `sigma_1/sqrt(n)` within 7.5%.
- **Limits:** no transverse dimensions, `u` is not recovered from `c`, and
  the side `u < v` is supplied.

## B3. How the criteria were set

- Exploratory runs with other seeds came first; exponent tolerances were
  chosen then.
- `s6` first printed two FAIL lines, both defects of the checks: a maximum of
  negative numbers initialised at 0.0, and a requirement that
  `lambda^2 > t_n` at scan entry, which the draft does not claim. The second
  is now reported as F8.
- `s4` T6, the free-intercept fits, `s6` L4b and `s2` D were added after the
  first complete run. None replaced a failing check.
- An early exact-cell integrator disagreed with the grid by a factor of order
  one. The cause was the toy's mirror window `u > v`, fixed by restricting
  probes to the block.
- No tolerance was loosened to turn a failure of the reviewed claim into a
  pass.
