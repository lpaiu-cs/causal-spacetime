# Review of local_null_challenge_r2.md (argument 2 revised and its steps written out)

- **Referee:** an independent AI referee in a fresh context, 2026-11-26.
  No human has reviewed this.
- **Scope:** the revision as committed at `cfe4b8e`, with
  [checks.py](checks.py) and [checks_output.txt](checks_output.txt).
- **Brief:** the referee was given the revision, its imports, the frozen
  original, the audit's items for argument 2 and the earlier referee's
  report, as claims to test. It was not given the author's assessment. It
  was told not to read the other reviews in `independent_audit/`.
- **Limits:** only short deterministic checks. No sampling of orders and no
  statistical campaign. CPU only.
- **Archiving:**
  - The referee returned the report below as text, and the author saved it
    with light formatting.
  - Its scripts were written in [referee/](referee/). They are archived
    with a three-line header, and their outputs `rN_output.txt` beside
    them. `r4_output_v1.txt` is the output of r4's first version, whose
    criterion failed (see r4's docstring).
  - The revision was corrected after the review; its Section 9 lists the
    corrections.

---

I re-derived every numbered result, re-ran `checks.py`, and wrote five
checks of my own. The proof holds, with ten minor fixes needed.

## 1. Verdict: CORRECT WITH MINOR FIXES

**What is proved.** Fix a known `0 < q <= 1` and the exact circular vacuum
plane wave, with `N` i.i.d. points from the normalised volume of `D_lam`
and the observer given only the unlabelled order.

- **Level, every `N >= 4`.** Under the background the test rejects with
  probability at most `(3/4) N^-2`. This holds for the exact test, and for
  any implementation whose interval enclosures and derivative bounds are
  sound (they contain the true values).
- **Power.** If `lam -> 0`, `N lam^3 / log N -> infinity`, and the
  implementation's derivative bounds and enclosure widths satisfy (H-enc)
  (exact evaluation does), the rejection probability under the wave tends
  to one, uniformly in phase. Hence `TV -> 1`.

The imports ([GFL] E, F, K, K2; [WCR] Theorem 2, Lemma 3, Section 7; [QOB]
Lemma 2.1; the calibration and its certificate) are each used within their
stated hypotheses. No step needs a new idea. The repairs are:
- a constant slip in Lemma 4.2(ii);
- two missing clauses in (H-enc);
- a lambda range in Lemma 4.1;
- a bias constant about 10x larger on `K` than the text quotes;
- one false status sentence (AUDIT.md has no Section 15);
- reporting and wording fixes.

## 2. Findings

| id | location | grade | what is wrong | fix | how checked |
|---|---|---|---|---|---|
| F1 | (H-enc); Thm 4.4(a); preamble "(level)" | MINOR | (H-enc) lists only bounded derivative bounds and an `o(lam^2)` width. It omits two things the implementation needs. (i) Soundness: enclosures contain the true values and the bounds dominate the true derivatives. Without it, neither (a) nor (b) holds for the code. (ii) No spurious skip: `derivative_boxes`/`difference_interval` also return `None` when an enclosure of `p0(q)` or of `s` fails, or when `calibrated` fails, so every pair with `R` inside `K_theta` must actually be evaluated. Prop 2.5 says "needs the enclosures to be valid"; Thm 4.4(a) and the preamble state the level without that condition. | Add (i) and (ii) to (H-enc). Say (a) is for the exact test or for a sound implementation. | Read `order_wave_local_null.py` lines 104-318. |
| F2 | Lemma 4.2(ii) and its proof's last sentence | MINOR (constant slip) | A box of half-width `r_0` containing a point within `r_0` of `theta_lam(xi_X)` is within `3 r_0` of `theta_lam(xi_X)`. That puts it within `4 r_0 = r_*` of `theta_0(K)`, not inside `K_theta` (the `3 r_0`-neighbourhood), and at `r_*` it can touch the boundary of `D_F`. | Take `r_0 = r_*/5` and `K_theta` the closed `4 r_0`-neighbourhood, or restate (ii) for boxes containing `theta_lam(xi_X)`. On the coverage event `R` contains `theta_X`, and Lemma 4.3 already places `R` within `r_0` of it. Lemma 4.3 and Thm 4.4 are unaffected. | By hand. |
| F3 | Lemma 4.1 constant; Section 4 "up to about 2.6e4" | MINOR (practicality) | The audit report also gives about 4e4 at `u = -0.2`, and on `K` (`|u| <= 0.3`) the constant is about 10x larger. `|F bias| / (lam^4 log(1/lam))` is 1.37e5 (q = 0.25) and 1.42e5 (q = 1) at `lam = 1/2560`, at `(u, eta, zeta) = (-0.3, 0.6, 0.4)`. The single-log law holds, with slope about 2.4e5 to 2.6e5 per unit log. So Thm 4.4's two-point bias term falls below `g_1 lam^2` only at roughly `lam < 3e-5`. | Quote the `K`-wide constant, or say the quoted one is for interior `u`. | r4 |
| F4 | Lemma 4.1 statement "for `lam <= lam_1`" | MINOR (scope) | The `F` part uses `F` Lipschitz on a neighbourhood of `theta_0(K)`. That needs `theta_lam(xi)` inside the neighbourhood, which the proof gives only for small `lam`. No counterexample: the wave `c` stays in the window at `u = +/-0.3` for `lam = 1/10` (q = 0.25) and `lam = 1/40` (q = 1). | State `lam < lam_b`. Only small `lam` is used. | By hand; a one-off computation with r4's code. |
| F5 | Section 6, "its new Section 15 records this" | MINOR (false status sentence) | At `cfe4b8e`, AUDIT.md ends at Section 14. No Section 15 exists. | Delete the clause, or mark it as to-do. | `git show cfe4b8e:.../AUDIT.md` |
| F6 | Section 7, N8 | MINOR | "The least kappa/lam^2 is 0.0101 (q = 0.25)" sits in a sentence covering `1/512`, where the output shows 0.009844. N8's samples also miss the binding corner. The adversarial minima are 0.00557 (q = 0.25) and 0.00974 (q = 1). Those are still above the floors 0.00252 and 0.00441, and above the leading-order floor `(3/32) c_L rho_1^2`. | Report 0.0098. Add the corner case: `u_1` near 0, `tau = pi/2`, `|z_1| = rho_1`, offset +/-pi/6, worst direction of `w`. | `checks_output.txt`; r3 |
| F7 | Remark 2.6, Section 6: "skipped at every reference count" | MINOR (wording) | True on the coverage event of Lemma 2.3, where `R` contains `theta_X` and `theta_Y`, hence the limit box. Not true surely. | Add "on the coverage event". | By hand. |
| F8 | Section 1.3 | MINOR (wording) | "Sections 3 and 4 do not use [q <= 1]" is false: Lemma 4.2's first bullet uses `p0` decreasing on `[0, q]`, which is the certificate, which needs `q <= 1`. Also, "strictly decreasing, so `p0^-1` is smooth" needs `p0' != 0`. The certificate does give that (`d log p0/dz <= -1/50`). | Fix the wording. | Read null_coordinate_challenge.md Section 2. |
| F9 | Notation; Section 0 | MINOR | Audit finding 8 is only partly addressed. `f` is still both `cos^-2` and the future fraction. `c` is overloaded. `B` is now both the integer and the wave matrix `B(s,t)`, a new clash. `s` means both the past-segment length and the tip distance `s(X)`. Separately: "Z_lam within a factor 5/4" (the actual range is [3/4, 5/4]; `d_+-` are right), and Section 0 says Lemma 4.1 "proves" the uniform bound when it imports K2(ii). | Rename; fix the two sentences. | Reading. |
| F10 | The checks (coverage) | MINOR | Section 6's central claim is INFO-only in N2, never gated. No check touches Lemmas 3.4, 3.5, 4.2(ii), 4.3 or Thm 4.4. N10 measures `|bias_X - bias_Y|` on the frozen rows, not the bias bound of Lemma 4.1. N9 is background-only. | r1, r4 and r5 cover part of this. | Reading `checks.py`. |

## 3. Verdict on each numbered result

- **Lemma 1.1 (`a = p^2 c`):** correct.
- **Lemma 1.2 (nested kernels):** correct, by transitivity and antisymmetry
  (strictness holds almost surely).
- **Section 1.3 (calibration):** correct as imports, including the domain
  `D_F` and the background identity (re-derived). Wording issue F8.
- **Lemma 2.1:** correct. The Bernstein root and its upper bound were
  re-derived.
- **Lemma 2.2:** correct (Hoeffding 1963, Section 5). Verified by the exact
  law of complete U-statistics in r2.
- **Lemma 2.3 and the choice of `B`:** correct.
  `6P + 6C(P,2) <= 6N + 3N^2`, and the failure bound is `(3/4) N^-2` for
  `N >= 1`. Both choices of `B` satisfy `B <= 6 + 4 log2 N`.
- **Lemma 2.4 (propagation, (4)):** correct. Both segments lie in `R`; the
  `hull(0, box_d)` part is needed.
- **Proposition 2.5 (level):** correct for every `N >= 4`, given sound
  enclosures (F1).
- **Remark 2.6:** correct. The numbers were reproduced independently (r1).
  Wording F7.
- **Lemma 3.1:** correct. The integration by parts, the three-term split,
  `|l|` and `c_L = 16 k s_T/(1.008 M^2)` all check. (5') follows from [GFL]
  E(iv), with the right sign.
- **Lemma 3.2:** correct. The measured ratio of `A` to its bound is 0.75,
  as `M ~ f tau lam` predicts.
- **Lemma 3.3:** correct. Constants reproduced; `-1/8 + 1/4 - 1/32 = 3/32`.
  Witness points are in `K` (`W_0 >= 1/8`, `zeta_0 >= 0.44`). On
  adversarial pairs the minimum slack is about 1.1x the leading-order floor
  (r3).
- **Lemma 3.4:** correct. Density bounds `d_+-`, the Fubini volume (arcs of
  total angle `2pi/3`) and the conditional sections all check.
- **Lemma 3.5:** correct. At most `4P^3` index-sharing terms;
  `I_ij I_ji = 0`; Chebyshev bound as stated.
- **Remark 3.6:** correct.
- **Lemma 4.1:** correct for small `lam` (F4), and within K2's hypotheses.
  `K` lies in `K_(0.2)`, `s_0 = 0.2 <= 1/4`, every interval used has length
  at least 0.2, and `p`, `f`, `a` are exactly K2(i)'s forms. The constant is
  large (F3).
- **Lemma 4.2:** (i) and (iii) correct; (ii) has the slip in F2.
- **(H-enc):** the right kind of hypothesis, but incomplete (F1).
- **Lemma 4.3:** correct given (H-enc) and the fixed Lemma 4.2(ii).
  `kk >= N/5` actually holds for every `N >= 12`. The three terms over
  `lam^2` are as stated.
- **Theorem 4.4:** (a) correct for every `N >= 4`. (b) correct under
  `lam -> 0`, `N lam^3/log N -> infinity` and (H-enc) with soundness. One
  power of `log N`, entering only through the Bernstein term. Phase
  invariance follows from [WCR] Lemma 3(vi).
- **Section 0 (change map):** accurate. The "Executable derivative
  certification" text and Section 5 are verbatim from the original, apart
  from the italic notes. The original's SHA matches and it was touched only
  by `7dec69a`.
- **Audit items 1 to 5:** all addressed. Item 3 is met by allowing a range
  of `B` rather than "the least B".
- **Section 6:** correct. All six q = 0.25 rows, and the q = 1 rows with
  `r <= 512`, are skipped. AUDIT item 4 does omit q = 0.25 at `r <= 512`.
  Of the ten rows with a positive decrease, at most two can be evaluated.

## 4. The checks

I re-ran `checks.py` once with system Python 3.11: 4.5 s, 26 PASS, 0 FAIL.
The output is byte-identical to `checks_output.txt` apart from its trailing
"exit 0" line.

**Can fail on the mathematics they test:**
- N1 window values.
- N2 nestedness and the X4 reproduction.
- N4: exact binomial test of Lemma 2.1, but for means only.
- N6, on a synthetic `F`; there is no negative control.
- N7a: the residual-scaling criterion would catch an `O(lam^2)` error.
- N7b sign of (5').
- N8: exact slack, but the sampling misses the worst corner (F6).
- N9(i), (iii).

**Weak or near-vacuous:**
- N1 "fast series" only tests code consistency.
- N3 is arithmetic.
- N5(b) checks Jensen's inequality on one kernel, so it can only catch code
  errors.
- N9(ii) follows from monotonicity of `p0`.
- N2 verdicts and N10 are INFO only and cannot fail.

## 5. The referee's own scripts

All are in [referee/](referee/), each with its `rN_output.txt` next to it.

- **`r1_window_rows.py`** (1.5 s, 9/9 PASS). Computes `p0` from the [GFL]
  K1 interval kernel by quadrature, not from the 0F1 series. Gets
  `p0(1) = 0.09591918` and `p0(0.25) = 0.09903149`, and `Lambda_0 = alpha`
  to 2e-15. Reproduces the Section 6 table to 4.8e-7; ten rows skipped,
  only q = 1 at r = 1024 and 2048 evaluable.
- **`r2_confidence_bounds.py`** (2 s, 8/8 PASS). Exact law of complete
  U-statistics, 384 cases: the `r_D` failure probability is at most 0.080
  of `4e^-B`. The negative control (`kk = C(n,2)`) violates by 29.7x.
  Union-bound counting holds for `N` up to 20000. `kk >= N/5` fails only at
  `N = 6, 7, 11`.
- **`r3_witness_slack.py`** (5.5 s, 10/10 PASS). Independent canonical map;
  adversarial witness pairs down to `lam = 1/32768`. All pairs related and
  in `K`; slack is preserved to 4.5e-16. Minimum `kappa/lam^2` is 0.00557
  (q = 0.25) and 0.00974 (q = 1). `|kappa - leading| / lam^3` stays at most
  5.5.
- **`r4_bias_on_K.py`** (7 s, 6/6 PASS). Exact interval volumes, wave
  against background, over `K`. Bias is `lam^4 (A + B log(1/lam))` with
  stable slopes. Constants as in F3. Its first version's criterion failed
  for `F` at both q (factor 2.04, pre-asymptotic). That output is kept as
  `r4_output_v1.txt`, and the change is explained in the docstring.
- **`r5_n9_exact_counts.py`** (information; imports `checks.py` to reuse
  N9's seeded draws). Exact counts: q = 1 at `2^-17` is 2000/2000, so
  Section 6's claim holds. q = 0.25 at `2^-17` is 1994/2000.

## 6. The step I trust least

Lemma 4.1 feeding Theorem 4.4. I took [GFL] K2(ii) as given. The theorem
then subtracts the bias bound at the two points separately
(`2 C_b lam^4 log(1/lam)`) rather than bounding the difference of biases.
With the measured `C_b` of about 2.5e5 on `K`, that term is larger than the
witness signal `g_1 lam^2` until `lam` is about 3e-5. My numerics confirm
the single-log law only down to `lam = 1/2560` at one `u`. The asymptotic
logic is sound; the gap between "eventually" and anything checkable is
widest here.

## 7. What I did not verify

- The proofs of [GFL] E, F, K, K2, [WCR] Theorem 2 and Lemma 3, and [QOB]
  Lemma 2.1; I used their statements only.
- The `p0` certificate itself; I checked monotonicity on a grid only.
- The implementation's rational-interval soundness, rounding, `bits`
  precision, and absence of spurious skips.
- Lemma 3.4's volume constants and Lemma 3.5 numerically; both were checked
  by hand only.
- Phase uniformity numerically; I used phase 0.37 only, plus the analytic
  argument.
- Geometry beyond `lam = 1/32768` and bias beyond `1/2560`.
- `local_null_results.md`, Section 5's provenance, and the audit referee's
  scripts. I relied on its report's quoted output line for "checked only
  r = 1024, 2048".
