# Review of wave_causal_relation.md (the wave's and the background's causal relations written out)

- **Referee:** an independent AI referee in a fresh context, 2026-11-26.
  No human has reviewed this.
- **Scope:** the note as committed at `1ba2e6e`, with [checks.py](checks.py)
  and both runs ([checks_output_run1.txt](checks_output_run1.txt),
  [checks_output.txt](checks_output.txt)).
- **Brief:** the referee was given the code and the documents, not the
  author's conclusions.
- **Archiving:** the referee returned the report below as text. The author
  saved it unchanged, with two exceptions in its Section 7. The scripts'
  folder now points to [referee/](referee/). The re-run row is marked "not
  archived".
- **The referee's scripts** are archived there with a three-line header,
  and its output `r2_output.txt` beside them.
- **Not archived:** the referee's re-run of `checks.py`. It matches
  `checks_output.txt` line for line.

---

# Referee report: `wave_causal_relation.md` and `wave_relation/checks.py` at `1ba2e6e`

## 1. Verdict for the mathematics: CORRECT WITH MINOR FIXES

I re-derived every step and found no false step. The findings are unproved side claims, an incomplete sufficient condition in a remark the note does not use, an incomplete justification sentence, and notation clashes.

- **Lemma 1.** Correct. It is consistent with the tips: under `T`, `p_- <= p_+` because `2h >= 0`. Under `-T` the relation reverses and the diamond is empty.
- **Theorem 2, both directions.**
  - The class of causal curves is adequate. "Only if" is proved for AC curves that may have zero velocity, and "if" produces a C^1 curve with nonzero velocity. So `J^+` comes out the same for piecewise-C^1, Lipschitz and AC curves.
  - Absolute continuity of F holds, as does the chain rule. `M(u1, .)` is C^1 on all of `I`, so the classical chain rule applies wherever `u` is differentiable, including where `u(tau)` sits at `u1` or `u2`.
  - The a.e. bound holds, including where `u' = 0`, and `p = M^-1 Delta x` is the maximiser.
  - The explicit curve is correct, and it is timelike when `eps > 0`.
- **Remark 2.1, Lemmas 3, 4 and 5.** Correct.
- **Lemma 5, the conformal argument.** It is complete. Φ is a diffeomorphism onto ℝ⁴, the conformal factor is positive, and `dT(d_u + d_v) > 0`. With the standard Minkowski relation for AC curves, the relation follows.

## 2. Verdict for the script and runs: CORRECT WITH MINOR FIXES

- I reran `checks.py` at HEAD and its output is byte-identical to `checks_output.txt`. The einsum assembly is correct.
- Fixing the step rule after run 1 was legitimate:
  - no tolerance was changed;
  - the rejected `d/100` also passes;
  - `d/400` converges further.
- Problems:
  - The note's account of why run 1 failed is wrong in its mechanism.
  - Two numbers in the Results paragraph are misquoted.
  - Four more checks cannot fail on the mathematics, and the note does not disclose this.
  - The script contains no check of the "only if" direction for the wave that could fail.

## 3. Findings

| id | location | grade | what is wrong or missing | how checked |
|---|---|---|---|---|
| F1 | §1, "for sampled events J+ and the chronological relation agree almost surely" | MINOR | Asserted, not proved. It needs two facts:<br>(i) points with strict inequality in (a) lie in `I^+`, because the eps>0 curve has `g(c',c') = -2eps/Delta u < 0`, or by the standard `I^+ = int J^+`;<br>(ii) `∂J^+(X)` is null: a graph `v = v(u, x)` for `u > u1`, plus the ray at `u = u1`.<br>It also needs the pair law to be absolutely continuous. No user document relies on `I^+` (all use `>=`), so the conclusion is safe. The sentence still needs its two lines. | Derived by hand. |
| F2 | §2 Lemma 1, §7 | MINOR | The time orientation is "declared". The note never says that S has exactly two orientations, or that the tips force `T`: under `-T`, `{p_- <= X <= p_+}` is empty. Section 7 does not list it among the definitions. | Theorem 2 applied to the tips. |
| F3 | Remark 2.2 | MINOR | "u extends to a non-decreasing function" is not sufficient for an extension. One also needs the extension's slab `u^-1([u1,u2])` to equal S's slab. This holds for the Brinkmann form, since Rosen→Brinkmann maps S onto `{U in I}`. Disclaimed and unused. | Argument. |
| F4 | Remark 2.1 | MINOR | `∂J^+(X)` also contains the ray `{u = u1, x = x1, v >= v1}`, a `d_v` null geodesic. The displayed equation omits it. | Closure computation. |
| F5 | Lemma 4 proof, "every product … contains one of these zero factors" | MINOR | The named zero factors are `Gamma^a_vb` and `Gamma^a_uu`. Some terms vanish only because every `Gamma^u` is zero, for example `Gamma^k_ju Gamma^u_ki` and `d_u Gamma^u_ij` in `R_ij`. The conclusion is right; the stated reason is incomplete. | Enumerated all terms. R1 confirms the result. |
| F6 | notation | MINOR | Several symbols carry two or three meanings:<br>- `T`: vector field (§1), matrix (§4), coordinate (§6);<br>- `Gamma`: Christoffel symbols and the unimodular matrix (Lemma 4);<br>- `c`: the curve and `cosh 2beta`;<br>- `Q`: here `Gamma^-1 Gamma'`, but the scalar `tr Q^2/8` in design.md. | Reading. |
| F7 | Lemma 5 / §7 vs r2 §1 | MINOR | Two gaps:<br>- r2 also asserts the density is positive, bounded and C^1 with bounded `d rho/dT` near the closed diamond. That is trivial (`d rho/dT ∝ 2√2 a^5 a'`) but not written.<br>- The note never says the coupling identifies the background's transverse coordinate with the wave's `y`. That identification is what makes Lemma 5, stated in `x`, the same thing as r2's (2), stated in `y`. | Comparison of r2 §1–2, one_point_strip_bound §1 and coordinate_oracle_barrier §1. |
| F8 | §8 and preamble, run-1 cause ("stencil reached u = 0.55, past the pole") | MINOR | The mechanism is wrong:<br>- `gamma = cos^2(√q u) Gamma` is entire, so sampling the metric at 0.55 is harmless.<br>- `gamma^-1` is evaluated only at Christoffel centres, which are `<= 0.50 < 0.5236`.<br>- Applying the same outer stencil to the exact Christoffel symbols reproduces the failure: 1.058 against 1.061.<br>The failure is truncation error of the outer derivative. Its reach `2h = 0.05` against a distance to the pole of 0.074. It has nothing to do with crossing the pole. | R3 |
| F9 | §8, "halving the step cut the error about 16-fold" | MINOR | At the failing point the successive halvings give 28.8×, then 17.9×, then 16.4×. 16× holds only asymptotically. | R3 |
| F10 | §8, step-rule choice | MINOR | The note does not record that `d/100` also passes: vacuum 3.8e-7, control 3.5e-7, background 8.2e-7 against 1e-6. The trial outputs were not kept. Recording these would show the fix was not tuned to pass. | R3 |
| F11 | §8 Results | MINOR | Two misquotes:<br>- "W1 below 2.4e-8", but the output is 2.44e-8.<br>- "W2 to 1e-13 or better", but the eigenvalue check gives 1.5e-13.<br>Also "seeded draws only test inequalities": they test identities too. | Output files. |
| F12 | §8, "What the checks do not test" | MINOR | Four checks cannot fail on the mathematics, and the note lists none of them:<br>- **W4 first part.** The ratio equals `1 + ∫δ'ᵀγδ'/act` exactly. The cross term `2pᵀ∫δ'` is at most 1.3e-13. It is Cauchy–Schwarz with the remainder built in, so the "1.0106" only reflects the perturbation size. It can fail only through a code inconsistency.<br>- **W2 "M = c m I − s J".** Both sides use the same quadrature nodes, so it is W2-inv summed with weights.<br>- **W3 "u stays affine".** `Gamma^u` comes from finite differences of constant entries.<br>- **W4 Lemma-1 sign test.** A tautology: about 11% of draws are in the tested class, and the case `w^u = 0` is never tested.<br>So no check in the script that could fail tests the wave's "only if" direction. | R4 |
| F13 | §7 | MINOR | Two omissions:<br>- **Background facts used but not listed:** the chain rule for C^1∘AC; AC preserved under a diffeomorphism; Minkowski's `J^+` for AC curves; `I^+ = int J^+` (behind F1).<br>- **Definitions not listed:** the time orientation (F2), and that the model's order is `J^+`. | Reading. |

## 4. My own numbers

**R1 — exact derivatives, explicit loops, full Riemann tensor, sign calibrated on S² (Ricci − g = 2e-16):**
- Circular family, 400 random (q, λ, φ, β, u) including β ≠ vacuum: `R_uu` matches `2q − s²/(2λ²)` to 1.5e-12 relative. All other components are exactly 0.
- Lemma 4's Christoffel symbols match the generic ones to 8e-15.
- Vacuum family: max |Ric|/(4q+1) = 3.3e-12.
- General Rosen γ(u), random and polynomial: `R_uu` matches `−½ tr P' − ¼ tr P²` to 8e-16, and matches design.md's form to 8e-16. Other components are 0.
- `tr Q²` equals `2s²/λ²` to 2.6e-13.

**R2 — direct test of the relation:**
- (A) 600 random piecewise-linear causal curves, each segment certified causal by a Lipschitz bound. They include 1007 segments that pause in u and 373 with zero velocity, a non-vacuum β = 1.3, and the background. None violates Theorem 2, none has u decreasing, and the minimum of `2Δv/(ΔxᵀM⁻¹Δx)` is 1.0020.
- (B) Minimum of `2Δv` over piecewise-linear curves with K = 1, 2, 4, 8, 16 segments, divided by the bound. The sequence is always ≥ 1 and decreases monotonically: 1.576 → 1.050 at (q=1, λ=1/4); 7.34 → 1.29 at (q=9, λ=1/8); 78.7 → 5.71 at (q=4, λ=1/16, β=1.3). Convergence is slow at high frequency, so I rely on the explicit curve for sharpness.
- (C) Lemma 3 by independent quadratures: M by adaptive `quad_vec` of the numerical inverse, and J and G by QUADPACK QAWO. M matches `c m I − s J` to 2.6e-15, the eigenvalues match `±t|G|` to 2.9e-15, and the trace of E is 4.7e-15.

**R3 — the W1 grid under several step rules (vacuum / control / background maxima):**

| step rule | vacuum | control | background |
|---|---|---|---|
| `λ/40` (run 1) | 1.06 | 0.967 | 2.8e-6 |
| `λ/80` only | 3.7e-2 | 3.4e-2 | — |
| `min(λ/80, d/100)` | 3.8e-7 | 3.5e-7 | 8.2e-7 |
| `min(λ/80, d/200)` (run 2) | 2.44e-8 | 2.20e-8 | 5.1e-8 |
| `min(λ/80, d/400)` | 1.3e-9 | 1.5e-9 | 3.3e-9 |

**R5 — checks.py's einsum Ricci on a generic 4-coordinate Lorentzian metric:** it matches the exact loops to 2.5e-10.

**Rerun:** `checks.py` at HEAD reproduces `checks_output.txt` byte for byte. Lines 4–21 are identical between run 1 and run 2, as expected, because W1 draws no random numbers.

## 5. The step I trust least

Section 1's sentence identifying the model's order with `J^+` and asserting almost-sure agreement with the chronological relation (F1). It is correct, but no proof is written, and the model documents never state which relation they use. A close second is the unstated coupling identification of the background's `x` with the wave's `y` (F7). It is the only link between Lemma 5 and r2's (2) in the `y` chart.

## 6. What I did not check

- Whether the model's code anywhere samples orders using `I^+` or `J^+`. The experiments I read use interval kernels, not pairwise relations.
- Any extension beyond the Rosen chart, beyond the Brinkmann argument in F3.
- The rest of `one_point_strip_bound.md` and `quantized_order_bound_r2.md` beyond the imports.
- `check_a_exact_family.py`.
- W3's accuracy, which I did not reproduce independently. Its conclusion is covered by R2(A) and (B) and by the analytic geodesic derivation.
- Whether the convergence in R2(B) reaches 1 at high frequency. The sequences are still above 1 at K = 16.

## 7. Scripts I wrote

Archived in [referee/](referee/):
- `r1_ricci_jets.py`
- `r2_relation_pl_curves.py` (output in `r2_output.txt`)
- `r3_run1_diagnosis.py`, which copies checks.py's finite-difference code to audit it
- `r4_vacuity.py`
- `r5_einsum_generic.py`
- `rerun_checks_output.txt`, the rerun of the author's script (not archived; see the header)

The repository was not modified.

---

## How each finding was handled

All changes are in [wave_causal_relation.md](../../wave_causal_relation.md);
its Section 9 has the full table.
- **Unchanged:** no statement of a theorem or lemma.
- **`checks.py`:** only a correction line was added to its docstring (F8).
  Its outputs are unaffected.

| finding | handling |
|---|---|
| F1 | applied: Section 1 writes out the explicit timelike curve, the measure-zero boundary and the absolutely continuous pair law |
| F2 | applied: Remark 1.1 |
| F3 | applied: Remark 2.2 states the slab condition |
| F4 | applied: Remark 2.1 adds the ray at `u = u_1` |
| F5 | applied: Lemma 4's proof lists every zero factor |
| F6 | applied in part: the vector field is now `Xi`, the curve `C` and the unimodular matrix `Sigma`. `T` (matrix in Section 4, coordinate in Section 6) and `Q` are flagged in Section 9 rather than renamed |
| F7 | applied: Remark 5.1 |
| F8 | applied: Section 8 and the preamble state the truncation mechanism. The script's docstring has a correction line. Commit `3d3f9ca`'s message cannot be changed and is corrected here |
| F9, F10, F11 | applied: Section 8 |
| F12 | applied: Section 8 names every check that cannot fail on the mathematics, and credits R2(A) for the "only if" direction |
| F13 | applied: Section 7 |
