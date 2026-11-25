# Review of argument 3, revision 2 (`gram_fourier_noise_r2.md`)

2026-11-26. **AI-level review; no human review.** An AI referee in a fresh
context, which had not seen this work, was given the revision as committed
at `c91f3ea` and these sources:
- the frozen original;
- `gram_fourier_challenge.md`;
- AUDIT.md;
- both audit referee reports on argument 3;
- the anchor-guard certificate;
- [checks.py](checks.py).

It was told to find what is wrong, and to run only short deterministic
checks: no sampling, nothing long. Its own checks are archived as
[ref_checks.py](ref_checks.py), verbatim apart from a three-line header.
The report below is the referee's, lightly formatted. The revision was
corrected after it; its Section 9 lists the corrections.

## Overall verdict: CORRECT WITH MINOR FIXES

`A_reg` (1a) is one well-defined object, and (5), (10) and (11) hold for
it. Three things need fixing:
- one statement is literally false (E1): harmless, with a one-line fix;
- Section 0 does not list every change (M1);
- the Section 8 headline overstates what is now written (item 10).

## Item verdicts

1. **`A_reg` well defined: yes, apart from E1.**
   - **Smoothness near the support.** Everything `H` uses is smooth on a
     neighbourhood of the support:
     - the square roots, since all six `b >= b_flat/2`;
     - `alpha(gap)`, since gaps are at least 0.096 or 0.248;
     - `G_A^-1`, since `det >= 0.00638^2 - (2e-5)^2`;
     - `Sigma`.
   - **The calibration.** `p0` is strictly decreasing, with `p0' != 0` on
     the anchor and probe s-ranges, for q in {1e-4, 0.01, 0.25, 1}.
   - The gap bounds 0.096 and 0.248 follow from the cutoffs.
   - **`supp w` lies in `Sigma(D_Sigma)`.**
     - For q in {1e-6, 0.25, 1}, probes have `eta >= 0.6167`,
       `zeta >= 0.1500` and `eta + zeta <= 0.995`.
     - Anchor states on the outer box have `eta >= 0.126`,
       `zeta >= 0.638` and `eta + zeta <= 0.985`.
   - **The anchor numbers.** `iota` is 12 functions of the 12 anchor raw
     moments, and the preimage of the closed box is compact inside its
     domain.
   - **The cutoffs** reach all 18 coordinates, but only from below in
     `b_iX`.
   - **`alpha`.** `alpha(L) = (1 - x cot x)/(2q)` with
     `x = sqrt(q) L <= 1 < pi`, so `alpha` is increasing. Checked against
     `wavegeom` to 1e-7.
   - **The phase `u~` is enough:** `F_omega` does not depend on the choice
     of `u~`, because `A_reg = 0` wherever `u~` may differ from `Sigma_u`.
2. **(R2): correct.**
   - `b_AX > 0` exactly when A is chronologically before X. A null pair
     has `b = 0`, so the word should be "chronologically related" (M5).
   - The identity also holds for null pairs; only spacelike pairs break
     it, so nothing turns on strict versus null.
   - Only facts available on the support are used. "No good-anchor
     condition" is right.
3. **(R3): correct for the accepted triples,** with one step only cited
   (M6).
4. **Section 3, (5): correct.**
   - For `K = 1` there is no `k = 2` term; `1/n` is absorbed by the
     remainder `omega/n` instead (M4).
   - The elementary-check variance is exactly right. Derived, including
     the diagonal `1/(6nm)`, and checked by quadrature.
5. **Section 5: the identity T1 + T2 + T3 + T4 is exact, and the four
   bounds are right.**
   - "Globally" is false; the bounds hold on any bounded set (E1).
   - T2 also needs `omega ~ 1/lambda`; the text says so only in T3.
   - The out-of-active-set remark is correct.
6. **Section 6: correct.**
   - `a_2x a_2y = G12 sqrt(G22/G11 - G12^2/G11^2)`. Its supremum on the
     box is 2.00625e-5, reached at the corner, since the product is
     monotone in `|G12|`.
   - The slack identity is exact. It reproduces `0.0289166 =
     0.096 - 0.01288/0.192` exactly, using `dm >= du`. It needs only the
     box, so it is not circular.
   - "The acceptance event is not needed under the null" is correct.
   - There is one numerical inconsistency (M3).
7. **Section 7: correct.**
   - The null vote bound is `(K C_K / log n)^2` on the scan range.
   - A block with no accepted triple votes zero, so no acceptance event is
     needed.
   - `n_P >= n` holds by construction.
   - Hoeffding gives `exp(-L/8) <= N^-2`, and the union over fewer than
     `4 sqrt(N)` points vanishes.
   - The scan-cap remark is right: a grid up to `sqrt(n)` gives
     `alpha > 2 + 2/(K+1)`.
8. **Section 4: correct.**
   - `|beta| <= C lambda min(du, lambda)` follows from `|f| <= f_max` and
     one integration by parts. The sup ratio is 5.06, 4.95 and 4.84 at
     r = 64, 256 and 1024.
   - `M_left = Mt(-1/2, u)` and `M_right = Mt(u, 1/2)`, because
     `m(-1/2) = 0`.
   - The eigenvalues are `±|beta|`.
9. **Completeness: AUDIT items 1 to 7 are all addressed.** What remains:
   - practical notes (3a F8 and F9; 3b F4 and F5), acknowledged in general
     terms;
   - 3b F6, which belongs to `gram_fourier_results.md`;
   - Section 0, which is incomplete (M1).
10. **Headline: not fair as worded** (see below).

## Errors and gaps

**E1: an error as stated; local; nothing downstream changes.**
- (R1) says `A_reg` is in `C_c^infinity(R^18)` with derivatives of
  `F_omega` bounded by `C_j (1 + omega)^j`. Section 5 cites this as
  "globally".
- But `chi_b(b_iX)` cuts off only small b. `H` depends on `sqrt(b_iX)`, so
  `A_reg` and its probe-direction derivatives grow without bound as
  `b_iX -> infinity`.
- At the frozen anchors with a background probe, `H = 2.7e-17` at the
  population value `b_3X = 0.0127`. It is 0.93 at `b_3X = 1`, 3.18 at 10
  and 10.3 at 100.
- **Corrected statement:** "`A_reg` is smooth on `R^18`. For every bounded
  `Q`, `supp A_reg ∩ Q` is compact and all derivatives are bounded on `Q`.
  Take `Q = [-1, 2]^18`." Every point the proof evaluates lies in `Q`.
  Alternatively, multiply by a smooth upper cutoff in each `b_iX`.

**M1. Section 0 is incomplete.** The diff shows changes that Section 0 does
not list:
- the whole preamble, which removed the status sentences;
- the new sentence after (1);
- the parenthetical on `B(u)` in Section 2;
- the reworded simplified-bound sentence in Section 3;
- in Section 6:
  - "outer good set" became "three-quarter set, where `chi_anc = 1`" (a
    substantive and correct strengthening);
  - the left side of (11) and "together with (R3)" were added;
- in Section 7:
  - "and regularized calibration being realized" was dropped;
  - the guard paragraph was rewritten, deleting "It must retain a failure
    without enlarging the neighbourhood or weakening its target
    inequalities".

Otherwise unmarked text is word-identical. The SHA-256 of the original
matches the certificate's recorded `f932f3f6...`, and `f1072f0` is
reachable from HEAD.

**M2.** `psi` and `chi_b` are existence choices, while (R4) calls them
fixed. Name one concrete construction, for example the `exp(-1/x)` smooth
step.

**M3.** `|ell_2| < 1.02` gives 0.0106972; 0.0106978 needs the
certificate's 1.009434. Both clear 0.01, and the grid coefficient is at
least 0.0101673 either way.

**M4.** "Absorbed into the `k = 2` term" should read "into the `k = 2`
term (`K >= 2`) or the remainder `omega/n` (`K = 1`)".

**M5.** In (R2), "chronologically related, in that order".

**M6.** The (R3) second bullet and the Section 5 out-of-active remark.
- "For every probe" needs `sup_X |mu_lam(X) - mu_0(X)| -> 0`, tips
  included, and anchors in a compact interior set.
- The anchor hypothesis is unstated: Section 5 should say it conditions on
  an accepted triple in the three-quarter set.
- The tip closeness is cited from the 3b report, not derived. (8) is stated
  only for active X.

**M7.** Notation clashes:
- `c` means both `a/p^2` and `sqrt(1 + 4 q lambda^2)`;
- `alpha` means both the sampling exponent and the interval coefficient;
- `K` means both the Taylor order and the kernel in (8).

## What was checked numerically (deterministic)

- `checks.py`: 11 PASS. Its `Im D_A` check rebuilds the first term from
  the certificate's own output, so it is not independent. Its elementary
  check covers only the cross term.
- `a_2x a_2y`:
  - supremum 2.00625e-5, monotone in `|G12|`;
  - the corrected bound is 0.0106978 with 1.0094, and 0.0106972 with 1.02;
  - the grid coefficient is 0.0101679.
- `alpha`: closed form and monotonicity at q = 1e-6, 0.01, 0.25 and 1.
- `b_low`: 7.18e-6, 7.58e-6 and 7.70e-6. The probe-to-anchor ratio is 80.1
  to 80.4.
- The anchor-pair slack, 0.0289166, as an exact rational.
- Probe and anchor membership in `D_Sigma`.
- The calibration `p0` is strictly decreasing.
- The growth of `H` behind E1.
- The bound on `|beta|`.
- The full elementary variance at n = 1, m = 2, by exact quadrature, to
  1e-10 at three frequencies.
- Scan-range domination and the `K = 1` absorption.

## Least-trusted remaining step

The coupling uniform in X.
- T2 and T3 sit exactly at `lambda/sqrt(n)` by (8), with no slack against
  the exponent 2.
- The revision's "every probe" extension rests on a one-line citation.
- (8) at the tips is still the original sketch; its proof is in the 3b
  report, Section 3.2.

## Item 10: is the headline fair?

Not as worded.

**Now fully written:** the statistic; (R1) to (R4); the Section 5 chain;
the corrected constants.

**Still outline or citation:**
- the index-counting proof of (5), complete as an outline, with its
  routine pieces in the 3a report;
- the tip estimates (6) to (9), re-derived in the 3b report;
- (11), which rests on the challenge's predicted limit (derived in 3b) and
  on two assertions: that anchor-phase terms vanish, and that calibration
  errors are `o(lambda^2)`;
- the geometric constants, which come from the author-code certificate.
  One of its steps was wrong; it is corrected by hand, not re-certified;
- the U-statistic moment and Hoeffding bounds, cited as standard.

**The bracket.** "The bracket becomes `4/3 <= alpha_low <= alpha_high <= 2`"
is not new: AUDIT.md Section 2 already states it. "Both ends written out"
is inaccurate:
- `quantized_order_bound.md` has not been revised since the audit;
- AUDIT argument 1, items 1 to 6, are still open in its text. They include
  the constant link, without which the claim fails at L = 1;
- its strict-separation sentence relies on the unaudited
  `coordinate_oracle_barrier.md`.

A fair wording: "Argument 3's regularised statistic and wave-side chain are
now written in full; the remaining steps are at outline level in the text,
and the referee reports supply them."

## Not checked

- the re-derivation of (6) to (9) at the tips;
- the Section 3 index count;
- the population expansion behind (11);
- the certificate's interval arithmetic, beyond the `Im D_A` term and the
  pair slack;
- values of q in (0, 1] beyond those sampled;
- the lower-bound document.

No simulation was run.
