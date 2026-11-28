# Review of gram_fourier_lemmas.md (argument 3's outline steps written as proofs)

- **Referee:** an independent AI referee in a fresh context, 2026-11-26.
  No human has reviewed this.
- **Scope:** the note as committed at `6257276`, with
  [checks.py](checks.py) and [checks_output.txt](checks_output.txt).
- **Brief:** the referee was given the note, its sources and the checks, as
  claims to test. It was not given the author's assessment, and it was told
  to ignore status sentences in the sources.
- **Limits:** only short deterministic checks. No sampling of orders and no
  Monte Carlo of the statistic. CPU only.
- **Archiving:**
  - The referee returned the report below as text, and the author saved it
    with light formatting.
  - Its scripts were written outside the repository. They are archived in
    [referee/](referee/) with a three-line header, and their outputs
    `rN_output.txt` beside them.
  - The note was corrected after the review; its Section 9 lists the
    corrections.

---

## Overall verdict: CORRECT WITH MINOR FIXES

I re-derived every numbered result and checked the delicate ones with my
own code. Every result the note relies on is true:
- the noise bound (Proposition D), including uniformity over anchor
  triples, degree at most `k-1` in `omega`, and the lone-index count;
- the canonical-map bounds at the tips with the stated constants (Lemmas E
  and F);
- the coupling and the closeness of moments (Lemmas G to J);
- the exact volumes (Lemma K) and the `lambda^4 log(1/lambda)` comparison
  (Corollary K2);
- the signal (Proposition L);
- the anchor selection (Lemma M).

Two items are graded above MINOR:
- **F1, a GAP.** Proposition L, and through Section 5 also (R3), applies
  Lemma I(d) to the anchors before showing that the anchors satisfy its
  hypothesis. As written this is circular. It closes in two lines with the
  note's own Lemma I(e).
- **F2, false as stated.** The `b` formula is missing a positive part
  `(kappa)_+`. Nothing downstream changes.

Everything else is wording, scope or a constant. The Section 0 map covers
every item in r2 Section 8. There are small overclaims in the note's
account of itself (F7 to F11).

## Findings

| id | location | grade | what is wrong | how checked | fix |
|---|---|---|---|---|---|
| F1 | Prop. L, bullet "anchors and probes are in `K_(0.14)`"; Section 5 "How (10) now stands", (R3) | GAP | The bullet uses Lemma I(d) to get `iota^0 = iota_lam + O(lambda^2)`, but I(d) assumes the anchors are already in `K_(s_0)`, which is what the bullet is proving. (R3) bullet 1 has the same gap: it moves from "wave state in `supp w`" to "canonical point in `N`" with no argument. | re-derivation | Lemma I(e) holds at every point, so each anchor's `(p, f, a)` differs between wave and canonical background by at most `C lambda^1.5`. The wave states lie in the compact set `Sigma^-1(outer box)` inside `D_Sigma`, and `Sigma` is Lipschitz on a neighbourhood of it. Background calibration is exact at canonical points. Hence `|u_i - u_i*| <= 0.002 + C lambda^1.5`, the anchors lie in `K_(0.14)`, and then I(d) applies. The same argument handles probes in (R3). |
| F2 | Cor. K2(i), last bullet; Remark K1 (`b_0`); K2(iii) (`eta^(0)`) | ERROR as stated; no downstream effect | `b = kappa^2 Lambda(u_A, u_X) / Lambda(-1/2, 1/2)` is false when `kappa < 0`. For a spacelike pair with `u_X - u_A >= l_0`, `b = 0` but the right side is positive. Lemma K needs `P_1 <= P_2`. | direct | Write `(kappa)_+^2`, or assume `xi_A <= xi_X`. Every use has slack at least 0.1. |
| F3 | Cor. K2(iii) | MINOR | The claim is made uniformly on `K_(s_0)`, but the proof only says "`Sigma` is Lipschitz near the compact probe region". `K_(s_0)` restricts only `u`, so it contains points with `eta = 0` or `zeta = 0`, where `Sigma` is undefined. The conclusion is still true where `p, f > 0`, by a multiplicative argument (`c_lam` depends only on `u`, and `eta_hat = eta (1 + O(lambda^4 log))`), which the note does not give. | re-derivation | Restrict to the active set and the anchor outer set (all that is used), or write the multiplicative argument. |
| F4 | Lemma M, pilot accuracy | MINOR | A triple is accepted "with a labelling". To certify whichever labelling is accepted, `b` must be accurate for all 6 ordered pairs: 9 + 6 = 15 moments, so the factor is 30, not 24. | r7 | Use 30 (or 6 x 24). The limits are unchanged. |
| F5 | Lemma H(1), proof | MINOR | "`gt'` and `gt` lie on opposite sides of 0" leaves out the case where the orientation flips. That case is impossible off `u_X = u_R`, because `Psi` preserves `u`, but the note does not say so. Case 1 also needs the step `|gt| <= |gt'| + sqrt 2 (d_X + d_R)`. | r6: no flips; max ratio 0.287 | Add both sentences. |
| F6 | Lemma I(a), via H(2) | MINOR | H(2) needs `lambda <= 2 l_0 = s_0`. Lemma I is stated for `lambda <= lambda_1`, and `lambda_1 > s_0` for small `q` (`lambda_1 = 0.62` at `q = 0.01`). The statement survives trivially with `C >= s_0^-2`. | constants | State `lambda <= min(lambda_1, s_0)`. |
| F7 | Section 0 map | MINOR (scope) | The row "closeness for every probe -> I(c), (d)" should read I(c), (e) plus Remark I1. Lemma J proves (9) only for `p = 2, 4` and only with `X` and the anchors in `K_(s_0)`, while r2's (9) says "every `p >= 2`". The cases r2 actually uses are covered. | comparison | Say so in the map. |
| F8 | Preamble | MINOR (overclaim) | "Reviewed by a fresh-context AI referee (Section 9)", while Section 9 reads "Pending". | reading | Change to "to be reviewed". |
| F9 | Section 8, G5 row | MINOR | G5's "TV" is the L1 upper bound `(pi/Z_0) int a^2 W_0 |det T - 1|`, not TV. The exact `TV/(lambda^4 log(1/lambda))` at `lambda = 1/512` is 8.07 (`q = 0.25`) and 30.9 (`q = 1`), not 10.7 and 41. | r2 | Relabel. |
| F10 | Section 8, "agree with the author's certificate" | MINOR | `a_3y` differs (0.0793745 against the certificate's 0.0793741). `Im D_A` and the grid coefficient also differ, because they use different inputs. Only `ell` and 0.0106978 are discussed. | r5 | Say which values agree and why the others differ. |
| F11 | Section 0, "(8), uniform in the fixed endpoint" | MINOR (wording) | It is uniform only over endpoints in `K_(s_0)`. For an arbitrary endpoint the note proves `lambda^1.5` (I(c)). | reading | Add "in `K_(s_0)`". |
| F12 | Prop. D, Step 4 | MINOR | Uniformity over anchor triples is argued from compactness alone. It also needs joint continuity in (anchors, `x`) of `mu_0` and `d_u mu_0`, and a uniform lower bound on `u_x - u_A` (from `b_0 >= b_flat/2` and `kappa_0 <= 1`). Both are routine. | re-derivation | Add one sentence. |

## Item verdicts

- **Lemma A.** True; the proof is correct. The Bell-number count and the
  degenerate-kernel count both hold.
- **Lemma B.** True; correct.
- **Lemma C.** True; correct (Stieltjes integration by parts, no boundary
  terms).
- **Proposition D.** True, and the proof proves it.
  - **Degree in `omega`.** The `omega^k` coefficient of `d^alpha F` is
    `A (-i)^k prod d u~`, which is zero at `mu_0` by (R2). Derivatives in
    anchor directions cannot raise the degree.
  - **Lone-index count.** Each surviving index occurs at least twice in
    `2L` slots, so there are at most `L` distinct indices. The
    normalisation is at most `2^L n^-L`, giving
    `C n^(-k/2) omega^(k-2)`.
  - **Remainder.** The `J(t)` factorisation gives
    `omega^K |eps|^(K+1)`.
  - **Uniformity.** `T_k = 0` identically unless `mu_0,anc` is in
    `supp chi_anc` (F12 applies).
- **Lemma E.** True with the stated constants, tips included.
  - Most constants are generous: the proof actually gives
    `|det T - 1| <= 1.29 e^2`, and `|delta v| <= 6 C_E lambda^2 / M`
    holds without `f_max`.
  - The one tight bound is `eps <= e(s)`, which reaches 0.993 of `e(s)` at
    `q = 0.01`. It comes from QOB Lemma 2.1.
- **Lemma F.** True; correct. The logarithm is genuine.
- **Lemma G.** True; correct.
- **Lemma H.** True. The proof is correct up to F5 and F6.
- **Lemma I.** True.
  - The tip-zone argument in (a) is correct: the window `w(s_R)` does not
    depend on `v~_R`, `kappa_0` has slope `+/- 1` in `v~_R`, and the mass
    is `int s lambda^2 (1 + s^(-1/2)) ds = O(lambda^2)`.
  - (b) to (e) are correct (F6 applies).
- **Remark I1.** True; correct.
- **Lemma J.** True for `p = 2, 4` (F7).
- **Lemma K.** True; correct.
- **Remark K1.** True, apart from F2.
- **Corollary K2.**
  - (i) True for related pairs (F2).
  - (ii) True with the stated constant, including intervals shorter than
    `lambda`. The `lambda^4 log` rate is attained.
  - (iii) True where it is used (F3).
- **Proposition L.** True. The proof is correct up to F1.
  - The expansion, the angular average and the Fourier step are right.
  - Every dropped term is of the stated order: `O(lambda^6 log)`,
    `O(lambda^5 log)`, `O(lambda^4)`, and an `O(lambda^3)` remainder.
  - The `delta ell` term vanishes exactly under the angular average.
  - `c_q` is valid. The actual coefficient is 0.0184, against the certified
    floor 0.0102.
- **Lemma M.** Valid even though the pilot references are shared. The union
  bound needs no independence, and the triples' positions are independent
  (F4 applies).

**Section 0 map against r2 Section 8.** Every outline step is covered:
- (5);
- (6) to (9);
- closeness of moments for every probe;
- the expansion behind (11);
- the guard constants;
- the U-statistic bounds.

None is left unproved, apart from the text gap F1 and the scope note F7.

## What I checked numerically

All scripts are deterministic and each runs in under 35 s. They are
archived as [referee/](referee/)`rN_*.py`, with `rN_output.txt`.

| script | what it checks | result |
|---|---|---|
| r1 | Lemma E, `q` in {0.01, 0.25, 1}, `lambda` up to `min(lambda_1, 1/2)`, `s` from `1e-10` | `eps/e <= 0.993`; `|det T - 1|/(4 e^2) <= 0.245`; `||T - I||/(4e) <= 0.121`; `sup |delta v| / (3 sigma (1 - sigma)(eps_L + eps_R)) <= 0.334`; `sup |delta v| / bound <= 0.082`; `sup |delta y| / bound <= 0.120`; identities to `1e-15` |
| r2 | Lemma F, Lemma G, K2(ii), Remark K1 | exact TV below (7) everywhere (ratio at most 0.026); `density * Z_0 <= 1.0001`; the K2(ii) bound holds, short intervals included (ratio at most 0.10); the deviation over `lambda^4` grows linearly in `log(1/lambda)`; the pair fraction equals `p_0` to `1.5e-15` |
| r3 | Lemma K by direct quadrature | (i) agrees to `1.5e-10`; (ii) agrees to `4e-14` |
| r4 | Prop. L, with the exact map and slacks, `delta ell` included, calibration taken as exact | relative error against the Step 5 main term: 0.19, 0.015, 0.0047, 0.0009 at `1/lambda` = 64, 128, 256, 512, also in a rotated and a reflected frame; reproduces challenge (4), 0.017294 + 0.008246i |
| r5 | the G9 constants in rationals | all reproduced; r2's 0.0106978 still holds with `|ell_2| <= 1.009434` (0.01069785) |
| r6 | Lemma H, seeded near-cone pairs | 24,000 pairs, 2,981 mismatches, all within the bound |
| r7 | Lemmas A, J, M | exact enumerations agree; Lemma J's bound ratio at most 0.69 |
| r8 | the `Lambda_2` deviation rate | slope per unit `log(1/lambda)` of -1.065, -0.995, -0.988, so the rate is `lambda^4 log` |

## The step I trust least

F1: getting from "`iota_lam` is in the three-quarter box" to "the canonical
anchors are in `K_(0.14)`". Every later use of I(d) depends on it:
Proposition L, (R3), and T2 and T3 through (9). As committed it is not
written. The repair is sound, but it has to be added.

Second, F12: Step 4's uniformity over anchor triples is only sketched.

## What I did not check

- r2's own T1 to T4 bounds, (R1) to (R4) and the scan, beyond reading them.
- The imported WCR and QOB results, and the certificate that `Sigma`
  inverts the background map.
- The full calibrated population signal. r4 takes the calibration as exact,
  and K2 is checked separately.
- Lemma G's constant and the mismatch masses of Lemma I, numerically.
- Proposition D and Lemma M on data (no sampling, as instructed).
- `q` outside {0.01, 0.25, 1}, or `lambda` below 1/512.

I made no changes in the worktree and ran no state-changing git command.
