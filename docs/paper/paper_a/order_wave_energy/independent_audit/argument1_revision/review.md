# Review of argument 1, revision 2 (`quantized_order_bound_r2.md`)

2026-11-26. **AI-level review; no human review.**

**Who reviewed it.** An AI referee in a fresh context, which had not seen
this work. It was given the revision as committed at `a032f0e` and these
sources:
- the frozen original `quantized_order_bound.md`;
- AUDIT.md;
- the audit referee's report `proof1_lower_bound/report.md`;
- `coordinate_oracle_barrier.md`;
- the author's [checks.py](checks.py).

**What it was told.** To find what is wrong, and to run only short
deterministic checks: no sampling of causal orders, nothing long, CPU only.
It could not write to the repository.

**Its scripts.** They are archived here verbatim apart from a three-line
header:
- [ref_checks.py](ref_checks.py): 28 PASS, about 8 s. Run it from this
  folder with the worktree root as its argument.
- [section_diff.py](section_diff.py).

The report below is the referee's, lightly formatted. The revision was
corrected after it. Its Section 9 lists the corrections, and Section 2
below says how each finding was handled.

## 1. The referee's report

### Overall verdict: ERROR (E1, E2)

No false step was found in the mathematics of Sections 2 to 5. Every
inequality, constant and case holds as written. The errors are in the
revision's account of itself: two statements about what the document
contains are false as written. Neither touches an inequality, so the proof
of (8) does not change. Once E1 and E2 are fixed, the grade would be
CORRECT WITH MINOR FIXES (M1 to M12).

### Item verdicts

1. **The steps of the revision.**
   - **Lemma 2.1: correct.**
     - Integration by parts gives `|I| <= lam (4 f_max - 2)`.
     - The trivial bound gives `|I| <= b <= f_max l`, and
       `f_max <= 4 f_max - 2`.
     - With `t <= 2 lam sqrt q` and `t < 1`, this gives
       `e <= C_E lam min(l, lam)` and `e < b` for every `lam`.
   - **Lemma 2.2: correct in all four cases.**
     - The identity `g (dT + |dz|) = 2 b dv - |dy|^2` holds.
     - First inclusion: `dT + |dz| >= dT >= dv/sqrt 2`, and
       `c0 delta >= 2 sqrt 2 e`.
     - Second inclusion: `g >= -2 sqrt 2 e dv/(b + dv)`.
     - `u1 = u2`: `|dZ| = |dT|`, so `g <= 0`.
     - `u1 > u2`: `g > 0` would force `b > 0`.
     - There is a notation slip (M7).
   - **Corollary 2.3: correct.** Lemma 2.4 is not even needed there, since
     the reverse margin already gives `dT + |dz| <= c0 delta`. Nothing
     downstream uses the corollary.
   - **Lemma 2.4: correct.**
   - **Section 3: correct throughout.**
     - (3a) is the right link. The `(1 + sqrt 3) h` shift is attained
       exactly (ratio 1.0000), so `L > 1 + sqrt 3` is exactly what this
       argument needs.
     - The non-edge argument: true margins are strictly larger than
       `c0 delta`, with the sign of the centre.
     - Orientation: correct, using the newly stated half-open aligned
       grid. Same-T-interval edges cost `C_amb N^2 h^4`.
     - The uniform bound: an edge implies `| |dT_c| - |dz_c| | <= L h`.
       This gives the `(L + (1 + sqrt 3)/2) h` strip and two intervals of
       total length at most `4w`, and hence `C_unc`.
     - The tree count and (5): the last vertex in parent-first order is a
       leaf; conditioning and induction; Cayley's formula with
       `C(N, K+1) <= N^(K+1)/(K+1)!`.
   - **Lemma 4.1: correct.** The referee re-derived (i), (ii) and (iii)
     independently.
     - In (ii), `[lo, hi]` does not depend on `s_i`, because the pairs are
       distinct. The rows through `i` and the cube factor in `s_i` do not
       depend on `tau_e`.
     - The Euclidean operator norm of the symmetric non-negative matrix
       `(c_ef)` is at most its largest row sum, which is at most `Lambda`.
     - The scaling `p_h(theta) = p(theta/h)` is right, and so are the bounds
       `R <= rho_max/rho_min` and `R' <= rho'_max/rho_min` for `h <= 1`.
     - The monotonicity behind `Lambda_K` holds, using `R >= 1`.
   - **Lemma 4.2: correct.** The descent inequality, then a split at the
     first and last zero.
   - **Lemma 4.3: correct.** The following all hold:
     - the reduction to signed rows;
     - the squeeze;
     - the interval length `2 sqrt(Lambda_K/2) c0 delta sqrt(m)/h`;
     - the constant `2^(M+1) M Lambda_K c0^2`;
     - joint convexity.
   - **Section 5: correct.**
     - The Hellinger convention, the component count, and (7) with
       `sqrt(C_H C_unc/2)` all hold.
     - The three identities at the chosen h hold, and so does
       `N^2 h^4 <= (N h)^3`, which needs `h <= N`.
     - (8) holds for every `N >= 1` and `0 < lambda <= 1`. The exponent
       `-(2K+1)(K-1)/(2K)` is right.
     - The case "h > 1" never occurs, because `h <= 1` whenever `N >= 1`
       and `delta <= 1`. This is harmless (M12).
     - "Not the minimiser" and "exactly (1)" are true: multiply
       `N^(K+1) h^K` by `(N delta/sqrt h)^(2K)`.
     - K = 1 is right, apart from a notation clash (M3).
     - The exponent coverage is right, and the stated K is minimal.
2. **AUDIT.md's items 1 to 6 for argument 1 are all addressed in the text.**
   - Item 1: (3a), the non-edge bullet, and `c0` and `C_E` in Lemmas 2.1
     and 2.2.
   - Item 2: Corollary 2.3, Lemma 2.4 and the Orientation bullet.
   - Item 3: Lemma 4.1 with the scaling paragraph. Its constant matches the
     audit's.
   - Item 4: the Convention paragraph and the component bullet.
   - Item 5: the K = 1 and "not the minimiser" paragraphs.
   - Item 6: the preamble and Section 8.

   F1 to F9 are all reflected.
3. **Against the frozen original.**
   - The SHA-256 `39217a16...564b` is verified, at `a032f0e` and in the
     worktree.
   - The original's Sections 1 to 5 are byte-identical to the `5dda03d`
     freeze.
   - r2's Sections 6 and 7 equal the original's, apart from the italic note
     and one trailing blank line.
   - Section 0 does not list every change (E2). "Later commits added
     records only" is loose (M4).
4. **Section 8 and the preamble: mostly fair.** They do not overclaim the
   theorem. The problems:
   - one false claim (E1);
   - loose statements of the constants (M1, M2);
   - imports that are not listed (M6);
   - a forward reference (M5);
   - wording that understates (M10);
   - stale text in Sections 6 and 7 (M11).

### Errors and gaps

**E1 (Section 8).** "The direction used here, that (1) implies (8) tends
to zero, is written above." It is not.
- The revision deleted the original's "Under (1), h tends to zero ... This
  proves (1)". No argument for general sequences replaces it.
- The Exponents paragraph covers only `N = lambda^-alpha`.
- Nothing downstream changes. The fix is two lines: with `N >= 1` and
  `lambda <= 1`, `(N lambda^2)^(2K) <= N^(3K+1) lambda^(4K) -> 0`, and
  `N lambda^4 <= N lambda^2`.

**E2 (preamble: "Section 0 maps every change").** False at the level of
individual changes. The following are not listed:
- In Section 1:
  - the definitions of `rho_min`, `rho_max` and `rho'_max`;
  - C renamed `C_1`.
- In Section 3:
  - the hypothesis "half-open cubes on a grid; T-intervals coincide or are
    disjoint". This is substantive: the orientation argument needs it;
  - the rewritten boundary-layer reason, and the explicit density
    constants;
  - the discard rule, changed from "orientation not fixed throughout its
    cubes" to "all same-T-interval edges";
  - the deleted "the null-margin functions are uniformly Lipschitz".
- In Section 5:
  - the deleted "Under (1) ... This proves (1)", which causes E1;
  - the new Exponents paragraph with the formula for K.

Nothing mathematical changes. The fix is to add these to the rows for
Sections 1, 3 and 5.

**M1.** "C_K contains ... 2^(M+1)" is loose. `C_H` enters only through
`sqrt(C_H C_unc/2)`, so the factor is `2^((M+1)/2)`: 7.0e13 at K = 14, not
4.9e27. The tree coefficient is `(K+1)^(K-1)/(K+1)!`.

**M2.** "alpha = 1.3 needs K >= 14, with a margin of at most 0.0035" holds
only at K = 14, where the margin is 0.00345. The margin grows with K: 0.0171
at K = 20 and 0.0366 at K = 50, tending to 0.05.

**M3.** The K = 1 display reuses `C_1`, which is already the one-point
constant.

**M4.** "Later commits added records only." In fact `9205b3a` also edited
an existing preamble sentence and replaced Section 6's bracket paragraph.
The claim that matters, that Sections 1 to 5 are unchanged, is true.

**M5.** "This revision is reviewed in argument1_revision/": at `a032f0e`
that folder held only checks.py.

**M6.** Section 8's import list omits two things:
- the model definitions: `P_{N,lambda}`, `P_{N,0}`, the effective
  background, the tips and slack, and the unlabelled order
  (`full_order_resolution_probe.md`);
- the background relation (2) and the `a^4` density. Section 1 asserts
  them; only the referee report, Section 3.1, derives them.

**M7.** In Lemma 2.2, `b` is defined as positive for `u1 < u2` only, but
the other two cases use it as the signed `Delta m`.

**M8.** "The smallest valid c0 is sup |g|/delta" needs "over pairs in D and
lambda in (0, 1]". The quoted referee values are sampled lower bounds of
that supremum.

**M9.** The Section 6 note says that checks.py covers Lemma 4.1.
- checks.py tests step (ii) only, at k = 2, plus Lemma 4.2's ratio on the
  uniform triangle.
- Its Lemma 2.4 check is vacuous.
- "The original's sketch reached the same constant structure" overstates:
  the original gave no constant.

**M10.** "Asymptotic only / no finite-N statement" sits awkwardly with (8),
which is stated for every N with constants that are explicit or computable.
Better: "C_K is not evaluated; no content at reachable sizes".

**M11.** Sections 6 and 7 are now stale.
- Section 6's bracket (`N lambda^3/log N`) disagrees with Section 8's
  `alpha_high <= 2`.
- Section 7's "Sections 1-5 here" now points at the rewritten sections.

**M12.** "If h > 1" in Section 5 is a vacuous case. Harmless.

### What was checked numerically

`ref_checks.py`: 28 PASS, 0 FAIL, about 8 s. The author's `checks.py`:
17 PASS.

- **R0.** The SHA-256. The original's Sections 1 to 5 are identical to the
  freeze, and r2's Sections 6 and 7 are identical to the original's.
- **R1, Lemma 2.1.** 1000 random segments: q up to 9.8, `lambda` up to 1,
  `l` down to 1e-7. The largest `e/bound` is 0.895 and the largest `e/b` is
  0.9875.
- **R2, Lemma 2.2, with the exact phase-dependent `E`.**
  - 151,200 pairs, with 0 violations both with `c0` and with the per-pair
    constant `2 sqrt 2 e/delta`.
  - On the wave null cone, `|g|/(c0 delta)` is at most 0.407.
  - The cases `u1 = u2` and `u1 > u2` pass on 20,000 draws.
- **R3, end to end at `h = h_min`.** q = 1, `lambda = 1/12`, L = 4, exact
  kernel, 1,619,100 pairs.
  - All 6,260 disagreements fall on edges and none on non-edges.
  - Every disagreement has `|margin| <= 0.222 c0 delta`.
  - The control with L = 1 gives 497 disagreements on non-edges, so the
    test has power.
- **R4.** The `(1 + sqrt 3) h` shift is tight (ratio 1.0000). The ratio for
  the uniform bound is 0.9996.
- **R5, Lemma 4.1 at k = 3,** with exact piecewise-polynomial `p` and
  `grad p`.
  - Coverage: all 8 triangle orientations and 8 two-row paths, with
    non-uniform densities.
  - Validated against brute-force grids (0.080401 exact against 0.080937
    on a 300^3 grid) and against finite differences (3e-11).
  - The largest Lipschitz ratio is 0.235; the largest row-wise `c_ef` ratio
    is 0.590.
- **R6.**
  - `|grad p|^2/(2 Lambda p)` is at most 0.222, and 0.149 near the cyclic
    zero set.
  - The ratio for the `sqrt(Lambda/2)`-Lipschitz bound on `sqrt p` is
    0.425.
  - The ratio for the width of the squeeze is 0.439.
- **R7.**
  - (8) holds on 124,025 grid points: N up to 1e40, `delta` from 1e-14 to
    1, and K in {1, 2, 3, 5, 14}.
  - An exact-rational check confirms that the stated K is minimal; the
    margins are as in M2.
- **R8, the imported one-point bound, computed deterministically.**
  TV/lambda^2 for `lambda` from 1/8 to 1/64:
  - 1.644 to 1.623 at q = 0.25;
  - 3.440 to 3.283 at q = 1;
  - 7.08 to 8.05 at q = 4.

  These agree with the audit referee's 1.61, 3.28 and 8.05.

### Least-trusted remaining step

**The imported one-point bound** `TV(mu_lam, mu_0) <= C_1 lambda^2`,
uniform in `lambda` up to 1 and in phase.
- The revision neither writes it nor gives its constant.
- Every `alpha` above 1 depends on it.
- The referee confirmed it only numerically (R8): three values of q, one
  phase, and `lambda` from 1/64 to 1/8.

**Among the steps that are written out:** Lemma 4.1 for general k. The
referee re-derived it in full, but tested it numerically only at k = 3. The
audit referee tested k up to 5, with uniform densities.

### Not checked

- An analytic re-derivation of the exact causal relation of the wave (R2
  and R3 assume it).
- The proof of the strip bound itself.
- `kappa(q) > 0`, and behaviour as q approaches `pi^2`.
- Lemma 4.1 for k >= 4.
- Explicit values of `C_bd` and `C_amb`.
- The audit referee's figures 1840 and 154, which were not re-run.
- The provenance hashes of Section 7.

## 2. How each finding was handled

All findings were applied to `quantized_order_bound_r2.md`; its Section 9
lists them.
- **E1:** Section 5 now has "(1) follows from (8)", the referee's two-line
  argument. It also has "(1) is exactly what (7) delivers", the referee's
  product identity, which replaces a citation.
- **E2:** the missing rows are added to Section 0.
- **M1 to M8 and M10 to M12:** applied as the report suggests.
- **M9:** the note at the head of Section 6 now says what each test of
  Lemma 4.1 covers. The author's Lemma 2.4 check in checks.py now ranges
  over `Delta v` in `[-1, 1]` (15,735 related pairs, 0 violations).
- **The least-trusted step** is now named in r2 Section 8, as an import
  that every `alpha > 1` depends on.

No inequality, constant or exponent of the proof changed.
