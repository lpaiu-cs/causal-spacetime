# Review of one_point_strip_bound.md (the one-point bound written out)

- **Referee:** an independent AI referee in a fresh context, 2026-11-26.
  No human has reviewed this.
- **Scope:** the note as committed at `3fc1ffc`, with
  [checks.py](checks.py) and [checks_output.txt](checks_output.txt).
- **Brief:** the referee was given the code and the documents, not the
  author's conclusions.
- **Archiving:** the referee returned the report below as text. The author
  saved it unchanged, with two exceptions in its Section 7. The scripts'
  folder now points to [referee/](referee/). The re-run row is marked "not
  archived". The referee's scripts are archived there with a
  three-line header. Their outputs are archived beside them.
- **Not archived:** the referee's re-run of `checks.py`. It matches
  `checks_output.txt` line for line, apart from the referee's added timing
  line `exit=0 elapsed=170s`.

---

## Referee report: `one_point_strip_bound.md` and `independent_audit/one_point_bound/` at `3fc1ffc`

### 1. Verdict for the proof: CORRECT WITH MINOR FIXES

I re-derived Steps 1–7 and found every one correct. The fixes concern Remark 4, the import list and notation. None of them touches an inequality or a constant.

### 2. Verdict for the script and results: CORRECT WITH MINOR FIXES

The numbers are right:
- I re-ran `checks.py` and its output is byte-identical to `checks_output.txt` (170 s).
- My own computation, by a different method, agrees to the printed digits.

Several checks are weaker than the note says, though. K5, "the main check", follows from K1 plus algebra on the same discretisation. One Results sentence is inverted, and the preamble misdescribes the earlier support.

### 3. Findings

| id | location | grade | what is wrong or missing | how I checked it |
|---|---|---|---|---|
| F1 | §4 Results, K5 bullet | MINOR | "The constant is at most about 5 to 20 times too large" is the wrong way round. On the grid `TV/(C_1(q,t)λ²)` runs from 4.2e-7 (q=9.8, λ=1) to 0.187. So `C_1(q,t)` is **at least** 5.3x too large, and about 2.4e6x at q=9.8, λ=1. As λ→0 the overestimate tends to `C_1^0/L` = 5.33, 6.48, 10.2, 12.5, 12.6. For the uniform `C_1(q)` that r2 actually uses, the small-λ overestimate is about 31x (q=0.25) to about 2000x (q=9.8). | output lines; `r1_constants.py` |
| F2 | `checks.py` K5 | MINOR | K5 cannot fail on its own. On the same u-nodes and θ-nodes, Step 5 holds at every node exactly once `e ≤ t m` (checked by K1) and `\|nEn\| ≤ e` (Cauchy–Schwarz on the computed G). Then `TV ≤ B ≤ C_1(q,t)λ² ≤ C_1(q)λ²` follows by summation. K5 confirms that Steps 4–6 are implemented consistently; it gives no independent evidence that the computed number is the TV. That evidence comes only from K2, at 10 points. | read the code path |
| F3 | K1, last two sub-checks | MINOR | These checks cannot fail. `dA` is evaluated from Step 1's own closed form, so `\|dA\| ≤ e/(m(m−e))` and `A = 1/m + dA ≥ 1/(m+e)` hold identically, up to rounding. They never test the inverse formula. Only K2 uses `np.linalg.inv`. | code reading; `r4`(i): 200k random matrices with a true inverse, 0 violations |
| F4 | phase grid; "840 cases" | MINOR | Phase `kπ/4` rotates y by `kπ/8`, which is exactly 32k spacings of the 256-point θ-midpoint grid. The discretised TV is therefore exactly phase-invariant. A wrong phase convention (`e^{2iφ}` or `e^{iφ/2}`) would also pass. The check is vacuous even as a test of the code, and the 840 cases are 105 distinct (q, λ). | arithmetic: π/256 spacing |
| F5 | K3 | MINOR | In the code, Δ_λ − Δ_0 = (h²/2)∬a²(1/S_λ − 1/S_0) at each node: the S_∩ terms cancel exactly. So K3 is blind to the `max()` kinks, the only place the midpoint rule is weak. At λ=2⁻¹⁰, \|V_λ−V_0\|/V_0 ≈ 7e-11 to 1e-9, so the closed-form side carries about 1e-6 cancellation error. The recorded 1.03e-6 is rounding noise there. | `r2 small`: (V_λ−V_0)/V_0 = −7.0e-11 (q=0.25), −2.7e-10 (q=1), −1.1e-9 (q=4) |
| F6 | K6 coverage | MINOR | Doubling was tested only at λ ∈ {1/8, 1/256}. At q=9.8, λ=1/64: <br>• the 256-point θ rule is off by 8.8e-6 relative (120.05579 against a converged 120.05473); <br>• doubling moves the value by 1.6e-5, more than the reported "largest change 3.9e-6"; <br>• the printed 120.06 should be 120.05. <br>This does not affect the bound. | `r3_convergence.py`: kink-aware GL20/GL40 against mid256/512/2048 |
| F7 | preamble; §4 "earlier computations"; K7 | MINOR | The note calls the earlier support "two independent deterministic computations … one phase … λ to 1/128". That is wrong on three counts: <br>• the audit referee's `check_c_one_point_law.py` is Monte Carlo in (u, y), with a random phase per case, reaching λ = 1/128; <br>• R8 is deterministic at phase 0.3, but only reaches λ = 1/64; <br>• R8 (`ref_checks.py one_point_tv`) uses the same polar-reduction algorithm, so K7 reproduces a second implementation of one method and is not an independent check. | read both scripts and the run log |
| F8 | §4 table, q=0.25, λ=1/8 | MINOR | 1.645 should be 1.644 (output: 1.64447). | output |
| F9 | Remark 4; import list | MINOR | Remark 4 uses `T(θ)` and the matrix form `E_λ = −t∫ f T(w/λ+φ) dw`, neither of which the note states. §1 gives only the eigenvalues ±t\|G\|, and eigenvalues alone do not give rotation covariance. The import bullet names `E_λ` without its form. Proposition 1 does not need Remark 4. | the rotation identity itself checks (`r4`(ii): 1e-15) |
| F10 | notation | MINOR | Here `s = √q`. The probe, the oracle note and both referee reports use `s = sinh 2β = 2λ√q` and `t = s/c`, which is a hazard when reading them side by side. | text |
| F11 | Step 4 | MINOR | The TV convention (sup\|μ(A)−ν(A)\| = 1 − ∫min) is used but not stated. It agrees with r2's [0,1] convention. | text |
| F12 | Step 5 (not an error) | — | The constant can be tightened. `(B−A)^+ ≤ e/(m(m+e)) ≤ e/m²`, and this is attained (ratio 1.000000000 in `r4`). So Δ_λ ≤ πh²(1+t)J and `max(Δ) ≤ πh²J/(1−t)`, which lowers `C_1(q)` by a factor `1+t_1` (1.71 to 1.99). At leading order `TV ≈ J/(πV̂)`, so the θ step loses a further factor π. | derivation, `r4` |

**Items checked and correct:**
- **(a) The proof.**
  - The 2x2 identity and (2): `A−B` runs over `[−e/(m(m+e)), e/(m(m−e))]`.
  - The v-sections and their intersection, including the strict `−1/2 < u < 1/2` and the u = ±1/2 slices being null sets.
  - The polar integral `h²/(2S)`, also for `S_∩`.
  - Identity (4).
  - `S_λ ≥ S_0/(1+t)`, `S_∩ ≥ S_0`, and the tip cancellation `m_∓²/M² ≤ 1`.
  - `V̂` and `W`.
  - `t = tanh 2β` is increasing in λ.
  - `R T(θ) Rᵀ = T(θ+φ)` with R the rotation by φ/2.
- **(b) Fit to r2.** Proposition 1 is what r2 needs. It uses:
  - the same y chart, in which `c` cancels and the tips sit at y = 0;
  - D_λ defined by the wave relation and D_0 by the background relation;
  - the same tips, (−1/2, 0, 0) and (1/2, 1, 0), consistent with the oracle's `tip < X < Y < tip` and `A+b+B=M`;
  - density a², which is r2's a⁴ in (T, Z, y);
  - uniformity in phase and in 0 < λ ≤ 1;
  - r2's first term of (7), `C_1 N δ`.

  I found no mismatch.

### 4. My own numbers

**Constants table.** I recomputed it with scipy `quad` on the defining integrals, with m also from `quad`. Every printed digit matches: f_max, C_E, W/V̂, C_1^0, t_1, (1+t_1)/(1−t_1) and C_1(q) for all five q. The closed forms for V̂ and W agree with quadrature to 12 digits.

**Independent TV** (`r2_indep_tv.py`). Method:
- m and G from cumulative GL-16 between nodes; the G spot-check against QAWO agrees to 4.9e-16;
- the 2x2 inverse by adjugate;
- V_λ from the ellipse formula πh²/√det(K₋+K₊), with no θ grid;
- V_∩ by kink-aware GL-20 in θ;
- TV taken as ½∫\|p_λ−p_0\|, not from (4).

Results:
- 25 (q, λ) points × 2 phases (0 and 1.234, the second off the grid): agreement with the printed TV/λ² to within 5-significant-figure rounding everywhere. The one exception is q=9.8 at λ=1/64 and 1/256, about 1e-5 (F6).
- Phase 1.234 gives the same value as phase 0 to 7 digits.

**Limit λ→0.** I derived the leading order: `TV/λ² → L(q) = (2√q/(πV̂)) ∫ a² ⟨|f(u) − f_max e^{iχ}|⟩_χ (m₊²+m₋²)/M² du`.

| q | L(q), my limit | note at λ=1/1024 |
|---|---|---|
| 0.25 | 1.61982 | 1.6198 |
| 1 | 3.26922 | 3.2693 |
| 4 | 8.06357 | 8.0637 |
| 9 | 81.96 | 81.44, still rising |
| 9.8 | 1124.5 | 765.7, only 68% of the limit |

This confirms λ² is the true order for q ≤ 4, which bears on Remark 3. It is not a proof.

**Further checks:**
- 24 seeded random cases (q up to 9.8, λ from 2⁻⁸ to 1, random phase): Proposition 1's chain holds in all of them, with max `TV/(C_1(q,t)λ²)` = 0.116.
- Cartesian 3001² brute force against the polar formula, at near-tip u (u = −0.4995, 0.499, 0.4999 at q=9.8; u = −0.49, 0, 0.4999 at q=9): volumes agree to 1e-8, and Δ_λ, Δ_0 to 7 digits.
- Refining the u-panels by 2x and 4x changes nothing at the 1e-10 level.

### 5. Step I trust least

Step 2's reduction depends on the imported exact wave relation, in the form `M̃ = mI + E` with `e = t|G|` and `t = tanh 2β`. Every constant scales through `t` and `e ≤ t m`, and every check, K2 and mine included, assumes that form. Among the note's own steps there is nothing I distrust.

### 6. What I did not check

- The exact wave causal relation in the Rosen/y chart, re-derived only in report §3.1, and the vacuum condition.
- Lemma 2.1 beyond re-deriving its integration-by-parts bound.
- Any step of r2 other than its use of C_1.
- The validity of the earlier sketch in report §3.3.
- q in (9.8, π²) numerically.
- λ > 1.
- κ(q) and the oracle note's other statements.

### 7. My scripts

Archived in [referee/](referee/):

| file | what it does |
|---|---|
| `r1_constants.py` | constants table and the limit L(q) |
| `r2_indep_tv.py` | independent TV; modes `grid`, `random`, `small`; outputs `r2_grid_out.txt`, `r2_random_out.txt`, `r2_small_out.txt` |
| `r3_convergence.py` | u and θ convergence study; output `r3_out.txt` |
| `r4_step1_and_cartesian.py` | Step 1 random test, Remark 4 identity, Cartesian brute force; output `r4_out.txt` |
| `rerun_author_checks.txt` | re-run of the author's `checks.py`, byte-identical to the recorded output (not archived; see the header) |

---

## How each finding was handled

All changes are in [one_point_strip_bound.md](../../one_point_strip_bound.md);
its Section 5 has the full table.
- **Unchanged:** no inequality, constant or statement of Proposition 1.
- **Not edited:** `checks.py` and its recorded output stay as they were
  run. The note now says what each check does and does not test.

| finding | handling |
|---|---|
| F1 | applied: "at least 5.3 times too large"; the overestimates of `C_1(q, t)` and `C_1(q)` are stated |
| F2, F3, F5 | applied: a new paragraph, "What each check actually tests", in the note's Section 4. The criteria fixed before the run are left as written |
| F4 | applied: Remark 4 and the Results say the phase grid is invariant by construction; the off-grid phase 1.234 is credited to the referee |
| F6 | applied: K6's coverage is stated; q = 9.8 at `lambda = 1/64` is 120.0547 |
| F7 | applied: preamble and Section 4 |
| F8 | applied |
| F9 | applied: Section 1 states the matrix form of `E_lambda` and lists it among the imports |
| F10 | applied: a notation warning in Section 1 |
| F11 | applied: the convention is stated in Step 4 |
| F12 | recorded as Remark 5, with the sharper uniform constants. Proposition 1 keeps the reviewed constant |

- **The referee's limit `L(q)`** is quoted in Remark 3 as the referee's,
  not re-derived.
- **The step the referee trusts least** is the imported exact wave
  relation in the form `mI + E`. The note's Section 5 records it, and it
  stays on the import list of quantized_order_bound_r2.md Section 8.
