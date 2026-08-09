# P14 §8 P3-E — discriminability, exploration stage (exploratory)

**EXPLORATORY.** No gate, no threshold, no verdict, no confirmation
seed window (design §8.2). P3-E explores and CHOOSES — the primary
operating point, the primary discriminator, the fit ranges — and its
choices are frozen inputs to **P3-C**, which runs on a fresh seed
block and is the only stage whose composite equivalence gate may
close anything. Component seeds `20260821–24`; burned seeds
(`20260808`, `777–781`, `20260811–14`) verified absent by the script
and a test.

**The termination sentence** (P3-C's, never P3-E's), frozen in the
design review:

> "동결된 운영점에서, 3a의 (i) 표준화 평균이동, (ii) 순위 판별(AUC),
> (iii) 중점-역치 분류 — 이 세 규칙으로는 분리가 해상되지 않았다."

It closes the chosen operating point under those three rules only —
not 3a's full distribution, not any other pure-order statistic (the
`N(0,1)` vs `N(0,4)` counterexample: mean shift, rank discrimination,
and a midpoint threshold are all blind to a pure variance
separation).

## Protocol (frozen in the design review, five rounds, in-session)

- **Ladder-E (exponent).** One box (slice-a1.0, square transverse),
  one density (`E[N] = 300`); only `w` varies, `A = w²`, nine rungs
  `A ∈ [2⁻⁶, 4]`. Regression axis is `log A` (a slope in `w` would
  read 2 where the heuristic says `O(|A|)`). Shared point sets across
  rungs → the ladder is one sprinkling cluster; every CI is a
  sprinkling-level cluster bootstrap. `D = gained + lost` and
  `Δr = gained − lost` are regressed SEPARATELY — the review's
  structure is `D = O(A)`, `Δr = O(A²)`, with the pilot's gained/lost
  ratio falling toward 1 as `A → 0`. Fits: `D` over the full ladder
  and `A ≤ 0.5`; `Δr` asymptotic over `A ≤ 0.0625` only, plus an
  `A ≤ 0.5` *finite-range effective exponent*. The three lowest
  rungs carry **total `n = 800`** (400 shared + 400 low-A-only),
  sized for a 95% slope-CI half-width ≤ 0.15. Bootstrap replicates
  with a non-positive rung mean are recorded as **fit failures**,
  never discarded; the drop-highest-rung refit is a *refit
  diagnostic* (a two-point slope), not an independent validation.
  `w = 0` is a pipeline null only. The `−A` contrast is the `x↔y`
  swap on the square box (the geometry admits no negative `w`);
  its prediction is distributional equality of `D` and of
  gained/lost — the axis swap itself is statistic #4's business.
- **Ladder-G (operating-point map).** P1's seven points, `n = 250`
  each (CV ≤ 0.31 → RSE(mean D) ≤ 2%) — never an exponent ladder
  (box and density co-vary). `D` as the §5.1 interval, gained/lost,
  per-arm 3a with the per-sprinkling samples stored (the P3-C power
  model refits on them), paired difference with absolute SE. AUC and
  `s` are REALIZED precision (the `√(2/n)` SE holds only near
  `s = 0`), never sizing bases.
- **Class C characterization (descriptive only).** aniso-a1.0,
  `E[N] ∈ {1200, 2400}`, fixed `n = 72`, both arms under the ONE
  eligibility rule (the curved guard's). Reported: causal-interval
  counts, cardinality profiles, opportunity counts `ΣC(m_i,k)`, raw
  chains, opportunity-weighted `Ĉ_k = Σchains_k / ΣC(m_i,k)` (the
  frozen primary normalization; no cardinality cutoff), per-arm AND
  union zero-denominator rates with exact CP intervals, and the
  full-pipeline wall-clock (median and p95, CPU recorded). **No
  pass/fail, no density freeze**: a `0/72` rule was rejected in
  review (true rate 1% passes it with probability ~48.5%, and
  expanding `n` on a failure is optional stopping) — the freeze gate
  is a later, power-sized decision. `C₂, C₃` are primary-candidate,
  `C₄` exploratory.
- **P3-C preliminary calibration (code only).** Margins frozen at
  `ε_s = 0.0806, ε_AUC = 0.0233, ε_BA = 0.0285` — the accuracy
  standard never moves to meet a power result; joint 90% power is
  met by SAMPLE SIZE. Preliminary candidate `n = 4800/팔` under the
  normal design model; certification is the exact CP 95% LOWER bound
  of the joint pass rate ≥ 0.90, never the point estimate. After
  P3-E, the power model is refit on the measured 3a distribution at
  the SAME margins, the final `n` is recomputed and committed, and
  P3-C is separately approved.

Everything below is RENDERED from the committed artifact
`p14_probe_p3e_results.json` by the script's table functions, and a
test asserts the doc embeds the renderings verbatim.

## 1. Ladder-E: exponents

| fit | range | slope | boot 95% CI | failures | drop-highest refit |
|---|---|---|---|---|---|
| d_full | A ≤ 4 | 0.932 | [0.928, 0.936] | 0/4000 | 0.941 |
| d_low | A ≤ 0.5 | 0.961 | [0.955, 0.968] | 0/4000 | 0.972 |
| dr_asymptotic | A ≤ 0.0625 | 1.977 | [1.827, 2.173] | 0/4000 | 2.057 |
| dr_finite_range | A ≤ 0.5 (effective) | 1.622 | [1.564, 1.693] | 0/4000 | 1.741 |

The **−A contrast** (x↔y swap at `w = 1`): paired `D` difference `+0.000146 ± 0.000196` (0.75σ); gained `0.01491` vs `0.01475`, lost `0.00082` vs `0.00083` — gained matches gained and lost matches lost across the polarization rotation, the corrected prediction (the axis of the gained pairs swaps, which only statistic #4 can see).

## 2. Ladder-G: the operating-point map

| point | mean D | RSE | D interval | g/l | Δf̄ ± SE | s (realized ± SE) | AUC [boot 95%] |
|---|---|---|---|---|---|---|---|
| slice-a0.3 | 0.00171 | 1.2% | [0.00171, 0.00171] | 6.6 | +0.00126 ± 0.00002 | 4.50 ± 0.17 | 1.000 [0.999, 1.000] |
| slice-a0.6 | 0.00617 | 1.0% | [0.00617, 0.00617] | 11.7 | +0.00520 ± 0.00007 | 6.70 ± 0.23 | 1.000 [1.000, 1.000] |
| slice-a1.0 | 0.01563 | 0.9% | [0.01563, 0.01563] | 18.3 | +0.01401 ± 0.00015 | 8.09 ± 0.27 | 1.000 [1.000, 1.000] |
| aniso-a1.0 | 0.05122 | 0.6% | [0.05122, 0.05122] | 98.3 | +0.05019 ± 0.00031 | 11.70 ± 0.38 | 1.000 [1.000, 1.000] |
| roomy-a0.2 | 0.00055 | 1.3% | [0.00055, 0.00055] | 1.0 | -0.00001 ± 0.00001 | -0.00 ± 0.09 | 0.499 [0.449, 0.550] |
| high-a2.0 | 0.04317 | 0.5% | [0.04317, 0.04317] | 2.4 | +0.01805 ± 0.00031 | 2.58 ± 0.12 | 0.966 [0.952, 0.978] |
| edge-a2.4 | 0.04393 | 0.5% | [0.04393, 0.04393] | 0.5 | -0.01536 ± 0.00030 | -1.57 ± 0.10 | 0.132 [0.101, 0.164] |

## 3. Class C characterization (descriptive)

| E[N] | arm | intervals/spr | k | Σopp | Σchains | Ĉ_k | zero-denom | union-zero |
|---|---|---|---|---|---|---|---|---|
| 1200 | curved | 108.1 | 2 | 12980 | 1223 | 0.0942 | 3/72 | 4/72 |
| 1200 | curved | 108.1 | 3 | 31763 | 62 | 0.0020 | 7/72 | 9/72 |
| 1200 | curved | 108.1 | 4 | 77437 | 0 | 0.0000 | 12/72 | 16/72 |
| 1200 | flat | 105.1 | 2 | 11493 | 1103 | 0.0960 | 4/72 | 4/72 |
| 1200 | flat | 105.1 | 3 | 26794 | 63 | 0.0024 | 9/72 | 9/72 |
| 1200 | flat | 105.1 | 4 | 63053 | 1 | 0.0000 | 16/72 | 16/72 |
| 2400 | curved | 410.9 | 2 | 159609 | 16186 | 0.1014 | 0/72 | 0/72 |
| 2400 | curved | 410.9 | 3 | 634905 | 1812 | 0.0029 | 0/72 | 0/72 |
| 2400 | curved | 410.9 | 4 | 2493153 | 64 | 0.0000 | 0/72 | 0/72 |
| 2400 | flat | 399.9 | 2 | 138559 | 14163 | 0.1022 | 0/72 | 0/72 |
| 2400 | flat | 399.9 | 3 | 512354 | 1575 | 0.0031 | 0/72 | 0/72 |
| 2400 | flat | 399.9 | 4 | 1859132 | 59 | 0.0000 | 0/72 | 0/72 |

## 4. P3-C preliminary power record

Preliminary candidate `n = 4800/팔` (normal design model): joint pass `18273/20000 = 0.9136`, exact CP 95% `[0.9097, 0.9175]` — the certification (CP lower bound ≥ 0.90) holds. Margins frozen at `ε_s = 0.0806, ε_AUC = 0.0233, ε_BA = 0.0285`; the final `n` is refit on the ladder-G 3a samples at these SAME margins and separately approved before P3-C runs.

## Reading

**One asymptotic structure is measured, one is deferred by the
probe's own rule, and the map has a sign reversal nobody ordered.**

- **Exponents (Ladder-E).** `Δr`'s asymptotic fit gives
  `1.977 [1.827, 2.173]` with the refit diagnostic (2.057) INSIDE
  the CI: the `O(A²)` scalar suppression is measured, not assumed.
  `D`'s low-A slope is `0.961 [0.955, 0.968]`, approaching 1 from
  below — **directionally consistent with the `O(|A|)` heuristic,
  but its refit diagnostic (0.972) falls OUTSIDE the CI, so by the
  frozen rule the asymptotic-recovery claim for `D` is DEFERRED:
  not yet asymptotic at measurable amplitudes.** The `A ≤ 0.5`
  effective `Δr` exponent (1.622) documents the same pre-asymptotic
  reality from the other side. Zero fit failures in 4000
  cluster-bootstrap replicates; the −A contrast matches
  distributionally.
- **The map (Ladder-G).** The unsigned probe `D` is strong wherever
  `a` is large (0.043–0.051 at aniso/high/edge). The SIGNED 3a
  shift REVERSES across operating points: `+0.050` at aniso-a1.0
  (g/l = 98), `+0.018` at high-a2.0 (g/l = 2.4), **`−0.015` at
  edge-a2.4 (g/l = 0.5)**. What this establishes is that the sign
  DEPENDS ON THE OPERATING POINT / window shape — high-a2.0 and
  edge-a2.4 differ in `dv`, transverse aspect, and density, not in
  `a` alone (edge's defocusing-`x` window is wider than its
  focusing-`y`, favouring losses; aniso's `dy ≫ dx` amplifies
  gains). A fixed-geometry zero crossing in `a ∈ (2.0, 2.4)` is NOT
  established; showing one needs a fixed-box, fixed-density ladder
  in `w` alone — an interesting follow-up diagnostic, not blocking
  P3-C. roomy-a0.2 is a LOW-SIGNAL control (`A ≠ 0`; `s = −0.00`,
  AUC 0.499) — the true null is `w = 0`, the pipeline check.
- **Primary-candidate reading** (P3-E's choice, frozen for P3-C):
  **aniso-a1.0** — the largest `|Δf|`, `s = 11.70 ±
  0.38`, AUC pinned at 1.000, and Class C characterization already
  co-located there. The shape-dependent reversal makes edge-a2.4 a
  poor 3a operating point (its mean channel is nearly dead) but a
  fascinating physics target for statistic #4.
- **Class C.** At `E[N] = 2400` every zero-rate is 0/72 and pooled
  `Ĉ_2`/`Ĉ_3` differ between arms at the ~1e-2 relative
  level — descriptive only; the density freeze is the later
  power-sized decision, and these numbers are its inputs.

**What P3-E hands to P3-C:** the frozen operating point
(aniso-a1.0), the three rules, the stored 3a samples for the power
refit at fixed margins, and a recomputed `n` to be committed and
separately approved.

## Scope and limits

- Ladder-E slopes are effective exponents over their stated fit
  ranges; the asymptotic `Δr` claim stands only if the refit
  diagnostic stays inside the bootstrap CI.
- Ladder-G's AUC and `s` are realized precision on paired-design
  marginals — classification-grade claims belong to P3-C's unpaired
  design.
- The Class C table is characterization: it informs, and does not
  make, the later power-sized density decision.
- Nothing in P3-E closes anything: the termination sentence is
  P3-C's, at its frozen operating point, under its three rules.

## Changelog

Initial record.

Post-merge review corrected two INTERPRETATIONS, no numbers: the
`D` asymptotic-recovery claim violated the probe's own frozen refit
rule (drop-highest 0.972 sits outside [0.955, 0.968]) and is
deferred -- directionally consistent with `O(|A|)`, "not yet
asymptotic at measurable amplitudes"; and the "zero crossing in
`a ∈ (2.0, 2.4)`" over-read the operating-point map -- high-a2.0
and edge-a2.4 differ in window shape and density, not `a` alone, so
what is established is a SHAPE-DEPENDENT sign reversal (a
fixed-geometry `w`-ladder would be the follow-up diagnostic).
roomy-a0.2 is relabelled a low-signal control (`A != 0`; the true
null is `w = 0`).

Review R1 corrected two instrument defects, neither moving a number
the campaign measured: the ladder bootstrap resampled each rung
independently, destroying the shared-sprinkling cluster the protocol
itself specified -- one index vector per stratum per replicate now,
and the CIs recomputed from the unchanged raw record (dr_asymptotic
tightened to [1.827, 2.173]; conclusions unchanged); and
`pair_flips` counted a pair ambiguous in one arm as a definite flip,
overstating §5.1's D_lower and double-counting in the upper bound
-- the ambiguous union is masked first (campaign ambiguity was
exactly zero, so every published count stands; the flaw is closed in
code and pinned by a one-arm-ambiguous regression test).
