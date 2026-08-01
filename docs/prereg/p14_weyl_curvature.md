# P14: which pure-order statistic reads a pure-Weyl deformation, and at what density? — design draft

Status: **DESIGN DRAFT v0.7 (2026-08-01). NO STOCHASTIC PROBE OR CAMPAIGN
HAS RUN. NOTHING IS FROZEN.** Two deterministic design checks HAVE run and
are committed beside this document (`docs/prereg/p14_checks/`), because an
output pasted into prose is not a check; their frozen, stamped copies move
to `docs/prereg/frozen/p14/` if and when this design freezes. The document
still states no gates, thresholds, sample sizes, or seed windows: those
belong in a freeze that follows Section 8's probe.

Dates are local (UTC+9, Asia/Seoul). **Dimension notation, fixed once:**
`n` is SPACETIME dimension throughout this document. The frozen P12/P13
records' handover phrase "`d >= 3`" reads consistently as SPATIAL
dimension and is not edited there; this document does not use `d`.

Claims are labelled **[VERIFIED]** (checked by a committed script or a
frozen artifact of this repository) or **[TO VERIFY]** (must be
established before any freeze, on the Section 3 rule that killed Stage C
v1 — an expansion coefficient nobody checked).

**Decision record.**

- 2026-08-01: both probes run (§8, §8.1) — the plane wave carries the
  detectability question; a separate, narrow Schwarzschild measurement
  prices the generalization path.
- 2026-08-01: the interval rule for accuracy gates is adopted house-wide
  (§6, `AGENTS.md`), with the computation showing it would have made P12
  itself infeasible under P12's own cap (§6.1).

**Review record. R1 (2026-08-01, in-session), seven findings, all
accepted, none cosmetic:**

- **R1.1** The question was mis-posed. Causal order determines the
  conformal class (Malament 1977; Hawking–King–McCarthy 1976), so at the
  continuum level a Weyl-deformed order DIFFERS from flat by theorem. The
  design's question is finite-density detectability, not existence. §1
  rewritten.
- **R1.2** The draft's control C0 rested on conflating the volume FORM
  with the volume of a causal INTERVAL. Refuted by computation: at equal
  proper time the diamond volume sits 0.4% (`wT = 1`) to 7.3% (`wT = 2`)
  above flat, inside the conjugate bound. C0 deleted; §4.4 added; design
  check 2 committed. The four-layer split of §3 replaces the old
  two-channel story.
- **R1.3** Statistic priorities reordered: paired relation disagreement
  first (unsigned, uses the same-points pairing directly, survives the
  quadrupole cancellation); chain vectors second; Myrheim–Meyer demoted
  to a relation-fraction discriminator; coordinate anisotropy demoted to
  diagnostic. §5 rewritten.
- **R1.4** The profile `A(u)(x^2 - y^2)` stays: trace-free is FORCED by
  vacuum, asymmetric coefficients would break Ricci-flatness, and the
  quadrupole symmetry pushes scalar responses to `O(A^2)` rather than to
  zero. A rotating-polarization `K(u)` is the named second profile. §4.5.
- **R1.5** The conjugate-point constant is analytic, not a guess: the
  focusing Jacobi propagator is `sin(w Du)/w`, first zero exactly
  `pi/sqrt(A)` (Harte & Drivas 2012). Upgraded from [TO VERIFY] to
  analytic-with-implementation-check; near-null ambiguity band and
  box-containment requirements added. §4.3, §8.
- **R1.6** The sandwich-profile escape was wrong and is deleted: compact
  support in `u` does not restore global hyperbolicity while `H` stays
  quadratic in the transverse plane (Flores & Sánchez 2003). The patch is
  a constant-`A` slab, stated explicitly; subquadratic pp-waves are the
  escape and a different family. §4.3, §10.
- **R1.7** Dimension notation was self-contradictory ("Weyl == 0 for
  d <= 3" next to "d >= 3 is the precondition"). Fixed: Weyl is
  identically zero for spacetime `n <= 3` and first non-trivial at
  `n = 4`; in `n = 3` the obstruction to conformal flatness is the
  Cotton tensor, a different channel and a different design. §2. Status
  line and §9 reproducibility fixed (scripts committed).

**R2 (2026-08-01, in-session): physics design PASSED; two P1 defects,
both fixed here. No redesign.**

- **R2.1** `D`'s undecided pairs had no construction, only an
  acknowledged tension — and a fixed exclusion band would have removed
  the near-cone population that carries the entire signal, then reported
  the removal as a small `D`. Replaced by an interval construction:
  verified error interval for `sigma` per arm, either arm's interval
  containing zero marks the pair `ambiguous`, denominator stays
  `C(N, 2)`, and `[D_lower, D_upper]` is reported with
  `ambiguous_fraction`. A band that eats the signal now shows as a wide
  bracket, never as a null. Conditional `D` is diagnostic only;
  escalation is adaptive precision, never a wider tolerance. §5.1.
- **R2.2** `D` is a pair-dependent **U-statistic** over one point set,
  so treating `C(N, 2)` pairs as independent samples would understate
  its standard error by a factor growing with `N`. The replication unit
  is the **sprinkling**, everywhere, for every quantity derived from
  `D`. §5.2, and §8 P3/P4 rewritten to match.
- **R2.3** §6.1 quoted the **break-even** `n` (`1538`) as the price a
  design must pay. Break-even is 50% power by construction. At P12's own
  90% convention the number is **`4204`**, `3.76x` rather than `1.38x`,
  and the frozen cap of `1300` is exceeded 3.2 times over rather than
  marginally. Corrected here and in `AGENTS.md`, with a simulation check
  (`0.499 / 0.801 / 0.901` at `n = 1538 / 3141 / 4206`). §6.1.
- **R2.4 (accepted, upgrade not defect)** The axis-volume `A^2` response
  has a closed coefficient: `V_A/V_0 = 1 + (wT)^4/252 + O((wT)^8)`,
  confirmed independently here and now asserted by the committed check.
  §4.4. The `D`-versus-scalar hierarchy stays `[TO VERIFY]`.

**R3 (2026-08-01, in-session): physics design PASSED again; one P1
protocol conflict and three smaller items, all fixed here.**

- **R3.1 (P1)** §4.2's pair rules contradicted §5.1's `D`. It allowed a
  fixed exclusion band for near-null pairs, and demanded diamond
  containment "for every pair entering any statistic" — but not every
  pair in a finite box has a contained diamond, so either would drop
  pairs or move `D`'s denominator, and either invalidates §5.1. Fixed by
  separating the rules by statistic class (§4.6): **Class R** (`D`, its
  split, global relation fractions) keeps every pair and `C(N, 2)` with
  no containment condition, since a predicate between two points does
  not care what the box holds; **Class C** (interval cardinality, chain
  vectors, diamond-based ordering fractions, the §8 P2 check) uses a
  frozen guard inset with its own reported denominator. The band is
  withdrawn from §4.2 outright, and §4.3 now makes the conjugate bound a
  property of the SLAB rather than a per-pair filter, which is what lets
  `C(N, 2)` be unconditional.
- **R3.2 (P2)** §8 P2 said "within Poisson error" for two counts taken
  over overlapping regions of the SAME realization, which are
  correlated. Fixed with both halves stated: a marginal per-arm check
  against known `rho V`, and a paired ratio whose uncertainty is
  estimated between independent sprinklings, with the exact
  `Cov(N_A, N_0) = rho Vol(I_A ∩ I_0)` and
  `Var(N_A - N_0) = rho Vol(I_A △ I_0)` as the cross-check.
- **R3.3 (P2)** `D` consumes two orders plus the point correspondence,
  so it is a paired counterfactual and not a property of a single causal
  set — §7's "a finite causal set carries..." did not follow from it.
  `D` is now scoped as a **paired-order sensitivity diagnostic**, and
  the single-poset claim is assigned explicitly to §5's items 2 and 3.
  §5 and §7.
- **R3.4 (P2, before freeze)** The exploration/confirmation boundary is
  now stated as a seed rule rather than left implicit: confirmation runs
  on a fresh block disjoint from everything the probe consumed. §8.2.
- **R3.5 (P3)** Two stale `d >= 3` mentions in
  `docs/paper/results_source_wip.md` corrected.

**R4 (2026-08-01, in-session): R3.1's P1 resolved; one new P1 where the
guard rule reached back into the claim. Four items, all fixed here.**

- **R4.1 (P1)** §4.6's guard inset decides Class C eligibility from
  COORDINATES, so a guarded chain vector needs the poset *and* an
  external interiority label — which means R3.3's assignment of the
  single-poset claim to "items 2 and 3" did not hold. Fixed by
  separating the two axes that v0.5 conflated: what a statistic COUNTS
  fixes its admissibility rule, what a statistic is GIVEN fixes the
  claim it can carry (§4.6.2's table). The strictly single-poset claim
  now rests on **item 3a alone** — the global relation fraction over all
  elements — and the guarded statistics are demoted to single-ARM
  diagnostics. An order-invariant eligibility rule that would promote
  them back is named and assigned to the probe, with the warning that
  selecting on interval cardinality is forbidden because truncation
  biases that very quantity. §4.6.2, §5, §7.
- **R4.2 (P2)** The inset's freeze time was stated two ways: §4.6 froze
  it before the probe, §8.2 let the probe pick it. Resolved by what it
  is — a closed-form function of the operating point, hence DETERMINED
  rather than chosen (§4.6.1). Definitions and derivation rule frozen
  now; the number follows the operating point and freezes with it before
  confirmation. §5's "may not be chosen after seeing its curves" is
  withdrawn as written — it contradicted §8.2 and would have forbidden
  the probe from doing its job.
- **R4.3 (P2)** "Independent per arm" was the wrong phrase: the
  marginals are Poisson, the counts are not independent. And
  `Var(N_A - N_0) = rho Vol(I_A △ I_0)` is exact for the raw difference
  but is not the cross-check for a RATIO. Replaced with the residual the
  ratio actually tests, `N_A - r N_0` at `r = V_A/V_0`, whose variance
  is `rho (V_A + r^2 V_0 - 2 r Vol(I_A ∩ I_0))`. §8 P2.
- **R4.4 (P3)** `assert_windows_disjoint_and_fresh` exists in **P12's**
  runner only; P13 enforces the same property via its own frozen window
  constants and tests. "P12/P13 runners" corrected. §8.2. Also: §4.6 had
  been inserted into the middle of §4.5, splitting the quadrupole
  discussion; moved below it.

**R5 (2026-08-01, in-session): P1s all cleared; one P2 estimand and two
stale sentences. No redesign.**

- **R5.1 (P2)** §8 P2 still estimated the uncertainty of the per-sample
  ratio `N_A / N_0`. Wrong estimand: with a finite Poisson denominator
  `E[N_A/N_0] != V_A/V_0` and `N_0 = 0` has positive probability, so the
  per-sample ratio is biased and adding sprinklings does not fix it. The
  right quantity was already in the document as a cross-check and is now
  the test: per sprinkling `Z = N_A - r N_0` with `r = V_A/V_0`,
  `E[Z] = 0` exactly, empirical uncertainty from the between-sprinkling
  variance of `Z`, analytic `Var(Z)` as the cross-check **of the same
  quantity**, and the ratio reported secondarily as a pooled
  `sum N_A / sum N_0`. §8 P2.
- **R5.2 (P3)** §10 q5 still said items 2 and 3 carry the single-poset
  claim — v0.5 wording that contradicts §4.6.2, §5, §7 and its own
  neighbour q8. Rewritten: a thin retained fraction costs the
  single-arm diagnostics, not the claim, whose lever is q9.
- **R5.3 (P3)** §4.4 still said the volume prediction is checked "within
  Poisson error", the phrase R3.2 replaced. It now points at §8 P2's
  marginal-plus-residual protocol directly.

**The pattern R5 exposes, recorded because it is the fourth instance:**
every round of this review found the corrected claim consistent at the
site of the correction and stale at a downstream consumer — §7's
assignment after R3, §10 q5 after R4, §4.4 after R3.2. Naming a risk is
not a mechanism (R2), a mechanism must be verified across the
statistic's whole input boundary (R4), and **a changed claim boundary
must be searched for in the open questions and every downstream
consumer, not only fixed where it was found (R5).**

---

## 1. The question, and the theorem that sharpens it

P12 and P13 closed the curvature track on 2026-07-31 and handed over the
same lever from two directions.

- **P12 Stage B-12** (`9.3`): the unchanged P11 estimator **recovers the
  curvature itself** — `Delta_B = -0.2637`, CI `[-0.3037, -0.2233]`,
  top-rung dimensionless recovery `0.1708 <= 0.25`. RECOVERS-CURVATURE.
- **P13 Stage A-13C** (`13.3`): the flat-normalized reading **does not
  notice** curvature out to `tau/ell = 1.5`, i.e. `R tau^2 = 4.5`.
  CURVATURE-ROBUST, control CLEAN.

Both refuse to promote the programme's general thesis, for the structural
reason each states: in 1+1D every metric is conformally flat, so the
causal order there is flat by construction and curvature reaches the
order only through the volume — the count.

The first draft of this document then posed P14's question as "can
curvature reach the causal order without passing through the measure?"
**That question is answered by a theorem, not an experiment** (R1.1). For
distinguishing spacetimes, the causal order determines the metric up to a
conformal factor (Malament 1977, J. Math. Phys. 18 1399, with
Hawking–King–McCarthy 1976). A pure-Weyl deformation changes the
conformal class, so it changes the continuum causal order. Existence of
the order-level signal: **yes, by theorem.**

What no theorem answers, and what this programme's own slogan makes the
right question: the P11 record fixed the division of labour as **"order
fixes conformal structure, the count fixes scale — order + number =
geometry."** P12 and P13 tested the *number* half in the only place 1+1D
offers. The *order-fixes-conformal-structure* half has never been tested
at finite density anywhere in this programme, because in every ensemble
so far the conformal structure was flat.

> **The question.** At finite sprinkling density, which pure-order
> statistic detects the conformal-order deformation produced by pure-Weyl
> curvature, at what effect size per sample, and at what cost?

P14 stands to Malament's theorem as P11 stood to the continuum question:
the theorem guarantees the signal exists in the limit; the experiment
asks whether an instrument at the densities this programme can afford
reads it. A null here is a statement about instruments and densities —
the same scope discipline P13 §13.3 applied — never a counterexample to
the theorem.

## 2. Dimension, stated precisely

**[VERIFIED — standard]** The Weyl tensor vanishes identically for
spacetime dimension `n <= 3` and is first non-trivial at `n = 4`. So the
construction below is 3+1D, and this is a floor, not a choice.

Two clarifications the first draft got wrong or blurred (R1.7):

- In `n = 3`, conformal flatness is NOT automatic: its obstruction is
  the Cotton tensor (see e.g. García, Hehl, Heinicke & Macías 2004,
  gr-qc/0309008). A 2+1D experiment against a Cotton-curved background
  is therefore conceivable — but it tests a different channel than the
  Weyl one the frozen handovers point at, and it is out of scope here.
- Raising dimension alone buys nothing. Every maximally symmetric
  spacetime — Minkowski, de Sitter, anti-de Sitter — has vanishing Weyl
  in every dimension, and FRW is conformally flat too. `dS_4` sits under
  exactly the same ceiling as `dS_2` at several times the cost. The
  requirement is **Weyl != 0**; `n >= 4` is only its precondition.

## 3. Four layers, separated — what actually made P12 affordable

The first draft told a two-channel story ("measure vs order") and
credited P12's affordability wholesale to conformal flatness. Review
R1.1/R1.2 showed that conflates four things that move independently.
Separating them is not pedantry; C0 died of exactly this conflation.

| layer | what it is | dS_2 (P12) | plane wave (P14) |
|---|---|---|---|
| (i) causal predicate | `p prec q`, pairwise | **flat** (conformal flatness) | **curved** |
| (ii) volume form / intensity | `sqrt(-g)`, the sprinkling measure | curved (`ell^2/eta^2`) | **exactly flat** (`sqrt(-g) = 1`) |
| (iii) interval domain | the diamond `J^+(p) ∩ J^-(q)`, a REGION defined by (i) | flat cones, curved weight | curved cones, flat weight |
| (iv) finite-poset statistics | anything computed from the sprinkled order | moved through (ii) only | move through (i) and (iii) |

Two corrections to the record the first draft distorted:

- P12's algebraic causality was a gift of layer (i) being flat. Its
  CLOSED-FORM diamond volume was **not** — that came from `dS_2`'s
  maximal symmetry (position drops out; `Vol = 4 ell^2 ln cosh(tau/2ell)`,
  P12 §10.1). A generic conformally flat metric gives no such formula,
  and §4.4 below exhibits a conformally NON-flat case whose axis-diamond
  volume reduces to a one-dimensional quadrature. Symmetry pays for
  closed forms; flatness pays for cheap predicates.
- "Curvature can only reach the order through the count" was correct in
  P12's setting as a statement about layers: with (i) flat, only (ii)
  carried curvature into (iv). In P14 the situation is the mirror image:
  (ii) is exactly flat, so **everything** that differs between the arms
  enters through (i) — including, via (iii), the interval volume itself.
  That last clause is R1.2, and §4.4 quantifies it.

The standing rule survives unchanged and binds §4.4's formula too: P12
§10.1 refused the small-diamond expansion because its coefficient is
"exactly the kind of unverified constant that killed Stage C v1".
Whatever this design uses must be exact or verified against an
independent computation to a pinned tolerance.

## 4. The construction: a vacuum plane wave, where the sprinkling measure is exactly flat

Take a pp-wave in Brinkmann form,

    ds^2 = 2 du dv + H(u, x, y) du^2 + dx^2 + dy^2,     H = A(u) (x^2 - y^2).

Three properties, all **[VERIFIED]** by the committed check
`p14_checks/p14_brinkmann_vacuum_check.py` (symbolic, exact, profile
`A(u)` left undetermined; output in §9):

1. **`det g = -1`, hence `sqrt(-g) = 1` exactly.** The volume form is
   Minkowski's, identically, for every profile. A uniform-coordinate
   Poisson sprinkling IS the covariant Poisson process of the true
   4-volume — exactly, at any size, with nothing to calibrate.
2. **`R_{mu nu} = 0` exactly.** Vacuum; `H` harmonic in the transverse
   plane is what buys it.
3. **`R^a_{bcd} != 0`** — eight components, all `+/-A(u)`. With Ricci
   zero, the curvature is **pure Weyl** (Petrov type N).

What property 1 licenses — and, after R1.2, only this: the two arms can
share **one** sprinkling. What it does not license is §4.4's subject.

### 4.1 The control is the same points, read with `A = 0`

Since the intensity is identical for every `A`, the curved arm and the
flat arm are **the same sprinkled points**; only the causal relation
differs. Consequences:

- No twin intensity calibration (P12 §10.6). Nothing to converge.
- No cross-arm `m` gate and no spurious-trip disclosure, because a
  sampling-density mismatch between the arms is not representable.
- Perfect pairing: every statistic is a within-point-set difference,
  the strongest form of the paired design P12 §10.3 argued for.

**The corrected statement of what this buys (R1.2):** any difference
between the arms is a functional of the relation change alone — it
cannot be a sampling artifact, an intensity mismatch, or a calibration
residual. It is NOT the statement that order statistics with volume
meaning (interval cardinality above all) stay flat. They do not; §4.4.

### 4.2 Geodesics, and therefore causality, stay tractable

**[VERIFIED — standard, re-derived from the Lagrangian]** `v` is cyclic,
`u` is affine along geodesics, and with `u` as parameter the transverse
equations are linear:

    x'' = A(u) x        (defocusing),
    y'' = -A(u) y       (focusing).

Hill's equation; elementary for constant `A`, a linear ODE for any
profile. This family is the one place `Weyl != 0` and computable
causality are simultaneously available.

**[TO VERIFY]** The step to a *decidable* predicate at campaign cost:
the plane-wave world function is explicit in the literature (Harte &
Drivas 2012, Phys. Rev. D 85 124039) and must be reproduced here and
checked against direct numerical null-geodesic integration to a pinned
tolerance. Two requirements R1.5 adds, both absolute — and **R3.1
corrects how the second one is scoped**:

- **No pair is ever dropped for being near-null.** "Machine precision
  for arbitrary pairs" is too strong a promise: pairs with `sigma ~ 0`
  sit inside floating-point ambiguity. The v0.4 draft offered a fixed
  published band as an alternative to resolving them, which contradicts
  §5.1 — a band is an exclusion, and an exclusion moves `D`'s
  denominator. **Withdrawn.** Undecided pairs stay in the sample and are
  classified `ambiguous`; the response to too many of them is adaptive
  precision or interval arithmetic (§5.1), never a band.
- **Containment is a requirement of COUNTING statistics only.** The
  v0.4 draft demanded that both arms' causal intervals lie inside the
  sprinkled box "for every pair entering any statistic". That is wrong
  as stated and it is the conflict R3.1 names: not every pair in a
  finite box has a contained diamond, so applied to `D` the rule would
  force dropping pairs or changing the denominator, and either
  invalidates §5.1. A relation predicate reads the geometry between two
  points and does not care what the box contains; only a statistic that
  COUNTS SPRINKLED POINTS INSIDE A CAUSALLY DEFINED REGION is biased by
  truncation. §4.6 sets the two rules separately.

### 4.3 Conjugate points: the analytic bound, and what does not escape it

**[VERIFIED — analytic]** For constant `A = w^2` the focusing direction's
Jacobi propagator is `sin(w Du)/w`, so the first conjugate point sits at
exactly

    |Du| = pi / w = pi / sqrt(A)

(Harte & Drivas 2012 give the same). The patch is therefore an explicit
**slab in `u`**, and R3.1 makes its role precise: **the slab's `u`-extent
is frozen strictly below `pi/sqrt(A)`, so patch validity is a property of
the BOX and never a per-pair filter.** Every pair of sprinkled points
then satisfies `|Du| < pi/sqrt(A)` automatically, by construction, and no
pair can be excluded on conjugate-point grounds — which is what lets
`D`'s denominator be `C(N, 2)` unconditionally (§5.1).

What remains [TO VERIFY] is only the implementation against this known
answer, not the constant itself: an assertion that the frozen slab extent
is under the bound, and a run-time check that the sprinkled `u`-range
respects it.

**The escape the first draft offered is wrong and is withdrawn (R1.6).**
Making `A(u)` compactly supported does not restore global hyperbolicity:
`H` remains quadratic in the transverse coordinates, and the classical
quadratic family is strongly causal but NOT globally hyperbolic
(Flores & Sánchez 2003, gr-qc/0211086). Within this family the patch
constraint is the honest device. If a future design needs global
hyperbolicity outright, the move is to a SUBQUADRATIC pp-wave — a
different family, a different design, and not this document's problem.

### 4.4 The interval volume is not flat, and the size of the shift is known

**[VERIFIED]** by the committed check
`p14_checks/p14_interval_volume_constant_a.py`, which reproduces review
R1.2's computation independently (derivation in the script; output in
§9). Constant `A = w^2`, central-axis pair. On the axis `H = 0`, so the
pair's proper time is IDENTICAL in both arms — the comparison is at
equal clock reading by construction. The diamond volume still moves:

    V_A = (pi T^2 / 4) * int_0^T ds / sqrt(a_x(s) a_y(s)),
    a_x = w [coth(ws) + coth(w(T-s))],   a_y = w [cot(ws) + cot(w(T-s))],

against flat `V_0 = pi T^4 / 24`, giving

| `wT` | `V_A / V_0` |
|---|---|
| 1.0 | `1.00400047` |
| 2.0 | `1.07300802` |

both inside the conjugate bound `wT < pi`. The volume form is flat; the
**boundary** of the interval is not; the interval volume follows the
boundary. Three consequences:

- **C0 is dead.** The first draft's negative control — "the P12/P13
  volume-based reading must return flat" — is factually wrong. The
  volume channel is not silent, and a control built on its silence
  would have failed on a correct implementation, which is the worst
  kind of control.
- **What replaces it is stronger: a quantitative prediction.** Interval
  cardinalities over axis-pair diamonds must reproduce the ratio
  `r = V_A/V_0` at the sprinkled density. **The test is §8 P2's
  marginal-plus-residual protocol** — each arm against its own marginal
  Poisson mean, and the paired residual `Z = N_A - r N_0` against zero —
  **not "within Poisson error"**, which is what v0.5 said here and is
  wrong for two counts taken over overlapping regions of one
  realization (R3.2, R5.1). An implementation check with a number on it,
  not a silence test. The flat arm needs no control at all: it IS the
  reference, on the same points.
- **The excess is quadratic in `A` with a closed coefficient** (R2,
  confirmed by the committed check):

      V_A / V_0 = 1 + (wT)^4 / 252 + O((wT)^8)

  The check asserts `(V_A/V_0 - 1)/(wT)^4 -> 1/252 = 0.003968253968`
  over `wT in [0.15, 0.5]` and that the residual scales as `(wT)^4`
  (coefficient `~3.2e-5`), ruling out a `(wT)^6` term. Below
  `wT ~ 0.1` the excess falls under `1e-9` and double-precision
  cancellation against 1 dominates, so the window is chosen where
  physics rather than arithmetic is the limit.

  This is consistent with the small-diamond expansion's leading
  corrections being `R` and `R_00` (Roy, Sinha & Surya 2013,
  arXiv:1212.0631) — both identically zero here — so the volume channel
  responds only at second order. It is exactly why §5 does not lead with
  a volume-type statistic, and it is now a number rather than a trend.

### 4.5 The profile stays, and why (R1.4)

`H = A(u)(x^2 - y^2)` was not an aesthetic choice the design is free to
soften: vacuum forces the quadratic profile's matrix to be trace-free,
so the eigenvalues come in `(+A, -A)` pairs — a constant plus/cross
polarization rotates back to the same form. Making the coefficients
unequal breaks Ricci-flatness and with it the entire §4 case.

The quadrupole symmetry has a cost and a mitigation, both stated now:
scalar expectations lose their `O(A)` term to the `x`/`y` cancellation,
so generic scalar responses start at `O(A^2)` — pushed down, not zero
(§4.4's volume shift is the worked example). The statistic that keeps
`O(|A|)` sensitivity is one that does not let the two transverse
directions cancel, which is §5's first candidate. If a second profile
is ever needed, it is a `u`-rotating trace-free polarization `K(u)` —
transverse equations still linear — not an asymmetric static one.

### 4.6 Which pairs enter which statistic (R3.1), and what that costs the claim (R4.1)

There are **two** admissibility rules, not one, because there are two
kinds of statistic. Collapsing them was the v0.4 defect. Neither rule
may depend on realized data — eligibility is a deterministic function of
the frozen geometry, or it is selection.

**Class R — relation statistics.** `D`, its gained/lost split, and any
global relation fraction over the box. These read the causal predicate
between two sprinkled points and nothing else.

- **Eligible: every pair. Denominator: `C(N, 2)`, unconditionally.**
- No containment requirement: the predicate is a fact about two points
  and the metric, not about what the box contains, so a diamond
  extending outside the box does not bias it.
- No conjugate-point filter either — §4.3's slab makes the bound a box
  property, so it is satisfied by every pair by construction.
- Undecided pairs are `ambiguous` and stay in the denominator (§5.1).

**Class C — counting statistics.** Interval cardinality, the chain
vectors `C_k`, any Myrheim–Meyer or ordering-fraction estimator computed
**on a diamond**, and §8 P2's volume check. These count sprinkled points
inside a causally defined region, so a diamond leaving the box is a
truncated count.

- **Eligible: pairs whose diamond is contained in the box, decided by a
  GUARD INSET rather than per-pair inspection.** Both endpoints must lie
  in an inner sub-box inset from the sprinkled box by the maximum
  transverse excursion a null geodesic can make over the slab's
  `u`-extent — computable in closed form from the same Jacobi
  propagators as §4.3 and §4.4.
- **The inset is the SAME in both arms** and is taken from the binding
  arm — the curved one, whose defocusing direction gives the larger
  excursion. Arm-dependent eligibility would give the two arms different
  eligible sets and destroy the pairing that §4.1 exists for.
- **Its denominator is the eligible-pair count, reported alongside**,
  and never silently substituted for `C(N, 2)`. A Class C number and a
  Class R number are never quoted against the same denominator.

#### 4.6.1 When the inset is frozen — one policy, not two (R4.2)

v0.5 said both "frozen before the probe" (here) and "the probe selects
it" (§8.2). Fixed, and the resolution is forced by what the inset is:

> **The class definitions and the derivation rule are frozen now. The
> numeric inset is DETERMINED — not chosen — by the operating point,
> and is frozen with it before confirmation.**

The inset is a closed-form function of the slab geometry and `A`. Once
the probe fixes those, the number follows; there is no free parameter to
tune. So it is not on the list of things the probe selects, and it is
not a pre-freezable constant either, because it cannot be computed
before the geometry it depends on exists. What §8.2's fresh-seed rule
protects is the operating point; the inset rides along with it.

#### 4.6.2 The guard costs Class C its pure-order standing (R4.1)

This is the correction R4 forced, and it narrows a claim rather than a
method. **The Class C eligibility label is coordinate-derived.** A chain
vector is computed from order alone, but deciding *which* intervals to
compute it on requires knowing which elements sit in the inner sub-box —
a fact about coordinates that no poset carries. So:

| statistic | admissibility | inputs it needs | claim it can carry |
|---|---|---|---|
| `D`, gained/lost split | Class R | two orders **+ the point correspondence** | paired-order sensitivity diagnostic (§5, §7) |
| global relation fraction over ALL elements | Class R | **one poset, nothing else** | **strictly single-poset, pure-order** |
| chain vectors `C_k` on guarded intervals | Class C | one poset **+ coordinate eligibility labels** | single-ARM, NOT pure-poset — diagnostic |
| Myrheim–Meyer / ordering fraction on a diamond | Class C | one poset **+ coordinate eligibility labels** | single-ARM, NOT pure-poset — diagnostic |
| interval cardinality (§8 P2) | Class C | one poset + coordinate labels | implementation check, carries no claim |

Two axes, and v0.5 conflated them: **what a statistic counts** decides
its admissibility rule, **what a statistic is given** decides the claim
it can carry. `D` is Class R and still not single-poset. The guarded
chain vectors are single-arm and still not pure-poset.

The consequence for §7 is real and is not softened there: **the strictly
single-poset claim now rests on the global relation fraction alone**
until an order-invariant eligibility rule exists.

**[TO VERIFY — and it would restore the chain vectors]** Whether Class C
eligibility can be defined order-invariantly, e.g. by requiring each
endpoint to have at least `k` elements below/above it inside the
sprinkled set, as an order-theoretic proxy for interiority. Such a rule
would make the guarded chain vectors genuinely pure-poset and hand §7's
claim a far more sensitive instrument. It is a proxy, not a
certificate — its agreement with true containment has to be MEASURED
against the coordinate guard before anything rests on it, and selecting
on interval cardinality directly is forbidden, since truncation biases
that very quantity. **Assigned:** §8 P1.

**[TO VERIFY]** The closed form for the maximum transverse excursion,
and the fraction of pairs Class C eligibility retains at working slab
geometries. If the inset is so large that few pairs survive, that is a
cost of the volume-type statistics and a reason §5 does not lead with
them — but it must be measured, not assumed. **Assigned:** §8 P1.

## 5. What the instrument reads, in probe order (R1.3)

**No gate is proposed for any of these here, and the probe is allowed to
choose (R4.2).** Which statistic, at what operating point, is exactly
what §8 P3 exists to decide, and deciding it from P3's curves is what a
probe is for. The v0.5 sentence "may not be chosen after seeing its
curves" said the opposite of §8.2 and is withdrawn. What may not happen
is choosing a gate **inside a frozen design** after seeing that design's
own data; the boundary that keeps those apart is the fresh-seed rule of
§8.2, not a prohibition on the probe learning from itself.

**What each statistic can and cannot be claimed for (R3.3, narrowed by
R4.1).** The classes of §4.6 split a second way, and the split runs
opposite to sensitivity. `D` takes **two orders plus the point
correspondence**, so it is a paired counterfactual and **not computable
from a single causal set**. But the guarded chain vectors are not
strictly single-poset either: their eligibility labels are
coordinate-derived (§4.6.2), so an observer holding one poset cannot
even tell which intervals to evaluate. The full picture is §4.6.2's
table; the operative consequences here:

> `D` is the most sensitive probe and the least claim-bearing. It is a
> **paired-order sensitivity diagnostic**: it answers "does pure-Weyl
> curvature move the order at all, at this density", which is what a
> probe needs. It does not answer "can a finite causal set be told from
> a flat one", and no single-poset claim may rest on it.
>
> That claim rests on **item 3 computed over ALL elements** — the only
> candidate here needing nothing but the poset. Items 2 and 3-on-a-
> diamond are single-ARM diagnostics while their eligibility is
> coordinate-guarded, and are promoted only if §4.6.2's order-invariant
> rule is validated.

**1. Paired relation disagreement (primary probe statistic).** With
`prec_A` and `prec_0` the two relations on the same `N` points,

    D = | prec_A  triangle  prec_0 | / C(N, 2),

reported together with its signed decomposition: relations GAINED
(flat-spacelike pairs pulled inside the curved cone — the focusing
direction's work) and relations LOST (flat-causal pairs pushed out —
the defocusing direction's), separately. Why it leads:

- it uses the same-points pairing directly — it is not even definable
  in a two-sprinkling design;
- it is unsigned at the pair level, so the quadrupole cancellation of
  §4.5 cannot reach it: a flip is a flip whichever transverse direction
  caused it. Expected sensitivity `O(|A|)` where scalars get `O(A^2)`
  — **[TO VERIFY]**, this is the heuristic P3 exists to measure;
- the gained/lost split is the quadrupole made visible: gains should
  concentrate among `y`-separated pairs, losses among `x`-separated.

**2. Chain vectors (secondary) — Class C, single-ARM.** Abundances
`C_2, C_3, C_4, ...` at fixed or conditioned cardinality, normalized,
over §4.6's eligible intervals. The `C_k` are where a small-interval
expansion places curvature information beyond the volume — but the
leading terms it places there are `R` and `R_00` (RSS 2013), both zero
here, so Weyl sensitivity, if any, sits at higher order. Measured second
for exactly that reason. **While its eligibility is coordinate-guarded
this is a single-arm diagnostic, not a pure-poset observable** (§4.6.2);
promotion depends on the order-invariant rule being validated.

**3. Relation fraction, as a discriminator only — and read carefully,
because two different statistics hide under this heading.**

- **3a. Global relation fraction over ALL sprinkled elements.** Class R,
  no eligibility rule, no coordinates: `(# related pairs) / C(N, 2)` on
  one arm's poset. **This is the only strictly single-poset, pure-order
  candidate in the list**, and §7's single-causal-set claim rests on it.
  Its sensitivity is unknown and may well be the weakest here — that is
  the price of being the one thing an observer holding a poset can
  actually compute.
- **3b. Myrheim–Meyer or ordering fraction ON A DIAMOND.** Class C,
  coordinate-guarded, therefore single-arm and not pure-poset, exactly
  as item 2. The MM machinery is characterized in this repository (P6b)
  and cheap, but in a generic curved spacetime its output may NOT be
  interpreted as a flat-ensemble dimension estimate. It runs as a pure
  discriminator — does the number move between arms — with the word
  "dimension" kept out of every sentence that quotes it.

**4. Transverse anisotropy, diagnostic only.** Any statistic that uses
the coordinates `x, y` directly (e.g. conditioning pairs on their
transverse separation axis) sees the quadrupole at `O(A)` but is not a
pure-order quantity. It is a DIAGNOSTIC for understanding the primary
statistics' behaviour, never a primary claim, unless a future design
reconstructs the polarization axes from order data alone — which would
be its own result and its own preregistration.

### 5.1 Undecided pairs, and why `D` is reported as an interval (R2.1)

`D` lives exactly where the predicate is hardest: a flip happens at the
cone, and the cone is where `sigma ~ 0` and floating point runs out. A
fixed exclusion band would therefore remove the signal-carrying
population and could report the removal as a small `D` — a null
manufactured by the tolerance. The construction below makes that failure
mode unrepresentable.

**Definition.** For each pair and **each arm**, compute not a value of
`sigma` but a **verified error interval** for it. If either arm's
interval contains zero, that pair is `ambiguous`: its relation is
undecided at the achieved precision, and no verdict is imputed to it.

**Reporting, and the denominator does not move.** `C(N, 2)` stays the
denominator throughout. Report three numbers together, never one:

    D_lower              = (definite flips) / C(N,2)
    D_upper              = (definite flips + ambiguous) / C(N,2)
    ambiguous_fraction   = (ambiguous) / C(N,2)

`D_lower` assumes every undecided pair agrees between the arms;
`D_upper` assumes every one of them flips. The truth is bracketed by
construction. **If the band eats the signal, that shows up as a wide
`[D_lower, D_upper]`, never as a small `D`.**

**A conditional `D` restricted to decided pairs may be reported and may
not be primary.** It conditions on an outcome correlated with the
quantity being measured — precision fails preferentially near the cone,
which is where flips are — so it is a diagnostic, in the sense item 4
above uses the word.

**When ambiguity is high, buy precision, not tolerance.** The response
to a large `ambiguous_fraction` is adaptive precision or interval
arithmetic on the offending pairs — never a wider fixed band, which
trades a measurable quantity for an invisible one. The escalation policy
and its cost are part of what §8 P1 must report. §4.2 previously offered
a fixed band as an alternative here; it is withdrawn (R3.1), because a
band is an exclusion and an exclusion moves this denominator.

**Nothing else moves the denominator either.** `D` is a Class R
statistic in §4.6's sense: every pair of sprinkled points is eligible,
there is no containment condition, and the conjugate-point bound is a
property of the slab rather than a per-pair test. `C(N, 2)` is therefore
the count of all pairs, without qualification.

The same construction applies to the gained/lost split, which inherits
`ambiguous` from the pair-level classification rather than resolving it
by fiat.

### 5.2 `D` is a U-statistic: the resampling unit is the sprinkling (R2.2)

`D` is an average over `C(N, 2)` pairs drawn from **one** point set, and
those pairs share points. They are not independent samples, and treating
them as such would understate `D`'s standard error by a factor that grows
with `N` — millions of pairs presenting as millions of degrees of
freedom when the actual replication is one sprinkling.

**Therefore: every uncertainty attached to `D`, `D_lower`, `D_upper`,
`ambiguous_fraction`, or the gained/lost split is computed across
independent sprinklings, and never across pairs within a sprinkling.**
The sprinkling is the sampling unit; `N` buys resolution inside one
sample, not replication. Effect size per sample in §8 P3 means per
sprinkling, and the power section that follows sizes the number of
sprinklings.

This is the same defect class as P11's C10 (a statistic that dropped the
twin's variance) and P12's fixed-denominator bootstrap, caught here
before it can be implemented rather than after.

## 6. The interval-versus-point rule, settled before it can be convenient

P12 §9.3 handed over a second, separable lever "available now at no
experimental cost", and this is the design that has to spend it.

The facts: P12's co-requirement (ii) is a **point** comparison — the top
rung's recovery `0.1708` against the threshold `0.25`. The realized 95%
interval reached `0.2644`, i.e. **across** the threshold. The frozen rule
passed; an interval rule would not have. §9.3 refused to re-score:

> Converting a point rule into an interval rule after seeing the interval
> would be the same move as loosening a threshold after seeing the data, in
> the opposite direction and equally forbidden.

That refusal was right, and it leaves a decision owed to the next design.

**ADOPTED 2026-08-01: an accuracy claim gates on the interval, and the
choice is made once, as a house rule, not per design.** The rule now lives
in `AGENTS.md` so it binds every subsequent design, not only P14.

Rationale:

- The threshold in a recovery gate is a claim about accuracy — "this
  instrument recovers `R tau^2` to 25%". A point estimate that clears 0.25
  while its interval reaches 0.2644 does not support that sentence at the
  confidence the rest of the record is quoted at. Every other number in
  §9.3 is reported with an interval; the one number that decides the
  verdict should not be the exception.
- It is directional and must be applied in the direction that hurts: an
  interval rule is **stricter**, so adopting it cannot be a way to rescue a
  marginal result later.

### 6.1 What it would have cost P12 — break-even is not the design number (R2.3)

The rule was adopted with its price on the table, read out of
`p12_stage_b_summary.json`. **The v0.3 draft then quoted the wrong price**:
it gave the break-even `n` and called it what a design must size for.
Review R2 caught the conflation, and the correction is a factor of three,
not a rounding.

Realized inputs: recovery `0.1708`, 95% CI `[0.0784, 0.2644]`, effective
`SE = 0.0474` (the record's `1.67 sigma` margin against `t = 0.25`), at
`n = 1117` per rung per arm.

**Break-even** is the `n` at which the realized point estimate's interval
just touches the threshold: `SE <= (t - theta)/1.96 = 0.0404`, giving
`n = 1538`, i.e. `1.38x`. A design sized there passes **half the time** —
by construction, since the interval lands exactly on the line. It is a
description of what happened, never a sizing target.

**Power-sized** is what a preregistration owes. P12's own convention is
90% (`N_SUP_COEFF` carries `1.960 + 1.282`, and the equivalence slot's
`1.645` is already the two-sided quantile for `beta = 0.10`), so applying
the house convention to the recovery gate:

| target power | `SE` needed | `n` per rung per arm | vs the `1117` run |
|---|---|---|---|
| break-even (**50%**) | `0.0404` | `1538` | `1.38x` |
| 80% | `0.02827` | `3141` | `2.81x` |
| **90% (house convention)** | `0.02443` | **`4204`** | **`3.76x`** |
| 95% | `0.02197` | `5200` | `4.66x` |

Verified by direct simulation, not by the closed form alone: pass rates
`0.499 / 0.801 / 0.901` at `n = 1538 / 3141 / 4206`.

**Under this rule P12 would have declared itself INFEASIBLE**, and by a
wide margin rather than a narrow one: its frozen cap was `N_CAP = 1300`,
which break-even already exceeds and the 90%-power number exceeds **3.2
times over**. Not because the result was wrong, but because it could not
afford the sentence it wanted to write. That is the rule working, and it
is stated here — with the right number — so no later design adopts it
believing it is cheap.

Two caveats on the arithmetic, since this is a design-stage estimate:
`SE ~ 1/sqrt(n)` and normal-theory intervals stand in for what is
actually a bootstrap CI, and the power figures are conditioned on the
true value sitting at P12's realized `0.1708`. A real design substitutes
its own design value and its own interval construction. Neither caveat
touches the finding: break-even is 50% power, and the gap to a house
90% is a factor of `2.7` in `n` on top of it.

Note also where P12's `n` came from: `n_sup = 101`, `n_eq = 1117`, both
from the **rate** gate. The recovery gate's precision was never a sizing
constraint and simply came along. Under the interval rule it becomes one,
and on these numbers it becomes the **binding** one by a large factor.
Every design from here sizes its accuracy gate for **power**, states the
target explicitly, and may not quote a break-even `n` as if it were one.

### 6.2 Why the cost is smallest exactly where the rule first binds

P12's two arms were **separate sprinklings on disjoint seed blocks**, so
the ratio `Q_hat` carries both arms' fluctuations independently.

In this design the two arms are the **same points** (§4.1). The common
fluctuation cancels in the paired difference, so the interval should be
materially tighter at equal `n` — which is to say the design that first has
to pay for this rule is the one where it is cheapest.

**[TO VERIFY]** That is a reasoned expectation, not a measurement. §8 P4
must report the paired variance directly, **between sprinklings** per
§5.2; if the pairing does not deliver, the power section pays the full
`3.76x` of §6.1 — the 90%-power figure, not the `1.38x` break-even — and
says so.

**This is not retroactive.** P12's record stands as written, on the point
rule it froze. §9.3 already says so and this document does not reopen it.

## 7. What this design does NOT claim, written before there is anything to protect

Per §10.9's discipline, stated now rather than after a verdict:

- **A positive `D` supports less than the v0.4 draft said (R3.3).** The
  draft claimed "a finite causal set carries a detectable pure-Weyl
  signature", which `D` cannot support: `D` consumes two orders and the
  point correspondence between them, so it is a paired counterfactual
  and not a property of any single causal set. What a positive `D`
  supports is: **pure-Weyl curvature moves the causal order measurably
  at the tested density, patch and profile** — the existence and size of
  the effect, not its accessibility to an observer holding one poset.
- **The single-causal-set claim rests on §5 item 3a alone (R4.1)** — the
  global relation fraction over all elements, the one candidate needing
  nothing but the poset. v0.5 assigned this claim to items 2 and 3
  together; that was wrong, because their eligibility labels are
  coordinate-derived (§4.6.2), so an observer holding one poset cannot
  even tell which intervals to evaluate, let alone evaluate them. Only
  3a may be quoted as "a finite causal set can be told from a flat one",
  and only at the effect size it itself achieves — which may be the
  weakest in the list. **That narrowing is a real cost of this design
  and is not hidden**: the most sensitive statistic (`D`) carries no
  single-poset claim, and the statistic that carries it may be too blunt
  to see anything.
- **The guarded chain vectors and diamond-based estimators (items 2 and
  3b) are single-ARM diagnostics** while their eligibility is
  coordinate-guarded. They may be reported, compared between arms, and
  used to understand `D` — they may not be quoted as what a causal-set
  observer could do. §4.6.2 names the order-invariant eligibility rule
  that would promote them, and makes its validation a probe task rather
  than an assumption.
- A design that finds `D` positive, 3a null, and the guarded statistics
  positive has learned three separate things and must report them in
  three sentences, not one. That combination is not a failure; it is the
  finding that the effect exists, that it is not reachable by the one
  unguarded pure-order statistic tried, and that an order-invariant
  eligibility rule is where the next design should look.
- Neither establishes quantitative recovery of the Weyl tensor —
  detection and recovery are different claims, and P12 needed a separate
  stage for the second.
- A **null** result is a statement about the tested statistics and
  densities, never about the continuum: Malament's theorem stands either
  way, and a null must not be phrased as "curvature does not enter the
  order". The honest phrasing is P10's — the question may be undecidable
  by this instrument family at affordable density — and §8.1's S1 price
  attaches so the closure names what the alternative would have cost.
- Plane waves are special twice over: Petrov type N, and `sqrt(-g) = 1`.
  The exactness of the measure is what makes the experiment clean and
  what makes it unrepresentative of generic `Weyl != 0` spacetimes.
  Generalizing beyond type N is a different design (§8.1).
- The `O(|A|)`-vs-`O(A^2)` hierarchy of §5 is a heuristic for `D` versus
  generic scalars until P3 measures it, and stays labelled that way. The
  axis-volume case is the exception and is no longer a heuristic: its
  `A^2` response has a closed coefficient, `V_A/V_0 = 1 + (wT)^4/252 +
  O((wT)^8)` (§4.4). If the primary statistic underperforms the
  heuristic, that finding is reported, not absorbed.
- A wide `[D_lower, D_upper]` is a **precision** result, not a physics
  result, and may not be reported as either a detection or a null. §5.1
  exists so the two cannot be confused, and §8 P1 must show the
  ambiguous fraction is small enough for the bracket to decide anything
  before P3's numbers mean what they appear to.

## 8. What comes next, and it is a probe, not a stage

This repository has a precedent for exactly this situation:
`docs/p7_fss_rescope_recon.md`, where a reconnaissance closed the FSS
rescope path in about thirty seconds of compute by asking whether the
instrument worked at the sizes the sampler could reach. The decision
record notes what killed it: "**죽인 것은 컴퓨트가 아니라 계측기 물리**".

P14 has the same shape and is answered the same way, **before any prereg
freeze and with nothing frozen**:

- **P0. Measure.** Confirm the uniform-coordinate sprinkling reproduces
  the true 4-volume in the slab (it must, given `sqrt(-g) = 1`; this is
  an implementation check, not a physics one).
- **P1. Causality, and the two admissibility rules.** Implement
  `p prec q` in the constant-`A` slab from the world function and check
  against direct null-geodesic integration to a pinned tolerance. Return
  a **verified error interval for `sigma` per arm**, not a value, so
  §5.1's undecided set is well defined; report `ambiguous_fraction` and
  the precision-escalation policy with its cost. Assert the frozen slab
  extent is under `pi/sqrt(A)` so no pair needs a conjugate-point filter
  (§4.3). Derive and report §4.6's **Class C guard inset** in closed
  form, and the eligible-pair fraction it leaves — Class R keeps every
  pair and `C(N, 2)`. Also measure §4.6.2's candidate **order-invariant
  eligibility rule** against that coordinate guard: report the agreement
  rate and the pairs each admits that the other rejects, since that
  number decides whether §7's single-poset claim has anything better
  than item 3a available to it.
- **P2. The volume prediction, with the right error model (R3.2).**
  Axis-pair interval cardinalities in the curved arm must match the
  `V_A/V_0` quadrature of §4.4 (`1.00400047` at `wT = 1`, `1.07300802`
  at `wT = 2`), over Class C eligible pairs only. This replaces the dead
  C0: a number, not a silence.

  **"Within Poisson error" was underspecified and, read naively, wrong.**
  `N_A` and `N_0` are counts over two OVERLAPPING regions of the SAME
  realization, so they are correlated. For a Poisson process the
  covariance is exact:

      Cov(N_A, N_0) = rho * Vol(I_A  intersect  I_0)

  Two things are reported, and R4.3 fixes how each is stated:

  1. **Marginal check per arm.** Each arm's count is compared against
     its own known mean, `rho V_A` and `rho V_0`, under its **marginal**
     Poisson law. Note the word: the marginals are Poisson, the two
     counts are **not independent of each other**, and v0.5's "independent
     per arm" was the wrong phrase for the right check. Each comparison
     stands on its own marginal and neither borrows the other's; this is
     the check that the sprinkling and the volumes are right.
  2. **Paired residual — the primary check, and the ratio is derived
     from it (R5.1).** v0.6 estimated the uncertainty of the per-sample
     ratio `N_A / N_0` between sprinklings. That is the wrong estimand:
     for a finite Poisson denominator `E[N_A / N_0] != V_A / V_0`, and
     `N_0 = 0` has positive probability, so the per-sample ratio is
     biased and more sprinklings do not remove it. The design already
     contained the right quantity and used it only as a cross-check.
     Promoted:

     - **Primary.** Per sprinkling `i`, with the predicted ratio
       `r = V_A / V_0`, form

           Z_i = N_A,i - r * N_0,i

       and test `E[Z] = 0`. `Z` is exactly zero-mean under the
       prediction, is defined whatever `N_0` is, and needs no
       denominator.
     - **Empirical uncertainty.** The between-sprinkling variance of
       `Z` (§5.2's replication unit, unchanged).
     - **Analytic cross-check, now of the same quantity the test uses:**

           Var(Z) = rho ( V_A + r^2 V_0 - 2 r Vol(I_A ∩ I_0) )

       Because the primary test is `Z` rather than a ratio, the analytic
       and empirical numbers are two estimates of ONE variance instead
       of two descriptions of different things. `Var(N_A - N_0) =
       rho Vol(I_A triangle I_0)` is exact for the raw difference and is
       NOT this check, which is what v0.5 got wrong.
     - **Ratio, reported secondarily.** As a pooled `sum N_A / sum N_0`
       or an explicit ratio-of-means, never as a per-sample average of
       ratios, and never as the primary test.
- **P3. Discriminability.** For a ladder of `A` and slab sizes, measure
  the §5 statistics in order — `D` first, reported as
  `[D_lower, D_upper]` with `ambiguous_fraction` and the gained/lost
  split (§5.1) — on the same point sets, over **independent sprinklings**
  as the replication unit (§5.2). Effect size per sample means per
  sprinkling. **No gates, no verdicts** — characterization, in the sense
  P10 §6 uses the word. This is also where the `O(|A|)` heuristic for
  `D` and the `O(A^2)` suppression of scalars get measured.
- **P4. Paired variance.** Report the between-sprinkling variance of the
  paired difference directly, so §6.2's expectation is measured rather
  than assumed and the power section — which sizes a NUMBER OF
  SPRINKLINGS, not a number of pairs — knows what the interval rule
  actually costs here.

### 8.1 The second probe: what a general `Weyl != 0` spacetime would cost

Decided 2026-08-01: the plane-wave probe runs, **and so does a separate,
narrower measurement on Schwarzschild.** The reason is §7's honest limit —
plane waves are Petrov type N, and a null result there cannot distinguish
"undetectable at this density" from "this family is too special". Leaving
the generalization path unpriced until after the first probe would repeat
the mistake the Wang-Landau record names: porting a kernel before
measuring the exponent that decides whether the port is worth anything.

This probe does **not** run a campaign and does not need one. It answers
one question:

- **S1. Cost of causality.** Implement `p prec q` in a Schwarzschild
  patch by null-geodesic integration and measure the wall-clock per pair
  and its scaling. Compare against the sprinkling sizes P11–P13 operated
  at (`n` in the low thousands per rung, so `O(n^2)` relations per
  sample). Report a price, not a verdict.

Two things follow from S1 regardless of P3's outcome. If the price is
affordable, Schwarzschild becomes the generalization stage and the
type-N limitation is a temporary one. If it is not, that number is what
closes the general-`Weyl` path for now, in the same quantified way the
tunneling exponent closed Wang-Landau — and a plane-wave null result then
has to be reported with a limitation that is **priced**, not merely named.

The Section 3 obstacle stands and S1 does not resolve it: Schwarzschild
has no closed-form diamond volume, so a campaign there would still have
to confront the small-diamond expansion that P12 §10.1 refused. S1
measures the causality cost only; the volume question is separate and
unaddressed.

### 8.2 Exploration and confirmation do not share seeds (R3.4)

P0 through P4 and S1 are **exploratory**: they may look at every curve,
compare candidate statistics, and pick the primary statistic, the
operating point, and the ambiguity policy. The freeze that follows may
set gates from what they found — that is what a probe is for, and §5 now
says so rather than contradicting it (R4.2).

**What the probe does NOT choose:** the guard inset. Per §4.6.1 it is a
closed-form function of the operating point, so fixing the point fixes
the inset; it is determined, not tuned. v0.5 listed it here as a probe
choice, which is why §4.6 and this section disagreed.

**The line between exploration and confirmation is the seed block, and
it is fixed here rather than left to the freeze:**

> **Confirmation runs on a fresh seed block disjoint from every seed the
> probe consumed.** No sprinkling that informed the choice of statistic,
> operating point, ambiguity policy, or threshold may appear in the
> confirmatory sample.

The repository already has machinery of exactly this shape, though it is
not one shared helper (R4.4): **P12's runner** calls
`assert_windows_disjoint_and_fresh`, which aborts on a window collision;
**P13** enforces the same property through its own frozen window
constants and tests, and uses `supports_disjoint` for the separate job
of packing disjoint boxes within a sample. The eventual preregistration
inherits the property and should inherit P12's helper rather than
writing a third variant. Recording the rule now is what keeps the
probe's freedom to explore from quietly becoming a licence to tune and
confirm on the same data.

### 8.3 Exit conditions

If P3 finds no separation at achievable sizes, P14 closes there, cheaply
and honestly — phrased as instrument-and-density scope per §7, with S1's
price attached so the closure records what the alternative would have
cost rather than implying none existed. If it does separate, P3's effect
sizes and P4's paired variance are what the power section of the actual
preregistration is built from — power-first, as P11 through P13 all
were, and now sized for an interval-rule accuracy gate per §6.

## 9. Design checks, committed and run

Both scripts live in `docs/prereg/p14_checks/`, run by hand at draft
stage, abort on failure rather than recording it, and pin every number
this document quotes so prose and computation cannot drift (the C25/C28
lesson applied at draft stage). Frozen, stamped copies move to
`docs/prereg/frozen/p14/` if and when this design freezes.

**Check 1 — `p14_brinkmann_vacuum_check.py`** (sympy, exact, `A(u)`
undetermined):

```
det g              = -1                     [exact]
sqrt(-g)           =  1                     [exact]
Ricci              =  0 (all 16 components) [exact]
Riemann != 0       : 8 components, all +/-A(u)
=> vacuum, curvature pure Weyl, volume form exactly flat: PASS
```

**Check 2 — `p14_interval_volume_constant_a.py`** (numpy quadrature;
derivation in the script; asserts the flat limit, both pinned values,
and the super-quartic growth):

```
flat-limit check (wT=1e-3): 1.000000000000
wT=1: V_A/V_0 = 1.00400047   (pinned 1.00400047)
wT=2: V_A/V_0 = 1.07300802   (pinned 1.07300802)
excess ratio (quartic => 16): 18.25
leading coefficient -> 1/252 = 0.003968253968, next term (wT)^8 with coefficient ~3.22e-05
=> the volume channel is NOT silent; C0 as drafted is dead: PASS
```

**Not a committed check, and named so it is not mistaken for one:** the
break-even-versus-power figures of §6.1 were computed and simulated at
draft stage but are arithmetic on P12's frozen artifact, not a property
of this design. They belong to the power section of the eventual
preregistration, which must derive them for its own design value and its
own interval construction rather than inheriting P12's.

## 10. Open questions this draft does not answer

1. Which pure-order statistic is Weyl-sensitive at affordable density,
   and whether `D`'s expected `O(|A|)` advantage over scalars is real.
   **Assigned:** §8 P3.
2. Whether the gained/lost decomposition of `D` cleanly tracks the
   polarization quadrupole, and whether the axes can ever be
   reconstructed from order data alone (which would promote §5's
   diagnostic to a claim).
3. The cost of deciding causality per pair. **Assigned:** §8 P1 for the
   plane wave, §8.1 S1 for Schwarzschild.
4. What `ambiguous_fraction` actually is at campaign densities. §5.1
   makes a large one visible rather than silent, but a bracket
   `[D_lower, D_upper]` too wide to decide anything is still a dead
   probe, and the escalation to interval arithmetic has a cost nobody
   has measured. **Assigned:** §8 P1.
5. The Class C guard inset of §4.6 in closed form, and what fraction of
   pairs it leaves eligible. A thin retained population is a cost of the
   **single-arm diagnostics** (§5 items 2 and 3b) and of §8 P2's
   implementation check — **not** of §7's single-poset claim, which
   rests on item 3a and needs no eligibility rule at all (§4.6.2, R4.1).
   The claim-side lever is q9, not this one. **Assigned:** §8 P1.
6. Whether the same-points pairing tightens the interval enough to
   offset the interval rule's cost, measured **between sprinklings**
   (§5.2). The benchmark it must beat is §6.1's `3.76x` at 90% power,
   not the `1.38x` break-even. **Assigned:** §8 P4.
7. What the diamond volume is in a Schwarzschild patch, if S1 says the
   causality cost is affordable. §3's standing refusal of unverified
   expansion coefficients applies, and this draft has no answer.
8. Whether **item 3a** — the global relation fraction over all elements,
   now the only strictly single-poset candidate (§4.6.2) — reaches any
   useful sensitivity. §7 rests the only observer-relevant claim on it
   alone, so a result where `D` separates and 3a does not is a live and
   reportable outcome rather than a failure.
9. Whether Class C eligibility can be defined **order-invariantly**, and
   how closely such a rule agrees with the coordinate guard. This is the
   lever that would return the chain vectors and diamond-based
   estimators to pure-poset standing and hand §7 a far more sensitive
   instrument than 3a. Validation is measurement against the coordinate
   guard, not argument. **Assigned:** §8 P1.

## References named by reviews R1 and R2

- D. Malament, J. Math. Phys. 18, 1399 (1977); S. W. Hawking, A. R. King,
  P. J. McCarthy, J. Math. Phys. 17, 174 (1976) — causal structure
  determines the conformal class for distinguishing spacetimes.
- M. Roy, D. Sinha, S. Surya, "Discrete geometry of a small causal
  diamond", arXiv:1212.0631 — small-diamond expansions; leading
  curvature corrections are `R` and `R_00`.
- A. I. Harte, T. D. Drivas, Phys. Rev. D 85, 124039 (2012) — plane-wave
  world function and Jacobi propagators; constant-`A` conjugate points.
- J. L. Flores, M. Sánchez, "Causality and conjugate points in general
  plane waves", gr-qc/0211086 — the quadratic family is strongly causal
  and not globally hyperbolic.
- A. A. García, F. W. Hehl, C. Heinicke, A. Macías, "The Cotton tensor
  in Riemannian spacetimes", gr-qc/0309008 — the `n = 3` obstruction.
