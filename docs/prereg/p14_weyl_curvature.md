# P14: which pure-order statistic reads a pure-Weyl deformation, and at what density? — design draft

Status: **DESIGN DRAFT v0.3 (2026-08-01). NO STOCHASTIC PROBE OR CAMPAIGN
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
tolerance. Two requirements R1.5 adds, both absolute:

- **Near-null ambiguity band.** "Machine precision for arbitrary pairs"
  is too strong a promise: pairs with `sigma ~ 0` sit inside floating-
  point ambiguity, and the implementation must either resolve them with
  interval arithmetic or EXCLUDE them by an explicit, published band —
  never silently misclassify. The band's width is a frozen constant of
  the eventual design, and both arms use the same band.
- **Containment.** For every pair entering any statistic, both arms'
  causal intervals must lie inside the sprinkled box, verified per pair
  — the analogue of P12's "causal box inside the patch" check. A
  truncated diamond is a biased count with no flag otherwise.

### 4.3 Conjugate points: the analytic bound, and what does not escape it

**[VERIFIED — analytic]** For constant `A = w^2` the focusing direction's
Jacobi propagator is `sin(w Du)/w`, so the first conjugate point sits at
exactly

    |Du| = pi / w = pi / sqrt(A)

(Harte & Drivas 2012 give the same). The patch is therefore an explicit
**slab in `u`**: all sprinkled pairs satisfy `|Du| < pi/sqrt(A)` with a
margin to be frozen, and the containment check of §4.2 enforces it
per pair. What remains [TO VERIFY] is only the implementation against
this known answer, not the constant itself.

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
- **What replaces it is stronger: a quantitative prediction.** The
  curved arm's mean interval cardinality over axis-pair diamonds must
  match `V_A/V_0` at the sprinkled density within Poisson error. That
  is §8 P2 — an implementation check with a number on it, not a silence
  test. The flat arm needs no control at all: it IS the reference, on
  the same points.
- **The excess is quadratic in `A`** (quartic in `wT`; the committed
  check asserts the growth), consistent with the small-diamond
  expansion's leading corrections being `R` and `R_00` (Roy, Sinha &
  Surya 2013, arXiv:1212.0631) — both identically zero here. The volume
  channel responds at second order. This is exactly why §5 does not
  lead with a volume-type statistic.

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

## 5. What the instrument reads, in probe order (R1.3)

No gate is proposed for any of these here. Which one, at what size, with
what threshold, is what §8 P3 measures and what may not be chosen after
seeing its curves.

**1. Paired relation disagreement (primary).** With `prec_A` and
`prec_0` the two relations on the same `N` points,

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

**2. Chain vectors (secondary).** Abundances `C_2, C_3, C_4, ...` at
fixed or conditioned cardinality, normalized. The `C_k` are where a
small-interval expansion places curvature information beyond the volume
— but the leading terms it places there are `R` and `R_00` (RSS 2013),
both zero here, so Weyl sensitivity, if any, sits at higher order.
Measured second for exactly that reason.

**3. Relation fraction / ordering fraction, as a discriminator only.**
The Myrheim–Meyer machinery is characterized in this repository (P6b)
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

### 6.1 What it would have cost P12, computed rather than asserted

The rule was adopted with its price on the table, read out of
`p12_stage_b_summary.json` rather than estimated:

| | value |
|---|---|
| realized top-rung recovery | `0.1708` |
| realized 95% CI | `[0.0784, 0.2644]` |
| effective standard error | `0.0474` (record: margin `1.67 sigma`) |
| SE needed for the interval to clear `0.25` | `0.0404` |
| implied `n` multiplier (`SE ~ 1/sqrt n`) | `1.38x` |
| implied `n` per rung per arm | **`1538`** |
| P12's frozen cap `N_CAP` | `1300` |

**Under this rule P12 would have declared itself INFEASIBLE** — not because
the result was wrong, but because it could not afford the sentence it
wanted to write. That is the rule working, and it is stated here so no
later design can adopt it believing it is free.

Note also where P12's `n` came from: `n_sup = 101`, `n_eq = 1117`, both
from the **rate** gate. The recovery gate's precision was never a sizing
constraint and simply came along. Under the interval rule it becomes one,
and on these numbers it becomes the **binding** one. Every design from here
sizes for its accuracy gate explicitly.

### 6.2 Why the cost is smallest exactly where the rule first binds

P12's two arms were **separate sprinklings on disjoint seed blocks**, so
the ratio `Q_hat` carries both arms' fluctuations independently.

In this design the two arms are the **same points** (§4.1). The common
fluctuation cancels in the paired difference, so the interval should be
materially tighter at equal `n` — which is to say the design that first has
to pay for this rule is the one where it is cheapest.

**[TO VERIFY]** That is a reasoned expectation, not a measurement. §8's
probe must report the paired variance directly; if the pairing does not
deliver, the power section pays the full `1.38x` and says so.

**This is not retroactive.** P12's record stands as written, on the point
rule it froze. §9.3 already says so and this document does not reopen it.

## 7. What this design does NOT claim, written before there is anything to protect

Per §10.9's discipline, stated now rather than after a verdict:

- A positive result supports: **a finite causal set at the tested density
  carries a detectable pure-Weyl signature in the tested pure-order
  statistic, over this patch and profile.** It does not establish
  quantitative recovery of the Weyl tensor — detection and recovery are
  different claims, and P12 needed a separate stage for the second.
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
- The `O(|A|)`-vs-`O(A^2)` hierarchy of §5 is a heuristic until P3
  measures it. If the primary statistic underperforms it, that finding
  is reported, not absorbed.

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
- **P1. Causality.** Implement `p prec q` in the constant-`A` slab from
  the world function; check against direct null-geodesic integration to
  a pinned tolerance; publish the near-null ambiguity band and apply it
  identically in both arms; verify per pair that both arms' intervals
  are box-contained and inside `|Du| < pi/sqrt(A)` (§4.2, §4.3).
- **P2. The volume prediction.** Axis-pair interval cardinalities in the
  curved arm must match the `V_A/V_0` quadrature of §4.4 within Poisson
  error (`1.00400047` at `wT = 1`, `1.07300802` at `wT = 2`). This
  replaces the dead C0: a number, not a silence.
- **P3. Discriminability.** For a ladder of `A` and slab sizes, measure
  the §5 statistics in order — `D` with its gained/lost split first —
  on the same point sets, and report effect size per sample with
  uncertainties. **No gates, no verdicts** — characterization, in the
  sense P10 §6 uses the word. This is also where the `O(|A|)` heuristic
  for `D` and the `O(A^2)` suppression of scalars get measured.
- **P4. Paired variance.** Report the variance of the paired difference
  directly, so §6.2's expectation is measured rather than assumed and
  the power section knows what the interval rule actually costs here.

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

### 8.2 Exit conditions

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
=> the volume channel is NOT silent; C0 as drafted is dead: PASS
```

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
4. The near-null ambiguity band: how wide it must be for correctness,
   and what fraction of pairs it excludes at campaign densities — a
   band that eats the near-cone population would eat `D`'s signal with
   it. **Assigned:** §8 P1.
5. Whether the same-points pairing tightens the interval enough to
   absorb the interval rule's cost, or whether the power section pays
   the full `1.38x` of §6.1. **Assigned:** §8 P4.
6. What the diamond volume is in a Schwarzschild patch, if S1 says the
   causality cost is affordable. §3's standing refusal of unverified
   expansion coefficients applies, and this draft has no answer.

## References named by review R1

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
