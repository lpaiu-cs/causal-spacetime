# P14: can curvature reach the causal order without passing through the measure? — design draft

Status: **DESIGN DRAFT v0.2 (2026-08-01). NOTHING IS FROZEN AND NOTHING HAS
BEEN RUN.** This document argues a spacetime choice and a discriminability
probe. It deliberately does not state gates, thresholds, sample sizes, or
seed windows: those belong in a freeze that follows Section 8's probe, and
writing them now would fix numbers against a construction whose
discriminability is still unmeasured.

Two decisions taken 2026-08-01, both recorded here with their price:

- **Both probes run** (§8, §8.1). The plane wave carries the
  discriminability question; a separate, narrow Schwarzschild measurement
  prices the generalization path so a type-N null result is not reported
  as if no alternative existed.
- **The interval rule is adopted** (§6), house-wide in `AGENTS.md`, with
  the computation showing it would have made P12 itself infeasible under
  P12's own cap (§6.1).

Dates are local (UTC+9, Asia/Seoul).

Claims below are labelled. **[VERIFIED]** means checked in this document by
symbolic computation or by citation to a frozen artifact in this
repository. **[TO VERIFY]** means it must be established numerically before
any freeze, on the Section 3 rule that killed Stage C v1 — an expansion
coefficient nobody checked.

---

## 1. The question, and why it is the binding one

P12 and P13 closed the curvature track on 2026-07-31 and handed over the
same lever from two directions.

- **P12 Stage B-12** (`9.3`): the unchanged P11 estimator **recovers the
  curvature itself** — `Delta_B = -0.2637`, CI `[-0.3037, -0.2233]`,
  top-rung dimensionless recovery `0.1708 <= 0.25`. RECOVERS-CURVATURE.
- **P13 Stage A-13C** (`13.3`): the flat-normalized reading **does not
  notice** curvature out to `tau/ell = 1.5`, i.e. `R tau^2 = 4.5`.
  CURVATURE-ROBUST, control CLEAN.

Both records then refuse to promote the programme's general thesis, and for
the same structural reason, which each states in its own words:

> **in 1+1D every metric is conformally flat**, so `dS_2`'s causal order is
> conformally equivalent to the flat twin's by construction, and the
> curvature can only reach the order through the conformal factor — that
> is, through the volume and hence the count. (P13 §13.3)

The thesis under test is not "curvature is recoverable". P12 already
settled that, in a setting where it could hardly have failed. The thesis is
the **exclusivity** clause:

> Does curvature enter the causal order ONLY through the counting measure?

In 1+1D that is close to a theorem about the geometry, not a measurement
about the instrument. It becomes a question — becomes falsifiable — exactly
where the Weyl tensor is non-trivial.

## 2. `d >= 3` is necessary and NOT sufficient

The handover sentence in P12 §9.3 says "the `d >= 3` question". Read as an
instruction to raise the dimension, it is a trap, and the trap is expensive
because it looks like progress.

**[VERIFIED — by definition]** Every maximally symmetric spacetime has
vanishing Weyl tensor in every dimension. Minkowski, de Sitter and
anti-de Sitter are maximally symmetric. FRW spacetimes are likewise
conformally flat. So `dS_4` sits under exactly the same ceiling as `dS_2`:
moving P12's construction to four dimensions buys nothing but cost.

The requirement is not dimension. It is:

> **Weyl != 0.**

Dimension enters only because `Weyl == 0` identically for `d <= 3`. `d >= 3`
is the precondition for the requirement to be satisfiable, not the
requirement.

## 3. What P12's affordability actually rested on

This is the part that makes P14 hard, and it is the same fact seen twice.

`p12_curved.py` opens by stating why the sprinkling was cheap:

> the conformal factor cannot alter causal structure

Because `dS_2` is conformally flat, its causal order **is** Minkowski's.
P12 decided every causal relation algebraically in null coordinates, and
got the diamond volume in closed form (§10.1, `Vol = 4 ell^2 ln cosh(tau /
2 ell)`, verified against quadrature at `1.0e-14`).

Both of those were gifts of conformal flatness. Giving up `Weyl == 0` gives
up both at once:

| | P12 (`dS_2`) | naive P14 (e.g. Schwarzschild) |
|---|---|---|
| causal relation | algebraic, flat null cones | integrate null geodesics per pair |
| diamond volume | exact closed form | no closed form known |
| twin arm | separate flat sprinkling, calibrated to matched `m` | ? |

And §10.1 records the standing rule about the obvious escape:

> Section 5 refused the small-diamond expansion because its coefficient is
> "exactly the kind of unverified constant that killed Stage C v1".

So P14 may not simply expand the volume in powers of `tau` and read a Weyl
coefficient off the series. Whatever it uses must be exact, or verified
against an independent computation to a pinned tolerance.

**A design is only affordable here if it recovers, in a Weyl != 0
spacetime, the two things conformal flatness used to give away free.**
Section 4 is a construction that does.

## 4. The construction: a vacuum plane wave, where the measure is exactly flat

Take a pp-wave in Brinkmann form,

    ds^2 = 2 du dv + H(u, x, y) du^2 + dx^2 + dy^2,     H = A(u) (x^2 - y^2).

Three properties, all **[VERIFIED]** symbolically in this document
(`sympy`, exact — see Section 9 for the check and its output):

1. **`det g = -1`, hence `sqrt(-g) = 1` exactly.** The volume element is
   Minkowski's, identically, everywhere, for every profile `A(u)`.
2. **`R_{mu nu} = 0` and `R = 0` exactly.** Vacuum. `H` harmonic in the
   transverse plane is what buys this; `x^2 - y^2` is the simplest choice.
3. **`R^a_{bcd} != 0`** — eight non-vanishing components, all proportional
   to `A(u)`. With Ricci zero, the curvature is **pure Weyl**.

Property 1 is the one that makes the experiment clean, and it is stronger
than the Ricci-flatness argument it supersedes. The usual reasoning is
"vacuum, so the small-diamond volume agrees with flat to leading order in
`tau`". Here we do not need leading order and we do not need a small
diamond:

> **A uniform-coordinate Poisson sprinkling IS a Poisson sprinkling of the
> true 4-volume, exactly, at any size.** The counting measure is
> Minkowski's by construction, not approximately and not asymptotically.

Therefore **any** difference the causal order shows between `A != 0` and
`A = 0` is, necessarily, curvature that did not pass through the measure.
That is the thesis's exclusivity clause, stated as a two-arm comparison
with nothing left to calibrate.

### 4.1 The control is the same points, read with `A = 0`

This deserves its own line because it removes the machinery that cost P12
three review rounds.

Since `sqrt(-g) = 1` for every `A`, the curved arm and the flat arm can be
**the same sprinkled points**. The two ensembles are not merely matched in
density, or calibrated to matched mean interval cardinality — they are
*identical as point sets*. Only the light cones differ.

Consequences, and they are large:

- No twin intensity calibration (P12 §10.6). Nothing to converge.
- No cross-arm `m` gate, and no cross-arm offset to disclose, because a
  systematic offset between the arms is not representable: the arms differ
  by the causal relation alone.
- Perfect pairing. Every statistic can be computed as a within-point-set
  difference, which is the strongest form of the paired design P12 §10.3
  argued for on much weaker footing.

The residual bias mismatch that P12 §10.6 had to price, and that §9.3
qualification (3) had to report as a limit of the construction, **does not
arise here.**

### 4.2 Geodesics, and therefore causality, stay tractable

**[VERIFIED — standard, and re-derived from the Lagrangian]** In Brinkmann
coordinates `v` is cyclic, so `u` is affine along geodesics. With `u` as
parameter the transverse equations are **linear**:

    x'' = A(u) x,        y'' = -A(u) y

Hill's equation. For constant `A` they are elementary (`cosh`/`cos`), and
for any profile they are a linear ODE, not a geodesic shooting problem in
four dimensions. This is why the plane-wave family is the one place the
`Weyl != 0` requirement and computable causality are simultaneously
available.

**[TO VERIFY]** The step from "geodesics are solvable" to "the causal
relation `p precedes q` is decidable to machine precision, for arbitrary
pairs, at the cost of a sprinkling campaign" is the load-bearing
implementation question. Plane waves have an explicitly known world
function in the literature; that must be reproduced and checked here
against direct numerical geodesic integration before anything freezes,
exactly as §10.1 checked the closed-form volume against quadrature.

### 4.3 The known defect of this family, and the constraint it imposes

**Plane waves are not globally hyperbolic** (Penrose). Conjugate points
form: for constant `A > 0` the `y` equation oscillates with period
`2 pi / sqrt(A)`, and focusing appears at `u`-separations of order
`pi / sqrt(A)`.

This is not fatal, it is a patch constraint, and it is quantitative:

    |Delta u| < pi / sqrt(A)      (first conjugate point)

**[TO VERIFY]** the exact constant and the correct condition for the chosen
profile. The design must state the patch and demonstrate that no sprinkled
pair straddles a conjugate point — the analogue of P12's "causal box inside
the patch" check, which `p12_curved.py` already enforces per pair.

The tension this creates is the design's central trade: `A` large enough
that the order is measurably deformed, small enough that the patch holds
enough points. Section 8 measures it rather than guessing it.

## 5. What the instrument must read, and the built-in negative control

P12's instrument reads `g_hat = m_open / (rho_hat * tau_hat^2)` — interval
cardinality over count-calibrated volume — and inverts
`G(s) = ln cosh(s)/s^2` to get the curvature.

In this construction **that instrument must read zero curvature**, because
the volume it normalizes by is exactly flat. That is not a problem; it is
the experiment's negative control, and it is free:

> **Control C0.** The P12/P13 volume-based reading, applied unchanged in
> the plane-wave patch, must return the flat answer. If it does not, either
> the sprinkling or the instrument is wrong, and no positive result may be
> reported until that is resolved.

The positive quantity must therefore be an order statistic that is **not** a
function of the interval volume. This is the genuinely open part of the
design and the document will not pretend otherwise. Candidates, in the
order Section 8 should probe them:

1. **Ordering fraction / Myrheim–Meyer dimension estimate** on the interval.
   Already implemented in this repository, already characterized on flat
   and near-critical ensembles (P6b), and it is a pure order statistic.
2. **Chain and antichain abundances** beyond the 2-chain, which is the
   volume. If the causal order carries Weyl information at all, the higher
   abundances are where a small-interval expansion would place it.
3. **Transverse anisotropy of the interval**, which is the one that matches
   the physics: `H = A(u)(x^2 - y^2)` focuses in one transverse direction
   and defocuses in the other, so the deformation is a quadrupole, not a
   scalar. A scalar statistic may be blind to it by symmetry — a real risk
   that Section 8 must test before a scalar gate is written.

**No gate is proposed for any of these here.** Which one, at what size, with
what threshold, is precisely what may not be chosen after seeing the probe.

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

- A positive result here supports **the exclusivity clause failing** — that
  curvature can reach the causal order without the measure — in a vacuum
  plane wave, over the tested patch and profile. It does not establish that
  the effect is recoverable, quantitatively, as curvature. Detection and
  recovery are different claims and P12 needed a separate stage for the
  second.
- A **null** result is a statement about the tested statistic and the tested
  patch, not about the thesis. The transverse quadrupole structure named
  in §5 candidate 3 is a way a scalar statistic could be blind by
  symmetry — `H` focuses in `x` and defocuses in `y` — and a
  null from a scalar statistic must be reported with that limitation
  attached in the same sentence.
- Plane waves are a special family: Petrov type N, `sqrt(-g) = 1`, not
  globally hyperbolic. The exactness of the measure is what makes the
  experiment clean and is also what makes it unrepresentative of general
  `Weyl != 0` spacetimes. Generalizing beyond type N is a different design.

## 8. What comes next, and it is a probe, not a stage

This repository has a precedent for exactly this situation:
`docs/p7_fss_rescope_recon.md`, where a reconnaissance closed the FSS
rescope path in about thirty seconds of compute by asking whether the
instrument worked at the sizes the sampler could reach. The decision record
notes what killed it: "**죽인 것은 컴퓨트가 아니라 계측기 물리**".

P14 has the same shape and should be answered the same way, **before any
prereg freeze and with nothing frozen**:

- **P0. Measure.** Confirm numerically that a uniform-coordinate sprinkling
  reproduces the true 4-volume in the patch (it must, given `sqrt(-g) = 1`,
  so this is an implementation check, not a physics one).
- **P1. Causality.** Implement `p precedes q` in the plane-wave patch and
  check it against direct null-geodesic integration to a pinned tolerance.
  Verify the conjugate-point constraint of §4.3 is respected by every pair.
- **P2. Control C0.** Show the P12/P13 volume reading returns flat.
- **P3. Discriminability.** For a ladder of profile amplitudes `A` and
  patch sizes, measure whether ANY of §5's order statistics separates
  `A != 0` from `A = 0` **on the same point set**, and with what effect
  size per sample. Report curves with uncertainties. **No gates, no
  verdicts** — this is characterization, in the sense P10 §6 uses the word.
- **P4. Paired variance.** Report the variance of the paired difference
  directly, so §6.2's expectation is measured rather than assumed and the
  power section knows what the interval rule actually costs here.

### 8.1 The second probe: what a general `Weyl != 0` spacetime would cost

Decided 2026-08-01: the plane-wave probe runs, **and so does a separate,
narrower measurement on Schwarzschild.** The reason is §7's honest limit —
plane waves are Petrov type N, and a null result there cannot distinguish
"the thesis holds" from "this family is too special". Leaving the
generalization path unpriced until after the first probe would repeat the
mistake the Wang-Landau record names: porting a kernel before measuring the
exponent that decides whether the port is worth anything.

This probe does **not** run a campaign and does not need one. It answers
one question:

- **S1. Cost of causality.** Implement `p precedes q` in a Schwarzschild
  patch by null-geodesic integration and measure the wall-clock per pair
  and its scaling. Compare against the sprinkling sizes P11–P13 operated
  at (`n` in the low thousands per rung, so `O(n^2)` relations per sample).
  Report a price, not a verdict.

Two things follow from S1 regardless of P3's outcome. If the price is
affordable, Schwarzschild becomes the generalization stage and the
type-N limitation is a temporary one. If it is not, that number is what
closes the general-`Weyl` path for now, in the same quantified way the
tunneling exponent closed Wang-Landau — and a plane-wave null result then
has to be reported with a limitation that is **priced**, not merely named.

The Section 3 obstacle stands and S1 does not resolve it: Schwarzschild has
no closed-form diamond volume, so a campaign there would still have to
confront the small-diamond expansion that §10.1 refused. S1 measures the
causality cost only; the volume question is separate and unaddressed.

### 8.2 Exit conditions

If P3 finds no separation at achievable sizes, P14 closes there, cheaply
and honestly, and the programme's thesis remains a hypothesis with a
documented reason — with S1's price attached, so the closure records what
the alternative would have cost rather than implying none existed. If it
does separate, P3's effect sizes and P4's paired variance are what the
power section of the actual preregistration is built from — power-first, as
P11 through P13 all were, and now sized for an interval-rule accuracy gate
per §6.

## 9. The Section 4 verification, in full

Run with `sympy 1.14.0`, exact rational/symbolic arithmetic, metric
`ds^2 = 2 du dv + A(u)(x^2 - y^2) du^2 + dx^2 + dy^2` with `A(u)` an
undetermined function:

```
metric determinant  det(g) = -1
volume element sqrt(-g)    = 1

Ricci tensor:  [0 0 0 0]
               [0 0 0 0]
               [0 0 0 0]
               [0 0 0 0]
Ricci scalar R = 0

non-vanishing Riemann components: 8
  R^v_{xux} = -A(u)      R^v_{xxu} =  A(u)
  R^v_{yuy} =  A(u)      R^v_{yyu} = -A(u)
  R^x_{uux} =  A(u)      R^x_{uxu} = -A(u)
  (and the two y-components)

Ricci == 0 ?  True
Riemann == 0 ?  False
=> vacuum with non-zero curvature (pure Weyl) ?  True
```

The check itself must be committed as a frozen artifact under
`docs/prereg/frozen/p14/` before any freeze, per the repository rule that
every frozen artifact records the `code_version` that produced it. It is
not committed by this draft, because this draft freezes nothing.

## 10. Open questions this draft does not answer

1. Which order statistic, if any, is Weyl-sensitive. §5 lists candidates and
   §8 P3 is the experiment; nothing here predicts the answer.
2. Whether the quadrupole structure of `H = A(u)(x^2 - y^2)` makes scalar
   statistics blind by symmetry, and whether a profile without that symmetry
   would be a better or a less standard choice.
3. The cost of deciding causality per pair, which sets whether the campaign
   sizes P11–P13 operated at are reachable here at all. **Assigned:** §8 P1
   for the plane wave, §8.1 S1 for Schwarzschild.
4. Whether the non-global-hyperbolicity of the family can be confined to a
   patch constraint cleanly enough to survive review, or whether it must be
   escaped by a sandwich profile (`A(u)` compactly supported), which
   restores global hyperbolicity at the price of a `u`-dependent profile.
5. Whether the same-points pairing of §4.1 tightens the interval enough to
   absorb the interval rule's cost, or whether the power section pays the
   full `1.38x` of §6.1. **Assigned:** §8 P4.
6. What the diamond volume is in a Schwarzschild patch, if S1 says the
   causality cost is affordable. §3's standing refusal of unverified
   expansion coefficients applies, and this draft has no answer.
