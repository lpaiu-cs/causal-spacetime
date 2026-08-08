# P14 §8 P2 — the volume prediction, measured (exploratory)

**EXPLORATORY.** No gate, no threshold, no verdict, no confirmation
seed window. The `equivalent` / `discriminates` labels below are
properties of confidence intervals at the stated coverage, never
program verdicts (design §8.2). Campaign seeds `20260811/12/13`, MC
sub-streams from `20260814`; burned seeds `20260808` (P1), `777`
(design checks) verified absent by the script and a test.

**The question.** Do the curved arm's anchor-diamond cardinalities
reproduce the §4.4 quadrature `r = V_A/V_0`? This replaces the dead
C0 with a number. Scope is **axis anchors only** — the one diamond
with an independent analytic prediction; off-axis curved volumes are
known only through MC, and a check against them would be MC-vs-MC
circular.

## Protocol (frozen in the design review, 2026-08-08, three rounds)

- **Error model.** Anchors are EXTERNAL — the fattest eligible axis
  diamond of the frozen geometry, never among the sprinkled points
  (the count's strict `u` window excludes them by construction,
  test-pinned) — and `sprinkle` draws `N ~ Poisson(ρV_box)`, so
  `N_A`, `N_0` are exactly Poisson over deterministic regions and
  `Cov(N_A, N_0) = ρ·V_∩` holds as written. Post-selecting anchors
  from a realization would void the whole model (§8 P2, R6.5).
- **Primary — equivalence, per the accuracy house rule** (interval,
  95%, two-sided): normalized error `θ = E[Z]/λ_0` with
  `Z_i = N_A,i − r·N_0,i`; the **Student-t** 95% two-sided CI on `θ`
  must lie inside `[−τ, +τ]`, with `τ = δ = r − 1` the chosen
  effect-scale margin. The licensed sentence is *"the error from the
  quadrature prediction is smaller than the flat-to-prediction gap"*
  — no coefficient-grade (`τ ≪ δ`) recovery is claimed. The t-CI was
  verified by compound-Poisson simulation at the frozen `n`:
  coverage `0.9495/0.9492/0.9473`, equivalence power `1.0000`.
- **Secondary — discrimination**: the same CI excludes the
  flat-truth value `θ = −δ`.
- **Marginal GOF per arm** — rate ratio (observed pooled count /
  predicted pooled mean) with its exact **Garwood** Poisson 95% CI
  inside `[1−τ_m, 1+τ_m]`, `τ_m = δ`. This catches common-mode
  error, which `Z` is exactly blind to (P1 review R7.1), and it is
  the constraint that BINDS `n` at the first two points.
- **Var(Z) diagnostic** — empirical between-sprinkling variance
  (bootstrap percentile 95% CI, `B = 4000`) against the analytic
  `r·ρ·V_dis` (CP bracket, 2e6-sample MC on a separate sub-seed).
  Reported, never judged: a cross-check of the error model, not an
  estimand.
- **Coverage** — pointwise 95% per operating point. Any sentence
  joining the three points is Bonferroni-adjusted and separate.
- **Sizing** — `n = 21 000 / 400 / 100` =
  `max(primary TOST at the V_dis 95% UCB, marginal at τ_m = δ,
  variance floor 100)`; exact minima `20 939 / 349 / 13`. Inputs:
  P1's committed `p14_probe_p1_sizing.json` (cross-checked at run
  time) and the review's 2e6-sample MC (120/79/737 hits, design seed
  777, burned).
- **Statistics implementation** — the repo carries no scipy; the
  Student-t quantile (incomplete beta) and the Garwood interval
  (incomplete gamma) are built in the script and test-pinned against
  published table values, the same way P1 hand-built
  Clopper–Pearson.

Everything below this line is RENDERED from the committed artifact
`p14_probe_p2_results.json` by the script's table functions, and a
test asserts the doc embeds the renderings verbatim — nothing here is
hand-transcribed (the P1 hits-column erratum, closed structurally,
PR #41).

## 1. Primary: θ equivalence, τ = δ

| point | n | θ 95% CI (t) | τ = δ | equivalent | discriminates |
|---|---|---|---|---|---|
| slice-a1.0 | 21000 | [-0.00084, +0.00005] | 4.000e-03 | yes | yes |
| high-a2.0 | 400 | [-0.01394, +0.01903] | 7.301e-02 | yes | yes |
| edge-a2.4 | 100 | [-0.02409, +0.02737] | 1.815e-01 | yes | yes |

Joint sentence over the three points (Bonferroni, per-point t at 99.17%): equivalent yes, discriminates yes. The primary record stays pointwise.

## 2. Marginal rate-ratio equivalence, τ_m = δ

| point | arm | rate ratio 95% CI (Garwood) | 1 ± τ_m | equivalent |
|---|---|---|---|---|
| slice-a1.0 | curved | [0.99653, 1.00086] | 1 ± 4.000e-03 | yes |
| slice-a1.0 | flat | [0.99691, 1.00126] | 1 ± 4.000e-03 | yes |
| high-a2.0 | curved | [0.99945, 1.07227] | 1 ± 7.301e-02 | yes |
| high-a2.0 | flat | [0.99584, 1.07120] | 1 ± 7.301e-02 | yes |
| edge-a2.4 | curved | [0.94322, 1.00667] | 1 ± 1.815e-01 | yes |
| edge-a2.4 | flat | [0.93916, 1.00810] | 1 ± 1.815e-01 | yes |

## 3. Var(Z) diagnostic (reported, not judged)

| point | Var(Z) emp [boot 95%] | analytic r·ρ·V_dis [CP 95%] | ratio |
|---|---|---|---|
| slice-a1.0 | 1.645 [1.608, 1.682] | [1.280, 1.898] | 1.04 |
| high-a2.0 | 1.399 [1.179, 1.646] | [1.645, 2.362] | 0.70 |
| edge-a2.4 | 16.955 [12.077, 22.366] | [11.913, 13.800] | 1.32 |

## Reading

**The volume prediction holds at every operating point, at the
effect scale.** All three primary CIs sit inside `[−τ, +τ]` and
exclude the flat value, pointwise 95% — and the joint Bonferroni
sentence holds too. The two §4.4 pinned values (`1.00400047` at
`wT = 1`, `1.07300802` at `wT = 2`) are now *measured* facts about
sprinkled counts, not only quadrature outputs, and C0's silence is
replaced by six passing interval checks.

Three observations the tables carry:

- **The `a = 1` pin is the strong one.** The marginal-driven
  `n = 21 000` makes slice-a1.0's primary CI (`±0.00045` on `θ`)
  about 9× tighter than its `τ = 0.004` requires — the measured
  ratio agrees with the quadrature well inside the equivalence
  margin, at the operating point where the effect itself is only
  0.4%.
- **All six marginals pass at `τ_m = δ`**: the absolute calibration
  (`ρ ×` volume, both arms, all points) is verified at the same
  scale as the effect. At high-a2.0 both arms sit ~3.5% above their
  predicted means *together* — the same-realization correlation
  (`Cov = ρV_∩`) moves the two marginals in lockstep while the
  contrast `Z` cancels it, which is precisely why the marginal check
  exists as a separate contract.
- **The Var(Z) diagnostic is consistent with the error model where
  it is precise and noisy where it is not.** slice: ratio 1.04 with
  the tightest brackets. high and edge land on opposite sides
  (0.70, 1.32) of their analytic brackets — single-comparison
  fluctuations of the size six 95% diagnostics are expected to
  produce, in opposite directions, with no common pattern. Reported,
  as the protocol requires; judged by nothing.

**What this feeds.** P2's contract is settled for the external-anchor
path: the counting instrument reproduces the analytic volume ratio
under the exact error model. P3 takes the split (external anchors
admissible at large `a`; sprinkled pairs starved there) and asks
whether the §5 statistics separate; P4 measures the paired variance
the power section will be built from.

## Scope and limits

- Axis anchors only; the claim covers the on-axis diamond's volume
  ratio, nothing else.
- `τ = δ` licenses "error smaller than the flat-to-prediction gap"
  per point, pointwise 95%. It does not establish the `a⁴/252`
  coefficient — that is a later campaign's question, needing
  `τ ≪ δ`.
- The marginal check verifies the absolute calibration (`ρ × volume`)
  at the effect scale; a common-mode miscalibration smaller than
  `τ_m` would pass both it and the primary.
- S1's eventual price closes (or opens) the Schwarzschild solver /
  domain / budget path it measures — not the general-Weyl question.

## Changelog

Initial record.
