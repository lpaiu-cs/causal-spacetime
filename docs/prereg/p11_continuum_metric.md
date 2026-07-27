# P11: a metric instrument for the continuum limit — power-first preregistration

Status: **DESIGN v1.1** (2026-07-27; v1.0 same day). Nothing below has
been run. Dates in this record are local dates (UTC+9, Asia/Seoul);
commit timestamps carry their +09:00 offset.

v1.1 review corrections (pre-freeze, pre-data): (i) the coordinate
normalization is now frozen explicitly — the unit-diamond `/N`
convention gives `E|I| = N tau^2 / 4`, so the chain estimator is
`L / sqrt(N)`, not the v1.0 `L / sqrt(2N)` whose assumed
`E|I| = N tau^2 / 2` matched no convention in the repository; (ii)
the verdict table is evaluated by precedence, making the categories
mutually exclusive (v1.0's IMPROVES and FLAT could both match a small
all-negative CI).

v1.2 review corrections (pre-freeze, pre-data): (iii) the hypothesis
no longer says "pure order" — the count supplies the density
calibration, so a pass supports order-plus-density reconstruction
(Sections 2, 3); (iv) Stage B's pair-disjointness support is now
defined (the realized meet-join bounding box — v1.1's causal-interval
rule was undefined for spacelike pairs); (v) Stage B gets its own
quarantined variance pilot, sample-size formula, and frozen seed
windows instead of inheriting Stage A's `n` computed from a different
estimator's variance.

v1.3 review corrections (pre-freeze, pre-data): (vi) the chain
estimator subtracts the closed-interval endpoints — v1.2's
`L / sqrt(N)` carried a `+2/(sqrt(N) tau)` relative bias that FALLS
with `N` (0.20 to 0.10 across the ladder), i.e. it would have faked
part of the improvement the gate exists to test; the estimator now
goes through the repository's shared definition
(`estimate_tau_from_longest_chain_1p1`, which subtracts the endpoint
contributions), per the shared-definition rule. (vii) Stage B's
estimator and support are DE-FROZEN to an addendum: review checked
100 samples per rung and found the v1.2 bounding-box rule cannot
yield 6 disjoint supports (the meet/join frontiers are Pareto
staircases spanning the null arms), and in-house derivation found the
v1.2 median-over-`M x J` estimator biased upward (far-from-corner
candidates contribute `tau > d`); constraints for the addendum are
recorded in Section 3. (viii) the feasibility projection weights
rungs by cost (`1 + 2 + 4`), not `3x` the bottom rung. (ix) the
ensemble language is corrected: the uniform permutation is a
sprinkling CONDITIONED on `N`, and `(i, pi(i))/N` is its rank grid —
the permutation-based chain theorems apply verbatim, but count
fluctuations are sub-Poisson and are labelled as rank-ensemble
statements.

Lineage. P10 closed with: "What a continuum-limit test now requires is
a reconstruction whose continuum accuracy grows with density — a
different instrument." This document designs that instrument. P10's
record is closed and is not reopened by anything here; P10 also taught
this design its two hardest lessons, and they are structural below:
(i) 20 samples per rung could not separate effects of the size the
question produces — so the POWER SECTION COMES FIRST and every gate is
powered before it exists; (ii) consecutive seed allocation let derived
streams collide across nominal replicates — so every sample owns a
stride-200 private seed window from day one.

The scientific bet, stated before any data: the P10 discriminator was
a fixed-relative-resolution device (fixed comparison budget, fixed
constraint count, rank-valued inputs), so its continuum-unit accuracy
COULD NOT grow with `N`. The estimators below are counting-measure
estimators whose convergence as density grows is not a hope but a
theorem — the longest-chain law. If the continuum limit is real for
emergent geometry, THIS instrument family is one through which it can
appear; if accuracy fails to improve even here, that is a finding
about the limit, not about an instrument's resolution.

---

## 1. Power design (first, by directive)

### 1.1 Primary contrast and predicted effect size

Ladder: `N in {600, 1200, 2400}` (geometric, factor 4 end to end; the
middle rung exists to expose non-monotonicity — the unexplained 900
dip of P10 — and takes no part in the gate).

Per-sample statistic: `y = log10( median over K pairs of the relative
proper-time error |tau_hat - tau_true| / tau_true )`, with the pair
protocol of Section 4. Rung statistic: median of `y` over the `n`
samples. Primary contrast:

    Delta = median_y(N = 2400) - median_y(N = 600)

Theory (Section 7) predicts relative error `~ m^(-1/3)` with
`m = E|I| = N tau^2 / 4` the interval cardinality (Section 3's frozen
normalization), so at fixed continuum separation:

    Delta* = -(1/3) * log10(4) = -0.2007 dex   (primary, chain estimator)
    Delta* = -(1/2) * log10(4) = -0.3010 dex   (secondary, volume estimator)

### 1.2 Sample-size formula (frozen before the pilot)

With `sigma` the per-sample SD of `y` at the bottom rung and medians
carrying the normal-approximation efficiency factor `pi/2`, the
two-rung contrast needs, for two-sided alpha = 0.05 and power 0.90
(`z = 1.960 + 1.282 = 3.242`):

    n_per_rung = ceil( pi * sigma^2 * 3.242^2 / Delta*^2 )
               = ceil( 819.8 * sigma_hat^2 )

bounded by **floor 12** and **cap 60**. The pilot (1.3) supplies
`sigma_hat`. If the formula demands more than the cap, Stage A is
declared INFEASIBLE-AS-DESIGNED and stops before running — that
refusal is this preregistration operating, exactly as P10's gates
were.

Worked anchors (illustrative only, the pilot decides): sigma_hat =
0.10 -> n = 12 (floor); 0.15 -> 19; 0.20 -> 33; 0.27 -> cap.

### 1.3 Stage P: the pilot, variance-only, quarantined

- Bottom rung only (`N = 600`), `n_pilot = 12` samples, seed window
  Section 5.
- Records: per-sample `y`, its SD `sigma_hat`, and wall time.
- **Quarantine**: the pilot computes NO cross-rung contrast (it runs
  one rung, so none exists), and its samples are never reused in any
  stage. Its two outputs are `n_per_rung` (via 1.2) and feasibility.
- Feasibility rule (frozen; v1.3 corrected the rung weighting):
  per-pair work scales linearly in `N` (interval extraction `O(N)`,
  interior LIS near-linear in `m ∝ N`), so the projection weights
  each rung by `N_r / 600`: projected Stage A wall time
  `= n_per_rung * (pilot wall time / 12) * (1 + 2 + 4)` — seven
  bottom-rung units per sample-set, not the v1.2 `3x` that priced
  every rung at bottom-rung cost. If it exceeds 12 hours,
  INFEASIBLE-AS-DESIGNED is recorded and the design returns to the
  drawing board without any outcome having been read.

**Stage P-B (power protocol frozen now; estimator by addendum).**
`d_hat` carries locator variability the timelike estimator does not,
so Stage B is powered from ITS OWN variance, never Stage A's (v1.2,
review): a second variance-only, quarantined pilot — bottom rung,
`n_pilot = 12` spacelike samples under the Stage B pair protocol —
feeds the SAME formula (1.2) with the SAME `Delta* = -0.2007`,
giving Stage B its own `n_per_rung`, floor 12 / cap 60, and the
rung-weighted feasibility stop above. Quarantine as in Stage P: one
rung, no contrast computable, samples never reused. The `d_hat`
DEFINITION this pilot runs is the Stage B addendum's (Section 3,
v1.3) — the power protocol here does not change when that addendum
lands, which is constraint (3) on it.

### 1.4 Verdict logic (frozen for every stage gate)

Each stage's gate reads the 95% bootstrap CI `[lo, hi]` of its
`Delta` (4000 resamples, `_stable_seed` labels fresh to this
experiment). Rows are evaluated IN ORDER and the first match is the
verdict — the categories are mutually exclusive by precedence:

| # | condition | verdict |
|---|---|---|
| 1 | `hi < 0` | **IMPROVES** — gate passes (a small established improvement still passes; its size against theory is the labelled slope/constant check, not the gate) |
| 2 | `lo > 0` | **DEGRADES** — gate fails, informative |
| 3 | `-0.067 < lo` and `hi < +0.067` (CI straddles 0 inside the margin) | **FLAT-WITHIN-MARGIN** — gate fails, informative |
| 4 | otherwise | **INCONCLUSIVE** — gate fails; indicts the pilot variance estimate; record and stop |

The equivalence margin `delta_eq = |Delta*|/3 = 0.067 dex` is frozen
NOW, before any data, so that a flat result is a verdict rather than a
shrug — the vocabulary P10 lacked until its post-hoc rounds. A
DEGRADES verdict from THIS instrument (a consistent estimator family)
would be a major anomaly and is left representable.

Completeness is part of every gate: a sample that cannot supply all
`K` pairs under the frozen pair rule is recorded incomplete, and gates
require full completeness at every rung (no
completeness-conditioning; the P10 Stage B lesson).

---

## 2. The question

Does emergent geometry deepen toward a continuum? Operationally here:
does there exist an instrument — inputs: relabeling-invariant order
quantities (chain lengths, interval cardinalities) PLUS the event
count as the density calibration — whose continuum-unit metric
accuracy IMPROVES as sprinkling density grows, at the rate theory
predicts?

Scope of a pass, stated before any data (v1.2, review): causal order
alone fixes only the conformal structure; the count supplies the
volume scale — the repository's own convention
(`src/causal_spacetime_lab/metrics.py`: order supplies the interval
count, density supplies the metric scale), and the causal-set thesis
"order + number = geometry". A passing experiment therefore supports
**reconstruction from order plus density calibration**, never a
pure-order metric claim: the same order with a rescaled diamond would
carry different proper times, and it is `N` in
`tau_hat = (L-2)/sqrt(N)` that fixes the scale. Dimensional scope: 1+1D first (the dictionary
regime where the longest-chain law is a sharp theorem); higher
dimensions are a listed lever, not part of this preregistration.

## 3. The instrument (inputs: order invariants + the count as scale)

Setting: uniform random order of size `N`, via the frozen
null-coordinate dictionary (`order_inputs` — imported, never
re-implemented; the shared-definition rule). Ensemble, stated
precisely (v1.3): this equals a sprinkling of the unit diamond
CONDITIONED on the count `N`, and the `(u, v) = (i, pi(i))/N`
coordinates form its rank grid — exactly one event per null-coordinate
rank, so interval counts have the same means as Poisson
(`E|I| = N du dv`) but sub-Poisson, negatively associated
fluctuations. The chain law needs no translation — the LIS theorems
of Section 7 are permutation theorems, i.e. theorems about THIS
ensemble; Poisson-flavoured variance intuitions are used nowhere as
gates and are labelled where they appear (the volume estimator's
noise reads conservative under sub-Poisson counts).

**Frozen coordinate convention (v1.1).** `order_inputs` returns
order-unit values `t = i + pi(i)`, `x = i - pi(i)`. The continuum
frame is the repository's `/N` convention: null coordinates
`(u, v) = (i, pi(i)) / N`, uniform on the unit square at density `N`.
The metric is `ds^2 = dt^2 - dx^2 = 4 du dv`, so for a related pair
`tau_true = 2 sqrt(du dv)` and for a spacelike pair
`d_true = 2 sqrt(|du dv|)`; a causal interval of proper time `tau`
occupies `(u, v)`-area `tau^2 / 4` and therefore has expected
cardinality **`E|I| = N tau^2 / 4`**. Every formula below uses this
and only this convention (v1.0 assumed `N tau^2 / 2`, which matches
no convention in the repository — review finding, corrected before
any implementation).

- **tau_hat (timelike, PRIMARY — the chain estimator).** For related
  events `x < y`: `L(x, y)` = longest chain in the closed interval
  `I(x, y)` (endpoints included, the convention of the shared
  `longest_chain_length`); then, through the repository's shared
  definition `estimate_tau_from_longest_chain_1p1(L, rho = N/2)`:

      tau_hat_chain = (L - 2) / sqrt(N)

  The endpoint subtraction is load-bearing (v1.3): the closed-interval
  `L` exceeds the interior LIS by exactly the two endpoints, and
  without the correction the estimator carries a `+2/(sqrt(N) tau)`
  relative bias — 0.20 at the bottom rung falling to 0.10 at the top,
  a systematic that would have faked improvement in the primary gate's
  own direction. Interior expected cardinality `m = N tau^2 / 4`; the
  interior longest chain converges to `2 sqrt(m) = sqrt(N) tau`
  (Section 7). Computed exactly: 2D dominance order makes it a
  longest increasing subsequence.
- **tau_hat_vol (timelike, secondary).** `2 sqrt(|I_open(x, y)| / N)`
  with `|I_open|` the interior count, via the shared
  `metrics.py` cardinality estimator mapped to this convention;
  faster predicted convergence (`m^(-1/2)`), noise labelled
  rank-ensemble (sub-Poisson).
- **d_hat (spacelike, Stage B) — NOT YET FROZEN (v1.3).** Two defects
  killed the v1.2 candidate before implementation: (a) review checked
  100 samples per rung under the frozen ensemble and the
  bounding-box-of-`{x, y} ∪ M ∪ J` supports cannot yield 6 disjoint
  pairs (the meet/join frontiers are Pareto staircases spanning the
  null arms, so each support covers most of the diamond); (b)
  in-house derivation: the median over `M x J` of interval proper
  times is biased upward — candidates far from the continuum corners
  contribute `tau > d`, so the estimator was not consistent as
  written. Spacelike distance from order data is a recognized hard
  problem; a candidate is frozen by ADDENDUM after Stage A reports
  and before any Stage B data, and must satisfy, in writing, all of:
  (1) a consistency argument for `d_hat -> d_true` as `N` grows;
  (2) a disjointness support that six pairs can realize at the
  bottom rung, demonstrated on synthetic checks that are quarantined
  from Stage B's seeds; (3) compatibility with the Stage P-B power
  protocol of 1.3 unchanged.
- **Coordinates (Stage C).** Declared here, frozen later (staged
  freezing, the B' precedent): reconstruction of per-event `(t, x)`
  from tau_hat / d_hat to anchor sets, gated on the continuum
  coordinate error falling. Its estimator addendum is frozen only
  after Stage B reports, and before any Stage C data.

Truth values `tau_true`, `d_true` come from the dictionary's continuum
coordinates — this is a positive control throughout, as the programme
has always run.

## 4. Pair protocol (scale-referenced, a priori)

- Separation band: `tau_true in [0.35, 0.45]` (timelike pairs,
  Stage A); `d_true in [0.35, 0.45]` (spacelike, Stage B). Continuum
  units, chosen a priori: with `m = N tau^2 / 4`, `m ~ 24` at
  `N = 600` and `~ 96` at `N = 2400` — clear of the discreteness
  floor and the diamond boundary, but deeper in the pre-asymptotic
  regime than v1.0's miscounted anchors suggested, which is one more
  reason the `-1/3` exponent is a labelled consistency check and
  never the gate.
- `K = 6` pairs per sample, drawn by the sample's own scoring stream
  (Section 5), accepted only if their SUPPORTS are pairwise disjoint
  (kills within-sample dependence); up to 200 rejection draws, else
  the sample is INCOMPLETE. The support: **Stage A**: the closed
  causal interval `I(x, y)`, which in `(u, v)` IS the pair's bounding
  box. **Stage B**: frozen with the Stage B addendum (v1.3 — the v1.2
  `{x, y} ∪ M ∪ J` box rule was checked infeasible, and its estimator
  inconsistent; Section 3), under constraint (2) there.
- Everything above is frozen before any run; the band never adapts to
  data (the B' scale-referencing lesson made structural).

## 5. Seed allocation (private windows from day one)

Stride 200; every sample owns `[s, s + 199]`; every derived stream an
implementation may ever take (pair selection at `s + 150`, any future
offset < 200) stays inside its own row's window. Windows:

| block | base | samples | window span |
|---|---|---|---|
| Stage P pilot (N=600) | 60000 | 12 | 60000-62399 |
| Stage A, N=600 | 64000 | <= 60 | 64000-75999 |
| Stage A, N=1200 | 76000 | <= 60 | 76000-87999 |
| Stage A, N=2400 | 88000 | <= 60 | 88000-99999 |
| Stage P-B pilot (N=600) | 100000 | 12 | 100000-102399 |
| Stage B, N=600 | 104000 | <= 60 | 104000-115999 |
| Stage B, N=1200 | 116000 | <= 60 | 116000-127999 |
| Stage B, N=2400 | 128000 | <= 60 | 128000-139999 |
| Stage C blocks | 140000+ | — | frozen with Stage C's addendum |

All spans sit above every range the programme has ever used (documented
maxima: 30000-30379, 40000-40168, 41000-41059, 43000-54999); a
regression test pins pairwise-disjoint windows and freshness against
the full documented list before the pilot may run.

## 6. Stages and gates

- **Stage P** (pilot): variance + feasibility only (1.3). No gate on
  outcomes; produces `n_per_rung`.
- **Stage A** (timelike): the primary gate of the experiment —
  `Delta` for tau_hat_chain per 1.4. Secondary, labelled consistency
  checks, never gates: (a) OLS slope of median `y` on `log10 N`
  across the three rungs, with bootstrap CI, against the predicted
  `-1/3` (chain) and `-1/2` (volume); (b) the constant-level band
  `median relerr ~ 0.89 m^(-1/3)`; (c) middle-rung monotonicity (the
  900-dip watch). Stage B runs only if Stage A passes.
- **Stage B addendum** (after Stage A reports, before any Stage B
  data): freezes `d_hat`, its support, and its Stage B seed use,
  meeting the three constraints of Section 3 in writing.
- **Stage P-B** (spacelike pilot): variance + feasibility only, per
  1.3 — supplies Stage B's own `n_per_rung`. Runs only after Stage A
  passes and the addendum is frozen.
- **Stage B** (spacelike): the 1.4 verdict logic applied to d_hat's
  `Delta`, at Stage P-B's `n`. Superiority gate only — locator noise
  makes an exponent gate overconfident; the slope is recorded,
  labelled.
- **Stage C** (coordinates): gate declared (continuum coordinate
  error falls, same verdict table), estimator frozen by addendum
  after Stage B, before any Stage C data.

Failure at any gate stops the programme at that stage with the verdict
recorded — no post-hoc rescue, no gate re-reading. Corrections, if
review finds defects, append as dated notes; artifacts regenerate only
from commits containing their implementation, with the dirty-aware
stamp (untracked files count) inherited from P10.

## 7. Theory appendix: why convergence is a theorem here

For a uniform random 2D order, the longest chain across an interval of
`m` points is the longest increasing subsequence of a uniform random
permutation: `E L_m = 2 sqrt(m)` (Vershik-Kerov; Logan-Shepp), with
fluctuations `m^(1/6)` and Tracy-Widom limit law, mean shift
`mu_1 ~ -1.77` (Baik-Deift-Johansson). With the Section 3 convention
(`m = N tau^2 / 4`, `2 sqrt(m) = sqrt(N) tau`),

    relative error of tau_hat_chain ~ ( |mu_1| + O_P(1) ) * m^(-1/3) / 2,

both bias and noise decaying as `m^(-1/3)`, `m` growing linearly in
`N` at fixed continuum separation. The proper-time reading of the
longest chain is the causal-set standard (Myrheim; Brightwell-Gregory)
on the sprinkling foundation of Bombelli-Lee-Meyer-Sorkin. The volume
estimator converges by the CLT at `m^(-1/2)`.

Worth recording: the G2 retraction (T1 positioning round) established
that this programme's convergence observations were known LPP/LIS
theorems rather than novelties. Those same theorems are now the
load-bearing guarantee of this design — the retraction did not shrink
the programme, it located the ground it stands on.

## 8. Implementation plan (after this document merges, before Stage P)

New module `experiments/positive_control/p11_metric.py`: interval
extraction by 2D dominance (`O(N)` per pair), the estimators, the
pair protocol, stage runners P/A. The shared-definition rule applies
from the start, three times over (v1.3): chain lengths through
`causal_spacetime_lab.chains.longest_chain_length` (whose
endpoints-included convention is exactly why the correction exists),
`tau_hat_chain` through
`causal_spacetime_lab.estimators.estimate_tau_from_longest_chain_1p1`
with `rho = N/2` (which performs the endpoint subtraction and the
`sqrt(2 rho)` normalization — equal to `(L-2)/sqrt(N)` under the
frozen convention, verified by a pin test), and `tau_hat_vol` through
the `metrics.py` cardinality estimator ("rho supplies metric scale;
the causal order alone provides the interval count") mapped to the
same convention. Nothing re-implemented. Tests
before any run: LIS against brute force on small posets; tau_hat
normalization pinned on constructed intervals; pair-disjointness and
incompleteness paths; seed-window privacy and freshness against the
full documented list; verdict-table logic on synthetic inputs. Stage P
runs only from a clean commit containing all of it.
