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

v1.9 review corrections (pre-data): (xxiii) the Bonett bound is
labelled what it is — an asymptotic approximation, not a
finite-sample distribution-free guarantee (no such variance bound
exists) — and is VALIDATED on the realized distribution instead of
trusted by construction: the pilot bootstrap-checks the bound's
coverage on its own samples (2000 resamples per rung) and, when the
nominal z falls short of 95% measured coverage, raises z until it is
met; the calibration record (coverage at nominal, z used) is
published in the pilot artifact, and the residual limitation — tail
mass absent from 200 pilot observations — is stated rather than
claimed away. (xxiv) provenance stamps reading "unknown" (git
unavailable: archives, copied trees) are refused by the preflight
exactly like dirty stamps — two unknowns would otherwise even pass
the prerequisite equality check.

v1.8 review corrections (pre-data): (xxii) the v1.7 chi-square
variance bound was itself a Gaussian-pivot result — the same class of
assumption v1.4 removed from the median factor, reintroduced one
level up (review). Replaced by a distribution-robust construction:
the pilot grows to 200 samples per endpoint rung (the 12-sample size
was priced for an MCMC-era cost model that the frozen verification
record falsified — 200 samples cost ~20 seconds), and each rung's
variance takes its kurtosis-adaptive Bonett upper bound at one-sided
95%, summed with Bonferroni so the pair-level guarantee is >= 90%
with no normal model anywhere. Seed windows renumbered for the
larger pilot blocks (nothing has run).

v1.7 review corrections (pre-data; the implementation exists but no
stage has run): (xix) the power formulas size from a conservative
variance bound, not the 12-sample point estimates — each pilot
variance is inflated to its one-sided 90% chi-square upper bound
(factor `11 / chi2_{0.10, 11} = 1.972` at nu = 11), so a downward
sampling error in a small pilot can no longer smuggle an infeasible
design past the cap or declare an unaffordable flat verdict
available. (xx) the first-`n`-complete estimand is named for what it
is — conditional on successful pair packing — with the frozen guard
that makes it innocuous: the verification record bounds per-sample
completion below one expected skip per campaign, zero realized skips
(the generic case) makes the conditional and unconditional estimands
coincide on the data, and any realized skip stamps the gate summary
with a SELECTION-CAVEAT label. (xxi) the introduction's "finding
about the limit" sentence is rewritten to defer to 1.4's
interpretation scope — one conservative voice, no contradiction.

v1.6 review corrections (pre-freeze, pre-data): (xvii) the contrast
variance uses BOTH endpoint rungs — Stage P now pilots `N = 600` and
`N = 2400` (variance and wall time only, cross-rung statistics
forbidden), and the formulas take `sigma_b^2 + sigma_t^2` in place
of `2 sigma_600^2`, which had assumed an equal-variance model the
pre-asymptotic ladder is owed nowhere. (xviii) completeness is
handled at campaign level: a per-sample 0.95 pin admits designs that
die almost surely (`0.95^180 ~ 1e-4`) under the all-complete gate;
the frozen rule is now seeds-in-order-first-`n`-complete with
per-window reserve slots and a skip cap, the verification pin is
raised to `>= 1998 of 2000` per rung (bounding expected skips below
one per campaign), and the no-completeness-conditioning concern is
answered structurally: completeness is decided by pair-packing
geometry BEFORE any estimator evaluation, `y` is never computed on
an incomplete sample, and skip counts are reported per rung. Seed
windows renumbered for the reserve slots (design phase; nothing has
run).

v1.5 review corrections (pre-freeze, pre-data): (xiv) the design is
sized for the flat verdict too — v1.4's superiority-only `n` made
FLAT-WITHIN-MARGIN unreachable (at `sigma = 0.20`, `n = 21` gives CI
half-width 0.121 against the 0.067 margin, so a truly flat world
would read INCONCLUSIVE and be mis-blamed on the pilot); `n` now
takes the max of the superiority and equivalence requirements, and
when the equivalence requirement exceeds the cap the pilot DECLARES
the flat verdict unavailable ex ante, which re-scopes INCONCLUSIVE.
(xv) the feasibility projection uses measured per-rung wall times
from the verification block instead of a scaling model — the v1.4
eligible-pool scan is `O(N^2)` per sample, so the linear `1+2+4`
weighting was wrong and any frozen exponent model would just invite
the next correction. (xvi) the endpoint-conditioned interval mean is
stated exactly — `(a-1)(b-1)/(N-2)` for rank gaps `a, b`, about 1.3%
below the asymptotic `N tau^2 / 4` at the bottom rung — and
`N tau^2 / 4` is labelled asymptotic wherever it appears.

v1.4 review corrections (pre-freeze, pre-data): (x) the rung
statistic is the MEAN of per-sample `y`, powered by the exact CLT
formula — v1.3's median-with-`pi/2` assumed a Gaussian parent the
statistic does not have; (xi) the pair-candidate pool is the
precomputed band-eligible set, with the 200-draw cap counting only
support-overlap rejections — review's 1000-sample check showed
uniform endpoint proposals complete only ~34% of samples; a synthetic
completeness pin (>= 95/100 per rung, verification-only seed block)
must pass before Stage P; (xii) the Stage B addendum must derive a
rate and a Stage-B-specific target effect `Delta*_B`, which Stage
P-B's formula uses instead of the chain contrast; (xiii) a frozen
interpretation-scope clause: only IMPROVES supports the existence
claim; FLAT / DEGRADES / INCONCLUSIVE are statements about this
estimator over this finite-density ladder (m ~ 24-96,
pre-asymptotic), never about the continuum limit itself.

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
about this estimator family over this finite-density ladder — a
sharper indictment than P10's, because here asymptotic convergence is
a theorem — but, per the interpretation scope frozen in 1.4, never
evidence that the continuum limit itself is absent. *(v1.7: this
sentence originally said "a finding about the limit"; review flagged
the contradiction with 1.4, and 1.4 is the voice of record.)*

---

## 1. Power design (first, by directive)

### 1.1 Primary contrast and predicted effect size

Ladder: `N in {600, 1200, 2400}` (geometric, factor 4 end to end; the
middle rung exists to expose non-monotonicity — the unexplained 900
dip of P10 — and takes no part in the gate).

Per-sample statistic: `y = log10( median over K pairs of the relative
proper-time error |tau_hat - tau_true| / tau_true )`, with the pair
protocol of Section 4 — the within-sample median keeps single-pair
outliers from owning a sample. Rung statistic (v1.4): the MEAN of `y`
over the `n` samples — the across-sample aggregator is the mean
precisely so the power formula rests on the CLT and not on a
distributional model of `y` (review: the `pi/2` median factor is a
Gaussian-parent result, and `y` — the log of a six-pair median of
integer-LIS errors — has no entitlement to it). Primary contrast:

    Delta = mean_y(N = 2400) - mean_y(N = 600)

Theory (Section 7) predicts relative error `~ m^(-1/3)` with
`m = E|I| = N tau^2 / 4` the interval cardinality (Section 3's frozen
normalization), so at fixed continuum separation:

    Delta* = -(1/3) * log10(4) = -0.2007 dex   (primary, chain estimator)
    Delta* = -(1/2) * log10(4) = -0.3010 dex   (secondary, volume estimator)

### 1.2 Sample-size formula (frozen before the pilot)

With `sigma_b`, `sigma_t` the per-sample SDs of `y` at the BOTTOM and
TOP rungs (both piloted — v1.6; the v1.5 `2 sigma_600^2` assumed an
equal-variance model this pre-asymptotic ladder is owed nowhere), by
the plain two-sample CLT with no distributional factor (v1.4), two
requirements are computed (v1.5 — the design must be able to AFFORD
every verdict it promises) from a CONSERVATIVE variance bound. The
bound is distribution-robust (v1.8; the v1.7 chi-square factor was a
Gaussian pivot, invalid for exactly the `y` this document says is
non-Gaussian): per endpoint rung, with `n_pilot = 200` samples,
sample variance `s^2`, and the kurtosis estimate
`g4 = n * sum((y - ybar)^4) / (sum((y - ybar)^2))^2`, the Bonett
upper bound

    se     = sqrt( ( g4 - (n-3)/n ) / (n-1) )
    bound  = s^2 * exp( 1.6449 * se )     # one-sided 95% per rung

    S^2_90 = bound_b + bound_t
             # Bonferroni at the NOMINAL levels; the pair-level 90%
             # is calibrated-approximate, not exact (v1.9)

    n_sup = ceil( S^2_90 * (1.960 + 1.282)^2 / Delta*^2 )
          = ceil( 260.9 * S^2_90 )            # 90% power at Delta*

    n_eq  = ceil( S^2_90 * (1.960 + 1.645)^2 / delta_eq^2 )
          = ceil( 2895.2 * S^2_90 )           # 90% P(FLAT) at Delta = 0

For Gaussian-kurtosis `y` the per-rung inflation is ~1.18 — cheaper
than the invalid 1.972 and honest; heavy tails raise `g4` and the
bound inflates adaptively, which is the point.

Validation instead of trust (v1.9): the Bonett interval is an
asymptotic approximation — no finite-sample distribution-free
variance bound exists — so the pilot validates it on its own
samples: 2000 bootstrap resamples per rung measure how often the
bound at the nominal `z = 1.6449` covers the full-pilot variance;
if measured coverage is below 0.95 the pilot raises `z`
(monotonically, bisection to the smallest sufficient value) and uses
THAT bound, publishing the coverage and `z` in its artifact. What no
construction can exclude — rare tail mass entirely absent from 200
observations — is recorded here as the bound's stated limitation,
not claimed away.

Frozen selection rule, floor 12 / cap 60:

- if `n_eq <= 60`: `n_per_rung = clamp(max(n_sup, n_eq))` and the
  FLAT verdict is AVAILABLE;
- else: `n_per_rung = clamp(n_sup)` and the pilot DECLARES, ex ante,
  that the flat verdict is unavailable at this variance — the
  vocabulary for that campaign is IMPROVES / DEGRADES / UNRESOLVED,
  and an UNRESOLVED there does not indict the pilot (1.4).

If `n_sup` alone exceeds the cap, Stage A is INFEASIBLE-AS-DESIGNED
and stops before running — that refusal is this preregistration
operating, exactly as P10's gates were.

Worked anchors (illustrative at GAUSSIAN kurtosis, factor ~1.18; the
realized `g4` decides, and heavier tails shift every threshold down):
equal sigma = 0.09 -> S^2_90 ~ 0.0191, n_eq 56 -> n = 56, flat
available; 0.095 -> n_eq > cap -> n = 12, flat declared unavailable;
0.15 -> n = 14; 0.20 -> n = 25; 0.31 -> n = 60 at the cap;
0.315 -> INFEASIBLE. The bound's price now scales with the measured
non-normality instead of a flat Gaussian-only factor — a 200-sample
pilot costs ~20 seconds by the frozen verification record, so the
added samples are free where the invalid assumption was not.

### 1.3 Stage P: the pilot, variance-only, quarantined

- BOTH endpoint rungs (v1.6): `n_pilot = 200` samples at `N = 600`
  and 200 at `N = 2400` (v1.8 — the 12-sample size was an MCMC-era
  cost assumption; the frozen verification record prices 200 samples
  at ~20 seconds, and the robust variance bound of 1.2 needs the
  kurtosis resolution), separate seed windows (Section 5).
- Records: per-sample `y`, the per-rung variances and kurtosis
  estimates feeding the 1.2 bound, and wall times. **Computing any cross-rung mean or
  contrast is forbidden by frozen rule** — the pilot script does not
  implement it, and the pilot artifact contains per-rung SDs and
  times only.
- **Quarantine**: pilot samples are never reused in any stage. The
  outputs are `n_per_rung` (via 1.2) and feasibility.
- Feasibility rule (frozen; v1.5 replaced the scaling model with
  measurement): the verification block (Section 4) already runs 100
  samples at EVERY rung for the completeness pin, and it records
  mean wall time per sample per rung, `t_600, t_1200, t_2400` —
  outcome-free quantities from discarded samples. Projected Stage A
  wall time `= n_per_rung * (t_600 + t_1200 + t_2400)`. No exponent
  model is frozen at all: v1.3 priced every rung at bottom-rung
  cost, v1.4's linear `1+2+4` was falsified by v1.4's own `O(N^2)`
  eligible-pool scan, and the correct lesson is to project from
  measured costs, not from the next model. If the projection exceeds
  12 hours, INFEASIBLE-AS-DESIGNED is recorded and the design
  returns to the drawing board without any outcome having been read.

**Stage P-B (power protocol frozen now; estimator by addendum).**
`d_hat` carries locator variability the timelike estimator does not,
so Stage B is powered from ITS OWN variance, never Stage A's (v1.2,
review): a second variance-only, quarantined pilot — BOTH endpoint
rungs, `n_pilot = 200` spacelike samples each (v1.8), under the Stage B pair
protocol (v1.6: the same both-variances rule as Stage P) — feeds the
formula of 1.2 with the ADDENDUM'S target effect
`Delta*_B` (v1.4, review: a merely-consistent `d_hat` may converge
slower than the chain rate, and reusing `-0.2007` would then
underpower the gate by construction; the addendum must derive its
estimator's rate and the implied `Delta*_B`, constraint (1) in
Section 3), giving Stage B its own `n_per_rung`, floor 12 / cap 60,
and the measured-cost feasibility stop above. Quarantine as in
Stage P: per-rung SDs and times only, cross-rung statistics
forbidden, samples never reused. The
`d_hat` definition and `Delta*_B` this pilot runs are the Stage B
addendum's — the power protocol SHAPE here does not change when that
addendum lands, which is constraint (3) on it.

### 1.4 Verdict logic (frozen for every stage gate)

Each stage's gate reads the 95% bootstrap CI `[lo, hi]` of its
`Delta` (4000 resamples, `_stable_seed` labels fresh to this
experiment). Rows are evaluated IN ORDER and the first match is the
verdict — the categories are mutually exclusive by precedence:

| # | condition | verdict |
|---|---|---|
| 1 | `hi < 0` | **IMPROVES** — gate passes (a small established improvement still passes; its size against theory is the labelled slope/constant check, not the gate) |
| 2 | `lo > 0` | **DEGRADES** — gate fails, informative |
| 3 | `-0.067 < lo` and `hi < +0.067` (CI straddles 0 inside the margin) | **FLAT-WITHIN-MARGIN** — gate fails, informative; row available only when the pilot declared the flat verdict affordable (1.2) |
| 4 | otherwise | **INCONCLUSIVE** — gate fails; when the flat verdict was declared AVAILABLE this indicts the pilot variance estimate; when it was declared unavailable ex ante, the verdict reads **UNRESOLVED** and indicts nothing (the margin was priced and not purchasable at the cap — 1.2). Record and stop either way |

The equivalence margin `delta_eq = |Delta*|/3 = 0.067 dex` is frozen
NOW, before any data, so that a flat result is a verdict rather than a
shrug — the vocabulary P10 lacked until its post-hoc rounds. A
DEGRADES verdict from THIS instrument (a consistent estimator family)
would be a major anomaly and is left representable.

**Interpretation scope (frozen, v1.4).** The verdicts bind to this
estimator family over this ladder, whose tested intervals
(`m ~ 24-96`) are pre-asymptotic. **IMPROVES** supports the
existence claim of Section 2 — an instrument through which
continuum-unit accuracy demonstrably grows over the tested densities,
with the theory rate as corroboration. **FLAT-WITHIN-MARGIN**,
**DEGRADES**, and **INCONCLUSIVE** are statements about this
estimator in this finite-density regime — grounds to extend the
ladder or change the estimator, and NEVER evidence that the
continuum limit itself is absent. P10's closure discipline, adopted
here before data instead of after review.

Completeness is part of every gate: a sample that cannot supply all
`K` pairs under the frozen pair rule is recorded incomplete and is
never scored; the campaign fills to `n` complete samples per rung by
the frozen first-`n`-complete rule of Section 4, whose skip counts
are published with the gate and whose skip cap stops the campaign
(v1.6). No verdict conditions silently on completeness — the P10
Stage B lesson, honoured by structure rather than by an
all-or-nothing rule the pair geometry cannot guarantee.

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
rank, with sub-Poisson, negatively associated count fluctuations.
Interval means, stated exactly (v1.5): conditioned on a pair of
endpoint events with integer rank gaps `a = Delta i`, `b = Delta pi`,
the OPEN interval's expected count is `(a-1)(b-1)/(N-2)` — about
1.3% below the asymptotic `N du dv = ab/N` at the bottom rung's
typical gaps — and every `N tau^2 / 4` in this document is that
asymptotic approximation, adequate for effect-size targets (the
`4^(-1/3)` ratio is unchanged) and labelled wherever it feeds a
finite-`N` quantity (the volume estimator's constant-level check
carries the conditioned mean). The chain law needs no translation — the LIS theorems
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
  own direction. Interior expected cardinality `m = N tau^2 / 4`
  (asymptotic; the endpoint-conditioned exact mean is Section 3's
  `(a-1)(b-1)/(N-2)`); the
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
  (1) a consistency argument for `d_hat -> d_true` as `N` grows,
  WITH a convergence-rate argument and the derived Stage-B target
  effect `Delta*_B` it implies over the ladder (v1.4 — consistency
  alone could hide a rate slower than the chain's, silently
  underpowering the gate); (2) a disjointness support that six pairs
  can realize at the bottom rung, demonstrated on synthetic checks
  that are quarantined from Stage B's seeds; (3) compatibility with
  the Stage P-B power protocol SHAPE of 1.3 unchanged, `Delta*_B`
  being the one number it supplies.
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
- `K = 6` pairs per sample. Candidate pool and draw semantics
  (v1.4, frozen — review's 1000-sample check showed uniform endpoint
  proposals with a 200-draw cap complete only ~34% of samples): the
  pool is the PRECOMPUTED set of all band-eligible pairs (relatedness
  plus the `tau_true` band, an `O(N^2)` scan); draws are uniform
  without replacement from that pool by the sample's own scoring
  stream (Section 5), accepted greedily if the support is disjoint
  from all previously accepted supports; **the 200-draw cap counts
  these eligible-pool draws only** (i.e. rejections for support
  overlap), else the sample is INCOMPLETE. Campaign-level handling
  (v1.6 — a per-sample 0.95 pin admits designs that die almost
  surely under an all-complete gate, `0.95^180 ~ 1e-4`): each block's
  samples are taken as **seeds in window order, first `n` complete**;
  the reserve slots of Section 5 supply the overflow, every skipped
  seed is recorded, and more than 20 skips in any block stops the
  campaign as INFEASIBLE-INCOMPLETE. Completeness is decided by
  pair-packing geometry BEFORE any estimator evaluation — `y` is
  never computed for an incomplete sample, by frozen order of
  operations — and per-block skip counts are published with the
  gate. The estimand, named for what it is (v1.7, review):
  first-`n`-complete estimates performance CONDITIONAL on successful
  pair packing, since completeness and the chain error share the
  permutation's geometry, and computing `y` afterward does not undo
  that selection. The frozen guard that keeps it innocuous: the
  verification pin bounds expected skips per campaign below one, and
  when the realized skip count is ZERO — the generic case — the
  conditional and unconditional estimands coincide on the data. Any
  realized skip stamps the stage summary with a **SELECTION-CAVEAT**
  label carried alongside the verdict, and the skip cap bounds how
  far the conditioning can reach before the campaign stops.
  Before Stage P may run, the synthetic completeness pin must pass:
  at least **1998 of 2000** synthetic samples complete at EVERY rung
  (bounding the expected skips per 180-sample campaign below one),
  generated from the verification-only seed block of Section 5 —
  quarantined from all experimental windows and carrying no outcome
  reading. The support: **Stage A**: the closed
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

Every block carries reserve slots for the first-`n`-complete rule of
Section 4: pilots hold 220 windows for 200 samples (v1.8), stage
rungs hold 80 windows for up to 60. Experimental windows begin at
200000 (v1.8 renumber for the larger pilots; nothing has run).

| block | base | fills | window span |
|---|---|---|---|
| design verification (non-experimental) | 190000 | 2000 per rung | 190000-195999, note below |
| Stage P pilot, N=600 | 200000 | 200 of 220 | 200000-243999 |
| Stage P pilot, N=2400 | 244000 | 200 of 220 | 244000-287999 |
| Stage A, N=600 | 288000 | <= 60 of 80 | 288000-303999 |
| Stage A, N=1200 | 304000 | <= 60 of 80 | 304000-319999 |
| Stage A, N=2400 | 320000 | <= 60 of 80 | 320000-335999 |
| Stage P-B pilot, N=600 | 336000 | 200 of 220 | 336000-379999 |
| Stage P-B pilot, N=2400 | 380000 | 200 of 220 | 380000-423999 |
| Stage B, N=600 | 424000 | <= 60 of 80 | 424000-439999 |
| Stage B, N=1200 | 440000 | <= 60 of 80 | 440000-455999 |
| Stage B, N=2400 | 456000 | <= 60 of 80 | 456000-471999 |
| Stage C blocks | 472000+ | — | frozen with Stage C's addendum |

All spans sit above every range the programme has ever used (documented
maxima: 30000-30379, 40000-40168, 41000-41059, 43000-54999); a
regression test pins pairwise-disjoint windows and freshness against
the full documented list before the pilot may run.

The verification block (v1.4; enlarged and relocated v1.6) is one
consecutive seed per sample, 2000 per rung: 190000-191999 (`N = 600`),
192000-193999 (`N = 1200`), 194000-195999 (`N = 2400`) — far above
every experimental window and free to grow if a later revision
enlarges the pin. Each sample uses a SINGLE generator for permutation
and pair draws — no derived offsets, hence no collision surface —
because these samples exist only to pass or fail the Section 4
completeness pin and are discarded; the stride-window discipline
binds experimental samples, whose streams feed frozen quantities.

## 6. Stages and gates

- **Stage P** (pilot): variance + feasibility only (1.3). No gate on
  outcomes; produces `n_per_rung`.
- **Stage A** (timelike): the primary gate of the experiment —
  `Delta` for tau_hat_chain per 1.4. Secondary, labelled consistency
  checks, never gates: (a) OLS slope of mean `y` on `log10 N`
  across the three rungs, with bootstrap CI, against the predicted
  `-1/3` (chain) and `-1/2` (volume); (b) the constant-level band
  `median relerr ~ 0.89 m^(-1/3)`, with `m` the endpoint-conditioned
  mean of Section 3 (v1.5); (c) middle-rung monotonicity (the
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
normalization pinned on constructed intervals (including the endpoint
subtraction); pair-disjointness and incompleteness paths; the
Section 4 synthetic completeness pin (>= 1998/2000 at every rung,
from the verification seed block); the first-n-complete fill and
skip-cap paths; seed-window privacy and freshness
against the full documented list; verdict-table precedence logic on
synthetic inputs, including CIs matching two rows. Stage P runs only
from a clean commit containing all of it.

## 9. Stage records (results)

### 9.1 Verification (2026-07-27)

2000/2000 complete at every rung against the 1998 pin — zero
incomplete samples anywhere, so the Section 4 conditional estimand
coincides with the unconditional one on everything that follows.
Measured wall times 0.006 / 0.022 / 0.083 s per sample. Artifact
`docs/prereg/frozen/p11/p11_verification_summary.json`, stamp
`99e7889` (regenerated at each implementation change under the
stamp-equality gate; completeness identical every time).

### 9.2 Stage P, the pilot (2026-07-27)

Both endpoint rungs, 200 samples each, zero skips. Variances
0.01187 / 0.01246 (per-sample `y` SD ~ 0.11 dex), kurtosis
`g4 = 4.01 / 4.24` — heavier-tailed than Gaussian, and the v1.9
calibration FIRED on its first live run: measured coverage of the
nominal Bonett bound was 0.914 / 0.918, short of the 0.95 target, so
`z` self-calibrated to 1.983 / 1.938. The invalid v1.7 chi-square
factor and the uncalibrated v1.8 nominal bound would both have
under-covered here; the review escalation that forced validation was
right on the data.

Power: `S^2_90 = 0.0311`, `n_sup = 9 -> floor 12`, `n_eq = 91 > cap`
— **the flat verdict was declared unavailable ex ante**; this
campaign's vocabulary is IMPROVES / DEGRADES / UNRESOLVED. Projected
Stage A wall time 1.4 s against the 12-hour stop. FEASIBLE. Artifact
`p11_pilot_summary.json` (with the 200 per-rung `y` values).

### 9.3 Stage A, the primary gate (2026-07-27): **IMPROVES**

12 complete samples per rung, zero skips, no selection caveat.
Mean `y` (log10 of the per-sample median relative proper-time
error): -0.618 / -0.637 / -0.805 — monotone in `N`, the middle rung
between the endpoints (no analogue of P10's 900 dip).

    Delta = -0.187 dex,  95% CI [-0.289, -0.091]  ->  IMPROVES

The interval sits entirely below zero, and the point estimate lands
within 7% of the frozen theoretical target `Delta* = -0.2007`.
Labelled consistency checks, none gating: the chain slope is
-0.311 with CI [-0.474, -0.150], containing the predicted `-1/3`;
the volume slope is -0.346 with CI [-0.567, -0.129], containing its
predicted `-1/2`; the constant-level check sits within ~25% of
`0.89 m_cond^(-1/3)` at every rung. Artifacts `p11_stage_a.csv`,
`p11_stage_a_summary.json`, stamp `99e7889`.

**What this establishes, exactly as scoped in 1.4**: the Section 2
existence claim is SUPPORTED over the tested densities — there is an
instrument, taking relabeling-invariant order quantities plus the
count as density calibration, whose continuum-unit proper-time
accuracy IMPROVES as the sprinkling densifies, at a rate consistent
with the longest-chain theorem. Through the P10 discriminator the
continuum limit was undecidable; through the chain estimator the
deepening is measured, with the theory exponent inside the interval.
Claims beyond the tested ladder and estimator remain unclaimed.
Next per Section 6: the Stage B addendum (spacelike), to be frozen
before any Stage B data.

## 10. Stage B addendum (frozen 2026-07-27, before any Stage B data)

Per Sections 3 and 6, this addendum freezes `d_hat`, its support, its
rate, and `Delta*_B`, meeting the three constraints in writing.

**The estimator: the dual-box chain.** For a spacelike pair `x, y`
(order them so `u_x < u_y`, hence `v_x > v_y`), the continuum meet
and join are the CORNERS `(u_x, v_y)` and `(u_y, v_x)` — determined
by the pair's own coordinates, so no discrete meet/join elements are
searched at all (the v1.2/v1.3 failure mode — Pareto frontiers
spanning the null arms — never arises). The meet-join interval is the
open rank box `B(x, y) = (u_x, u_y) x (v_y, v_x)`, and

    d_hat = L_box / sqrt(N)

with `L_box` the longest chain among the events interior to
`B(x, y)` (unanchored: the box has NO endpoint events — `x` and `y`
sit on its closed boundary — so there is no endpoint correction to
apply; through the shared definition with
`chain_counts_endpoints = False`, i.e. `L_box / sqrt(2 rho)` at
`rho = N/2`).

**Constraint (1) — consistency, rate, and the derived `Delta*_B`.**
In 1+1D the meet-join interval's proper time EQUALS the spacelike
separation: `d = 2 sqrt((u_y - u_x)(v_x - v_y))`. The box's interior
is a uniform sub-permutation with asymptotic expected count
`m = N (u_y - u_x)(v_x - v_y) = N d^2 / 4` (endpoint-conditioned
exact mean `(a - 1)(b - 1)/(N - 2)` for rank gaps
`a = N (u_y - u_x)`, `b = N (v_x - v_y)`, as in Section 3), so
`E L_box = 2 sqrt(m) = sqrt(N) d` by the same Vershik-Kerov law, with
the same BDJ `m^(-1/3)` relative bias and fluctuation. The estimator
is STRUCTURALLY IDENTICAL to Stage A's — the same LIS statistic on a
rank box of expected count `N (separation)^2 / 4` — so the rate is
the chain rate and

    Delta*_B = -(1/3) log10(4) = -0.2007 dex,

derived, not borrowed. The known asymmetry with Stage A: no endpoint
correction (no endpoints exist in the box), and the box is anchored
by construction rather than by relatedness.

**Constraint (2) — support and realizability.** The support IS the
closed dual box `[u_x, u_y] x [v_y, v_x]` — a bounding box of area
`d^2/4`, the same statistic as Stage A's timelike support, so
pairwise disjointness has the same packing behaviour that the Stage A
verification measured at 2000/2000. Pair protocol: the pool is all
spacelike pairs with `d_true` in the SAME frozen band
`[0.35, 0.45]`; `K = 6`; greedy without-replacement draws; the
200-cap counts support-overlap rejections; first-`n`-complete fill;
all unchanged from Section 4. The Stage B verification block (2000
single-stream samples per rung, completeness pin `>= 1998`, wall
times for the feasibility projection) uses consecutive seeds
600000-601999 / 602000-603999 / 604000-605999 — far above every
experimental window, quarantined, discarded.

**Constraint (3) — the power protocol.** Stage P-B runs the Section
1.3 protocol UNCHANGED (both endpoint rungs, 200 spacelike samples
each, calibrated Bonett bounds, cross-rung statistics forbidden) in
its frozen windows (Section 5), and feeds the Section 1.2 formulas
with the `Delta*_B` above. The verdict table, interpretation scope,
and stage gating of Sections 1.4 and 6 apply verbatim, with the
Stage B verification pin gating Stage P-B exactly as Section 4's pin
gated Stage P.

Implementation note: the runner lands only after this addendum
merges (design-first, as Stage A's did); the Stage B verification
pin must pass before Stage P-B may run.
