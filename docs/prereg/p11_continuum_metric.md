# P11: a metric instrument for the continuum limit — power-first preregistration

Status: **DESIGN v1.0** (2026-07-27). Nothing below has been run.
Dates in this record are local dates (UTC+9, Asia/Seoul); commit
timestamps carry their +09:00 offset.

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
`m = E|I| = N tau^2 / 2` the interval cardinality, so at fixed
continuum separation:

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
- Feasibility rule (frozen): projected Stage A wall time
  `= 3 * n_per_rung * (pilot wall time / 12)`, and if it exceeds
  12 hours, INFEASIBLE-AS-DESIGNED is recorded and the design returns
  to the drawing board without any outcome having been read.

### 1.4 Verdict logic (frozen for every stage gate)

Each stage's gate reads the 95% bootstrap CI of its `Delta`
(4000 resamples, `_stable_seed` labels fresh to this experiment):

| CI position | verdict |
|---|---|
| entirely < 0 | **IMPROVES** — gate passes |
| inside (-0.067, +0.067) | **FLAT-WITHIN-MARGIN** — gate fails, informative |
| entirely > 0 | **DEGRADES** — gate fails, informative |
| otherwise | **INCONCLUSIVE** — gate fails; indicts the pilot variance estimate; record and stop |

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
does there exist a pure-order instrument whose continuum-unit metric
accuracy IMPROVES as sprinkling density grows, at the rate theory
predicts? Scope: 1+1D first (the dictionary regime where the
longest-chain law is a sharp theorem); higher dimensions are a listed
lever, not part of this preregistration.

## 3. The instrument (all inputs pure order)

Setting: uniform random order of size `N` == Poisson sprinkling of the
unit causal diamond at density `N`, via the frozen null-coordinate
dictionary (`order_inputs` — imported, never re-implemented; the
shared-definition rule).

- **tau_hat (timelike, PRIMARY — the chain estimator).** For related
  events `x < y`: `L(x, y)` = longest chain in the closed interval
  `I(x, y)`; then

      tau_hat_chain = L / sqrt(2 N)

  (In dictionary units the interval of proper time `tau` has expected
  cardinality `N tau^2 / 2`, and the longest chain over a uniform
  2D-order interval of `m` points converges to `2 sqrt(m)`.) Computed
  exactly: 2D dominance order makes the longest chain a longest
  increasing subsequence — patience sorting, `O(m log m)`.
- **tau_hat_vol (timelike, secondary).** `sqrt(2 |I(x, y)| / N)` —
  the volume estimator; faster predicted convergence (`m^(-1/2)`),
  less theory-rich fluctuations.
- **d_hat (spacelike, Stage B).** For spacelike `x, y`: let `M` = the
  maximal elements of `past(x) ∩ past(y)` and `J` = the minimal
  elements of `future(x) ∩ future(y)`;

      d_hat = median over (a, b) in M x J of tau_hat_chain(a, b)

  In 1+1D continuum the meet-join interval's proper time EQUALS the
  spatial separation; the median over the discrete ambiguity set is
  the frozen tie-break. Samples where `M` or `J` is empty (pair too
  close to the boundary) are excluded by the pair rule, not by the
  estimator.
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
  units, chosen a priori: `m ~ 48` at `N = 600` and `~ 192` at
  `N = 2400` — inside the law's working range, away from both the
  discreteness floor and the diamond boundary.
- `K = 6` pairs per sample, drawn by the sample's own scoring stream
  (Section 5), accepted only if their closed intervals are pairwise
  disjoint (kills within-sample dependence); up to 200 rejection
  draws, else the sample is INCOMPLETE.
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
| Stage B blocks | 100000+ | — | frozen with Stage B's addendum |

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
- **Stage B** (spacelike): same contrast for d_hat. Superiority gate
  only — the meet-join construction adds an `O(N^(-1/2))` locator
  noise whose pre-asymptotic mixing makes an exponent gate
  overconfident; the slope is recorded, labelled.
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
`mu_1 ~ -1.77` (Baik-Deift-Johansson). Hence

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
extraction by 2D dominance (`O(N)` per pair), LIS by patience sorting,
the three estimators, the pair protocol, stage runners P/A. Tests
before any run: LIS against brute force on small posets; tau_hat
normalization pinned on constructed intervals; pair-disjointness and
incompleteness paths; seed-window privacy and freshness against the
full documented list; verdict-table logic on synthetic inputs. Stage P
runs only from a clean commit containing all of it.
