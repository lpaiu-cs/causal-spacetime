# Referee report: Route A hit instrument

Saved by the orchestrator on 2026-11-24 from the text the referee returned.
The referee is an AI agent in a fresh context. It was given
`flip_probe.py`, `hit_probe.py`, `analyse_hits.py`, `design.md` and
`design_v2.md` at commit `c1a6474`, no 3+1 result, and instructions not to
read or run the main stages. Its scripts were written outside the
repository and are not archived here. Everything below the line is the
referee's text.

---

The core code is right: the chain samples the conditional law, the partner
rule is exact, and the estimator is unbiased. The problems are in the
"exactly 4" claim, in the gate logic for the slope, and in what gets pooled.

**Incident.** One of my background commands ran in the worktree root and
created two untracked 223-byte files there (`t2_b_exact_n4_s102.txt`,
`t2_c_broken.txt`). They existed for about eight minutes before I deleted
them. `git status --porcelain` is empty now; no tracked file was touched and
no `__pycache__` was written. A clean-tree preflight started in that window
would have refused.

**Scope.** The four stamped files match the control manifest hashes, and
three stamped control orders reproduce bitwise. I did not read
`evidence_v2/main_a`, `main_b` or `RECORD.md`, and ran no 3+1 chain-hit
statistics.

## Verdicts

| # | Item | Verdict |
|---|---|---|
| 1 | Row lemma, `links_and_critical` | CORRECT |
| 2 | Chain law; `run_block` vs `chain` | CORRECT. GAP: not ergodic on R(C). ERROR: rounding can drop a real root pair (1+1 only in practice) |
| 3 | Chord-end partner is a wall | CORRECT; the check is weaker than it could be |
| 4a | Unbiased for independent blocks | CORRECT |
| 4b | "Correlation only lowers K" | Level: holds in expectation. Applied to the slope: ERROR |
| 4c | `1 / mean Q` | GAP |
| 4d | `mixing_ratio` | GAP |
| 5 | "Exactly 4" | ERROR in the claim; estimator unbiased |
| 6 | Hits vs wall density | CORRECT within a cell; GAP for the scaling claim |
| 7 | Constants | CORRECT arithmetic; slope 0.40 not reproduced |
| 8 | alpha = 8/(5 - 2 kappa) | CORRECT algebra |
| 9 | Other | list below |

## Evidence

**1. Row lemma.** Links plus transitivity give C in C'. If (a,b) is in C'
and not in C, antisymmetry makes a, b unrelated in C. Take a' minimal in
{x <= a : x not <= b}, then b' maximal in {y >= b : y not >= a'}; (a',b') is
critical, and a' <= a <' b <= b' forces a' <' b'.

- 480 configurations (1+1 and 3+1, N = 3 to 34, tips included) against
  brute force: 0 mismatches.
- Every unrelated pair had a critical pair above it.
- All 219^2 pairs of labelled posets on 4 points: 0 violations.

**2. Chain.**

- **Root formula:** the Minkowski norm of delta - t d is
  a t^2 - 2 b t + c, as coded. Every simple root is a status change, so the
  nearest negative and positive roots bound exactly the connected set.
- **Reversibility:** the segment is the same from either end, so the
  proposal is symmetric. Restricting to it costs irreducibility, not
  stationarity.
- **`chain` vs `run_block`:** same seed gives bitwise-identical
  trajectories (3 cases). The `(1 - 1e-9)` factor and the two `|a|` branches
  are O(1e-9) effects.
- **1+1 exact law, N = 4:** 50,000 exact-start replicas; rank means
  |z| <= 1.9, variances within 1.7%.
- **1+1 exact law, N = 8:** 62,000 checkpoints; sums of u and v within
  0.36% of exact (pooled z: u +0.1, v +2.0). An early +2.8 sigma shift in u
  did not persist.
- **3+1 density acceptance (N = 5, 12):** the unconditional law is
  preserved (|z| <= 1.7). The same test with the weight removed fails at
  15.8 sigma.
- **Gap:** in 1+1, R(C) is a disjoint union of realisers and the chain never
  changes realiser.
- **Error:** `disc >= 0` drops a genuine double root when the line passes
  within about 5e-9 rad of a partner's apex. I hit it once in 6.9M moves at
  N = 8: a corner event, a line through the past tip, tip disc = -1.7e-18,
  and `chain` raised `ValueError`. In `run_block` the same case records a
  spurious hit on index 0 and proposes NaN. `FLAT2` accepts it (its
  `log_weight` ignores the point); `BACKGROUND` rejects it. In 3+1 apex
  lines have codimension 3, so there is no practical effect on the campaign.

**3. Partner = wall.** At 72,000 chord ends on sampled configurations (1+1
N = 40; 3+1 N = 40, 150), stepping just beyond the end changed exactly one
relation, with the reported partner. The flipped pair was a link or critical
pair in the exact orientation every time. The 6,704 tip ends only removed a
tip link, so only minimal or maximal events hit tips. `hi % total` is right.
The symmetrised candidate set never gives false alarms but cannot detect an
orientation error.

**4. Estimator.**

- **(a)** Non-extremal events have exactly 4N hits per block, so there is no
  ratio bias. For extremal events it is small: 1+1 boundary events give
  Q = 0.33359 against 1/3.
- **(b)** E[Q_L] = sum p^2 + sum_i Cov(f_b(i), f_{b+L}(i)). The second term
  is a variance for frozen modes or unvisited components, and was positive
  at every lag I measured. That is a statement about the level, in
  expectation, at each N. Nothing bounds how the excess varies with N, so
  the slope can be biased either way. The sentence in design_v2 Section 5
  that the count-limited reading "does not need the mixing gate" is
  therefore wrong.
- **(c)** `1/mean Q` is a harmonic mean of per-event participations,
  dominated by the lowest-participation class.
- **(d)** `mixing_ratio` only sees relaxation between about 2N and 14N
  sweeps. In the control it is 1.005 to 1.008, while rms displacement had
  reached only 81%, 64%, 59% of equilibrium (0.573/sqrt(N)) and the hit
  statistic needs about N^2 sweeps.

**5. Known answer.** The exchangeability argument is right only for four
distinct neighbours. If the u-neighbour is also a v-neighbour,
p = (1/2, 1/4, 1/4) and the participation is 8/3. The `interior` flag does
not exclude these; 50 of the 1,486 control events are of this kind.

| | K |
|---|---|
| Stated target | 4 |
| Full conditional law (twin-swapped realisers mixed), same events | 3.956 |
| Within the sampled realiser (what a perfect chain reaches) | 3.933 |
| Stamped | 3.921 |
| Stamped, four-distinct events only (1,436) | 3.987 |
| Stamped, shared-neighbour events (50) | 2.659 |

- **No chain:** 20,000 independent exact draws through the instrument's own
  partner finder give Q = 0.250000 +/- 0.000002 (K = 4.0000) for
  four-distinct events and 0.3750 for shared-neighbour events.
- **Longer chain** (N = 64, 40,960 sweeps, 320 blocks of 128): the excess
  over 1/4 for four-distinct events falls with lag.

| Lag (sweeps) | Excess over 1/4 |
|---|---|
| 128 to 256 | +1.24% |
| 512 to 896 | +0.47% |
| 2048 to 3968 | +0.09% |
| beyond 4096 | +0.01% +/- 0.01% |

So 84% of the deficit against 4 is shared neighbours, 16% is finite lag, and
none is estimator bias. Attributing it to "slow collective modes" is not
supported.

**6. Hits vs wall density.** For a convex cell,
p_i = (omega_{d-1} / |S^{d-1}|) A_i Lbar_i / V, with Lbar_i the
cosine-weighted mean chord from facet i. Three random polygons agree to
0.002. The factor is logarithmic in aspect ratio: the short sides of a
1000:1 rectangle get 2.7 times their area share. What Section 3 omits is
that rho_e weights configurations by S/V while hits weight them equally. The
scalings agree if cell size and shape given C have N-independent relative
spread, and could differ otherwise. I could not test this.

**7. Constants.**

- **(a)** Vol(D) = (pi/6) TOTAL^2 = 0.625 (Monte Carlo 0.6247) and the
  interval volume is (pi/24)(2rs)^2, so rho V = N r^2 s^2 / TOTAL^2 holds.
  The 1.09 is the e-fold point: 98% to 94% of links lie below it at
  N = 256 to 2048, with the median at 0.26 to 0.34 of it. The density weight
  (mean 0.90, range 0.59 to 1.00) moves the constant to 1.03 to 1.34.
- **(b)** The arithmetic is right. Candidates per event per unit N|m| below
  4 are 1.00 to 1.08 for N = 64 to 2048, and rows below 16/N number
  (1.01 to 1.06) 8N. The Wendel count ignores that each row touches 8 of 4N
  unknowns, so w* is good to a factor of a few.
- **Candidate counts** (8 configurations each) against the quoted single
  seeds:

| N | Mine | Quoted |
|---|---|---|
| 64 | 39.4 +/- 2.9 | 44.6 |
| 128 | 53.5 +/- 1.1 | 62.0 |
| 256 | 75.3 +/- 0.8 | 78.8 |
| 512 | 106.6 +/- 1.6 | 103.0 |

The slope is 0.48 to 0.50, not 0.40.

**8. Algebra.** N (delta^2 N^2 / K)^2 = N^(5 - 2 kappa) lambda^8, so
alpha = 8/(5 - 2 kappa), giving 8/5 and 2. Hidden assumptions:

- independent, unit-variance per-event statistics, with the signal a mean
  shift;
- sum rho_e proportional to N exactly;
- one K for all events, a pure power law;
- no cross-event coherence;
- kappa <= 1/2 (beyond that alpha exceeds the proved bound of 2);
- K is the wall-density participation, not the hit participation.

**9. Other.**

- The FAR mean averages per-lag means, so the single pair (block 0,
  block 7) carries 25% of the weight. No bias, avoidable variance.
- Zero non-tip hits in a block gives NaN, which propagates to K and D. Not
  expected: the smallest F I saw was 7.
- The percentile bootstrap on 12 orders slightly under-covers for per-size
  intervals.
- F counts partners here and ordered pairs in design.md; they differ by
  under 1% from N = 128.
- Lag slices, lag indexing, seeds, integer types: fine.

## Defects, ranked by effect on "grows" vs "flat"

1. **Two populations are pooled.** Minimal or maximal events are a large and
   shrinking share of all events, and have 2.2 to 2.6 times the F of the
   rest (3+1 samples, no chain).

   | N | Extremal share | Mean F, extremal | Mean F, others |
   |---|---|---|---|
   | 256 | 48% | 102 | 47 |
   | 512 | 38% | 161 | 72 |
   | 1024 | 29% | 256 | 105 |
   | 2048 | 22% | 397 | 154 |

   - Pooled F slope is 0.50; within classes it is 0.65 and 0.57. K is a
     harmonic mean, F an arithmetic one.
   - Composition alone moves kappa_K by +0.03, +0.05, +0.09, +0.15 if
     Q_ext/Q_bulk is 1.25, 1.5, 2, 3 (negative if the ratio is below 1).
   - "Bulk flat, extremal growing" would give D of about -0.05, which reads
     as undecided. Pure count-limited loses 0.03.
   - The sign in 3+1 is unknown to me. Stratifying by extremal status needs
     no rerun: it is recomputable from the seeds.

2. **The count-limited verdict is exempt from mixing.** If under-mixing is
   worse at small N, kappa_K is inflated: a 10% differential across the grid
   is about +0.05. The mixing gate cannot certify the opposite case either
   (item 4d).

3. **No positive control.** The 1+1 control validates the partner rule and
   the estimator at a flat answer only. Whether the instrument can return
   growth, and whether R(C) is connected in 3+1, are untested.
   Non-ergodicity lowered K by 0.6% in 1+1.

4. **The "known value 4" is mis-stated** (item 5). This is harmless for the
   3+1 estimate.

5. **`disc >= 0`** (item 2): a rare, loud failure in 1+1.

**Direction of bias in 3+1.** For the level, only downward: finite lag and
unvisited components. For the slope, either sign, from defects 1 and 2.

## Not checked

- Anything needing 3+1 hit statistics: relaxation time, the sign of
  defect 1, connectivity of R(C).
- The LP search in `flip_probe.py` (`lp_step`, `realise_flip`).
- The main-stage branch of `analyse_hits.py`, which I read but did not
  execute.
- The scaling arguments of design.md Section 2, beyond their arithmetic.
