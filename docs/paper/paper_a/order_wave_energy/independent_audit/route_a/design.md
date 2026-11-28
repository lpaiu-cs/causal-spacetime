# Route A: is the order's null-margin resolution cell-limited or count-limited?

2026-11-24. Frozen design, written before any four-dimensional flip result
was computed (Section 9 lists everything that was seen first). Nothing in
this file is a result. Section 2 is heuristic analysis by one author and has
not been reviewed. Generator coordinates are used throughout: the probe is a
diagnostic of the background order law, not an observer and not a wave test.

This is Route A of [next_question.md](../next_question.md): decide the
structural unknown U1 on the background alone, before choosing which side of
the `4/3 <= alpha_low <= alpha_high <= 2` bracket to attack.

## 1. The quantity

Background `q = 1`, fixed-N iid volume sampling in the fixed diamond D,
Minkowski chart `(T, Z, y1, y2)` of
[quantized_order_bound.md](../../quantized_order_bound.md) Section 1, so the
background relation is `T_l - T_k > |z_l - z_k|` and the density is
proportional to `cos(u)^4`.

For a labelled order C let `R(C)` be the set of configurations in `D^N` whose
order is C. The conditional law of the configuration given C is the sampling
density restricted to `R(C)`. For an ordered pair `e = (k, l)` write
`m_e = (T_l - T_k) - |z_l - z_k|` for its oriented null margin.

- **Candidate.** An ordered pair whose single relation can be changed
  without violating transitivity: a link (remove it) or a critical pair (add
  it). Write `C^e` for the changed order. `F_j` is the number of candidates
  involving event j.
- **Wall.** A candidate e is a wall at configuration x when `C^e` is the
  order of a configuration reached from x by moving one endpoint of e only.
- **Realisable flip.** A candidate e for which `R(C^e)` is non-empty and
  open. Walls are realisable flips. `U_j` counts realisable flips of event j.
- **Wall density.** `rho_e(C)` is the conditional density of `m_e` at zero.
  A displacement `delta` of every causal boundary changes the order through
  pair e at first-order rate `delta * rho_e(C)`.

For any configuration the cones of the other events cut a neighbourhood of
`x_j` into cells of width of order `1/N` with a bounded number of facets
(Section 2.3), so `sum_e rho_e(C)` over the candidates of j is of order N in
every scenario. U1 is how that total is shared among the candidates of j.
Define the participation `K_j = (sum rho_e)^2 / sum rho_e^2`.

- **Cell-limited:** `K_j = O(1)`. The order fixes which O(1) pairs are
  walls; each has `rho_e` of order N.
- **Count-limited:** `K_j` of order `F_j`, about `N^(1/2)`. Every candidate
  is a wall at some typical configuration with order C; each has `rho_e` of
  order `N^(1/2)`.

Why this picks the branch (heuristic, from next_question.md Section 3): the
phase-free per-event even statistic has size `delta^2 N^2 / K` per event, so
with `K ~ N^kappa` its threshold is `N^(5 - 2 kappa) lambda^8 ~ 1`, that is
`alpha = 8 / (5 - 2 kappa)`: `8/5` at `kappa = 0` and `2` at `kappa = 1/2`.
`K_j <= U_j`, because a pair that is never a wall has `rho_e = 0`.

## 2. Local analysis (A1). Heuristic, unreviewed

### 2.1 Cones near an event are null hyperplanes

Put event j at the origin. An event at `(t, r n)` with `|t|` close to `r`
has, near the origin, the causal boundary `x^0 + n . x = s` (past events) or
`x^0 - n . x = -s` (future events), where `s = r - |t|` is its null margin.
The apex distance r drops out. The other events therefore give a Poisson
process of null hyperplanes with intensity `N Lambda(n) dOmega ds`, where
`Lambda(n)` is the density-weighted `r^2 dr` integral along the generator of
direction n inside D. The number of cones within null distance w of `x_j` is
of order `N w`.

### 2.2 Candidates per event grow like the square root of N

A related pair at null margin s and apex distance r spans an interval of
proper time squared `2 r |s|` and four-volume proportional to `r^2 s^2`, so
it is a link with probability `exp(-c N r^2 s^2)`: links have
`|s| <~ N^(-1/2)`. An unrelated pair at margin s leaves a sliver
`J^-(k) \ J^-(l)` of volume proportional to `s^2` as well (thickness s,
transverse widths `sqrt(s)` twice), so critical pairs have the same range.
Hence `F_j ~ c N^(1/2)`. Seen while timing (one seed per size, Section 9):
44.6, 62.0, 78.8, 103.0 candidates per event at N = 64, 128, 256, 512, a
log-log slope of 0.40 on this grid.

### 2.3 Moving one event: a bounded number of walls

With the others fixed, `x_j` lies in one cell of a Poisson hyperplane
tessellation. Its width is of order `1/N` and its number of facets is
bounded in N (the mean over typical cells of any stationary hyperplane
tessellation in dimension 4 is 8; the cell containing a fixed point is
larger). So the number `S_j` of walls at a given configuration is O(1), and
their margins are of order `1/N`. Both scenarios agree on this.

### 2.4 Moving all events: a counting argument for rigidity

Let `xi` be a displacement of all events. The order is unchanged to first
order iff every irredundant row c (links and critical pairs, Section 3)
satisfies `S_c g_c . xi >= -|m_c|`. Rows with `|m_c| < w` number about
`c_2 N^2 w` for `w < N^(-1/2)`, against `4N` unknowns. A cone cut out by P
half-spaces through the origin in generic position in dimension D has
non-empty interior with probability that collapses once `P > 2D` (Wendel).
So displacements of size w are available in a positive fraction of
directions only while `c_2 N^2 w <~ 8N`, that is `w <~ c/N`. A candidate at
margin far above `1/N` then cannot be brought to its cone, and

> **A1 predicts cell-limited:** `U_j = O(1)`, realisable flips have margins
> of order `1/N`, and the fraction `U_j / F_j` falls like `N^(-1/2)`.

The assumption is generic position of the rows. The continuum limit of the
row conditions, `g(Delta, xi(x + Delta) - xi(x)) = 0` for every null
`Delta`, is the conformal Killing equation, whose solution space is
15-dimensional in 3+1 dimensions. Section 2.5 shows what happens when the
assumption fails.

### 2.5 The solvable case: 1+1 dimensions

A two-dimensional order is the intersection of the two null-coordinate
rankings, and with a flat density every pair of rankings is equally likely.
A single relation changes exactly when two events adjacent in one ranking
are swapped. So `U_j <= 4` for every N, each realisable flip has
`P(C^e) = P(C)`, and every other candidate has `P(C^e) = 0` (up to orders
with more than one realiser). Wall participation is cell-limited exactly.

At the same time the position of an event given the order has spread
`N^(-1/2)`: the rows form two chains, far pairs are implied by adjacent
ones, the generic-position assumption of 2.4 fails, and the solutions of the
continuum condition are all monotone reparametrisations (the conformal group
is infinite-dimensional). **Position uncertainty and wall participation are
different questions.** This probe measures the second.

### 2.6 The competing argument

next_question.md Section 3 records why count-limited is plausible: the order
of cone crossings near an event is not witnessed by third events at scale
`1/N` (the lens between two cones holds about `1/N` events), and three
unrelated accountings return exponent 2 in that scenario. Neither argument
is a proof. That is why the question is put to a measurement.

## 3. Instrument

[flip_probe.py](flip_probe.py). Events are rows of an array; the two diamond
tips are two extra fixed rows, so "stays in D" is two more cone relations.

**Rows.** If every link of C is a relation of C' and no critical pair of C
is a relation of C', then `C' = C`: C is generated by its links, and every
unrelated ordered pair `(a, b)` of C has a critical pair `(a', b')` with
`a' <= a` and `b <= b'`, so `a < b` in C' would force `a' < b'`. So the
links and critical pairs of C and of `C^e`, minus e itself, pin every
relation except e.

**Searches.** Each ends in an exact comparison of the order matrix with
`C^e`, and requires every pinned margin and the flipped margin to exceed
`1e-9` in absolute value, so a success is an interior point of an open set.

- *single*: one endpoint of e moves; rows are all its relations.
- *joint*: all N events move. Sequential linear programming: minimise the
  linearised margin of e in a trust box, every other row kept above one
  tenth of its current margin; accept a step only if the exact order is
  still C and the margin of e fell; shrink the box otherwise; stop on
  success, after two consecutive accepted steps with relative decrease below
  `0.01`, when the box is below `1e-4 / N`, or after 40 programmes. The
  first box is `max(|m_e|, 1/N)` and may grow fourfold.
- *chain*: Metropolis line moves that keep the order exactly C, started at
  the sampled configuration. That configuration is an exact draw from the
  conditional law given C, so every state of the chain is one too. At each
  checkpoint the single search is repeated for every candidate of the tested
  events. After the last checkpoint the joint search is repeated from the
  final state ("second start").

A success is a certificate. A failed joint search proves nothing by itself.
The searches are local: a realisation of `C^e` far from every typical
configuration with order C would be missed, and would also carry no wall
density, which is the quantity U1 is about.

## 4. Campaign

Runner [run_route_a.py](run_route_a.py), CPU only, clean-tree preflight,
commit stamp, no overwriting. Seeds are
`SeedSequence([20261124, dimension, N, r])` for order `r = 0, 1, ...`.
Eight tested events per order, drawn uniformly from the N events. Every row
is retained.

| stage | geometry | N | orders per N | chain |
|---|---|---|---|---|
| control | flat 1+1 diamond | 32, 64, 128 | 8 | 2 checkpoints of 5 sweeps |
| pilot | background q = 1 | 96 | 8 | 8 checkpoints of 250 sweeps |
| main | background q = 1 | 64, 128, 256, 512 | from the pilot, Section 6 | 8 checkpoints of 250 sweeps |

The pilot size is off the main grid. Its rows are used to set the number of
main-stage orders by the rule in Section 6 and for nothing else. They are
reported. N = 1024 is left out for cost (about 40 CPU hours at the pilot's
settings); adding it later is an extension with its own freeze, not a
repair of this campaign.

## 5. Recorded quantities

Per candidate: kind, margin at the sampled configuration, which search
certified it from the sampled configuration (single-self, single-partner,
joint), programme counts, displacement of the certificate, the single-search
outcome at every chain checkpoint, and the second-start joint outcome.

Per tested event: `F_j`; `S_j` (walls at the sampled configuration); `R_j`
(certified from the sampled configuration by single or joint); `U_j`
(certified by any search); the wall sets at the nine states and their
overlaps by lag.

Per size: means of F, S, R, U; the fraction of walls still walls after a lag
of L checkpoints; the participation estimate `K_hat = M^2 / O`, with M the
mean number of walls per state and O the mean overlap at lags of at least 4
checkpoints (1000 sweeps); quartiles of `N |m_e|` for all candidates and for
certified ones; the chain's root-mean-square displacement times N.

`K_hat` is descriptive. It is a lower estimate of the participation whenever
the chain has not decorrelated, and the overlap by lag is reported so that
can be seen.

## 6. Gates

**Instrument gates.** Fixed thresholds; the decision gate is not read if
one fails.

1. *Control sensitivity.* In 1+1 dimensions the realisable flips are known
   (Section 2.5). The searches from the sampled configuration must find them:
   the lower end of the 95% two-sided Wilson interval of the found fraction
   must be at least `0.99`. With about 760 known flips this passes with up
   to two misses. Certified flips outside the rule are counted and reported.
2. *Completeness in 3+1.* Among candidates certified by the chain or by the
   second start, the fraction also certified from the sampled configuration
   must have a 95% two-sided Wilson lower end of at least `0.90`, pooled over
   the main grid. The property "`C^e` is realisable" does not depend on the
   starting point, so a deficit is a deficit of the search.
3. *The grid can discriminate.* The 95% interval of the log-log slope of the
   mean of `F_j` must lie above `0.25`.

**Decision gate.** Let `kappa_U` and `kappa_F` be the least-squares slopes
of `log mean U_j` and `log mean F_j` on `log N` over the main grid, and

`D = kappa_U - kappa_F / 2`.

Count-limited realisability predicts `kappa_U = kappa_F` and `D > 0`.
Cell-limited predicts `kappa_U = 0` and `D < 0`. The interval is the 95%
two-sided percentile interval of a cluster bootstrap that resamples orders
within each size, 10000 draws, seed 20261124.

- interval entirely below 0: **cell-limited realisability**;
- interval entirely above 0: **count-limited realisability**;
- otherwise: **undecided**.

This is a choice between two stated exponents, so it is decided on the
interval and sized for power as the house rule requires.

**Sizing for 90% power.** Allow a finite-size drift of `0.05` in `kappa_U`
against either hypothesis. At the design value `kappa_F = 0.40` the
distance from either hypothesis to the line is `0.20 - 0.05 = 0.15`, so the
requirement is `(1.960 + 1.282) SE(D) <= 0.15`, `SE(D) <= 0.0463`. On the
grid `sum (ln N - mean)^2 = 5 (ln 2)^2 = 2.402`, so the standard error of
`log mean U` per size must be at most `0.0463 * 1.550 = 0.0717`. With
`rel` the between-order standard deviation of the per-order mean of `U_j`
divided by its mean, measured in the pilot,

`orders per size = ceil((rel / 0.0717)^2)`, clipped to `[12, 48]`.

After the main stage the realised bootstrap `SE(D)` is compared with
`(kappa_F_hat / 2 - 0.05) / 3.2415`. If it is larger the verdict is
**unresolved: underpowered** whatever the interval says. Break-even
(`1.960 SE`) is not the price and is not used.

`D` computed from `R_j` instead of `U_j` is reported alongside. It is not
the gate.

## 7. What each verdict triggers, and what it does not show

- **Cell-limited realisability.** Then `K_j <= U_j = O(1)`: U1 is
  cell-limited. By the stop rule of next_question.md, Route B as stated is
  dropped (exponent 2 cannot be necessary on the heuristic accounting) and
  work moves to Routes C and D. It does not show that any exponent below 2
  is attained: that needs a statistic and a proof.
- **Count-limited realisability.** Realisable flips grow with the
  candidates. U1 is then decided by how the wall density is shared, and
  `K_hat` with its overlap-by-lag curve is read next. If `K_hat` also grows
  like `F`, the stop rule sends all effort to Route B. If `K_hat` is flat
  while U grows, the two disagree and a weights probe has to be designed
  before any route is dropped.
- **Undecided or unresolved.** Recorded as such. No route is dropped.

Heuristic consequence to keep in view, not a claim: in the cell-limited
case the 1+1 lesson of Section 2.5 does not transfer, because the conformal
group in 3+1 is finite-dimensional. If positions are then pinned to `O(1/N)`
up to that group, the "phase known" row of the accounting in
next_question.md applies and the heuristic exponent is `4/3`, the proved
lower bound. That is a target for Routes C and D, and the chain displacement
is recorded so the premise can be looked at.

## 8. Limits

- Sizes 64 to 512. The two scales `1/N` and `N^(-1/2)` differ by a factor
  8 to 23 on this grid. The gate is a slope, not a collapse.
- One background (`q = 1`), one diamond, tested events drawn uniformly, so
  boundary events are included.
- Realisability is searched locally. Certificates are exact; absences are
  not proofs. Gate 2 measures the searches against each other, not against
  the truth.
- The chain uses single-event moves. Collective modes relax slowly under
  it, and the four conformal modes that preserve the diamond do not change
  wall sets at all.
- No wave, no likelihood, no total-variation statement. The exponents in
  Sections 1 and 7 are heuristic accounting.

## 9. What was seen before this design was frozen

- The 1+1 control, in development: 563 of 563 known flips found at
  N = 32, 64, 128, none outside the rule, after one fix. The fix: an
  earlier version accepted two configurations whose flipped margin was
  exactly zero; the clearance condition of Section 3 was added.
- Self-checks of the row lemma (400 random perturbations, no mismatch) and
  of the chain keeping the order.
- A 3+1 timing run that printed, for one seed per size, the number of
  candidates per event (Section 2.2), row counts, and programme times of
  0.01 to 0.23 s. Its programme pushed one event toward the past tip and its
  solution was not printed or inspected.
- A 3+1 chain smoke test at N = 60 printed the largest coordinate change
  after 20 sweeps, 0.35. That is one loosely constrained event at small N.
  It is not a gated quantity.
- No single, joint or chain flip search was run in 3+1.
