# Route D, first measurement: how far do central events move at fixed order? DRAFT

2026-11-24. **Draft. Not frozen. Not implemented. Not run.** It becomes a
design only when the user approves it; freezing is then a separate commit.
Generator coordinates throughout: a diagnostic of the background order law,
not an observer and not a wave test. AI-written; not reviewed.

## 1. Why this measurement

Route A ([results_v3.md](../route_a/results_v3.md)) found that the walls
of a central event do not spread over more partners as N grows. It did not
measure how far a central event can move while the order stays the same.
An earlier version of that record read the chain's all-event displacement
as an answer and was corrected: that number is dominated by non-central
events.

This is the structural unknown U2 of [next_question.md](../next_question.md)
(phase transport): given the order, how well is the position of an event
fixed? Two things can move it.

- **The order-preserving family.** The conformal maps of the diamond that
  keep both tips fixed preserve every causal relation. They form a
  7-parameter group (Section 2). Only the rotation about the wave axis also
  preserves the sampling density; the other six are held by the density
  alone, at an amplitude of order `N^(-1/2)`.
- **Everything else,** the *residual*: what is left of the motion after the
  best member of that group is removed.

The question is the size of the residual for central events.

- **Cell-scale residual.** It stays a fixed multiple of the single-move
  scale, which is of order `1/N`. Then, up to seven global numbers, the
  order fixes central positions far below the wavelengths of the open range
  `N^(-3/4) < lambda < N^(-1/2)`.
- **Count-scale residual.** It grows relative to the single-move scale,
  like `N^(1/2)` if it sits at the count scale. Then the order does not fix
  relative positions below about `N^(-1/2)` and phase coherence is limited.

## 2. The family, exactly

In the Minkowski chart let the tips be p and q, and `I(y) = y / (y . y)`
the inversion. The map `phi(x) = I(x - q) - I(p - q)` sends q to infinity
and p to the origin, and the open diamond onto the interior of one null
cone at the origin. Events lie strictly inside the diamond, so `phi` is
regular on them. Every map

`g = phi^-1 o (lambda Lambda) o phi`, `lambda > 0`, `Lambda` an
orthochronous Lorentz transformation,

is a conformal map of the diamond onto itself that preserves the causal
order of any set of events. That is 1 + 6 = 7 parameters. The fit of
Section 3 is over this finite group, not over its linearisation: at these
sizes the family moves events by a tenth or more of the diamond, where the
second-order part of a finite map exceeds the single-move scale.

## 3. Instrument

The chain is that of Route A, unchanged in law (one sweep moves every event
once along a random line, uniformly on the segment where the order is
exact, accepted by the density ratio; started at the sampled configuration,
which is an exact draw from the conditional law given the order).

- **Chain length in units of `N^2`.** A collective coordinate shared by N
  events diffuses under single-event moves on a time of order `N^2` sweeps.
  Blocks are `N^2 / 8` sweeps; all positions are stored at every block end
  and the order is compared exactly there.
- **State pairs.** For two stored states at a lag of L blocks, the
  displacement field of the central events (central as in Route A, by the
  sampled position) is split into the best group element and a residual:
  minimise `sum_k |g(x_s,k) - x_t,k|^2` over the seven parameters by
  non-linear least squares, started at the identity.
- **Recorded per pair:** mean square displacement of central events;
  mean square residual, corrected for the seven fitted parameters; the
  fitted parameters; the residual of the non-central, non-extremal tracked
  events under the same fitted map, as an out-of-sample check.
- **Unit.** `s1(N)`: the root-mean-square displacement of central events
  over one sweep, measured at the start of every block. The gated quantity
  is `R(N) = rms residual at lag N^2/2, divided by s1(N)`.

## 4. Controls with known answers

1. **Fit recovery, no chain.** Apply a random group element of the
   expected size to sampled configurations, add independent noise of a
   known size to every event, and fit. The recovered residual must equal
   the injected noise: 95% interval of the ratio inside `1 +/- 0.05`, at
   N = 128, 256, 512. Also: the order matrix before and after a group
   element alone must be identical.
2. **The chain reaches count-scale collective motion.** In 1+1 dimensions
   with a flat density, given the two rankings the null coordinates are
   uniform order statistics, so the mean square difference between two
   independent draws, averaged over events, is exactly `1 / (3 (N + 1))`
   in the chart `(T, Z)`. Almost all of that is collective motion. Chains of `8 N^2` sweeps in 16 blocks
   at N = 64 (48 orders) and N = 128 (24 orders); the mean square
   displacement between states at least `N^2` sweeps apart, divided by the
   exact value, pooled over both sizes, must have a 95% interval inside
   `1 +/- 0.15`. This is the positive control Route A lacked for the chain:
   a case where large collective motion exists and must be found.
   Sizing: assumed relative standard deviation per order `0.26`, 72 orders,
   standard error `0.031`, against the requirement `0.15 / 3.2415 = 0.046`.
   If the realised standard error is larger the control is unresolved.

The 1+1 control validates the chain on collective modes, not the group
fit: the order-preserving family of 1+1 is infinite-dimensional. Control 1
validates the fit. Nothing validates both together where the answer is a
count-scale residual; Section 7 keeps the wording that follows from that.

## 5. Campaign, in two stages that can be approved separately

CPU only; twelve single-thread workers at below-normal priority, recorded,
as in Route A; clean-tree preflight, commit stamp, no overwriting. Twelve
orders per size. The verdict is not changed after results are seen.

| stage | what | sizes | sweeps | estimated wall time |
|---|---|---|---|---|
| 0 | controls of Section 4 | 1+1: 64, 128; fit: 128, 256, 512 | `8 N^2` | 45 min |
| 1 | first look, **no verdict** | 3+1: 256 | `2 N^2` | 30 min |
| 2 | scaling and the gate | 3+1: 128 (`2 N^2`), 512 (`N^2`) | | 2.5 h |

**Stage 1 reports and decides nothing about U2.** At N = 256: the
displacement of central events and its residual, in units of `s1`, at lags
`N^2/8` to `15 N^2/8`; the fitted parameters; the mixing ratios of
Section 6. A single size cannot separate "a fixed multiple of the cell
scale" from "a larger multiple", because the multiple is not known in
advance. Stage 1 exists to see whether the instrument works at the size
where it is cheap, before the hours of stage 2 are spent. Its numbers will
have been seen when stage 2 runs; the gate below is frozen before stage 1.

## 6. Gates (stage 2)

Instrument gates: both controls; exact order at every block end; the fit
converges for every state pair.

Mixing gate, every size: the ratio of the rms residual at lag `N^2/2` to
that at lag `N^2/4` has a 95% interval upper end of at most `1.10`; at
N = 128 and 256 the same for lag `N^2` against lag `N^2/2`. A mode slower
than the chain is invisible to this, as in Route A.

Decision gate: `kappa_R`, the least-squares slope of `log R(N)` on `log N`
over N = 128, 256, 512; `D = kappa_R - 1/4`; 95% two-sided percentile
interval, cluster bootstrap over orders within each size, 10000 draws.

- instrument gate failed, mixing not shown, or underpowered: **unresolved**,
  in that order of precedence;
- interval entirely above 0: **count-scale residual**;
- interval entirely below 0: **no residual growth seen by this instrument**
  (not "positions are fixed at the cell scale": see Section 7);
- otherwise: **undecided**.

Sizing for 90% power: the hypotheses are slopes 0 and 1/2, the line is 1/4,
drift allowance `0.05`, so `SE(D) <= 0.20 / 3.2415 = 0.0617`. On three
sizes a factor 2 apart that is a standard error of `log R` of at most
`0.0605` per size. Planning value for the between-order relative standard
deviation `0.174` (the Route A pilot's value; the realised values in Route
A were 0.03 to 0.05): twelve orders give `0.050`. After the run the
realised `SE(D)` is compared with `0.0617`.

## 7. Limits, stated before any data

- **Asymmetry.** If a count-scale residual exists and relaxes more slowly
  than `N^2 / 2` sweeps, the chain under-reports it at every size and the
  slope may or may not show it. A below-zero interval is therefore recorded
  as "not seen by this instrument". The 1+1 control shows the chain can
  equilibrate count-scale collective motion in `N^2` sweeps there; it does
  not show the same in 3+1.
- **Reachability.** Single-event moves may not reach every configuration
  with the given order.
- Central events only, one background (`q = 1`), one diamond, three sizes
  to 512. At these sizes the estimated family amplitude is about 6 to 13
  times the single-move cell radius, so the two scales are separated by
  one order of magnitude at most.
- The group fit uses central events only. With 18 central events at
  N = 128 it estimates seven parameters from 72 numbers.
- Whatever the verdict, it is about positions given the order on the
  background. It is not a statement about a wave, a statistic or the
  exponent. The heuristic consequence (a cell-scale residual would make the
  "phase known up to a few global parameters" row of next_question.md
  Section 3 the relevant one) stays heuristic.

## 8. What has been seen that bears on this

Everything in the Route A records, and from the independent check
([checker_report.md](../route_a/checker_report.md)): a central event's
single-move cell radius is about `9.3 / N` with the others held fixed; the
checker's estimate of the family's equilibrium amplitude is `5.2 / sqrt(N)`;
the all-event displacement of the Route A chains (about `3.5 / sqrt(N)`,
mostly reached within `2N` sweeps, dominated by non-central events). No
displacement of central events has been recorded or seen.
