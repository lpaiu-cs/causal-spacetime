# P7 FSS re-scope reconnaissance: the samplable and instrument-operable windows do not overlap

Date: 2026-07-16 KST. Status: **reconnaissance verdict — the N <= 120
re-scope route is closed as designed; nothing frozen, no G(beta)
measured, no MCMC run.** Probe: `experiments/positive_control/
p7_fss_rescope_probe.py`; full table with configuration in
`docs/p7_fss_rescope_probe_results.json` (deterministic seeds).

## Question

The P7 freeze (docs/p7_enhanced_sampling.md) left one escape route:
re-scope the FSS axis to the demonstrably samplable range — Wang-Landau
tunnels at N <= 80 and projects feasible to N ~ 120 — under a fresh
prereg. Before drafting any prereg, the prerequisite question: **does
the geometry instrument that defines G operate at those N at all?**

## Why it structurally cannot (envelope argument)

The frozen protocol requires six disjoint 25-tick chains plus >= 20
bracketed targets: at least ~170 elements before statistics are even
considered. Sharper: the expected longest increasing subsequence of a
uniform permutation is ~ `2 sqrt(N)`, so for `N <= (25/2)^2 ~ 156` even
**one** 25-tick chain typically does not exist. Across the entire
samplable range (N <= ~120, comfortably below the ~156 envelope) the
frozen instrument is therefore envelope-blocked in the
selector-independent sense: no extractor, greedy or optimal, can find a
chain that typically is not there. Above the envelope the story changes
— see the selector note under the table.

## Measured windows (uniform ensemble = the positive control, G must ~ 1)

Frozen protocol (6 x 25, >= 20 targets), 6 seeds per N:

| N | outcome |
|---|---|
| 200, 300, 400 | greedy chain extraction fails, 6/6 (selector-dependent — see note) |
| 500 | operates: 4/6 evaluate, all 4 pass, G median 0.858 |

Selector note. The N = 200–400 failures are properties of the frozen
protocol's *greedy* selector, not of the orders themselves: by Greene's
theorem a uniform permutation at these sizes can contain six disjoint
25-chains (the RSK shape's sixth row is still near `2 sqrt(N)` there),
and `select_disjoint_chains` consumes its longest chain whole, so an
optimal multi-chain selector could plausibly operate somewhere below
500. That would, however, be a *different instrument* — the greedy
selector is part of the frozen protocol — so N ~ 500 is the measured
operability boundary of the instrument **as frozen**, not a
selector-independent bound. None of this touches the conclusion: the
samplable window lies entirely below the ~156 envelope, where the
impossibility *is* selector-independent.

Re-scoped candidate specs (3-4 chains x 8-10 ticks, 20-30 targets,
scaled constraint counts), 8 seeds per cell, N = 100..160 — best cells:

| N | best spec | evaluates | passes G >= 0.5 | G median |
|---|---|---|---|---|
| 100 | 4x8_t30 | 4/8 | 2 | 0.475 |
| 120 | 3x8_t30 | 8/8 | 4 | 0.400 |
| 140 | 4x8_t30 | 8/8 | 3 | 0.075 |
| 160 | 4x10_t20 | 7/8 | 5 | 0.940 |

The bipartite crystal blocks structurally at every probed N (negative
side intact).

Verdict under two operability criteria, reported side by side so the
conclusion is not a threshold artifact:

- **Generous** (>= half of surviving runs pass): first soft-spec hit at
  N = 100 — but via cells where half of a handful of survivors pass.
  An order parameter whose beta = 0 baseline fails half the time cannot
  locate a geometry-loss crossing; the P7 design's frozen calibration
  requirement is "G ~ 1 on the uniform ensemble".
- **Calibration-grade** (G median >= 0.5 AND >= 3/4 pass — still far
  weaker than "G ~ 1"): **no re-scoped spec qualifies at any probed
  N <= 160.** The first spec to come close is 4x10_t20 at N = 160
  (median 0.94, 5/7 = 71% pass, just under the bar); the frozen
  instrument qualifies from N = 500.

## The samplable boundary, priced with measured throughput

Measured tunneling cost (docs/p7_tunneling_data, clean walks) and
benchmarked move throughput (numba path, single core; 13.9k/s at N = 60
falling to 9.6k/s at N = 120), scaling tau from the measured
tau(80) = 10.5M moves/round-trip by the fitted N^5.88, and assuming ~30
round trips for one converged WL + MUCA run:

| N | tau (moves/rt) | one run | 3 seeds, parallel |
|---|---|---|---|
| 60 | 1.3M (measured) | ~45 min | ~45 min |
| 80 | 10.5M (measured) | ~6.5 h | ~6.5 h |
| 100 | ~39M | ~1-1.5 days | ~1-1.5 days |
| 120 | ~114M | ~4 days | ~4 days |
| 160 | ~620M | ~3 weeks | ~3 weeks |

So the compute budget alone would have permitted an N <= 100 FSS
(roughly a desktop-day in parallel) — the re-scope does not die on
compute. It dies on the instrument: the smallest N where any probed
spec approaches calibration-grade operability (~160) costs ~3 weeks per
converged run, and the frozen instrument's window (>= 500) was already
priced at ~278,000 core-hours.

## Verdict

**Fourth quantified elimination for the P7 program, closing the
re-scope route:** after (1) local Metropolis, (2) the replica-ladder
arithmetic, and (3) the Wang-Landau tunneling exponent at N = 600, the
FSS-at-samplable-N alternative fails because the sampler-feasible
window (N <= ~120) and the instrument-operable window (N >= ~500 for
the protocol as frozen, greedy selector included — an optimal selector
could lower that boundary, but not below the ~156 envelope;
N >= ~160 for the most permissive re-scoped spec probed, and only
marginally) are disjoint. P7 remains frozen at its N = 600
characterization.

What would reopen the route (either suffices):

- a qualitatively different sampler reaching N >= 500 (the standing
  condition from docs/p7_enhanced_sampling.md), or
- a qualitatively different *small-N geometry instrument* — not a
  parameter re-scope of the gate-margin discriminator, but a new
  observable that passes its own positive control with ~10-30 targets;
  that is new instrument design with its own calibration and freeze,
  and it would answer a related but not identical question (G is
  defined by the frozen discriminator's gates).

## Deviations / scope notes

- Reconnaissance only: seeds fixed in the probe, 6-8 per cell; enough
  to establish window disjointness at the reported magnitudes, not to
  estimate pass rates precisely. Within the samplable range the blocks
  are envelope-typical (selector-independent); the N = 200-400 blocks
  are measured failures of the frozen greedy selector (see the selector
  note), and the marginal small-N cells fail the calibration bar by a
  wide margin, not a close call.
- The re-scoped specs alter constraint counts alongside chain specs;
  all such changes were made once, before any results were read, and
  are recorded in the probe's committed configuration block.
