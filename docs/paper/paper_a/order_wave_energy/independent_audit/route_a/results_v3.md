# Route A, third design: results

2026-11-24. Design frozen at `8ff268e` ([design_v3.md](design_v3.md)),
control and main stages run at `8ff268e` and `c88e371`. AI-run; the
instrument was reviewed by an AI referee before the design was amended and
has not been reviewed since; no human review. Generator coordinates: a
diagnostic of the background order law, not an observer.

## 1. Verdict of the frozen gate

> **No growth seen by this instrument.**

Between N = 256 and N = 2048, for events in the central region of the
diamond, the number of transitivity-safe flips per event grew like
`N^0.515` and the participation of the wall-hit distribution grew at most
like `N^0.03` (its slope interval just excludes zero; "flat" would overstate
it). By the user's condition this is not recorded as "cell-limited
participation": the instrument has no control in which the right answer
grows with N (Section 5).

| quantity, sizes 256 to 2048, central events | estimate | SE | 95% interval |
|---|---|---|---|
| slope of candidates per event, `kappa_F` | 0.515 | 0.008 | [0.499, 0.532] |
| slope of participation, `kappa_K` | 0.014 | 0.007 | [0.001, 0.027] |
| `D = kappa_K - kappa_F / 2` | -0.244 | 0.007 | [-0.257, -0.231] |

Gates, all passed: complete grid; no hit outside the exact-orientation
candidates and no skipped move; `kappa_F` interval above 0.25; mixing
(below); power (`SE(D) = 0.007` against the requirement `0.064`).

## 2. Numbers

Central events, twelve orders per size. K is the participation from the
gated window; the interval is the bootstrap over orders.

| N | events | candidates per event | K | 95% interval | K long / K window | upper end |
|---|---|---|---|---|---|---|
| 64 | 113 | 25.5 | 7.79 | [7.20, 8.48] | 1.018 | 1.035 |
| 128 | 212 | 33.4 | 8.98 | [8.60, 9.40] | 1.008 | 1.018 |
| 256 | 463 | 42.2 | 9.36 | [9.12, 9.61] | 1.003 | 1.006 |
| 512 | 761 | 59.2 | 9.38 | [9.19, 9.57] | 1.000 | 1.002 |
| 1024 | 768 | 86.9 | 9.64 | [9.47, 9.81] | 0.999 | 1.001 |
| 2048 | 768 | 122.3 | 9.58 | [9.45, 9.74] | (same window) | |

- Mixing gate: upper ends 1.006, 1.002, 1.001 against the ceiling 1.10;
  trend across the three sizes -0.001 against the allowance 0.05.
- The participation is the same at every block lag that was measured
  (at N = 256: 9.34 at a lag of one block, 9.36 at 63 blocks).
- Median per-event participation 9.6, 9.6, 9.9, 9.8 at the four gate
  sizes; its quartiles do not move (8.4 to 8.6 and 11.0 to 11.4).
- Distinct partners hit over a whole chain: 26.5, 27.4, 28.4, 28.3.
- Other strata, reported and not gated: non-central extremal events 9.18,
  9.64, 9.48, 9.86 with 103 to 383 candidates; non-central other events
  8.50, 9.13, 9.53, 9.52 with 47 to 161 candidates.

Control, 1+1 dimensions, stamp `8ff268e`: 3.983, 95% interval
[3.974, 3.990], for 1408 events with four distinct neighbours (exact value
4, tolerance 0.1); 2.656 [2.647, 2.678] for 68 events with a shared
neighbour (exact value 8/3). Passed.

## 3. What the chain did while the walls stayed the same

The root-mean-square displacement of the events from the sampled
configuration, with the order exact at every block end:

| N | after `2N` sweeps | at the end of the chain | end, times `sqrt(N)` | end, times N |
|---|---|---|---|---|
| 256 | 0.175 | 0.248 | 3.98 | 64 |
| 512 | 0.121 | 0.163 | 3.70 | 84 |
| 1024 | 0.086 | 0.109 | 3.50 | 112 |
| 2048 | 0.057 | 0.073 | 3.29 | 149 |

The displacement follows `N^(-1/2)`, not `1/N`: it is 60 to 150 times
`1/N`.

**Correction (after the independent check, Section 7).** The first version
of this section read that number as "configurations are loose at the count
scale while each event keeps its walls", compared it to the solvable 1+1
case, and suggested the conformal maps of the diamond as its source. The
check shows that reading is not supported:

- The number is a root mean square over **all** events and is dominated by
  non-central events. In the first sweeps central events carry 1 to 4% of
  the mean square. The displacement of central events alone was not
  recorded.
- Most of it (about 78%) is reached within the first block. A reversible
  chain cannot move along a collective family that fast: the conformal
  family can account for at most 1 to 16% of the first-block value. So the
  fast part is local motion of loosely caged events, not family motion.
  Only the slow drift afterwards has the size a freely diffusing family
  would give.
- For a central event with the others held fixed, the single-move cell has
  radius about `9.3 / N` at every size.

So this table says nothing settled about how far central events move, and
the 1+1 analogy is withdrawn as a description of these data. What the gate
measured stands on its own: the walls of central events and their weights.

## 4. What this shows, and what it does not

Shows, within the limits of Section 5:

- On the background, at these sizes, the growth of the candidates is not
  matched by any growth in how many partners share an event's walls. About
  ten partners carry the weight whether an event has 42 candidates or 122.
- The first design's counting argument (design.md Section 2.4) predicted
  this flat behaviour. That is a successful prediction of a heuristic, not
  a proof of the argument: its rigidity scale was not measured, and the
  displacement of central events was not recorded (Section 3).

Does not show:

- That the order's null-margin resolution is cell-limited. The reading
  needs a control with a growing right answer, and needs the hit
  participation to scale like the wall-density participation.
- Anything about a wave. No statistic, no likelihood, no total-variation
  statement. The exponents `8/5` and `4/3` attached to the cell-limited
  branch in next_question.md and design.md are heuristic accounting with
  the assumptions listed under item 8 of the referee report.
- That exponent 2 is not sharp. The stop rule of next_question.md would,
  on a cell-limited verdict, drop Route B and move to Routes C and D. With
  the verdict recorded as above, that step is a recommendation to the user
  and not an automatic consequence.

## 5. Limits

- **No positive control for growth.** Not repaired. The 1+1 control returns
  two different exact values, and the N = 64 to 128 runs return lower
  values than the larger sizes, so the instrument is not stuck at one
  number; neither is a case where participation is known to grow with N.
- **Reachability.** The chain uses single-event moves and may not reach
  every configuration with the given order. A component it never enters
  would lower K. The mixing gate cannot see a mode slower than the chain.
- **Hit participation is not wall-density participation** (design_v3.md
  Section 4). Their scalings agree only under an untested condition.
- One background (`q = 1`), one diamond, central events, sizes to 2048,
  where the ratio of the link scale to the counting scale is 5 to 14.
- The third design's code (stratum, exact-orientation check, longer
  chains, analysis) was written after the referee's review. Its 1+1 control
  passed, and its main-stage analysis ran for the first time on these data.
  The independent check of Section 7 was made afterwards.
- **Mixing is shown as stationarity only.** Every chain at N of 256 or more
  has 32,768 sweeps, which is `N^2 / 2` at 256 and `N^2 / 128` at 2048. A
  mode slower than that is invisible to the gate, and N = 2048 has no
  long-lag comparison at all.
- **The participation is dominated by the heaviest walls.** A growing tail
  of weak walls would not register in it. The number of distinct partners
  hit rises only from 26.5 to 28.3 over the gate sizes.

## 6. Run record

- Stamp `c88e371`, clean tree at start, 2026-11-24T09:55Z to 13:03Z.
- Twelve workers, one compute thread each. Runner and all workers at
  Windows priority class below normal (read back as 16384 by every worker
  and written to the manifest). Checked once from outside at the start:
  twelve workers below normal, using 11.2 of 16 logical processors (70%) in
  a five-second sample. The virtual-environment launcher stub of the runner
  showed class normal and uses no processor time.
- Seconds per order: 20, 88, 389, 947, 2469, 7351 for N = 64 to 2048. The
  largest size took 1.5 times the estimate.
- Evidence: [evidence_v3/main](evidence_v3/main/summary.json),
  `events.jsonl` SHA-256
  `5489b25a4e4f0edf90f51225ea4c0bbd3218b905958dac9df17dd28543a30dfb`.
- Nothing in the design was changed after the first 3+1 number was seen.

## 7. Independent check

[checker_report.md](checker_report.md). An AI agent in a fresh context,
given the frozen design, the instrument and the raw evidence, and told not
to open this record or any summary before it had its own numbers.

- **The gate reproduces.** Recomputed from the raw events with its own
  code: every point value agrees to all printed digits, intervals to the
  third or fourth decimal, and the verdict is the same.
- **No variant changes the sign.** Mean of per-event participations,
  event-weighted orders, median or geometric means, other lag windows,
  dropping either end size, a two-stage bootstrap and leave-one-order-out
  give `D` between `-0.251` and `-0.239`. `D = 0` would need the
  participation to rise by a factor 1.71 across the gate sizes; it rose by
  1.024.
- **No code defect found.** The second and third versions of the chain are
  bitwise identical in trajectory; the exact-orientation check was verified
  at 32,256 segment ends; stratum flags, candidate counts and tracked events
  were regenerated from seeds for nine orders.
- **A context for the number 9.5.** With the other events held fixed, the
  participation of the walls of a central event is about 6.3 at every size;
  along the chain it is about 9.5, with about 28 distinct partners against
  12. So the chain does turn walls over, within one block, and the result
  is as flat in N as the frozen-environment value.
- **Correction to Section 3** of this record, made above.
- **Its ranking of what could weaken the reading:** mixing shown only as
  stationarity; no positive control for growth; participation dominated by
  the heaviest walls; no long-lag check at N = 2048; the slope is slightly
  positive, so "at most `N^0.03`" and not "flat". These are in Sections 1
  and 5.
- **Its proposal for the cheapest next measurement:** rerun N = 256 only,
  recording at each block end the displacement of central events alone in
  units of the single-move cell radius, and the residual after fitting the
  seven conformal generators that fix the tips. That would separate motion
  along order-preserving families from an environment that is effectively
  frozen around central events. It is a proposal; nothing has been designed
  or run.

This check is AI review. It did not examine the linear-programme code of
the first design, the connectivity of the set of configurations with a
given order, or any chain longer than four sweeps.
