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
`N^0.515` and the participation of the wall-hit distribution did not grow.
By the user's condition this is not recorded as "cell-limited
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
`1/N`. So the configurations with a given order are loose at the count
scale, and over that range of motion each event keeps the same walls with
the same weights.

This is the pattern of the solvable 1+1 case in design.md Section 2.5:
position spread `N^(-1/2)`, wall set fixed by the order. In 3+1 the exact
order-preserving maps are the conformal maps that preserve the diamond, a
finite-dimensional family whose spread is set by the sampling density and
is also of order `N^(-1/2)`. Whether the measured displacement is that
family alone, with a remainder of order `1/N`, was not measured. The chains
at N = 64 to 256 had not stopped drifting at their end, and the chain at
2048 is the shortest in units of N.

## 4. What this shows, and what it does not

Shows, within the limits of Section 5:

- On the background, at these sizes, the growth of the candidates is not
  matched by any growth in how many partners share an event's walls. About
  ten partners carry the weight whether an event has 42 candidates or 122.
- The first design's counting argument (design.md Section 2.4) predicted
  this flat behaviour. That is a successful prediction of a heuristic, not
  a proof of the argument: its rigidity scale was not measured, and the
  displacement above shows that rigidity "of positions at `1/N`" is false
  as stated unless the conformal family accounts for it.

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
  An independent recomputation is recorded in Section 7 when it exists.

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
