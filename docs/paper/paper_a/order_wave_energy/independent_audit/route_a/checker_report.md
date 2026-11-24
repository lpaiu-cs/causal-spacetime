# Independent check of the third Route A campaign

Saved by the orchestrator on 2026-11-24 from the text the checker returned.
The checker is an AI agent in a fresh context. It was given the frozen
design, the instrument and the raw evidence at commit `dc23dc2`, and was
told not to open the results record or any summary before it had its own
numbers. Its scripts were written outside the repository and are not
archived here. Everything below the line is the checker's text.

---

The gate reproduces exactly and the verdict stands: **no growth seen by
this instrument**. I found no code defect. The main weakness is that the
mixing gate shows stationarity, not global mixing.

## A. Gate recomputed from `events.jsonl`

Own code, per-order means as the bootstrap unit, 10,000 cluster draws.
`summary.json` was opened only afterwards.

| Central events | orders / events | mean F | K, mine [95%] | K, summary [95%] |
|---|---|---|---|---|
| N=256 | 12 / 463 | 42.229 | 9.3598 [9.130, 9.600] | 9.3598 [9.123, 9.606] |
| N=512 | 12 / 761 | 59.247 | 9.3806 [9.187, 9.571] | 9.3806 [9.191, 9.570] |
| N=1024 | 12 / 768 | 86.895 | 9.6413 [9.476, 9.818] | 9.6413 [9.470, 9.810] |
| N=2048 | 12 / 768 | 122.276 | 9.5805 [9.452, 9.745] | 9.5805 [9.450, 9.735] |

| | mine | summary |
|---|---|---|
| kappa_K | +0.01404, SE 0.0068, [+0.0008, +0.0275] | +0.01404, SE 0.0068, [+0.0007, +0.0273] |
| kappa_F | +0.51541, SE 0.0084, [0.4989, 0.5318] | +0.51541, SE 0.0084, [0.4993, 0.5318] |
| D | -0.24366, SE 0.0066, [-0.2567, -0.2306] | -0.24366, SE 0.0066, [-0.2568, -0.2307] |
| K_long/K_window at 256 / 512 / 1024 (upper end) | 1.0026 (1.0061) / 1.0001 (1.0015) / 0.9993 (1.0008) | same |
| mixing trend | -0.0008 | -0.0008 |
| power: SE(D) <= bound | 0.0066 <= 0.0641 | same |
| hits outside candidates / skipped moves | 0 / 0 | 0 / 0 |

No discrepancy. Points agree to every digit; intervals differ in the third
or fourth decimal from the bootstrap seed. The control (my recomputation)
gives 3.983 [3.976, 3.989] for the four-neighbour class, inside 4 +/- 0.1.

Robustness:

- **Mean of 1/q:** K = 9.83, 9.85, 10.08, 10.00; kappa_K +0.011
  [-0.003, +0.024]; D -0.247 [-0.260, -0.233].
- **Event-weighted orders:** D -0.2439 [-0.2570, -0.2305].
- **Other variants:** median or geometric K, `q_long`, `q_near`, all lags,
  dropping either end size, two-stage bootstrap and leave-one-order-out all
  give D between -0.251 and -0.239.

No variant changes the sign or the verdict. D = 0 would need
kappa_K = 0.258; K(2048)/K(256) is 1.024 against the 1.71 needed.

## B. Raw fields

- **Integrity:** the events hash matches the manifest; 5,389 events, no
  duplicates, no NaN, every q in (0,1], F >= 15 for central events at gate
  sizes.
- **One zero:** the lag-63 overlap of event 256/5/205 is exactly 0. It is a
  non-central extremal event with tip fraction 0.53, whose walls turned over
  completely; it is not in the gate.
- **Reproduction:** `q_long` is reproduced from `overlap_by_lag` to 7e-16.
  At N=2048 `q_window` and `q_near` are reproduced exactly too, because
  `q_long` is identical to `q_window` there.
- **Stationarity:** the first-8-block window does not differ from the full
  chain at the same lags (|t| <= 1.6).
- **Tracking:** tracked central events equal min(64, central) in all 72
  orders. Central events per order average 38.6, 69.9, 137 and 277.
- **Stratum:** the gate population contains 4 extremal events at N=256 and
  1 at N=512; excluding them gives D = -0.247.
- **Regeneration:** nine configurations (N = 64 to 2048) regenerated from
  their seeds reproduce the central and extremal flags, F, the counts and
  the tracked ids, using my own region and link/critical code.

## C. Code audit

- **Chain law unchanged:** v2 and v3 `run_block` give bitwise-identical
  trajectories, counts and generator state (3+1 at N = 128, 256, 512; 1+1
  at N = 64).
- **Orientation check correct in all three cases:** over 32,256 segment
  ends (4 sweeps, N = 64 to 2048), stepping beyond the end flipped exactly
  one relation, with the reported partner, in the orientation the branch
  tests. Reversing the time test flags 99.9% of unrelated ends.
- **Guards:** both skips are null moves that preserve reversibility; none
  occurred.
- **Selection:** uniform, drawn before the chain; counting consumes no
  random numbers.
- **Pair sets match the design:** window has 10 pairs, near 13, far is
  lag >= blocks/2.
- **No cap or saturation:** `partners_hit` <= F everywhere, falling from
  66% to 24% of F.
- **Process gap, immaterial:** the main analysis never reads the control
  gate.

Flatness is also what any local statistic gives. In chain-free checks
outside the design (others pinned at the sample), all three levels are flat
in N:

| Statistic | Value | Slope over 256 to 2048 |
|---|---|---|
| Participation of walls visible from the sampled point | 3.9 at every N | +0.005 |
| Participation over the visible part of the event's own frozen cell | 6.2 to 6.4 | -0.003 |
| Chain, same events | 9.4 to 9.7 | +0.013 |

## D. Displacement

| N | rms times sqrt(N) at block 1 / block 8 / last |
|---|---|
| 64 | 2.88 / 3.80 / 4.67 |
| 128 | 2.99 / 3.78 / 4.17 |
| 256 | 2.80 / 3.41 / 3.98 |
| 512 | 2.74 / 3.49 / 3.70 |
| 1024 | 2.74 / 3.45 / 3.50 |
| 2048 | 2.56 / 3.29 / 3.29 (block 8 is the last) |

Exponents over the gate sizes are -0.54 +/- 0.035 (block 1), -0.52 +/- 0.03
(block 8) and -0.59 +/- 0.03 (last block, equal 32,768 sweeps). About 78% of
the final value is reached in block 1; after that the creep exponent is at
most 0.1.

Facts bearing on the question:

- **All-event number:** in the first 4 sweeps central events carry 1 to 4%
  of the mean square.
- **Frozen cell:** a central event's frozen-cell radius is 9.3/N at every N.
- **Exact families:** the tip-fixing conformal group has 7 parameters. Its
  equilibrium amplitude is 5.2/sqrt(N), plus one unconfined rotation.
- **Kinetic bound:** a reversible chain has MSD(t) <= t MSD(1 sweep), so
  these families can supply at most 1 to 2% (N=2048) to 16% (N=256) of the
  block-1 mean square.
- **Late creep:** 0.016 and 0.0029 at N = 256 and 512, the size one freely
  diffusing family would give (0.017 to 0.032, 0.0018 to 0.0035).

**Opinion.** Neither (a) nor (b) as posed. The fast part cannot be family
motion; it is local caged motion, dominated by non-central events. It is not
literally frozen for central events either: the chain shows 28 partners at
K 9.5, against 12 at K 6.3 for the frozen environment. That turnover
completes within one block; lag-1 excess overlap is 0.5% at N=256 and 0.02%
at N=1024. Only the slow creep fits (a).

**Cheapest discriminator.** Rerun N=256 only (about 6.5 CPU-minutes per
order). At each block end record the rms displacement of central events
alone, in units of 9.3/N, and the residual after a least-squares fit of the
7 conformal generators. Under (a) the ratio is far above 1 and the fit
removes it; under (b) the ratio stays near 1.

## What could weaken "participation does not grow", ranked

1. **Mixing is shown only as stationarity.** Every chain at N >= 256 has
   32,768 sweeps: N^2/2 at 256, N^2/128 at 2048. A slower mode is invisible,
   and local statistics are flat anyway. Mitigating: at N = 64 the chain is
   2N^2 sweeps and shows no slow component above about 0.5%.
2. **No positive control for growth** (R3, unrepaired).
3. **Sum p^2 is dominated by the heaviest walls.** A growing tail of weak
   walls would not register, though distinct partners at equal hit count
   rise only from 26.5 to 28.3.
4. **N=2048 has no long-lag check.** Without it, D = -0.239
   [-0.260, -0.218].
5. **kappa_K is slightly positive.** Its interval just excludes zero, so the
   accurate wording is "at most N^0.03", not "flat".
6. **Stratum set by initial position.** Central events' later displacement
   was not recorded.
7. **Percentile bootstrap on 12 orders under-covers slightly.** The margin
   is about 37 SE.

## Not checked

- The linear-programme code in `flip_probe.py`.
- `results_v3.md`, `RECORD.md` and the control `summary.json` were not
  opened.
- No chain beyond 4 sweeps was run.
- The frozen-cell value covers only the part of the cell visible from the
  sampled point.
- The ray and frozen-cell checks used 2 to 6 orders per size, not all 12.
- Connectivity of the set of configurations with a given order.
- An unexplained 0.3 to 0.4% feature at N=128: mid-lag overlap sits below
  the long-lag level.
