# Route A, third design. DRAFT

2026-11-24. **Draft. Not frozen. Not implemented. Not run.** It becomes a
design only when the user approves it; freezing is then a separate commit
that changes this status line and adds the code. It answers the
[referee report](referee_report.md) on [design_v2.md](design_v2.md), whose
3+1 campaign has not been run and should not be run as frozen: two of the
referee's findings (R1, R2 below) change what its decision gate would mean.

The question, definitions and hypotheses are those of
[design.md](design.md) Section 1. The instrument is that of design_v2.md
Section 3, with the changes listed here. What has been seen so far is
design.md Section 9 and design_v2.md Section 2, plus the referee's numbers
(candidate counts per event and the extremal shares in 3+1, all computed
without a chain). No 3+1 chain-hit statistic exists.

## 1. Changes, each tied to a finding

**R1. Two populations were pooled.** Minimal and maximal events are 48% of
all events at N = 256 and 22% at N = 2048, with 2.2 to 2.6 times the
candidates of the rest. A pooled slope mixes the scaling of interest with a
change of composition.

- *Change.* The gate is evaluated on **central events** only: events inside
  the concentric sub-diamond scaled by 0.6 about the midpoint of the tips,
  in the Minkowski chart. This is a fixed geometric region, so its
  population has the same definition at every N. Generator coordinates are
  already in use; the stratum is a diagnostic choice, not an observable.
  Expected central events per order: about 33 at N = 256, about 270 at
  N = 2048. Track all of them up to 64, plus 32 non-central events.
- *Reported, not gated:* the same statistics for non-central events, split
  into extremal and other.

**R2. The "count-limited" reading was exempt from the mixing gate.** That
was wrong for a slope: under-mixing that varies with N can bias the slope
in either direction. And the near/far ratio of design_v2 only sees
relaxation inside one chain.

- *Change.* Both verdicts need mixing shown.
- *Change.* Chains longer than the estimator window wherever it is
  affordable: `128 N` sweeps at N = 128 and 256, `64 N` at 512, `32 N` at
  1024, `16 N` at 2048, all in blocks of `2 N`. The estimator of the gate
  uses the same window at every size: the first 8 blocks, block pairs at
  lags 4 to 7, every pair weighted equally.
- *New mixing gate.* At N = 256, 512, 1024 compare, on the same orders and
  events, the window estimate with the estimate from block pairs at least
  half a chain apart. The 95% interval of `K_long / K_window` must have an
  upper end of at most `1.10` at each of the three sizes, and the three
  point values must not increase with N by more than `0.05` in total. At
  N = 2048 only the near/far ratio of design_v2 is available; that is a
  stated limit, not a gate that can pass.

**R3. No positive control for growth.** Not repaired. The 1+1 control has
two exact values (below), which shows the instrument can return different
values, not that it returns growth with N. Two things are reported to bear
on it: N = 64 is run as a soft-regime reference, where the pilot found most
candidates realisable; and the chain's displacement at every block end.
Connectivity of the set of configurations with a given order in 3+1 is
unknown and stays a limit.

**R4. The control's known value was mis-stated.** An event whose
u-neighbour is also its v-neighbour has participation 8/3, not 4.

- *Change.* Control classes: four distinct real rank neighbours, exact
  value 4; a shared neighbour, exact value 8/3. Gate on the first class:
  95% interval inside `4 +/- 0.1`. The second class is reported; it has
  about 50 events and is not sized as a gate. The tolerance is set knowing
  that the stamped control gave 3.987 for the first class.

**R5. Rounding can drop a double root.** *Change.* A move whose segment end
is not finite is skipped and counted; the count is reported.

**Smaller items.** The candidate check uses the exact orientation of the
crossed pair instead of the symmetrised set. `1 / mean Q` stays the gated
size-level estimate (it is what a pooled information sum would use) and the
median of the per-event participations is reported next to it. A block with
no non-tip hit is an error, not a NaN.

## 2. Campaign

CPU only, clean-tree preflight, commit stamp, no overwriting, new seed tag.
Twelve orders per size.

| stage | geometry | N | sweeps |
|---|---|---|---|
| control | flat 1+1 | 64, 128, 256 (8 orders each) | 16 N |
| main | background q = 1 | 64, 128, 256 | 128 N |
| | | 512 | 64 N |
| | | 1024 | 32 N |
| | | 2048 | 16 N |

Estimated cost: about 25 CPU hours, a little over 2 hours on 12 workers.
N = 2048 is two thirds of it.

## 3. Gates

Instrument gates: control (R4); zero hits outside the exact-orientation
candidates; the 95% interval of the slope of the mean candidate count of
central events over the gate sizes lies above `0.25`; mixing (R2).

Decision gate, sizes 256, 512, 1024, 2048, central events:
`D = kappa_K - kappa_F / 2`, 95% two-sided percentile interval, cluster
bootstrap over orders within each size, 10000 draws.

- mixing not shown: **unresolved: mixing not shown**, whatever the interval;
- interval entirely above 0: **count-limited participation**;
- interval entirely below 0: **cell-limited participation**;
- otherwise: **undecided**.

Sizing for 90% power: drift allowance `0.05`; design value
`kappa_F = 0.50` (the referee's candidate counts give 0.48 to 0.57
depending on the class); requirement
`SE(D) <= (0.25 - 0.05) / 3.2415 = 0.0617`, that is a standard error of
`log K` of at most `0.0956` per size on four sizes a factor 2 apart.
Planning value for the between-order relative standard deviation `0.174`
(the pilot's value for 8-event means; central-event means over 33 or more
events should vary less): twelve orders give `0.050`. After the run the
realised `SE(D)` is compared with `(kappa_F_hat / 2 - 0.05) / 3.2415`, and
a larger value makes the verdict **unresolved: underpowered**.

## 4. Limits that remain

- No positive control for growth with N (R3).
- Mixing at N = 2048 is not checked against a longer chain. A mode slower
  than the longest chain at any size is invisible.
- The chain may not reach every configuration with the given order. In 1+1
  that lowered the participation by 0.6%.
- Hit participation is not wall-density participation: within a cell they
  differ by a shape factor, and across configurations the density weights by
  surface over volume while hits weight equally. Their scalings agree only
  if the relative spread of cell size and shape given the order does not
  change with N. Untested.
- One background, one diamond, sizes to 2048, where the ratio of the link
  scale to the counting scale is 5 to 14.
- The exponents attached to each verdict in design.md Section 1 are
  heuristic accounting with the assumptions the referee lists under item 8.
