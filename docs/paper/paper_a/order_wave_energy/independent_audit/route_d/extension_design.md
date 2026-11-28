# Route D, first measurement: extension to N = 512

2026-11-25. **Frozen with the commit that adds this file under this name.**
The user approved freezing and running the draft in this session
("동결하고 실행"); the draft is in the history at `1c8fed6` and its text is
unchanged apart from this paragraph and Section 6. Implemented in
[drift_extension.py](drift_extension.py) and
[analyse_extension.py](analyse_extension.py). Generator coordinates
throughout: a diagnostic of the background order law, not an observer and
not a wave test. AI-written; this design was not reviewed before the
freeze. It adds one size to the measurement of [design.md](design.md) and
changes nothing in its estimator, which was checked independently
([checker_report.md](checker_report.md)).

## 1. Why

The first measurement found the residual of central events at `12.4 / N`
for N = 128, 192 and 256 and no size dependence. Its size-dependence fit
`R^2 = A + B N` had three sizes within a factor 2, which the design called
a weak lever: a count-scale term `c / sqrt(N)` was bounded only at about
`c < 0.4` to `0.5`. One more size at twice the largest widens the lever by
the most that is affordable.

It does not address the other doubt, that the chain may not reach all
configurations with the given order. Two attempts to design a test of that
did not produce a design that could be frozen ([RECORD.md](RECORD.md)).

## 2. What is added

N = 512, twelve orders, seeds `SeedSequence([20261127, 4, 512, r])` (the
frozen formula at the new size). Chain, central and outer sets, state
pairs, group fit, residual, unit `s1`, family thermometer and intervals:
exactly as design.md Section 3, calling the same functions.

**The one difference is the chain length.** `N^2` sweeps, not `4 N^2`: at
this size `4 N^2` is about eight hours per order. Stored states: after
`N 2^j` sweeps for `j = 0, ..., 7`, after `N^2/2` and `N^2`, and one sweep
after each. Lags: `N 2^j` up to `64 N`, and `N^2/2`. The lags `N^2` and
`2 N^2` of the first measurement are not available at N = 512.

**Control.** The fit-recovery control of design.md Section 4 at N = 512, 24
configurations, same requirement (95% interval inside `1 +/- 0.05`, orders
identical, all fits succeed). If it fails the chains are not run. The 1+1
chain control is not repeated; its result stands
([results.md](results.md) Section 2), including that a chain reports about
three quarters of collective motion at lag `N^2/2`.

## 3. What is reported, and how it will be described

For N = 512: the curves of design.md Section 6 at the available lags.

Described in advance, for the full residual and separately for the wave
coordinate (which the independent check found not shown flat):

- **Plateau at N = 512 up to `N^2/2`.** The per-order ratio of R at lag
  `N^2/2` to R at lag `16 N`, mean and 95% t interval: "plateau" if the
  upper end is at most 1.10, "creep" if the lower end is above 1.10,
  "neither" otherwise. The same ratio is recomputed from the existing
  records for N = 128, 192, 256 and shown beside it. This is a shorter long
  lag than the `2 N^2` of the first measurement, and is said so wherever it
  is quoted.
- **Size dependence over four sizes.** Ordinary least squares of the
  per-order `R^2` on N over the 48 orders (36 existing, 12 new), at lag
  `16 N` and at lag `N^2/2`: B with its 95% interval, the implied
  `c = s1 N sqrt(max(B, 0))` with the pooled `s1 N` of the four sizes, and
  the smallest c the data could have told from zero. "Count-scale term
  found" if the interval of B lies above zero at either lag; "none found"
  otherwise, quoted with the upper end of c.
- The residual in chart units times N at N = 512, beside 12.4, 12.4, 12.5.

Expected sensitivity, from the spread already seen: with the fourth size
the half-width of B falls from about 0.030 to about 0.009, so the smallest
c that could be told from zero falls from 0.50 to about 0.28. This is a
planning figure, not a gate.

## 4. What it cannot establish

Everything in results.md Section 4 still applies. In particular:

- no description here is a verdict on U2;
- "plateau" is "no residual growth seen up to `N^2/2` sweeps by this
  instrument" at this size, a shorter range than before;
- "none found" bounds a count-scale term that is present at these lags and
  sizes. It says nothing about a term the chain does not reach or one that
  appears later;
- one background, one diamond, central events.

## 5. Campaign

CPU only; twelve single-thread workers at below-normal priority, recorded;
clean-tree preflight, commit stamp, no overwriting; evidence in its own
folder beside the existing one, which is not modified.

| stage | what | estimated wall time |
|---|---|---|
| 0 | fit-recovery control at N = 512 | 1 min |
| 1 | twelve chains of `N^2` sweeps at N = 512 | about 2.2 hours |

The estimate uses the 29 ms per sweep seen at N = 512 under the same load
in Route A. A run that would take more than twice the estimate is stopped
and reported, not shortened.

## 6. What was seen before this freeze

Everything in the Route A and Route D records, and one development check
of the new code:

- the schedule at N = 512 by arithmetic: lags of 1, 2, 4, 8, 16, 32, 64 and
  256 times N, last stored state at `N^2`;
- fit recovery at N = 512 on three configurations: 0.999, order matrices
  identical, every fit succeeded;
- 3+1, to see that the code path runs: two orders at N = 40 with the
  extension schedule, analysed together with the existing records. Printed:
  32 state pairs, the list of sizes, 5 seconds per order. No displacement,
  residual or summary value from that run was printed or read.

No chain at N = 512 with stored states has been run.
