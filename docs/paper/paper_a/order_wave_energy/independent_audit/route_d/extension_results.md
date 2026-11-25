# Route D, first measurement, extension to N = 512: results

2026-11-25. Design frozen at `e302056`
([extension_design.md](extension_design.md)); control run at `e302056`,
chains at `458a5f2`. **No verdict on U2 was defined and none is recorded.**
AI-run; no human review. Checked independently by an AI checker in a fresh
context (Section 5); the wording of Sections 1 to 3 was corrected after
that check. Generator coordinates: a diagnostic of the background order
law, not an observer.

## 1. What was seen

At N = 512 the residual of central events, after the seven-parameter
order-preserving group is fitted out, is the same as at the smaller sizes:
`12.2 / N` at lag `16 N` and `12.5 / N` at lag `N^2/2` in chart units,
beside `12.1` to `12.6` at N = 128, 192, 256. Over the four sizes no
count-scale growth with N was found, for the full residual or for the wave
coordinate alone: a term above about `0.26 / sqrt(N)` (full residual) or
`0.16 / sqrt(N)` (wave coordinate) at these lags would have been seen.

By the descriptions fixed in the design:

| | full residual | wave coordinate |
|---|---|---|
| N = 512, R at `N^2/2` over R at `16 N` | 1.031 [0.999, 1.063]: **plateau** | 1.057 [0.971, 1.144]: **neither** |
| size dependence at lag `16 N` | **none found**, B = -0.0063 [-0.0142, 0.0017] | **none found**, B = -0.0008 [-0.0036, 0.0020] |
| size dependence at lag `N^2/2` | **none found**, B = -0.0037 [-0.0117, 0.0043] | **none found**, B = -0.0016 [-0.0045, 0.0013] |

"Plateau" here means no residual growth seen up to `N^2/2` sweeps by this
instrument, a shorter range than the `2 N^2` of the first measurement.

## 2. Numbers

Twelve orders per size; 95% t intervals; R in units of the one-sweep
displacement `s1`.

| N | `s1 N` | R at `16 N` | R at `N^2/2` | ratio | chart residual times N at `N^2/2` |
|---|---|---|---|---|---|
| 128 | 2.92 | 4.37 [4.02, 4.72] | 4.36 [4.04, 4.69] | 1.004 [0.930, 1.079] | 12.6 [11.6, 13.6] |
| 192 | 2.88 | 4.38 [4.08, 4.68] | 4.43 [4.14, 4.72] | 1.015 [0.966, 1.065] | 12.6 [11.8, 13.4] |
| 256 | 2.93 | 4.19 [3.83, 4.56] | 4.23 [3.83, 4.64] | 1.020 [0.921, 1.118] | 12.3 [11.0, 13.6] |
| 512 | 2.95 | 4.13 [4.01, 4.25] | 4.25 [4.14, 4.36] | 1.031 [0.999, 1.063] | 12.5 [12.1, 13.0] |

The first three rows are recomputed from the existing records at the same
lag. By the same rule N = 128 and 192 are "plateau" and N = 256 is
"neither" at this lag, its upper end being 1.118; at lag `2 N^2` the first
measurement had it as "plateau" with 1.0945.

Wave coordinate, R at `16 N` and at `N^2/2`: 1.91 and 2.13 (N = 128), 2.02
and 2.06 (192), 1.93 and 1.83 (256), 1.88 and 1.98 (512). The ratio is
"neither" at every size: 1.12 [0.99, 1.26], 1.02 [0.93, 1.11],
0.99 [0.87, 1.12], 1.06 [0.97, 1.14].

**Size dependence, `R^2 = A + B N`, 48 orders.** The implied count-scale
coefficient c (residual `c / sqrt(N)` in chart units) is zero in every
fit. Its upper end and the smallest value the data could have told from
zero:

| | upper end of c | smallest c distinguishable from zero |
|---|---|---|
| full, lag `16 N` | 0.12 | 0.26 |
| full, lag `N^2/2` | 0.19 | 0.26 |
| wave, lag `16 N` | 0.13 | 0.15 |
| wave, lag `N^2/2` | 0.11 | 0.16 |

With three sizes the corresponding figures were 0.39 to 0.43 and 0.50.

The upper ends lie below the smallest distinguishable values only because
every slope came out negative. The smallest distinguishable value is the
bound to quote. Without any fit, c cannot exceed 0.55 (upper end 0.58) at
N = 512 and lag `N^2/2`, because a count-scale part cannot exceed the whole
residual.

**At N = 512 along the lags.** Total displacement of central events grows
from 4.1 to 11.7 in units of `s1` between lag N and lag `N^2/2`, while R
stays between 3.9 and 4.3. Outer events: total from 8.1 to 18.4, residual
under the same fitted map between 8.0 and 8.5. The stiff family parameters
read 0.91 [0.63, 1.20] of their predicted spread at `N^2/2` (a thermometer
calibrated only to a factor of 1.5 to 2, see
[results.md](results.md) Section 3). About 71 central and 236 outer events
per order. No fit failed. One move was skipped by the rounding guard.

Control: fit recovery at N = 512, 0.9993 [0.9964, 1.0022], order matrices
identical, every fit succeeded.

## 3. What this does and does not establish

Adds to the first measurement:

- at lag `N^2/2` the residual follows `1 / N` over a factor 4 in N, not 2.
  At lag `16 N` the sizes are at different stages of relaxation (Section 5),
  so that lag is the weaker comparison;
- a count-scale term at these lags would have been seen above about a
  quarter of `1 / sqrt(N)` for the full residual and about a sixth for the
  wave coordinate, against a half before. The size U2's own estimate
  suggested was about 0.7, which also exceeds the fit-free bound of 0.55.

Does not establish, as before:

- that positions are fixed at the cell scale. Every limit of results.md
  Section 4 applies: the reading is one-sided, the chain may not reach all
  configurations with the given order, and two attempts to design a test of
  that did not yield one that could be frozen ([RECORD.md](RECORD.md));
- that the residual at N = 512 is mostly cell-scale. A count-scale part
  just below the sensitivity (c = 0.26) would be about 47% of its rms and
  about 22% of its mean square;
- anything beyond lag `N^2/2` at N = 512. The 1+1 control showed that a
  chain reports about three quarters of collective motion at that lag;
- that the wave coordinate is flat in the lag. It is "neither" at every
  size. What the size fit shows for it is no growth with N, which is a
  different statement;
- anything about a wave, a statistic or the exponent.

## 4. Run record

- Chains: stamp `458a5f2`, clean tree, 2026-11-25T09:12Z to 11:08Z (116
  minutes; estimate 132). Seconds per order: 6966.
- Twelve workers, one compute thread each, runner and workers at Windows
  priority class below normal (read back as 16384, in the manifest).
  Checked once from outside at the start: twelve workers below normal,
  11.2 of 16 logical processors (70%).
- `records.jsonl` SHA-256 begins `2020e15afdcde347`; the full value is in
  the manifest.
- Nothing in the design was changed after the first number was seen.

## 5. Independent check

[extension_checker_report.md](extension_checker_report.md). An AI checker
in a fresh context, given the designs, the instruments and the raw records,
and told not to open this record or any summary until it had its own
numbers. A first checker was stopped by a usage limit partway through; a
second, which did not see the first one's files, finished. The first one's
partial numbers agree with the second's. AI-level check; no human review.

- **Every number reproduces.** 305 values, to `6e-15` relative, with every
  description the same.
- **Same estimator.** The extension's task differs from the first
  measurement's only in its schedule. At N = 64 both reproduce each other's
  shared rows bit for bit. At N = 512 the first 4096 sweeps of two orders
  were replayed and the stored rows matched bit for bit. Seeds, workers,
  priority class and file hashes match the manifests.
- **The size fit errs on the conservative side.** The N = 512 orders spread
  about a third as much as the others. Weighted, robust and bootstrap
  intervals of B are 1.4 to 1.6 times narrower, and "none found" holds in
  all of them. At lag `16 N` the weighted interval lies below zero. That is
  a decrease, not a count-scale term, and it is not robust: it goes away
  when the per-size mean of `s1` is used.
- **What the record had to correct.** The upper ends of c rest on negative
  slopes; the sensitivity (0.21 to 0.26 for the full residual, depending on
  the fit) is the honest bound. Dividing by a per-order `s1` lowers B by
  about 0.0015; with the per-size mean of `s1` the upper ends become 0.17
  and 0.24. At lag `16 N` the stiff-family thermometer reads 0.84, 0.86 and
  0.64 at N = 128, 192 and 256 against 0.25 at N = 512, so that lag
  compares different stages of relaxation; at `N^2/2` it reads 0.9 to 1.6.
  Sections 1 to 3 were reworded to say these things. The table entries and
  the frozen descriptions are unchanged.
- **Wave coordinate.** "Neither" at every size is correct by the frozen
  rule. The ratio shows no trend with N (slope p = 0.73); pooled over the
  four sizes it is 1.048 [0.996, 1.100]. Its size means at `N^2/2` are not
  monotone (4.66, 4.47, 3.36, 3.95; quadratic p = 0.07), and weighting
  flips the sign of its slope.
- **Not checked.** The text of this record; the instrument's internals
  beyond the schedule; N = 512 rows beyond lag 2048; the pair-difference
  bins; chain reach.
