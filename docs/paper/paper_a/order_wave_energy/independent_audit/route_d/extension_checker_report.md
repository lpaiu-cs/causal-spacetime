# Independent check of the N = 512 extension

2026-11-25. **AI-level check; no human review.** Written by an AI checker
in a fresh context that had not seen this work, with its own code, from the
raw records. It did not read [extension_results.md](extension_results.md)
or run `analyse_extension.py`; it compared its numbers with
`evidence_extension/summary.json` only after producing them. It wrote only
to a scratch folder outside the repository.

Two checkers were started. The first was stopped by a usage limit partway
through part C. Its files were not shown to the second, which was told not
to read them. The second finished. Its report is below, lightly formatted.
Afterwards, the first checker's partial numbers were compared with the
second's: per-size ratios, all four size fits (B, its interval, upper end
of c, smallest c told from zero) and the control agree. The session that ran
the extension then recomputed, with a third script, the claims in part C
that change the wording of the results: the thermometer at `16 N` by size,
the weighted fit, the fit with the per-size mean of `s1`, and the fit-free
bound. All agree with the figures below.

## Summary

All 305 numbers recomputed agree with `summary.json` to within `6e-15`
relative, and every description is the same. The extension uses the same
estimator as the first measurement: rows reproduce bit for bit. The
extension's point is right: no count-scale term is found. But the quoted
upper ends of c rest on a negative slope B, and two descriptions need care.

## A. Recomputed against `summary.json`

| N | full ratio `R(N^2/2)/R(16 N)` | | wave ratio | | residual times N at `16 N` / `N^2/2` |
|---|---|---|---|---|---|
| 128 | 1.004 [0.930, 1.079] | plateau | 1.121 [0.985, 1.257] | neither | 12.64 / 12.61 |
| 192 | 1.015 [0.966, 1.065] | plateau | 1.022 [0.930, 1.114] | neither | 12.51 / 12.64 |
| 256 | 1.020 [0.921, 1.118] | **neither** | 0.992 [0.865, 1.119] | neither | 12.12 / 12.27 |
| 512 | 1.031 [0.999, 1.063] | plateau | 1.057 [0.971, 1.144] | neither | 12.17 / 12.54 [12.07, 13.00] |

Size fit, ordinary least squares as frozen, `sigma = 2.919`:

| | B [95%] | upper end of c | smallest c told from zero | description |
|---|---|---|---|---|
| full, `16 N` | -0.00627 [-0.01422, +0.00168] | 0.120 | 0.260 | none found |
| full, `N^2/2` | -0.00369 [-0.01165, +0.00427] | 0.191 | 0.260 | none found |
| wave, `16 N` | | 0.129 | 0.154 | none found |
| wave, `N^2/2` | | 0.106 | 0.157 | none found |

N = 512 curves: R from 3.98 to 4.25 between lag N and `N^2/2`; the family
thermometer rises from 0.25 at `16 N` to 0.91 [0.63, 1.20] at `N^2/2`;
no fit failures; one skipped move (seed 9). Control: 0.99934
[0.99644, 1.00224], all 24 order matrices identical, all fits succeeded.

Flags:

- N = 256 is "neither" at `N^2/2`, so "plateau at every size" would be false
  at this lag. (The results say so.)
- Cosmetic: the first measurement's 12.4, 12.4, 12.5 are at lag `2 N^2` and
  use `mean(R) * mean(s1 N)`; reproduced as 12.44, 12.45, 12.46. The
  per-order product the extension uses gives 12.26, 12.35, 12.30 at that
  lag. Like for like at `N^2/2` (12.61, 12.64, 12.27 against 12.54) the
  agreement is as close.

## B. Same estimator

- Code: `extension_task` differs from `main_task` only in its lags (`N^2/2`
  only, so its stored states differ), `dp.`-qualified names and the
  docstring. As a result `s1_squared` is averaged over a different set of
  states, and `sweeps` and `skipped_moves` cover a shorter chain.
- N = 64, seeds 3 and 8: all 14 shared rows bitwise identical in all 13
  fields; largest difference 0.0. The meta fields that should match do;
  `s1_squared` differs as expected and the checker's own chain reproduces
  both values exactly.
- N = 512, seeds 3 and 9: the first 4096 sweeps replayed. All six stored
  rows per seed (lags 512, 1024, 2048), the event counts and the Hessian
  eigenvalues match bit for bit, confirming the seed formula and the code.
- Manifests: twelve workers with one thread each, priority class 16384
  (below normal) everywhere. Stamps `e302056` and `458a5f2` are ancestors of
  HEAD. All five file hashes and both `records_sha256` match the files on
  disk. The first measurement's evidence and instruments are unchanged
  since committed.

## C. Is the size fit sound?

1. **Unequal spread.** The standard deviation of per-order `R^2` is 4.2 to
   5.5 at N <= 256 and about 1.5 at N = 512. Common-variance least squares
   gives `SE(B) = 0.0040`; HC3, weighted least squares and a 20 000-draw
   stratified bootstrap give 0.0025 to 0.0029, so the frozen interval is
   conservative. "None found" survives every method. At `16 N` the robust
   intervals lie entirely below zero (weighted: [-0.0113, -0.0007]): a
   resolved decrease, not a count-scale term. At `N^2/2`, weighted: upper
   end of c 0.10, smallest c 0.21.
2. **Per-order `s1`.** Dividing by a per-order `s1^2` inflates mean `R^2`
   by 3 to 5% at N <= 256 but 0.4 to 0.7% at N = 512, which pushes B down by
   about 0.0015. With the per-size mean of `s1` the upper ends of c are 0.17
   and 0.24; in chart units 0.19 and 0.26. With the per-size mean the
   decrease at `16 N` is no longer resolved, so it is fragile.
3. **Curvature.** None in the full residual (z = +0.4 and +0.3; quadratic
   p = 0.70 and 0.82). The wave coordinate at `N^2/2` is not monotone: size
   means 4.66, 4.47, 3.36, 3.95; quadratic p = 0.07; weighting flips the
   sign of B.
4. **How much could be count-scale at N = 512.**

   | c | basis | share of rms residual | share of mean square |
   |---|---|---|---|
   | 0.19 | frozen upper end at `N^2/2` | 34% | 12% |
   | 0.26 | sensitivity | 47% | about 22% |

   With no fit at all, c <= 0.55 (upper end 0.58), because the count-scale
   part cannot exceed the whole residual. A term of U2's size (about 0.7) is
   excluded at these lags.

## D. Opinion (the checker's)

- "About 12.5/N at N = 512 as at the smaller sizes": fair at `N^2/2`. At
  `16 N` the value is 12.17.
- "No size dependence over a factor 4": fair as "no count-scale increase",
  but loose as plain English. The quoted upper ends of c (0.12, 0.19) sit
  below the sensitivity (0.26) only because B came out negative. Better:
  "a count-scale term above about `0.26 / sqrt(N)` at these lags would have
  been seen".
- Wave "neither" at every size: correct by the frozen rule. No trend with N
  (slope p = 0.73); pooled 1.048 [0.996, 1.100]. The limit is precision,
  not drift.
- Most important thing still not shown: whether the chain reaches the
  configurations the order allows, and what happens beyond `N^2/2` at
  N = 512. The 1+1 control shows only about three quarters of collective
  motion at `N^2/2`.

## Ranked weaknesses

1. The longest lag at N = 512 is `N^2/2`, four times shorter than the first
   measurement's `2 N^2`, and chain reach is untested.
2. The upper ends of c rest on B < 0 and on the per-order `s1` bias; the
   sensitivity (0.21 to 0.26) is the honest bound.
3. N = 256 is "neither" at `N^2/2`, and the N = 512 plateau still allows a
   rise of 0 to 6%.
4. Lag `16 N` is at different relaxation stages across sizes (thermometer
   0.84, 0.86, 0.64 at N = 128, 192, 256 against 0.25 at N = 512), so the
   `16 N` fit mixes them. (The checker wrote "about 0.85 at N <= 256"; the
   figures here are from the recomputation.)
5. Least squares overstates `SE(B)` by 1.4 to 1.6 times; this errs on the
   conservative side.
6. Cosmetic: the 12.4 values mix lag and averaging definitions.

## Not checked

- The text of [extension_results.md](extension_results.md).
- The internals of `drift_probe.py` and `hit_probe_v3.py`, beyond whether
  segmenting the schedule changes the random stream.
- N = 512 rows beyond lag 2048 (only two orders replayed), and any N = 512
  `s1`.
- The pair-difference bins.
- Chain reach.
