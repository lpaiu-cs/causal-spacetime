# Independent check of the Route D first measurement

Saved by the orchestrator on 2026-11-25 from the text the checker returned.
The checker is an AI agent in a fresh context. It was given the frozen
design, the instrument and the raw records at commit `f145421`, and was
told not to open the results record or any summary before it had its own
numbers. Its scripts were written outside the repository and are not
archived here. Everything below the line is the checker's text.

---

**Outcome:** every number in both `summary.json` files reproduces from the
raw records, and the instrument is sound. The "plateau" reading is fair as
a statement about this chain, with caveats ranked at the end. Nothing in
the repository was touched (`git status --porcelain` is empty).

## A. Recomputation against the summaries

All entries in `evidence/main/summary.json` and
`evidence/control/summary.json` agree with my code to 1e-15 relative. No
discrepancies. Values are in units of s1, with 95% t intervals over 12
orders.

| N | R(16N) | R(N^2/2) | R(2N^2) | ratio 2N^2/16N | description |
|---|---|---|---|---|---|
| 128 | 4.371 [4.019, 4.723] | 4.362 | 4.268 [3.762, 4.773] | 0.9735 [0.9052, 1.0418] | plateau |
| 192 | 4.381 [4.077, 4.684] | 4.434 | 4.317 [3.982, 4.653] | 0.9899 [0.9165, 1.0633] | plateau |
| 256 | 4.190 [3.825, 4.556] | 4.232 | 4.257 [3.946, 4.568] | 1.0243 [0.9541, 1.0945] | plateau, margin 0.0055 |

| N | total(2N^2) | wave(2N^2) | outer residual(2N^2) | outer total(2N^2) | thermometer(2N^2) | s1 N |
|---|---|---|---|---|---|---|
| 128 | 10.66 | 2.147 | 6.550 | 14.93 | 1.18 [0.76, 1.60] | 2.915 |
| 192 | 11.38 | 1.991 | 7.499 | 18.05 | 1.05 [0.58, 1.52] | 2.883 |
| 256 | 11.89 | 1.991 | 7.680 | 19.25 | 1.10 [0.62, 1.58] | 2.926 |

- **Chart units:** the central residual times N is 11.8 to 12.7 at every
  lag from 16N up, at all three sizes.
- **Size fit (36 points):**
  - At 16N: B = -0.0119 [-0.0414, 0.0175], c = 0, upper end 0.385.
  - At N^2/2: B = -0.0077 [-0.0373, 0.0219], c = 0, upper end 0.431.
  - Smallest c distinguishable from zero is 0.499 and 0.500.
  - My own fit on the absolute residual (no s1) gives upper ends 0.40 to
    0.52.
- **Far/near pair difference:** the summary reports the rms ratio, which the
  design does not specify. At 2N^2 it is 1.181 [1.085, 1.278] (neither),
  1.197 [1.131, 1.264] (rising) and 1.180 [1.130, 1.231] (rising). In mean
  square that is about 1.40.
- **Control, fit:** 0.9984 [0.9889, 1.0080], 0.9944 [0.9859, 1.0029],
  0.9972 [0.9907, 1.0038]; orders identical and all fits succeeded.
- **Control, chain:** 0.574, 0.757, 0.873, 1.003, then 1.078
  [0.965, 1.192] at 4N^2; SE 0.0571 against 0.0617 required. It passes with
  0.008 to spare at the upper end.
- **Provenance:** file and record hashes match both manifests.

## B. Audit of `drift_probe.py`

1. **`group_map`: CORRECT.**
   - It agrees with my own point-by-point coding of
     `phi^-1 o (lambda Lambda) o phi` to 4e-13, and p is fixed exactly.
   - 120 elements at N = 64 to 256, with rms displacement up to 0.86, gave
     zero order mismatches including the tips.
   - The log-Jacobian matches a finite-difference determinant to 2.3e-7;
     conformality holds to 2e-9 and closure to 6e-13.

2. **`fit_group` and `pair_stats`: CORRECT, with gaps.**
   - The dof factor is right for iid noise: recovered over injected is
     0.998 [0.989, 1.007] at N = 128 and 1.001 at 256.
   - The optimiser does not return the identity and shows no local minima:
     0 of 96 multi-start fits on real pairs found a lower cost, and 144 of
     144 recorded long-lag elements were recovered from an identity start.
   - GAP: a common shift is absorbed 96 to 97%. If all interior events
     shift together, the outer residual keeps only 9 to 16%, so the design's
     out-of-sample check covers only a shift of the central cluster relative
     to the outer events (which keeps 0.96 to 1.09).
   - GAP: the pair-difference bins confound a variance gradient with
     correlation. In my chains the per-event residual variance rises
     2.4-fold from the centre to the edge of the central set. That alone
     gives a far/near ratio of about 1.3 in mean square, against 1.35 to
     1.46 observed and 0.95 for iid noise through the same fit.
   - Near pairs do carry a weak excess covariance, about +0.09 of the mean
     square above the fit-induced null.
   - The residual is anisotropic (the T component is 13 to 14% of the mean
     square), so the dof factor is approximate, but it is constant in lag.

3. **`central_mask` and `outer`: CORRECT.** The mask is identical to my own
   and to `hit_probe_v3`'s on 19,200 events. Central and outer counts
   reproduce for 36 of 36 orders. Extremal events are in neither set.

4. **Schedule and unit: CORRECT, one immaterial gap.**
   - At N = 192 the stored states 256N and 512N listed in the design are
     missing; no lag needs them.
   - s1 is a sensible unit: each event moves exactly once per sweep, so
     s1^2 is the mean squared single move, and R at one sweep is 0.89 to
     1.04 in my chains.
   - Per-order s1 is noisy (CV 13 to 16%); this cancels in lag ratios.
   - Using both pairs introduces no bias: the log rms of (t, 2t) over
     (0, t) is +0.013 [-0.010, +0.036].

5. **Thermometer: GAP.**
   - The quantity is what the design says. I reproduced it to 1.5e-6 with
     my own Hessian; the stiff four are the dilation and three boosts.
   - "About 1" is not the fully relaxed value. For exact independent family
     draws (group Metropolis, 6 configurations per size) it reads
     0.93 +/- 0.07 and 1.35 +/- 0.12 for (x0, draw) pairs at N = 128 and
     256, and 1.69 +/- 0.18 and 1.96 +/- 0.30 for (draw, draw) pairs.
   - The noise floor from the residual is 0.18 at N = 128 and 0.03 at
     N = 256.
   - The Hessian is taken at a non-stationary point: the exact rotation row
     is non-zero and the lowest eigenvalue is negative in every order. This
     touches only the soft block.
   - So it shows the stiff modes near their Laplace spread, to within a
     factor of about 1.5 to 2. Rotations are still growing at 2N^2, and the
     total there is 2.2 to 2.7/sqrt(N).

6. **RNG: CORRECT.** The chain stream is invariant to how the schedule
   segments it and is independent of masks and fits. The main seeds are as
   stated. The fit control uses 9 where the formula has the dimension,
   which is harmless.

## C. Direct tests

1. **Injection (40 configurations per size).** With noise only, the
   recovered residual times N is 11.95 and 11.92 against 12 injected. Adding
   a field of rms 0.5/sqrt(N):

   | field | surviving fraction of mean square, N = 128 / 256 |
   |---|---|
   | quadratic, projected off the family | 1.02 / 1.10 (dof factor included) |
   | linear shear, local dilation, local boost | 0.93 to 1.20 |
   | raw random quadratic | 0.69 / 0.81 |
   | sinusoid | 0.86 to 1.02 |
   | iid | 0.98 / 1.07 |
   | common translation | 0.035 / 0.043 |

   R rises by a factor 1.10 to 1.13 at N = 128 and 1.17 to 1.22 at N = 256.
   So c = 0.5 emerging after 16N would have broken the plateau at every
   size. From the upper ends, the smallest late-emerging c that would show
   is about 0.33 per size, or 0.23 pooled (pooled ratio 0.996
   [0.959, 1.033]).

2. **Short real chains.**
   - Same seeds as the stage (N = 128, six orders, 4097 sweeps), with my
     own mask, group map and optimiser: 60 recorded rows reproduce to
     1.9e-7.
   - Fresh seeds at N = 128 (six orders): R(N) = 4.32 [3.62, 5.02] and
     R(16N) = 4.39 [3.64, 5.14]; total/s1 at 16N is 5.50, against 5.59
     recorded.
   - Fresh seeds at N = 64 (twelve orders): R(16N) = 4.31 [3.66, 4.96],
     residual times N = 11.0.
   - Below lag N, R is 1.0, 1.4, 1.8, 2.2, 2.8, 3.2, 3.6 at 1 to 64 sweeps,
     the same in sweeps at both sizes.

## D. Interpretation (opinion)

Yes, "about 12/N from 16N to 2N^2" is fair. It is flat to within 4% over
three decades of lag, and the pooled log-log slope is -0.001
[-0.013, 0.012]. Total displacement meanwhile doubles and is tracked out of
sample.

The strongest alternative is non-ergodicity: count-scale freedom reachable
only by joint moves through necks that single-event line moves do not cross
in 2N^2 sweeps. The internal counter-evidence is that the stiff family
modes, also collective and count-scale, relax in about 0.06 N^2 sweeps.

The cheapest measurement to reduce that doubt is the chain-free joint
linear programme the referee withheld. With all events free and the seven
family parameters pinned, maximise a central event's displacement at
N = 128 to 512 and read whether it scales as 1/N or N^(-1/2).

## What weakens the description, ranked

1. Connectivity under single-event moves is untested.
2. A fast count-scale term already present at 16N is bounded only by the
   three-size fit, at c below about 0.4 to 0.5.
3. The wave coordinate alone is not a plateau by the same rule:
   [1.00, 1.23], [0.89, 1.10], [0.94, 1.20] per size, pooled 1.06
   [0.998, 1.122].
4. Saturation takes about 8N to 16N sweeps, not "a few": R(2N^2)/R(N) is
   1.13 [1.02, 1.23] at N = 128 and 1.09 [1.03, 1.14] at 256.
5. A common shift of the whole interior at c of about 0.5 is invisible to
   the residual, the out-of-sample check and the thermometer (which would
   rise by only 0.3 to 0.4).
6. The thermometer is calibrated to no better than a factor of about 1.5
   to 2.
7. "12/N" averages 9/N to 15/N across the central set, and "rising with
   separation" is mostly that gradient.

## Not checked

- `results.md` was not opened.
- The chain law was not re-verified.
- No chain above N = 128 or beyond 32N sweeps.
- The N-scaling of the near-pair covariance.
- Commit ordering and provenance beyond the hashes.
- Convergence of the group Metropolis beyond a lag comparison.
