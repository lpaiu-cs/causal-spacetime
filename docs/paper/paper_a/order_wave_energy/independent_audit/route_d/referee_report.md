# Referee report on the first Route D draft

Saved by the orchestrator on 2026-11-24 from the text the referee returned.
The referee is an AI agent in a fresh context. It reviewed the draft at
commit `764ef85` (in the history as `route_d/design_draft.md` at that
commit), ran no 3+1 chain, and wrote its scripts outside the repository;
they are not archived here. Everything below the line is the referee's
text.

---

The family and the fit are sound, but the decision gate is not: run on the
design's own 1+1 positive control, the mixing gate returns "unresolved"
every time. As drafted the design can return "no growth" or "unresolved",
never "count-scale residual". No 3+1 chain was run; the worktree is
untouched (`git status --porcelain` empty).

| # | Item | Verdict |
|---|---|---|
| 1 | Family | CORRECT (a to d). GAP: "held at N^(-1/2)" fails for two of the six modes at these sizes |
| 2 | Fit | CORRECT for family motion plus independent noise. GAP: blind to a common shift of the central cluster |
| 3 | Unit and R | Unit CORRECT. ERROR: count-scale does not predict slope 1/2 here, so the 1/4 line and the 90% power claim do not hold as stated |
| 4 | Time scale | CORRECT in 1+1; control achievable. GAP: in 3+1, N^2 is only a lower bound |
| 5 | Gate logic | ERROR: mixing gate rejects the positive control. GAP: creep tolerance, verdict wording |
| 6 | Question | GAP: one pooled number; the Route D stop rule needs dependence on separation |
| 7 | Other | GAP: the gated estimator is not specified |

## Evidence

**1. Family**

- **Order preservation:** 153 random group elements at N = 128, 256, 512,
  rms displacement 0.1 to 0.9 (diamond height 1.48), gave an identical order
  matrix with every event inside. A 0.02 translation flips 5 entries. `phi`
  alone reverses the order exactly.
- **Symmetries:** the Fisher matrix of the seven generators has exactly one
  zero eigenvalue, the wave-axis rotation. Transverse reflections are also
  exact but discrete. The tip-swapping reversal preserves the density and
  reverses the order. There is no boost-like symmetry.
- **Amplitude (Laplace, two independent draws):** 5.24/sqrt(N) over all
  events, matching the checker's 5.2, but 3.22/sqrt(N) for central events.
  Section 7's "6 to 13 cell radii" is therefore 4 to 8 for the gated events.
- **Two soft modes:** the tilts of the wave axis have Fisher eigenvalue
  0.021 against 1.1 to 2.4 for the rest. Their angle sd is 0.99, 0.85,
  0.43 rad at the three sizes, so they are nearly free.

**2. Fit**

| N | central events (mean) | recovered/injected, noise on target | noise on both states |
|---|---|---|---|
| 128 | 18 (minimum 8) | 0.993 +/- 0.008 | 0.983 +/- 0.008 |
| 256 | 35 | 0.994 +/- 0.005 | 0.997 +/- 0.005 |
| 512 | 70 | 1.002 +/- 0.004 | 0.997 +/- 0.003 |

- **Conditioning:** condition number about 2, and no dependence on the
  central count.
- **Convergence:** 0 failures in 480 fits started at the identity, for
  rotations up to pi.
- **Composition:** the group law holds to 1e-13.
- **Absorption:** smooth non-family fields at count scale keep 75 to 82% of
  their mean square at every N, so there is no N-dependent efficiency. A
  common translation of the central cluster keeps 3%.
- **Out-of-sample:** for a truly global map the residual on non-central
  events is 1.10, 1.04, 1.02 of the noise, sharp enough to gate.

**3. Unit and R**

- **Unit:** in 1+1, s1 N is 0.546, 0.556, 0.562 at N = 32, 64, 128, and R
  has slope 0.503 from 32 to 64. I estimate the acceptance effect on s1 in
  3+1 at under 2%.
- **Shape:** R^2 = A + B N, not a power law. Taking 9.3/N as the cell-scale
  residual, slope 1/4 corresponds to a count-scale residual of 0.6/sqrt(N).
- **Power:** the stated 90% power needs about 1.8/sqrt(N), over half the
  whole family amplitude. U2's own sqrt(d/N) at d of about 0.5 is
  0.7/sqrt(N), which sits on the line.

**4. Time scale (1+1)**

| lag / N^2 | 1/8 | 1/4 | 1/2 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|
| MSD / exact, N = 32 | 0.41 | 0.55 | 0.71 | 0.88 | 0.96 | 1.00 |
| MSD / exact, N = 64 | 0.43 | 0.57 | 0.72 | 0.86 | 0.97 | 0.90 +/- 0.04 |

- **Relaxation:** the first mode relaxes in 0.62 N^2 (N = 32) and 0.60 N^2
  (N = 64). The gated lag N^2/2 is under one relaxation time.
- **Exact value:** `1/(3(N+1))` is right within a realiser; stored states
  reproduce the exact variance at 0.990 +/- 0.024.
- **Control:** its statistic came out at 0.95 pooled, not 1, with per-order
  relative SD 0.15 against the assumed 0.26. It passes.
- **3+1:** Route A's recorded late creep gives 0.019 all-event mean square
  at lag N^2/2 for N = 256, against a family equilibrium of 0.107. The
  family is at most 18% relaxed there.

**5. Gate logic**

- **Mixing ratio in 1+1, lag N^2/2 against N^2/4:** 1.137 [1.125, 1.150] at
  N = 32 and 1.127 [1.112, 1.145] at N = 64, against a ceiling of 1.10.
- **Lag N^2 against N^2/2:** 1.118 [1.101, 1.135] and 1.091 [1.067, 1.115].
- **Full replay:** the stage-2 gate on 1+1 at N = 32, 64, 128 (chains 2N^2,
  2N^2, N^2, 12 orders) gives kappa_R = 0.46, SE 0.04, D interval
  [0.14, 0.30]. The verdict is "unresolved: mixing" in 300 of 300 resampled
  campaigns.
- **What the gate needs:** a single mode passes only if its relaxation time
  is below about 0.15 N^2.
- **(a) Creep:** a residual creeping like t^0.13 in sweeps passes the gate
  and adds up to 0.27 to the slope, because the lag grows like N^2. The
  drift allowance is 0.05.
- **(a) Slow modes:** a family-sized mode that is 18% relaxed can pass
  mixing and read as "no growth". Only the wording guards that.
- **(c) N = 512:** unbiased, because the start is stationary. It has 5
  overlapping pairs and no N^2 comparison.
- **(d) Stage 1:** N = 256 has zero weight in an equally spaced three-point
  slope, so stage 1 cannot leak into kappa_R. There is no rule for what
  happens if stage 1 fails a gate.

## Changes before freezing

1. **MUST. Run the identical gate on the 1+1 control and require it to
   return "count-scale".**
   - Either gate at lag 2N^2 or more with chains of 4N^2 or more. In 1+1
     the ratios there are 1.04 to 1.06, and 12-order campaigns still pass
     only 53 to 75% of the time. At N = 512 that costs about 8.4 hours per
     order, extrapolated from Route A's timing.
   - Or drop the plateau requirement for the count-scale verdict and state
     it as a lower bound.
2. **MUST. Replace slope against 1/4.** Fit R^2 = A + B N, or word the
   verdicts as "slope above/below 1/4" and report whether the interval
   contains 0 and 1/2. State the smallest detectable count-scale
   coefficient. A slope of 0.15 +/- 0.03 must not be recorded as "no
   residual growth".
3. **MUST. Tie the mixing tolerance to the drift allowance.** That means
   about 1.02 per lag doubling, or requiring kappa_R to agree at a fixed
   number of sweeps and at the N^2-scaled lag.
4. **MUST. Freeze the estimator.**
   - the pair set, the pooling across orders, and s1 per order or pooled;
   - the optimiser, the convergence criterion, and which lags "every state
     pair" covers;
   - the central and out-of-sample sets ("tracked" is undefined here);
   - the interval method for the mixing ratios, and the seeds.
5. **SHOULD. Add a family thermometer.** Compare the fitted parameters'
   variance at the gated lag with their exact conditional law, which is
   computable without a chain, and report the fraction relaxed with the
   verdict.
6. **SHOULD. Close the translation blind spot.** Gate the out-of-sample
   residual, or fit two disjoint regions and compare the elements.
7. **SHOULD. Record what U2 needs.** Residuals by component, with the wave
   coordinate separate, and pair differences against separation. That is
   what distinguishes sqrt(d/N) from flat.
8. **SHOULD. Add a stage 1 stop rule,** and report the two half-slopes as a
   curvature check.
9. **SHOULD. Store short lags** (N to 64N sweeps) so the cell-scale plateau
   is measured directly.
10. **SHOULD. Correct the text.** Central amplitude 3.2/sqrt(N); two
    near-free modes; control expectation below 1; "three sizes" is two for
    the slope.
11. **SHOULD. Fix the interval method.** A percentile bootstrap on 12
    orders covers about 91 to 93% (analytic estimate), so use a t-based
    interval or say so.

## Not checked

- Any 3+1 chain, s1 in 3+1, or whether a count-scale residual exists.
- Connectivity of the set of configurations with a given order.
- The family amplitude beyond Laplace: importance sampling failed
  (effective sample sizes 3 to 130). Only the two soft modes were done
  exactly, each alone.
- A chain-free linear-programme bound on joint slack. I did not run it
  because it would show the answer.
- Bootstrap coverage by simulation.
- N = 64 in 1+1 has 14 orders; N = 128 has 12 chains of N^2 only.
- Any Route D code, since none exists.
