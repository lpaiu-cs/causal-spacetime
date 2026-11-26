# Referee report on the U2 transport draft

2026-11-26. **AI-level review; no human review.** An AI referee in a fresh
context, which had not seen this work, was given the draft as committed at
`7131252` ([transport_design_draft.md](transport_design_draft.md)), the
Route D design, code, results and earlier reviews, Route C's record,
next_question.md and AGENTS.md. Its instructions:

- find what is wrong rather than confirm;
- run no 3+1 chain, and compute no transport readout on any 3+1 chain
  state;
- checks without a chain were allowed in 3+1;
- the 1+1 chain was allowed at N <= 32 with at most 16 orders, every
  number to be reported for disclosure.

Its scripts and their output are archived in
[transport_review/](transport_review/). Each script carries a three-line
header and is otherwise verbatim. The report below is the referee's,
lightly formatted.

## Outcome: do not freeze as drafted

- **(a)** Both positive descriptions can be produced without transport:
  - "growing with d" by cell-scale short-range correlation;
  - "count-scale transport found" by the N-dependent composition of the
    interior set.
- **(b)** C2's gate cannot detect loss of the field it injects.
- **(c)** The existing record already excludes the draft's motivating
  heuristic inside the central set. The draft should say so, and say what
  is left.

## Verdicts

**1. The question: a gap.**

- The posterior spread of a pair difference modulo the family is the right
  object.
- But `S_ij <= 2 I_ij`, so count-scale transport needs a count-scale
  per-event residual. The frozen central bound already limits that.
- Converted without a chain (central fit, post-fit mean square per
  `kappa^2 / N`):
  - `beta` = 0.18 for the full residual (isotropic Lévy field) and 0.090
    for the wave coordinate.
  - The extension's sensitivities, 0.26 and 0.157, then give
    `kappa <= 0.61` and `kappa_m <= 0.52` at break-even.
  - The heuristic, as a one-dimensional Brownian field in `u`
    (`Var du = |du| / N`), gives a central `c_m` of 0.29 to 0.30 at
    N = 128 to 512. That is 3.6 times the break-even B.
  - So it is excluded inside the central set at those lags. In variance,
    its constant would have to be at most 0.28 of the heuristic's.
- What is left: separations beyond the central set, and lag `2 N^2` at
  N = 512.
- **The far bin (`d >= 0.75`):**
  - an outer event is in 99.75 to 99.94% of its pairs;
  - it is 98.4 to 99.2% spacelike and only 0.5 to 0.6% near-null
    (`|tau^2| < 0.1 d^2`);
  - so it says little about pinning walls along cones (Route C).
- **The fit itself is benign.**
  - `Y` under independent noise: 0.976, 0.991 and 0.996 at N = 128, 256
    and 512.
  - With uneven noise (outer events 3 times, scale above 0.9 6 times):
    0.983, 0.993 and 0.997.
  - It absorbs 17, 13 and 12% of an isotropic Brownian far-bin field, and
    about 0% of a one-dimensional `u` field.

**2. The readout: an error in the descriptions.**

`X = 1` under independence whatever the variances: correct. But
`Y > 1.10` is neither necessary nor sufficient.

- **Cell-scale null.** A Gaussian correlation of length 0.2 (a third of
  the per-event variance) plus independent noise gives `Y` = 1.17, 1.20
  and 1.22, with `S_far N^2` flat (565, 581, 588). That is "growing with
  d" at every size, with no transport.
- **Count-like, noise `10 / N` per state.** `kappa = 1` gives `Y` = 1.085,
  1.225 and 1.444. `kappa = 0.3` gives 0.987, 1.015 and 1.047.
- **Transport shows as dependence on N.**
  - The near-bin covariance `C N^2 = (I - S) N^2 / 2` is 9, 35 and 82
    under `kappa = 1`, against 35, 44 and 49 under the correlated null.
  - The residual rise in the null is the fit's `1/n` artefact, which can
    be computed.
- **Opinion:** the central records already hint at such a near-pair
  covariance. `S_near N^2` is 244 to 280 and flat, against twice the
  central residual, about 310 to 320. The variance-gradient caveat
  applies.

**3. C3.**

- **The formula is correct.** Direct sampling gives 0.9998 ± 0.0010 at
  N = 32, and 0.9977 and 0.9993 ± 0.0010 at N = 64.
- **The premise is wrong, but harmless.**
  - The realizer is not unique: there is the reflection, and there are
    incomparable twins (0.95 per order, none in 38% of orders).
  - The chain never crosses a light cone, so it samples the starting
    realizer. That is what the formula with the starting ranks describes.
    The design should say so.
- **The pass rule is weakly informative.** At N = 32, with 16 orders and
  `8 N^2` sweeps (87 s), and bins scaled by `1/1.48`:
  - the recovered far fraction is 0.67, 0.76, 0.98, 1.00 and 1.05 at lags
    `N^2/4` to `4 N^2`, with a per-order spread of 0.44 at `N^2`;
  - the near fraction is 0.76, 0.77, 0.90, 1.06 and 1.08;
  - `Y` is 1.9 to 2.6.

  So "growing with d" is vacuous in 1+1, and `F >= 0.5` is near-certain.
- **With the bins as written (unscaled),** the far bin has 5.0 pairs per
  order among all events and 0.8 among interior events. It was empty in 5
  of 16 orders, and in 13 of 16 among interior events.
- **Power.** Under the draft's own assumptions it is 0.9997, not 0.9.
- **Reach.** C3 cannot speak to whether non-family count-scale modes mix
  in 3+1. The missing fit is not the main gap.

**4. C2.**

- The covariance is positive definite (smallest-to-largest eigenvalue
  ratio 0.0045, 0.0012, 0.0003).
- The linearisation is fine. Field-only recovery with and without the
  group element agrees within 0.03, and all 162 nonlinear fits with the
  full injection succeeded.
- **Error: the gate is on total S, which noise dominates.** The field's
  share of the expected far S is 2.1, 4.2 and 8.2% at `kappa = 0.3`, and
  19, 33 and 50% at `kappa = 1`.
  - Losing the whole field at `kappa = 0.3` still passes at every size.
  - Losing half of it at `kappa = 1` and N = 128 gives 0.905, which
    passes.
- "Expected" must include the noise, or C2 fails by construction.
- Per-draw relative spread: total S 0.04 to 0.15; field only 0.23 to
  0.43; the field part with noise at `kappa = 0.3`, 0.50 to 1.57.

**5. Power and cost.**

- **The draft's numbers reproduce under its own assumptions:** `SE(B)`
  0.076 to 0.125 for spreads of 90 to 150, giving `kappa` 0.455 and 0.21
  for the wave coordinate. But those are break-even (50% power). At 90%
  they are 0.58 and 0.27.
- **Projection (opinion).** With an interior far `S N^2` of about 700 to
  1050 and a relative spread of 0.25, `SE(B)` is about 0.2. That gives:
  - `kappa` 0.67 at break-even and 0.86 at 90% power;
  - for the wave coordinate, 0.32 and 0.41;
  - no gain over the central bound for the full residual.
- **The heuristic's scale.** U2's `d` is `|du|`. The far bin's mean
  `|dm|` is 0.29 to 0.31, against a mean `d` of 0.85 to 0.875. So the
  heuristic adds about `0.62 N` (80, 156 and 326), not `1.7 N`.
- **The size model is wrong.** The interior set drifts outward with N:

  | N | mean scale | share above scale 0.9 | outer events | mean far-bin d |
  |---|---|---|---|---|
  | 128 | 0.639 | 3.5% | 32 | 0.852 |
  | 256 | 0.682 | 4.5% | 95 | 0.863 |
  | 512 | 0.713 | 8.7% | 246 | 0.875 |

  - A cell-scale-only model, with variance proportional to
    `(1 - s)^-gamma` and calibrated to the central `155 / N^2`, gives
    `B` = 0.27, 0.93 and 2.8 for `gamma` = 0.5, 1 and 1.5. That is a
    spurious `kappa` of 0.56, 1.03 and 1.79.
  - `gamma` between 1 and 1.5 brackets the observed rise of the outer
    residual over the central one: observed 2.50 to 3.96; the model gives
    2.24 to 3.14 and 3.65 to 6.38.
  - `I_far` does not separate the two: under `kappa = 1` its B is 0.71,
    against 0.82 for S.
- **The cost table is consistent:**
  - 5.55, 8.64, 13.18 and 26.6 ms per sweep (`4.3` to `5.2e-5 N`);
  - 3.0 min, 29 min and 3.87 h per order;
  - 75.6 s per 1+1 order at N = 64;
  - the 1+1 stage at N = 128 needs about 64 min, not 55;
  - total about 6.3 to 6.5 h.

**6. Wording.**

- Section 7, "a flat result cannot trigger 'toward 4/3'": correct.
- "Signature of transport ... lower bound" (Sections 6 and 7) overclaims.
  The chain biases low, but composition and loose cells bias high.
- "Count-limited rows" names the wrong rows: the 10/7 row is cell-limited
  U1 with `d_c`.
- The flat wording is adequate. But "whatever C3 returns" contradicts
  "the main stage runs only if C3 passes".
- House rule: C2 is sized in the right form; C1 has no power line; the
  sensitivity is quoted at 50% power.

**7. Waste of CPU: yes, as drafted.** The likely outputs mislabel
cell-scale correlation and composition as transport.

## Must fix

- **M1.** Replace the `Y / 1.10` rule and "B > 0 on `S_far`". Use the
  N-dependence of the near and mid covariance `C(d) N^2`,
  projection-corrected, jointly with `S N^2`. Transport means `B_C > 0`
  with a consistent `kappa`. Power it without a chain.
- **M2.** Fix the event set geometrically, for example the sub-diamond at
  scale 0.85 or scale bands. Use the same set for the fit and the pairs,
  and report per band.
- **M3.** In C2, gate recovery of the field without noise, or as a paired
  difference at `kappa = 1`. Report `kappa = 0.3` only. State whether the
  interval is over draws or over configurations.
- **M4.** Rewrite Section 6's sensitivity: 90% power, interior spreads,
  `0.62 N`. State the existing central bound and what this run adds.
- **M5.** Correct the Section 6 and 7 claims listed under item 6.
- **M6.** Freeze the estimator.
  - Bins come from state 0: at N = 128 the `2 N^2` displacement is about
    0.25, a whole bin width.
  - `X` and `Y` as ratios of per-order means, with the interval method
    for `Y`.
  - C3's bin scaling and event set.

## Should fix

- **S1.** Bin the wave readout by `|dm|`. Add a near-null stratum for
  Route C.
- **S2.** In C3:
  - drop the "growing" clause;
  - add a 1+1 two-parameter fit that fixes the endpoints, with its exact
    projected expectation;
  - correct the premise and the power;
  - cut N = 128 to about 24 orders.
- **S3.** Add a power line for C1. Use the mean far-bin `d` per size, or
  a pair-level regression. Report a weighted fit alongside.

## Not checked

- any 3+1 chain, or the real interior covariance;
- whether the composition or cell-correlation models hold;
- C2 and C3 at full size, and C3 at N >= 64;
- chain reach;
- the code, which does not exist yet.

## Checked afterwards by the session that wrote the draft

- **The claim behind (c),** without a chain, with its own script:
  - A Brownian field in `u` with `Var du = |du| / N`, applied as a
    displacement of `m`, with the group fitted on central events, gives a
    post-fit wave-coordinate coefficient `c_m` of 0.24, 0.28 and 0.27 on
    average at N = 128, 256 and 512. The per-configuration spread is about
    0.1, so the root mean square is about 0.29.
  - The extension's wave-coordinate size fits put the upper end of that
    coefficient at 0.11 to 0.13, and could tell 0.16 from zero.
  - So a transport of the heuristic's size is not seen inside the central
    set at those lags. This agrees with the referee and is one-sided like
    every reading of this chain.
- **The near-pair covariance hint,** from the records: near-pair
  `S / s1^2` of central events is 29 to 32 at every size and lag, against
  twice the central residual, 34 to 39.
