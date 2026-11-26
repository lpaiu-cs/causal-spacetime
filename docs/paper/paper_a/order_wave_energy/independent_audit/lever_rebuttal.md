# Rebuttal of the lever recommendation in STATUS_AND_LEVERS.md

2026-11-26. **AI-level desk review; no human review; nothing was run.** A
fresh-context AI referee was asked to argue against the recommendation of
[STATUS_AND_LEVERS.md](STATUS_AND_LEVERS.md) as committed at `21a9892`
(lever L1, chain-free LP brackets on smooth-mode extents). It was given the
status document, next_question.md, AUDIT.md, the slack draft and its review,
and the Route A, C and D records. It was told to run no computation. The
text below is the referee's report, lightly formatted. One factual
correction by the session that wrote the status document is marked at the
end.

## Verdict: replace L1 as the top recommendation

In its own terms L1 can only ever return "rigid" or "stop". At
N = 2048 to 4096, even its rigid reading separates the cases by only a
factor of 2 to 4.

## The strongest case against L1

**1. It is one-sided, and its stop rule hides the interesting answer.**
The upper (outer) bound is valid. The lower (inner) bound cannot show
softness.

- **Maximum is not typical for smooth modes either.**
  - The ball (about 1 SD) and the cube (about 1.7) are the two most
    favourable shapes. For a convex body, extent over spread along a
    direction ranges from `sqrt 3` to about the dimension.
  - The only exactly known order cell, in 1+1, is a product of simplices.
    Along a smooth reparametrisation its extent is O(1) per event against
    a typical spread of `1/sqrt(N)`.
  - That cell is soft, so it does not refute the rigid-case claim. But
    order constraints generate exactly the shape the ball and cube
    calibration omits.
- **Curvature.**
  - The slack review found that straight steps along exactly free
    generators keep the order only to 1.0 to 5.7/N rms, which is 0.1 to
    0.5 of the 12.4/N plateau.
  - A convex inner region cannot follow a curved soft direction much past
    that chord.
  - The far-start referee's path-following repair reached only 1.5 to 2.0
    plateaus along smooth fields after 120 rounds.
- **The 1+1 figure does not transfer.** The 0.99 to 1.00 was single-event
  slack with polyhedral cones, where the inner model is exact. It says
  nothing about 3+1 curvature or about inscribed cones.
- **So in a soft world the reading would be:** inner at about one cell or
  less, outer at the box, ratio above 3, then "stop" and "recommend L6".
  "Either outcome steers" is not honest.

**Point 4: the inner region is inside the true cell,** by the row lemma,
provided:
- the row set is complete, tips included;
- each inscribed cone is built around its own link vector;
- certificates keep a clearance for solver tolerance;
- the opposite orientation of each critical pair cannot flip inside the
  box.

In a rigid world the inner side is probably tight enough: at about one
cell, the sagitta on long binding links is about 6% of a wall at N = 2048.

**2. Even the rigid side is marginal at the sizes L1 reaches.**

- **The extent per event in a rigid world.** For modes carried by central
  events, an isotropic cell gives `sqrt(N / n)`, about 2.7 typical spreads
  (`route_d/RECORD.md` Section 3), not 1. It is about 1 only for modes over
  all events, which outer and extremal events dominate and the box clips.
- **The loosening is underestimated.**
  - The document's "comparable at 1024, a third at 4096" needs pair
    separations L of about 1.
  - Central events' walls are their longest near-null pairs, which end on
    the boundary at L of about 0.35 to 0.5 (the half-height is 0.74).
  - Weighting walls by `L^2`, the loosening is about 1 to 2.3 walls at
    N = 2048 and 0.5 to 1.2 at 4096.
- **Short pairs are not the main issue.** Pairs shorter than `N^-1/4` are
  only about 4 to 7% of walls. The issue is that every central wall is
  short in absolute terms.
- **The written loosening is incomplete.**
  - `|delta y|^2 / Delta m` is only the transverse block. For the
    light-cone margin the exact remainder is
    `|delta zeta - omega delta m|^2 / (2 Delta m)`.
  - That is larger by up to `(1 + |omega|)^2` for the roughly two thirds of
    null directions with `|omega| > 1`.
  - A Minkowski-type margin avoids this. The document must say which margin
    it uses.
- **Net separation.** In a rigid world the upper bound comes out at about
  1.5 to 5.5 cells at N = 2048, against a 5-cell box and an 11-cell count
  scale. That is a factor of 2 to 4, the same "about a factor 3" that sank
  the far-start draft.
- **Enlarging the box does not help.** By the document's own constants,
  loosening over walls at a count-scale box is 2.4 to 3.2 at every N.
- **Prediction (opinion):** the N = 2048 feasibility check fires the stop
  in either world.

**3. Value even on success.**

- At N = 4096, five cells (0.015) is about the wavelength at `alpha = 2`
  (`N^-1/2 = 0.016`), and about `lambda/2` in the wave coordinate. So a
  rigid bracket at these sizes shows phase coherence for no `alpha` below
  2. Only a trend in N could, over a factor of 2 to 4, and loosening that
  falls like `1/N` biases that trend toward "rigid".
- The prior already points rigid. The proof work a rigid result would steer
  to is Route C's joint fit, which lacks a likelihood and an analytic U2.
  L1 supplies neither.
- **"Soft modes give a re-embedding coupling" is wrong twice.**
  - A coupling needs probability mass and a Jacobian, not one extremal
    configuration.
  - It needs freedom at the wave's wavelength `lambda` and amplitude
    `lambda^2`, which is below the cell scale for `alpha < 2`, not in
    low-order modes at the count scale.
- **Box locality.** The box clips non-central events: their residual is
  20 to 25/N, and the path-following LP threw them 74 to 156/N. So a
  "rigid" verdict can be an artifact of the box. Freedom reached only
  through a large non-smooth rearrangement is invisible to L1, the same
  blind spot the chain has. L1 removes the chain's slowness, not its
  locality.

## Recommended instead

The best next use of effort (opinion) is to **write the regularised
statistic into argument 3** and review it independently.

- AUDIT.md Section 4: "existence is therefore not in doubt; what is
  missing is the text".
- It is needed whatever U1 and U2 turn out to be.
- Route A's reason for deferring it is gone.
- It moves the fully written upper bound from 3 to 2.

Then **stop the U1/U2 numerics (L6)**, with L5 as a cheap option. L1
survives only as an optional rigidity check, labelled one-sided.

## Changes if L1 is kept

- **Drop the inner bound and the ratio rule.** The chain's 12.4/N already
  gives a one-sided lower side. Gate the upper bound against a stated
  threshold in cell units, with "inconclusive" as the only other outcome
  and no "soft" branch.
- **Run the feasibility check at N = 4096 as well as at 2048,** since that
  is where the claim lives.
- **Specify the margin function** and a bound on its full Hessian. Define
  "cell width".
- **Objective and modes.** Put the objective on central events,
  orthogonalise the modes against the family on the pin set, and report
  whether the result changes when the box on non-central events is
  doubled.
- **Calibrate first:** at N = 256, compare the upper bound's extent along
  each mode with the chain's spread along the same mode.
- **Positive control:** with one generator released, the upper bound must
  reach the box.
- **Re-cost.**
  - Interior-point time grew 5.5 to 7 times from N = 512 to 1024. That
    means about 5 to 10 minutes per solve at 4096 and 6 to 13 hours per
    configuration, not 2 hours.
  - Budget for a redraft: all three Route D drafts were withdrawn on
    review.
- **Delete** "the only candidate that can give a two-sided reading" and
  "either outcome steers".

## The status sections

- Section 3: the 1+1 curve is cited only as under-reporting. It is also the
  record's one positive control showing the chain reaches full collective
  motion by `2 N^2`, a lag the 3+1 runs at N <= 256 reached. Say both.
- Section 5:
  - The 0.99 to 1.00 figure omits that the cones were polyhedral, so the
    model was exact there.
  - The 9 to 12 s timings came from 3 to 4 programmes with a random or
    one-event objective.
- Section 2, Route A:
  - The review found five items, two of them gate defects, not "four
    defects".
  - The first design's pilot did run.
- "Separates rigid at a few cells from soft at the count scale": a 5-cell
  box sits below the count scale. It separates "bounded within the box"
  from "not bounded".
- The referee also wrote that "31 hours" is the cost of `N^2` sweeps and
  that `2 N^2` sweeps at N = 1024 is about 62 hours. **That is not right;
  see the correction below.**

## What would change the referee's mind

- A released-generator control at N >= 2048 in which the inner bound
  reaches the box along a curved free direction.
- An N = 256 calibration showing the upper-bound extent within about 2
  times the chain's typical spread along the same modes.
- Dual solutions with a simple, N-stable structure that could seed a
  rigidity argument. For example, row flows of O(1) hops at cell-scale
  cost.
- The user saying the goal is U2 steering rather than finishing the
  bracket. Then the reduced L1 is the best steering lever, still one-sided.

## Not settled without a run (estimates; nothing was run)

| Run | Estimated cost |
|---|---|
| Inner and upper bounds along a released generator, N = 2048, 3 generators, both signs (12 solves) | about 15 to 20 minutes |
| Loosening against margin on the rows active at the upper bound's optimum | none extra: from the duals of the run above |
| One upper-bound solve at N = 4096, for time and memory | about 5 to 10 minutes, 1 to 3 GB |
| Upper-bound extent against chain spread at N = 256 | about 1 hour on 12 workers if chain states must be regenerated; the LPs take seconds each |
| Sensitivity to doubling the box at N = 2048 | about 2 extra solves per mode, a few minutes each |

## Correction by the session that wrote the status document

The N = 512 extension's records show 262 145 sweeps per order (`N^2 + 1`)
for 6966 s, not `N^2/2` sweeps. The measured cost per sweep is
`5.2e-5 N` s at N = 512. So `2 N^2` sweeps at N = 1024 is about 31 hours
per order, as the status document says. The rebuttal's 62 hours rests on a
misreading of the sweep count. Nothing else in the rebuttal depends on it.
