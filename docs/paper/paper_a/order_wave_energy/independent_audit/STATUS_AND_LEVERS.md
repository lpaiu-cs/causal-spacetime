# Is exponent 2 sharp? Where the campaign stands, and which lever is left

2026-11-26. Written at the user's request ("정리 후 레버 점검", relayed by the
session that started this one): a short status of the whole campaign, then
an assessment of the means that could make its readings two-sided, with one
recommendation. **No new computation was run for this document.** It reads
the records and estimates costs from them. AI-written. A fresh-context AI
reviewer was asked to argue against the first recommendation, and it
succeeded ([lever_rebuttal.md](lever_rebuttal.md)). Sections 2, 3 and 5 were
corrected after it, and the recommendation was replaced (Section 6). The
version before the rebuttal is in the history at `21a9892`.

## 1. The question and what is proved

Exact circular vacuum plane wave, fixed known `0 < q <= 1`, N events
sampled i.i.d. by physical volume, order-only observation, `N = lambda^-alpha`.
The question ([next_question.md](next_question.md)): is
`alpha_low = alpha_high = 2`?

**Proved, as audited** ([AUDIT.md](AUDIT.md)): `4/3 <= alpha_low <= alpha_high <= 2`.

- 4/3 is a total-variation bound for an observer handed each event's cell,
  so it hands that observer the phase. Nothing is proved at `alpha = 4/3`
  itself.
- 2 rests on an audited draft (Gram-Fourier, `alpha > 2 + 2/K` for each
  fixed K). Its regularised statistic still has to be written into its
  text. The fully written upper bound is 3.
- All three are asymptotic statements. Their constants put any finite-N
  test out of reach: about `1e24` events for one witness of argument 2, and
  `lambda ~ 1e-40` for the scan of argument 3. The audit's words: this line
  produces proofs, not simulations.
- Scope: one family. Nothing is claimed for unknown q, other waveforms,
  Poisson N or generic spacetimes. The audit is AI-only.

Between 4/3 and 2 the open structural question is how finely the order
itself fixes null margins (U1) and how far it transports position (U2).
The heuristic table of next_question.md Section 3 turns assumptions about
these into exponents. Those exponents are not results.

## 2. What each route showed, and did not

**Route A, U1 on the background** ([route_a/RECORD.md](route_a/RECORD.md)).

- **Shown.** Central events, N = 256 to 2048, along the order-preserving
  chain:
  - flip candidates grow as `N^0.515` [0.499, 0.532];
  - the participation of wall hits stays at 9.4 to 9.6, slope 0.014
    [0.001, 0.027], so at most about `N^0.03`;
  - `D = -0.244` [-0.257, -0.231];
  - all gates passed, reproduced by an independent check.
  - Frozen verdict: "no growth seen by this instrument".
- **Not shown.** "Cell-limited": no control in which the right answer grows
  with N. The first two designs did not reach a main run. The first one's
  pilot ran: its LP search covered only 0.64 of its target. The second's
  review found five items, two of them defects in the gate.

**Route D, how far positions move given the order**
([route_d/RECORD.md](route_d/RECORD.md)).

- **Shown.**
  - The residual of central events after the seven-parameter family is
    fitted out stays at `12.4 / N` (N = 128 to 256, lags to `2 N^2`) and
    `12.5 / N` (N = 512, lags to `N^2/2`).
  - No size dependence. A count-scale term `c / sqrt(N)` above about 0.26
    (full residual) or 0.16 (wave coordinate) would have been seen at
    those lags. Without any fit, `c <= 0.55` at N = 512.
  - The review of the transport draft showed one consequence: inside the
    central set, a transport of the heuristic's size (variance of a
    `u`-difference equal to `|du| / N`) would leave a wave coefficient near
    0.29, against an upper end of 0.11 to 0.13. So it is not seen.
- **Not shown.**
  - That positions are pinned at the cell scale.
  - Anything beyond the central set, or beyond `N^2/2` at N = 512.
  - That the wave coordinate is flat in lag: it read "neither" at every
    size.
- **Three further drafts were reviewed and not frozen:**
  - a joint-slack LP: a maximum over a single event is not a typical
    spread;
  - a far-start test: it resolves only about a factor 3 at affordable
    sizes;
  - a transport draft: cell-scale correlation and an outward-drifting event
    set would read as transport.

**Route C, an explicit even statistic** ([route_c/RECORD.md](route_c/RECORD.md)).

- **Shown, analytically.**
  - The wave gives each event a spin-2 cone deformation of amplitude
    `2 sqrt(q) lambda^2 f(u)`, whose square no symmetry removes.
  - But an event's own relations absorb it exactly. The local sign
    statistic of next_question.md therefore has no response that grows
    with N: stop rule (i) fired, conditionally.
  - The two "per event only" rows of the heuristic table lose their
    premise. One of them is one of the three accountings behind the
    "coincidence at 2".
  - What the wave changes is the relations between partners in different
    directions (Weyl curvature).
- **Not shown.**
  - The zero was argued, not computed.
  - The joint per-event fit has no likelihood: free shears break
    transitivity.
  - Its resolution is a U2 question.

**Route B (raise the lower bound): not attempted.** D1 (a finite-n check of
the Gram-Fourier noise law): not run; optional.

## 3. The caveat on every numerical reading

- **The chain's reach.** Every reading of Routes A and D comes from a
  chain that moves one event at a time, started at the true configuration.
  - It is stationary for the conditional law given the order. But it can
    under-report: in 1+1 it showed 0.57, 0.76, 0.87 and 1.00 of the exact
    collective motion at lags `N^2/4`, `N^2/2`, `N^2` and `2 N^2`.
  - The same curve is also the record's one positive control. It shows the
    chain reaching the full collective motion by `2 N^2`, a lag the 3+1
    runs reached at N <= 256. How fast collective modes relax in 3+1 is
    not known.
  - Collective long-range modes are the slowest, and they are what U2 is
    about.
  - So every "not seen" is one-sided. A plateau is not evidence that
    positions are pinned, and nothing in Routes A or D can trigger the
    "toward 4/3" branch.
- **Resolution at affordable chain sizes.** At N = 128 to 512 the count
  scale (about `3 / sqrt(N)`) is only 3 to 6 times the cell scale. A chain
  costs about `N^3` per `N^2` sweeps; N = 1024 at `2 N^2` sweeps would take
  about 31 hours per order.
- **Common to everything else:**
  - generator coordinates, a diagnostic and not an observer;
  - one background (q = 1) and one diamond;
  - AI-level review throughout;
  - finite N, whereas the question is asymptotic.

## 4. The levers of next_question.md: pulled and left

| lever | state |
|---|---|
| A1, analytic count of realisable flips | not done |
| A2, numerical flip predictability | done: "no growth seen by this instrument" |
| B1, second-moment bound with symmetry averaging | not attempted |
| B2, sequential (cavity) replacement | not attempted |
| C, first step: analytic coefficient | done; stopped by rule (i) for the local form |
| C, frozen quadrature | superseded: check 3 and the cone-constraint record cover it |
| D, posterior spread of `u` differences | central set only, one-sided; transport draft not frozen |
| D1, Gram-Fourier noise law at finite n | not run, optional |
| Route C Section 6 finite-N test (referee's sketch) | not drafted |
| writing argument 3's regularised statistic (AUDIT.md Section 3) | not done |

## 5. Lever check

The aim is a reading of U1/U2 that is two-sided, or else a decision that no
lever is worth its cost. Each candidate is judged by what it can show, what
is left if it fails, its cost, and an honest estimate of its chance.

### L1. Chain-free brackets on the cell's extent along smooth collective modes

- **What it is.**
  - Take a sampled configuration and its order.
  - Pin the seven-parameter family by its exact linear stationarity
    condition (as in the slack review).
  - Restrict to a box of a few cell widths per coordinate.
  - For a fixed set of smooth displacement modes (low-order fields over the
    diamond, orthogonal to the family), solve two linear programmes per
    mode and direction. Only link and critical-pair rows are needed,
    because they determine the order.
- **Upper bound (rigorous inside the box).**
  - Link rows: tangent half-spaces. They relax, because the margin is
    concave.
  - Critical-pair rows: tangent half-spaces loosened by the largest
    curvature term over the box, about `|delta y|^2 / Delta m`. That makes
    them a valid outer approximation.
- **Lower bound (constructive).**
  - Link cones replaced by inscribed polyhedral cones, and critical rows by
    plain tangents. This is a convex region inside the true cell.
  - Every solution is an order-preserving configuration. In 1+1, where the
    cones are polyhedral and such a model is exact, the slack review found
    it certified 0.99 to 1.00 of the exact single-event slack. That says
    nothing about 3+1 curvature.
- **What it can show (corrected after the rebuttal).**
  - **Upper side.** A valid, rigorous bound inside the box: nothing in the
    box moves further along the mode with the order kept.
  - **Lower side: it cannot show softness.**
    - The cell is curved. Straight steps along exactly free generators
      keep the order only to 1.0 to 5.7 / N, so a convex inner region
      cannot follow a curved soft direction far.
    - Along a smooth mode, maximum and typical spread need not agree. The
      ball and cube are the favourable cases. The exactly known 1+1 cell
      (simplices) is the counterexample in shape.
  - **What it separates.** "Bounded within the box" from "not bounded".
    The box sits below the count scale, so this is not "rigid" against
    "soft at the count scale".
- **Why it can reach large N, and how far that helps.**
  - It needs no chain. Rows grow as about `2.6 sqrt(N)` per event: 76k at
    N = 1024. Scaling as `N^1.5` gives about 215k at 2048 and 610k at 4096.
  - The 9 to 12 s solves at N = 1024 were timed on 3 to 4 programmes with a
    random or one-event objective.
  - The rebuttal re-estimated the curvature loosening for central walls,
    whose pair separations are about 0.35 to 0.5, with the full Hessian of
    the light-cone margin. It comes to about 1 to 2.3 walls at N = 2048 and
    0.5 to 1.2 at 4096.
  - Its net separation in a rigid world is therefore a factor of 2 to 4,
    the same factor that sank the far-start draft. Enlarging the box does
    not help.
- **If it fails.** The bracket is too wide and it stops; nothing else is
  lost. Even a success says nothing about components of the cell far from
  the truth, freedom reached through large non-smooth rearrangements, or
  modes outside the chosen set. It removes the chain's slowness, not its
  locality.
- **Value even on success.**
  - At N = 4096 a five-cell bound (0.015) is about the wavelength at
    `alpha = 2`. It therefore bears on phase coherence for no `alpha`
    below 2.
  - Its steering would only confirm a prior that already points the
    cell-like way.
- **Cost.**
  - One to two sessions to build, then a review.
  - At N = 4096, interior-point solves of about 5 to 10 minutes put a
    configuration at 6 to 13 hours on one worker. Memory may limit the
    worker count.
  - No GPU. No new dependency: HiGHS is in SciPy 1.18.
- **Chance.** Low as a two-sided lever. As a one-sided rigidity check it
  is feasible, but its result would add little to what the record already
  says.

### L2. A collective-move sampler, validated where the answer is exact

- **What it is.** Add proposals along fixed smooth modes, with an exact
  order check and a density ratio, to the existing chain. Validate in 1+1,
  where the conditional law is given exactly by order statistics.
- **What it can show.**
  - Constructive in one direction: if it finds count-scale motion along a
    smooth mode, soft modes exist.
  - If it finds none, the reading is still one-sided, though the loophole
    of slow single-event mixing is narrowed.
- **If it fails:** validation in 1+1 shows its own limits; nothing new
  about 3+1.
- **Cost.** About one session to build and review. Hours of CPU. Limited to
  N <= 512 by the chain's cost, where the count scale is only 3 to 6
  cells.
- **Chance.** High as an instrument; low as a decider.

### L3. Route B proofs (second moment with symmetry averaging, or cavity)

- **What it can show:** a theorem raising 4/3.
- **Assessment.** Route C showed the information sits in the partners'
  mutual relations, the part the cone-constraint residual measures. A
  second-moment or cavity bound must control exactly the posterior
  localisation that U1/U2 ask about. The prior-art search found no worked
  example with two geometric hypotheses.
- **If it fails:** its own stop rule records a method barrier after two
  sessions.
- **Cost:** analytic sessions plus heavy review. No CPU.
- **Chance:** low without a known U1/U2.

### L4. A rigidity theorem for the order cell

- **What it can show:** cell-scale pinning modulo the family, with a rate.
- **Assessment.** It would need a rate better than `N^-1/2`. Ordinal
  embedding and rigidity results were outside this repository's prior-art
  search and are not assessed here. No route to such a rate is known in
  this programme.
- **Chance:** very low.

### L5. Route C Section 6 finite-N test

It would turn the conditional stop into a computed one. It is cheap, but it
cannot steer U1/U2 or move the bracket.

### L6. Stop the U1/U2 numerics

The record is coherent as it stands:
- the bracket holds;
- the steering evidence is one-sided and points the cell-like way;
- the per-event channel is closed;
- the open question is named.

Stopping costs nothing more. On the rebuttal's assessment it loses no
lever that could make the steering two-sided at reachable sizes.

### L7. Write the regularised statistic into argument 3 (not a U1/U2 lever)

- **What it is.** The upper bound 2 rests on a draft whose regularised
  statistic is described but not written down. AUDIT.md: "Existence is
  therefore not in doubt; what is missing is the text". Its Section 4
  supplies one explicit choice, which no referee has checked.
  - Both referees of argument 3 named this interface as the least-trusted
    step: the cutoff must act on all 18 raw moments, so that the noise
    bound, the wave-side bound and the signal bound hold for one object.
- **What it can show.** The fully written upper bound moves from 3 to 2,
  so the bracket rests on written proofs at both ends.
- **If it fails.** A referee finds that no single regularised object
  carries all three bounds. That would reopen the upper 2, and it is worth
  knowing.
- **Cost.** One analytic session to write it into the draft's text,
  following AUDIT.md Section 4, then an independent review. No CPU.
- **Chance.** High. Existence is not in doubt, and the work is to write and
  check it.

### Recommendation: L7, then stop the U1/U2 numerics (L6)

This replaces the first recommendation (L1); see Section 6.

- **First.** Write the regularised statistic into argument 3 and have it
  reviewed in a fresh context. It is needed whatever U1 and U2 turn out to
  be. It completes the one part of the proved bracket that is still only
  described.
- **Then.** Stop the numerical work on U1/U2.
  - The readings stay as recorded: one-sided, pointing the cell-like way.
  - No lever found here makes them two-sided at a cost worth paying:
    - L1 can only return "bounded in the box" or "stop", by a factor of 2
      to 4;
    - L2 stays one-sided and small-N;
    - L3 and L4 need the localisation they would assume.
- **Optional, cheap.** L5 (Route C's Section 6 finite-N test), to turn the
  conditional stop into a computed one.
- **If the user's goal is U2 steering** rather than finishing the bracket,
  a reduced L1 is the best steering lever, still one-sided. It would take
  the rebuttal's changes: upper bound only, a threshold in cell units, a
  calibration against the chain at N = 256, a released-generator control,
  and a feasibility check at N = 2048 and 4096.

Nothing here was run, drafted, frozen or approved.

## 6. Rebuttal

[lever_rebuttal.md](lever_rebuttal.md): a fresh-context AI desk review,
told to run nothing and to argue against L1. Its verdict was to replace L1.

- L1's inner bound cannot follow curved soft directions, so it cannot show
  softness.
- Its rigid-side separation at N = 2048 to 4096 is a factor of 2 to 4 once
  the loosening is taken over central walls with the full Hessian.
- A five-cell box at N = 4096 is about the wavelength at `alpha = 2`, so a
  rigid result bears on phase coherence for no `alpha` below 2.
- "Either outcome steers" overclaimed, and the claim about a re-embedding
  coupling was wrong: a coupling needs probability mass and freedom at
  amplitude `lambda^2`, below the cell scale.
- It proposed L7 first and L6 after, which this document adopts.

**Corrections it made to the status sections, now applied:**
- the 1+1 curve is also the positive control;
- the first Route A design's pilot did run;
- the second design's review found five items, not four defects;
- the 1+1 inner-model figure came from polyhedral cones;
- the LP timings came from a random objective.

**One claim of the rebuttal is not right.** The N = 512 extension ran
`N^2 + 1` sweeps per order in 6966 s, so `2 N^2` sweeps at N = 1024 is about
31 hours per order (Section 3), not 62.
