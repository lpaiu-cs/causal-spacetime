# Is exponent 2 sharp? Where the campaign stands, and which lever is left

2026-11-26. Written at the user's request ("정리 후 레버 점검", relayed by the
session that started this one): a short status of the whole campaign, then
an assessment of the means that could make its readings two-sided, with one
recommendation. **No new computation was run for this document.** It reads
the records and estimates costs from them. AI-written. A fresh-context AI
reviewer was asked to argue against the recommendation; that rebuttal is in
Section 6.

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
  with N. The first two designs did not reach a run: an LP search covered
  only 0.64 of its target, and a reviewed design had four defects.

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
  - Every solution is an order-preserving configuration. In 1+1 the
    slack review found that such an inner model certified 0.99 to 1.00 of
    the exact slack.
- **What it can show.** For each configuration and mode, a bracket on how
  far events can move together along that mode near the truth with the
  order unchanged.
  - Along a smooth mode, the per-event extent of a typically rigid cell is
    of the order of the typical spread. A ball gives about one SD; a cube
    about 1.7.
  - That is unlike a single event's extent, where a maximum and a typical
    spread differ by `sqrt(N)` (the slack review).
  - So a bracket separates "rigid at a few cells" from "soft at the count
    scale".
- **Why it can reach large N.**
  - It needs no chain. The rows grow only as about `2.6 sqrt(N)` per
    event: 76k at N = 1024, with interior-point solves of 9 to 12 s.
    Scaling as `N^1.5` gives about 215k at 2048 and 610k at 4096.
  - The curvature loosening relative to the walls falls like `1/N` for a
    box of fixed cell count. On rough estimates it is comparable to the
    walls at N = 1024 and about a third of them at N = 4096.
  - There the count scale is 8 to 16 cells. So the method is most
    informative exactly where chains cannot go.
- **If it fails.**
  - Most likely, the bracket is too wide: the loosening or the inner model
    is weak, and the method cannot separate the two cases at the sizes it
    reaches. That would be recorded, and nothing else is lost.
  - Failure does not touch the earlier readings.
  - A success says nothing about components of the cell far from the
    truth, nor about modes outside the chosen set.
- **Cost.**
  - One to two sessions to build: the rows, the loosening, the pin, the
    modes and both programmes. Then an independent review.
  - CPU: at N = 4096, about 20 modes × 2 directions × 2 bounds per
    configuration at minutes each, so about two hours per configuration
    on one worker. Twelve configurations per size over a few sizes is
    several hours. Memory, at perhaps 1 to 2 GB per interior-point solve,
    may limit the worker count.
  - No GPU. No new dependency: HiGHS ships with SciPy 1.18, which is
    installed. A cone solver would make the inner bound tighter but is not
    required.
- **Chance.** Moderate. The physics of the loosening favours large N. The
  open risks are the inner model's tightness in 3+1, which has not been
  tried, and solver behaviour at 600k rows.

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

### L6. Stop

The record is coherent as it stands:
- the bracket holds;
- the steering evidence is one-sided and points the cell-like way;
- the per-event channel is closed;
- the open question is named.

Stopping costs nothing more and loses the one lever that could make the
steering two-sided.

### Recommendation: L1, behind a small feasibility check

L1 is the only candidate that can give a two-sided reading, and it avoids
the caveat of Section 3 by not using a chain. It also reaches the sizes
where cell and count scales separate by a factor of 8 to 16.

**Proposed feasibility check (not run).**
- One configuration at N = 2048.
- Three smooth modes, a box of five cell widths.
- Both programmes, recording the solve times and the ratio of upper to
  lower bound.
- About ten minutes of CPU.

**Rule, stated now.** If the ratio exceeds 3 for most modes, the bracket
cannot separate the cases at that size: L1 stops and is recorded, and
stopping (L6) becomes the recommendation.

Either outcome of a full L1 would steer the proof work:
- soft modes found constructively would give a re-embedding coupling the
  freedom it needs for a lower bound toward 2;
- a rigid bracket would point at statistics below 2 that use the partners'
  mutual relations.

## 6. Rebuttal

Not yet written.
