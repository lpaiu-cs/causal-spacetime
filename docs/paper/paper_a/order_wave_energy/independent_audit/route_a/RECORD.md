# Route A: progress record

2026-11-24. **No verdict on U1.** The second design is frozen and its
control passed, but an independent review found two defects in its decision
gate, so its 3+1 campaign was not run. A third design answers the review
and is frozen after the user's approval (Sections 5 to 7); its results, when
they exist, are recorded in results_v3.md. AI-run work throughout; no human
review.

## 1. Why this route first

Route A of [next_question.md](../next_question.md) needs no wave, runs on a
CPU, and decides which side of the bracket to attack: whether to try to
raise the necessary exponent toward 2 (Route B) or to look for an order
statistic below 2 (Routes C and D). Completing the written proof of the
upper bound 2 is left until after it, because Route A decides whether that
bound is half of a sharp law or an intermediate bound.

## 2. Stamps

| commit | what |
|---|---|
| `741ddfc` | first design, instrument, runner, analysis frozen ([design.md](design.md)) |
| `228cc84` | first design, 1+1 control stage recorded (ran at `741ddfc`) |
| `510bc46` | pilot recorded (ran at `228cc84`); second design frozen ([design_v2.md](design_v2.md)) |
| `c1a6474` | second design, 1+1 control stage recorded (ran at `510bc46`) |

## 3. First design (linear-programming flip search)

- Control, 1+1 dimensions: 750 of 750 known realisable flips found, none
  certified outside the rule. Gate passed.
- Pilot, N = 96, 3+1 background: completeness 1398 of 2187, `0.639`,
  Wilson 95% `[0.619, 0.659]`, floor `0.90`. Joint searches ran out of
  their programme budget while still making progress. Cost 440 s per order.
- Main stage: **not run.** As frozen it would have returned "unresolved:
  instrument gate failed" at a cost above 30 CPU hours for N = 512 alone.

## 4. What the pilot showed, and what it does not

Per tested event at N = 96: 48.8 candidates, 15.3 walls at the sampled
configuration, 35.4 flips certified by some search. Of the walls at one
chain state, about 0.6 were still walls 250 to 2000 sweeps later. The chain
moved events by 0.32 to 0.38 (root mean square) in a diamond of proper
height 1.48, with the order exact at every checkpoint.

So at N = 96 most transitivity-safe flips are realisable and the embedding
is loose at the scale of the discreteness length. **This does not
discriminate the two hypotheses.** The link scale is about
`(2 to 5) / sqrt(N)` and the counting scale of the rigidity argument is
about `16 / N`; their ratio is of order one at N = 96 and reaches 4 to 10
only near N = 1000 to 2000 (design_v2.md Section 2). The pilot was declared
for sizing and is not read as a verdict.

The first design's analysis (design.md Section 2.4) predicted cell-limited
behaviour from a counting argument. That prediction is untested. Its
unstated constant is what put the first grid in the wrong range.

## 5. Second design (wall participation from chain chord ends)

- Instrument: the order-preserving chain itself. The partner whose cone
  bounds each line move is a wall found by exact geometry. Participation is
  estimated from products of partner frequencies in different blocks.
- Control, 1+1 dimensions (known value 4): `3.921`, 95% interval
  `[3.895, 3.942]`, tolerance `0.3`; no hit outside the candidates. Gate
  passed. The small deficit below 4 is explained in Section 6.
- Campaign, 3+1 background, N = 128 to 2048, 12 orders each: **not run**,
  and after the review of Section 6 not to be run as frozen.

**Authorisation.** The freezes and runs above were made after the user's
go-ahead to start Route A in this session. A condition relayed afterwards
by the coordinating session says designs are not to be frozen or run
without the user's approval. The stage `main_a` had been started at
`c1a6474` (2026-11-24T08:18Z). It was stopped before any order finished,
its output folder held only an in-progress manifest and was removed, and
nothing from it was seen. The campaign waits for approval. Until it runs,
design_v2.md can still be amended without having seen 3+1 data from its
instrument; the pilot numbers of Section 4 were seen and are disclosed
there.

Cost if approved: CPU only, about 20 minutes for N up to 1024 and about 85
minutes for N = 2048 on 12 workers.

## 6. Independent review of the second design

[referee_report.md](referee_report.md). An AI referee in a fresh context,
given the code and both designs and no 3+1 result.

Found correct: the row lemma and its implementation; the chain's stationary
law (checked against exact order-statistic laws in 1+1 and against the
density in 3+1); the rule that the partner at a segment end is a wall, with
the exact orientation, at 72,000 segment ends; unbiasedness of the
cross-block estimator (independent exact draws give 4.0000 in 1+1); the
constants of design_v2.md Section 2.

Found wrong or missing, in the referee's order of importance:

1. Extremal events are 48% of the population at N = 256 and 22% at 2048 and
   have more than twice the candidates. Pooling them biases the slope by an
   amount of unknown sign.
2. design_v2.md says the count-limited reading does not need the mixing
   gate. For a slope that is false, and the mixing ratio it defines cannot
   see slow relaxation.
3. There is no control in which the right answer grows with N.
4. The control's "known value 4" holds only for events with four distinct
   rank neighbours; with a shared neighbour it is 8/3. That, not slow modes,
   is 84% of the control's deficit (3.921 against 4): events with four
   distinct neighbours gave 3.987. design_v2.md Section 5 attributes the
   deficit to slow collective modes; that attribution is withdrawn here.
5. A rounding case can drop a root in 1+1. No practical effect in 3+1.

The design value 0.40 for the candidate slope was a single-seed reading;
the referee's eight configurations per size give 0.48 to 0.50.

The referee reported that one of its commands briefly created two untracked
files in the worktree root and removed them. The tree is clean; no stamped
run overlapped that interval's preflight.

## 7. Where this stands

- design_v2.md is left as frozen and its campaign is **not run**: findings
  1 and 2 change what its gate would mean.
- [design_v3.md](design_v3.md) answers the findings: a fixed central
  stratum for the gate, mixing required for both verdicts with longer
  chains to check it, corrected control values, the rounding guard. The
  user approved it with two conditions (its Section 5) and it is frozen
  with the commit that adds it under that name. Cost about 25 CPU hours.
- Finding 3 is not repaired, and the draft says so.
- Nothing has been done on Routes B, C, D or on writing the regularised
  statistic of the upper-bound draft.
