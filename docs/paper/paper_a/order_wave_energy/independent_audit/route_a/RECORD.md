# Route A: progress record

2026-11-24. **No verdict on U1.** The second design is frozen and its
control passed. Its 3+1 campaign has not produced a number and is waiting
for the user's approval (Section 5). AI-run work throughout; no human
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
  passed. The small deficit below 4 is under independent review.
- Campaign, 3+1 background, N = 128 to 2048, 12 orders each: **not run.**

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

## 6. Pending

- Independent review of the instrument and estimator by an AI referee that
  was given the code and designs and no 3+1 result. Its report is to be
  added to this folder, and any defect it finds is to be handled by an
  amendment before the campaign runs.
- The campaign and its gates.
- Nothing has been done on Routes B, C, D or on writing the regularised
  statistic of the upper-bound draft.
