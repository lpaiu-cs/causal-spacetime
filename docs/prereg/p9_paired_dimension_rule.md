# P9: dimension selection in 3+1D under a within-seed decision rule

Status: **PREREGISTERED; NOT RUN** (2026-07-25).
Depends on P8 (`p8_3plus1d.md`, PR #19) for the scene and the measurement
script; P9 changes **only the decision rule**.

No number below was derived from any data this rule will be applied to.
There is no calibration stage, because there is nothing to calibrate —
see Section 4, which is the point of the design.

## 1. Why this exists, stated plainly

P8 ran its calibration and **placed no gate**: the `d = 3` and `d = 2`
truth-error populations overlap at their extremes (0.1395 against
0.1328), as do the held-out populations (0.0990 against 0.0950). Its
Section 6 named that outcome in advance as admissible and forbade
repairing the scene afterwards.

But P8's own record shows the instrument was not indifferent to
dimension: `d = 3` beat `d = 2` in **20 of 20** seeds, without exception.
What failed was the *decision procedure* — inherited from P2 — which
places one absolute gate across seeds and therefore requires the two
*populations* to separate, not merely the within-seed comparison to come
out right.

**This rule was chosen after seeing that failure.** That is exactly why
P9 is a new preregistration on entirely fresh seeds rather than an
amendment to P8: a rule selected in light of some data may not then be
evaluated on it. P8's counts appear in Section 7 as motivation and are
**not** an input to any P9 threshold, count, or tolerance.

A possibility worth stating up front, because it changes how a P9 pass
should be read: the absolute-gate procedure may have been the wrong tool
all along and merely got away with it in low dimension. Going from three
spatial dimensions to two discards a third of them; going from two to one
discards half. The underfit is milder in higher dimension, so the
populations crowd. If that is the story, P9 is not a repair to P8 but a
correction to a decision rule the programme has been using since PC-V1,
and it should be reported that way.

## 2. Hypotheses

Let `t_d` be a seed's truth-order error at fit dimension `d`, `h_3` its
held-out violation at the gate dimension, and `c_3` the density-matched
geometry-free control's held-out violation on the same seed.

- **H-KNEE (confirmatory).** The truth error has a knee at the true
  dimension. Two legs, both required:
  - **underfit**: `t_3 < t_2` on at least 16 of 20 valid seeds;
  - **no-improvement**: `t_4 < t_3 - 0.05` on at most 4 of 20.
- **H-SPEC-PAIRED (confirmatory).** `h_3 < c_3` on at least 16 of 20.
- **H-CEILING (descriptive, not gated).** The fraction of seeds with
  `h_3 <= 0.10`, the programme's standing held-out ceiling, is reported
  as a diagnostic. It does not enter any decision.

## 3. Why two legs, and why neither alone

A within-seed test of `t_3 < t_2` **on its own would be close to
vacuous**: `d = 3` has more free parameters than `d = 2`, so it can only
fit the same constraints at least as well, and a fit that is merely
better is no evidence of a knee. The no-improvement leg is what supplies
the other side — if extra parameters were simply buying error reduction,
`d = 4` would beat `d = 3` in the same systematic way, and the rule
requires that it does not.

Together the two legs say: *error falls from 2 to 3 and then stops
falling.* That is the dimension-selection claim. Neither leg is the
claim, and neither is reported as if it were.

**Honest asymmetry.** The underfit leg is a positive finding; a sign test
can support it. The no-improvement leg is an *equivalence* claim, and a
sign test cannot establish equivalence — only fail to find a difference.
It is therefore stated as "no evidence that `d = 4` improves on `d = 3`
by more than the tolerance", and a pass on it is a non-rejection, not
support. Anyone reporting P9 must keep that distinction.

## 4. Why there is no calibration stage

Every constant in Section 2 is either a count or inherited:

| constant | value | provenance |
|---|---|---|
| pass count | 16 of 20 | P2's `pass_min`/`denominator`, frozen since P2 |
| saturation tolerance | 0.05 | P2's decision rule, frozen since P2 |
| seed count | 20 | P2's denominator |

**Nothing is fitted to data.** That removes the failure mode P8 hit
entirely — there is no gate to place, so there is no interval for
clusters to overlap in. It also removes the temptation P8's outcome
created, since a rule with no free parameter cannot be nudged.

The price is that P9 makes a weaker absolute statement than P2 did. P2's
gate says "held-out violation below 0.10", which is interpretable on its
own; P9's paired rule says only "better than the control on this seed".
That is why H-CEILING is reported alongside — it preserves the absolute
reading as a diagnostic without letting it gate anything.

## 5. Null model and what the counts mean

Under the null that the instrument is indifferent to fit dimension, the
sign of `t_3 - t_2` is exchangeable across seeds, so the underfit count is
`Binomial(20, 1/2)`. A threshold of 16 gives a one-sided
`p = 0.0059`. The same applies to H-SPEC-PAIRED.

The two confirmatory hypotheses are read as a **conjunction**: P9 is
supported only if H-KNEE (both legs) and H-SPEC-PAIRED all pass. No
partial credit, and no reporting of one leg without the others.

## 6. Design

- **Scene and measurement**: P8's, unchanged —
  `build_scene_3plus1d`, `R = 12` chains on a Fibonacci shell, 96 ticks,
  `n_events = 7000`, targets 30-44, fits at `d = 2, 3, 4` through the
  frozen PC-V1 pipeline. P8's script already emits every column P9 needs;
  only the decision function is new.
- **Seeds**: **500-519**, confirmatory, used once. Fresh with respect to
  P8's scene selection (0-3) and P8's Stage A (100-119), and with respect
  to every other arm in the programme.
- **No Stage A.** Section 4.
- **Invalid seeds** are recorded with their reason and excluded from the
  denominator, which is not topped up.

## 7. Motivation from P8, which is not evidence for P9

P8's Stage A table (`p8_stage_a_calibration.csv`, seeds 100-119) shows,
under the rule above:

| leg | P8 count |
|---|---|
| `t_3 < t_2` | 20 / 20 |
| `t_4 < t_3 - 0.05` | 0 / 20 |
| `h_3 < c_3` | 20 / 20 |

This is why P9 is worth running and is recorded so the reasoning is
visible. **It is not evidence for P9's hypotheses**, because the rule was
selected with these numbers in view. P9 stands or falls on seeds 500-519.

For the same reason, a P9 pass should be reported as "confirmed on fresh
seeds under a rule motivated by P8", never as "22 of 20 seeds" or any
pooled count across the two.

## 8. Stop rules

The rule is fixed by this document. There is no repair route that
preserves P9: if the counts fall short on seeds 500-519, P9 is reported
as not supported and any further attempt is a new preregistration with
new seeds. In particular the tolerance may not be widened, the pass count
may not be lowered, and the two legs may not be reported separately to
salvage one of them.

## 9. Claim boundary

- P9 concerns the *discriminator's dimension selection*, not the
  identifiability theory. It does not test G4c and does not run at the
  `R >= d + 2` threshold; it runs at 12 observers.
- Truth-order error is scored against coordinates the builder knows. It
  is a validation instrument, not something the pipeline sees.
- A pass licenses "the instrument selects `d = 3` on 3+1D sprinkled
  order, by a within-seed criterion". It does **not** license any
  absolute statement about held-out violation, which is what P8's gate
  would have provided and P9 deliberately does not.
- Nothing here speaks to continuum emergence, to quantum dynamics, or to
  causal sets not built by sprinkling.

## 10. Deviations from the programme's standing procedure

- **E1.** The decision rule is within-seed rather than an absolute gate.
  This is the substance of P9 and the reason it exists; see Sections 1
  and 4.
- **E2.** No calibration stage. Consequent on E1 — there is nothing to
  calibrate.
- **E3.** The rule was selected after seeing P8's outcome. Mitigated by
  entirely fresh seeds and by admitting P8's counts as motivation only.

## 11. Confirmatory outcome (post-run factual record)

Not yet run.
