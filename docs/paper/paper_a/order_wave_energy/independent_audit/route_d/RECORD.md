# Route D: progress record

2026-11-25. AI-run work throughout; no human review. No verdict on U2
(how well the order fixes positions) has been defined or recorded.

## 1. First measurement: done

[design.md](design.md), [results.md](results.md),
[checker_report.md](checker_report.md). Frozen at `4bdead6` after the
user's approval, run at `265d73a`.

Along an order-preserving chain at N = 128, 192, 256, the displacement of
central events keeps growing with the lag, and after the seven-parameter
order-preserving conformal group is fitted out the residual stays near
`12.4 / N` from `16 N` to `2 N^2` sweeps. In the design's words: no
residual growth seen up to `2 N^2` sweeps by this instrument. Not shown:
that positions are fixed at the cell scale. The wave coordinate alone is
not shown flat. The largest doubt left by the independent check is
non-ergodicity of the single-event chain.

## 2. Second measurement as first drafted: withdrawn

[slack_design_draft.md](slack_design_draft.md),
[slack_referee_report.md](slack_referee_report.md). The draft proposed a
linear programme in which all events move at once, with the family pinned,
to see how far a central event can be displaced, with the answer certified
against the exact order. It was reviewed before any freeze and is
withdrawn. Nothing in it was implemented or run. The review's reasons:

- **A maximum is not a typical spread, and in high dimension the gap is
  large.** A ball-like cell in `4N` dimensions whose typical per-event
  spread is the measured `12.4 / N` has an extent of about `8.8 / sqrt(N)`
  in any one direction. So "count-scale joint slack certified" is what a
  typically rigid cell would return, and the draft's stated consequence
  (the chain's plateau could not be read as rigidity) overclaimed.
- **The positive control cannot pass.** Releasing one family pin frees a
  motion that the final seven-parameter fit then removes, so the control's
  ratio is about 1 by construction.
- **The certification under-reports.** Against the exactly known slack in
  1+1 dimensions the drafted loop certified 19 to 65% of it, so "model
  slack not certified" would be returned whether or not slack exists.
- The default solver is too slow at N of 512 and above; an interior-point
  method is needed.

Kept from the review for later use: the pin is the exact stationarity
condition of the group fit; a links-only programme is a rigorous upper
bound on the pinned slack inside the box; exact link cones with tangent
critical rows give certificates without a step search but need a
second-order cone solver, which is not installed.

## 3. Far-start test: drafted, reviewed, not to be frozen as written

[return_design_draft.md](return_design_draft.md),
[return_referee_report.md](return_referee_report.md). The idea: push the
configuration along a collective non-family direction with its order kept
exactly, start the chain there, and see whether the push is forgotten.
Nothing in it was implemented or run by the orchestrator; the referee
built pushed starts chain-free and ran the two controls and a 1+1 version.

What the review found:

- **The path as written does not get far.** None of 36 paths reached twice
  the plateau distance (means 0.7 to 0.9), so the design could only return
  "not far enough". A repaired step rule (floors relative to the starting
  margin, a per-link trust region) reached 2 to 3 plateaus, of which only
  1.5 to 2.0 lie along the pushed field, and threw single outer events far.
  Whether such starts are sound was not tested.
- **The controls do not exercise the main statistic.** They were scored on
  fitted group parameters, not on the memory statistic. The rotation
  control could not even be built at N = 128 (the required angle exceeds a
  half turn). Scored through the main pipeline the free rotation is
  remembered and the stiff direction is forgotten, as they should be.
- **The descriptions lack power.** The ratio of two equilibrated distances
  has a per-order spread of 0.14 to 0.19, so with twelve orders "neither"
  is about as likely as "forgotten" under a perfect return. Twenty-four to
  forty-eight orders are needed.
- **"Forgotten" means narrow relative to the push, not cell scale.** In
  1+1 dimensions, where the conditional law is broad at the count scale,
  the pipeline still returns "forgotten" for a push of 2.8 plateaus.

Two things follow that do not depend on running it.

- **Resolution at these sizes is low.** At N = 128 to 256 the count scale
  (about `3 / sqrt(N)`) is only 3 to 4 plateau distances. A push of 2
  plateaus sits inside a count-scale spread, so the test could at best
  tell a width of one plateau from a width of three or more.
- **The reach of the repaired path is what an isotropic cell would give.**
  For a ball-like cell whose typical spread is the measured plateau, the
  extent along a collective direction carried by the central events is
  about `sqrt(N / n)` typical per-event spreads, n being the number of
  central events: roughly 2 plateaus. The path saturated at 1.5 to 2.0.
  This is an observation from the referee's chain-free runs, not a
  measurement of this programme, and by the draft's own rule it is not
  evidence.

## 4. Where Route D stands

- The first measurement stands as recorded: a plateau of the central
  residual at `12.4 / N`, no size dependence found, with the limits of
  [results.md](results.md) Section 4.
- Two attempts to reduce the non-ergodicity doubt numerically did not
  yield a design that can be frozen: one cannot discriminate in principle
  (a maximum is not a typical spread), the other needs a rebuilt path,
  rebuilt controls and 24 to 48 orders, and would still resolve only a
  factor of about three at the affordable sizes.
- A repaired far-start design is possible and its requirements are listed
  in the referee's report. It has not been written.
- The simpler way to strengthen the first measurement was to extend its
  frozen estimator to a larger size. That was done with the user's
  approval: [extension_design.md](extension_design.md) (frozen `e302056`),
  [extension_results.md](extension_results.md), checked independently in
  [extension_checker_report.md](extension_checker_report.md). At N = 512
  the central residual is about `12.5 / N`, as at the smaller sizes. No
  count-scale growth was found over the four sizes. A term above about
  `0.26 / sqrt(N)` (full residual) or `0.16 / sqrt(N)` (wave coordinate)
  would have been seen at these lags. Without any fit, c <= 0.55 at
  N = 512. The longest lag at N = 512 is `N^2/2`, and chain reach is still
  untested.

## 5. U2 transport measurement: drafted, reviewed, not to be frozen as written

[transport_design_draft.md](transport_design_draft.md),
[transport_referee_report.md](transport_referee_report.md). The user chose
it as the next draft after Route C's first step showed that both routes
turn on U2 ([../route_c/RECORD.md](../route_c/RECORD.md)). The idea: the
same chain and group fit, fitted over all interior events, with the
correlation of residuals of event pairs read against their separation, up
to the diamond's scale. Nothing was implemented or run.

What the review found:

- **The drafted descriptions would mislabel.** Cell-scale correlation
  between nearby events reads as "growing with d". The outward drift of
  the interior set with N reads as "count-scale transport found". Transport
  shows only in how the near and mid-range covariance changes with N.
- **C2's gate cannot fail when it should, and C3 passes almost surely.**
- **The existing record already answers part of the question.** Inside the
  central set, a transport of the heuristic's size (variance of a
  `u`-difference equal to `|du| / N`) would leave a wave-coordinate
  count-scale coefficient of about 0.29. The extension bounds that
  coefficient at 0.11 to 0.13 at lags up to `N^2/2` (N = 512) and `2 N^2`
  (N <= 256), and could tell 0.16 from zero. Confirmed by the referee and
  rechecked here without a chain.

So, by this instrument and one-sided: **no transport of the heuristic's
size is seen inside the central set.** This does not show that positions
are pinned against each other at the cell scale, and it does not trigger
Route D's "toward 4/3" branch, which needs a two-sided statement.

What a revised measurement could still add:

- separations beyond the central set, where loose outer cells and the
  composition of the event set confound;
- a near-null stratum of pairs for Route C;
- lag `2 N^2` at N = 512.

On the referee's projection, its sensitivity for the full residual would be
no better than the central bound. The six must-fix items are listed in the
report. No revision has been written.

Nothing is to be frozen or run without the user's approval.

## 6. Numerics stopped (2026-11-26)

After the lever check
([../STATUS_AND_LEVERS.md](../STATUS_AND_LEVERS.md), Sections 5 to 7), the
user approved stopping the numerical work on U1 and U2. For this route:
- **The readings stand exactly as recorded in Sections 1, 4 and 5:**
  - the central residual is about `12.5 / N`;
  - no count-scale growth was found, at about 0.26 (full) and 0.16 (wave)
    sensitivity, and `c <= 0.55` without a fit;
  - no transport of the heuristic's size is seen inside the central set.
- All of these are one-sided, from a single-event chain whose reach in 3+1
  is untested. They show neither that positions are pinned at the cell
  scale nor that they are not.
- Nothing was re-read or re-scored on stopping.
- The far-start repair and a transport revision are not written. Neither
  is to be revived without the user's approval.
