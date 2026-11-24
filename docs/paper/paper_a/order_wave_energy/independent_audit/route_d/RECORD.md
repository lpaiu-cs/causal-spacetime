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

## 3. What the review points to instead. Not drafted

The doubt is whether the chain's plateau reflects the whole set of
configurations with the given order or a part the chain is confined to.
The direct test of that is the one the review names as a discriminator:
build a configuration with the same order that is far, in residual, from
the sampled one (by joint moves, certified by the exact order), start the
chain there, and see whether the residual relative to the sampled
configuration comes back to the plateau value. A known-answer control
exists in 1+1 dimensions, where a monotone map of a null coordinate gives
an atypical configuration with the same order and the equilibrium is exact.

This is a different measurement from the one that was drafted. It has not
been designed, and nothing is to be frozen or run without the user's
approval.
