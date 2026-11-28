# Prospective three-band development pilot

2026-09-22. Freeze before the six seeds below are generated. This is a new
development diagnostic, not a replacement verdict for the first pilot or
an accuracy confirmation. Preserve its unfavorable magnetic estimates.

The exact two-slice certificate shows five magnetic directions can be
absorbed by inferred positions at first order. Its three-slice example has
rank ten. The present question is whether a third time band helps in actual
finite orders at the SAME event budget. Rank alone does not answer that.

Keep all mechanics, amplitude populations, sample sizes, covariance cutoff,
clock-selection threshold, optimizer starts/tolerances, reporting, and failure
handling from `joint_pilot_protocol.md`. Change only:

- Use independent seeds `2026092230..2026092235`, amplitude-major and
  repetition-minor over A=(0,0.25,0.75), two repeats per amplitude.
- Add a central time-proxy band `|t| <= 0.075`, retaining the original outer
  bands `0.15 <= |t| <= 0.32`. The proxies still use selection-order counts.
- All bands require the same `r^2 >= 0.01` and endpoint-clock constraints.
  Compute degree/radius scores over their combined candidate population.
  Select at most 21 past, 22 central and 21 future landmarks, then apply
  the same degree-eight peeling rule. Do not replenish or adapt after seeing
  measurement data. Fewer than sixteen landmarks is a recorded failure.
- Archives carry one extra algorithm field, `selection_design="three"`.
  It selects this frozen algorithm, not a generator population or amplitude.
  The estimator rejects any other value. Old five-field archives continue
  to select the unchanged two-band algorithm.
- Retain and inspect BOTH electric and magnetic spectra; do not decide this
  diagnostic from electric success alone. No significance or accuracy gate
  is invented for the comparison of two replicates per population.

The reconstruction model and ten unknown tensor components remain unchanged.
No coordinate, profile, true tensor or paired flat sample may enter the blind
reader. Verify a prior archive still replays unchanged and exercise the new
selector on deterministic flat clocks before data. Run from a clean pre-run
commit; save the new six archives and full result separately.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_joint_pilot.py --selection-design three --output outputs/intrinsic_weyl_joint_pilot_v2
```
