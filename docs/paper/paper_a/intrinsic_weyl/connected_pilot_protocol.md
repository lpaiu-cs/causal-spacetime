# Prospective connected-three-band development pilot

2026-09-22. Freeze before generating seeds `2026092240..2026092245`.
This is a selection/conditioning pilot, not an accuracy or coverage gate.
The two preceding six-case pilots and their failures remain unchanged.

## The specific defect

The exact two-slice obstruction absorbs all five magnetic first-order
responses into unknown positions. Three distinct time slices can remove it,
but the previous selector retained only 0--7 central landmarks. A subsequent
order-only diagnostic on all twelve previous selection samples found ZERO
central candidates with at least four count-qualified clocks to EACH of the
21 past and 21 future landmarks chosen using that three-band score. Merely
reranking the central candidates cannot meet this connectivity condition on
those selected outer populations. This was exploratory diagnosis, not new
independent validation data; no measurement spectra entered that diagnostic.

## Frozen change

Retain the blind quartic joint fit, all ten unconstrained Weyl components,
covariance/optimizer settings, physical-volume sampler, causal predicate,
amplitude populations, and measurement count from `joint_pilot_protocol.md`.
Change the selection design as follows:

1. Use 4096 independent selection events, with 262144 measurement events.
   Seeds are amplitude-major/repetition-minor over A=(0,0.25,0.75), two
   repetitions each. RNG streams remain `[seed,0]` and `[seed,1]`.
   Total events per case are 266240 plus the supplied outer tips. This is
   not an equal-cost comparison to the preceding pilots.
2. Obtain t and r-squared proxies exclusively from selection-order endpoint
   counts, using the same formulas as before. Set outer bands to
   `0.24 <= |t| <= 0.36`, central band to `|t| <= 0.075`, radial range to
   `0.075^2 <= r^2 <= 0.2^2`, and endpoint minimum counts to two.
   The selection-time floor associated with expected count two falls from
   about 0.2102 to 0.1487 in flat units; this is an expectation calculation,
   not a guarantee about noisy selected intervals.
3. Score all eligible candidates by count-qualified degree times
   `min(r^2,0.04)` as before. First retain up to 21 past and 21 future
   candidates by that score. Count-qualified means at least two selection
   events in the interval.
4. A central candidate must have at least four qualified clocks to EACH
   selected outer band. Rank these eligible central candidates by the
   smaller of their two degrees times `min(r^2,0.04)`. Keep up to 22.
   Exact ties retain input order; labels were independently randomized.
5. Apply the previous degree-eight peeling rule to the combined set. After
   peeling, require at least 16 past, 8 central, and 16 future landmarks,
   with every central point still connected by four qualified clocks to
   each retained side. Otherwise record selection failure before generating
   any measurement events. Do not refill, change bands, weaken thresholds,
   or replace seeds. These are deterministic instrument-construction
   requirements, not statistical accuracy claims. They do not prove good
   angular coverage, full rank, or adequate magnetic precision.
6. Archives select this algorithm with `selection_design="connected"`;
   no population or amplitude metadata enters the blind reader. On selection
   failure preserve the selection order in a separate `_selection_failed.npz`
   file, record its hash and failure, and proceed to the next fixed case.
   Successful cases save the same order-query fields used by the prior
   six-field archives. The reader recomputes all selection decisions.
7. Report retained bands and minimum central clock degrees, all previous fit
   diagnostics, both spectra, every failure, and actual run times. No claim
   about precision improvement is decided from two replicates per amplitude.

Before stochastic data, retain the old deterministic flat checks, check flat
recovery on a fully populated connected three-band configuration, and ensure
removing all central-to-future clocks causes selection failure. Replay one
old archive from each previous design without importing the generator and
verify all ten coefficients. Freeze source hashes and a clean reachable
commit. No new measurement population has been generated while designing
these changes.

If central retention is repaired but the magnetic estimates remain noisy,
report those as two separate findings. Selection feasibility does not achieve
the full reconstruction goal. Powered accuracy/coverage validation and
varying-profile tests would still be required for a stronger claim.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_joint_pilot.py --check
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_joint_pilot.py --selection-design connected --output outputs/intrinsic_weyl_joint_pilot_v3
```
