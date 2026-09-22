# Frozen nested-landmark comparison at a fixed event budget

2026-09-22. Commit this protocol and code before the new orders. The
[curvature-square pilot](square_pilot_results.md) left roughly 82--85% of
nominal tidal variance in the magnetic components. This development experiment
tests whether retaining additional relations from the SAME sampled events
reduces that variance and actual reconstruction error. There is no accuracy
or coverage certification gate and no result is yet implied by this design.

## Single fixed design change

The baseline is the connected selector with a maximum of 64 interior landmarks.
The candidate is the UNION of that entire baseline and the same connected
selector run with maximum 128. Preserve all baseline IDs in their original
order before appending new interior IDs; keep the outer tips at the endpoints.
The union has at most **192 interior landmarks**, not necessarily 128. Record
its actual size. Both proposals must pass their existing connectivity checks;
do not replace failed proposals or discard inconvenient landmarks.

Every old pair clock remains available, and the first two spatial frame
vectors are identical. The new landmark coordinates remain unknown fit
parameters. More correlated clocks and more nuisance coordinates do not
guarantee an observed improvement; measure the full marginalized covariance
and actual tensor error. Preserve all covariance cutoffs and optimizer rules.

Both arms use the same fifth-degree vacuum jet plus the fixed curvature-square
distance/volume correction. All ten Weyl and all 24 gradient coefficients
remain free, with no wave profile, polarization or algebraic-type constraint.
The six-field archive records only causal relations, selection IDs/design and
event count. The design name selects an algorithm, not generator information.
The new selector must replay its own IDs from the selection order alone.

## Paired data and unchanged error requirements

Use the existing exact rotating-wave generator with the fixed ordered cases

```text
(A,w) = (0,0), (0,.5), (.25,0), (.75,0),
        (.25,.5), (.75,.5), (.25,-.5), (.75,-.5).
```

Case i uses **2026100600+i**, with selection/measurement streams `[seed,0]`
and `[seed,1]`. Each pair shares 4096 selection events and the SAME 262144
measurement events, generated in chunks of 4096. Query their relations to the
union of landmarks once, then extract each arm's six-field archive. Do not
give either estimator coordinates, truth, profiles or a paired flat sample.
The fixed budget is **2,129,920 sampled events plus 16 supplied tips**.
The larger design requires more relation queries/storage/computation, which
must be reported rather than called free. No event, seed or case replacement.

Attempt both fits independently and retain generation, selection, zero-clock,
rank, optimizer and boundary failures in the assigned arm. After each pure fit
returns, its evaluator may read separate truth. Require the paired true
internal-frame coefficients to agree below 1e-10. Report full ten-component
error modulo reflection, electric and magnetic variance, largest variance,
nominal confidence radius, ranks and actual retained landmark/clock counts.

For each arm also fit continuum volume clocks on its OWN landmarks with its
selected mask and empirical covariance. Require 32/64 quadrature refinement
below 1e-10 and preserve the relative 2% bias condition (flat-zero tolerance
1e-10). Do not interpret lower sampling variance as success if bias or failure
rates worsen. Report every paired ratio, with no selected successful subset.

## Validation and decision boundary

Before data, the deterministic exact-flat check verifies baseline inclusion,
frame preservation, unique IDs, the 192-interior cap, full Weyl rank and zero
curvature under the corrected model. Replay a prior order through the unchanged
connected estimator and regenerate its baseline archive with the paired
generator, comparing all six fields exactly. These use already existing data.

The pre-data plumbing replay used old square-pilot case 5 (seed 2026100505).
All six regenerated connected fields were exactly equal to its existing
archive; nested membership bits on common IDs and the frame prefix were also
exactly equal. Its nested selector retained 120 interior landmarks. This
replay checked data routing, not the nested method's scientific performance.

Record a reachable clean generating SHA, all relevant source hashes and every
output/failure. After completion verify shared-event membership bits on common
IDs, source/artifact hashes and full-rank diagnostics; independently replay one
new nested archive through the pure estimator without truth/generator imports.

This is an eight-case prospective development comparison, not a powered
population claim. Adoption requires genuine uncertainty and finite-error
improvement without losing bias control, followed by fresh coverage and
power-sized accuracy validation. Neither rank nor additional observations alone
achieves the research objective. The finite-resolution obstruction and the
regular-convex-vacuum scope are unchanged.

```powershell
rtk proxy E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_design_pilot.py --design nested --model-order 5 --square-correction --data-seed 2026100600 --protocol nested_pilot_protocol.md --output E:/lab/causal-spacetime/outputs/intrinsic_weyl_nested_pilot_v1
```
