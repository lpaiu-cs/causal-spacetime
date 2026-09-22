# Three-band pilot: the selection did not retain the intended geometry

2026-09-22. This diagnostic is completed; the full reconstruction goal is not.
The prospective [protocol](three_band_pilot_protocol.md) and code were committed
at `cc924bb` before its six new seeds. All six cases ran, with 1024 selection
events and 262144 measurement events each, without changes or removed seeds.
The complete records are [the result](three_band_pilot_result.json) and
[the order-derived retention audit](three_band_retention_audit.json).

The intended maximum was 21 past, 22 central and 21 future landmarks. After
degree-eight peeling, the actual retained populations were:

| Case | A, evaluation only | Past | Central | Future |
|---|---:|---:|---:|---:|
| 0 | 0 | 21 | 2 | 21 |
| 1 | 0 | 21 | 4 | 21 |
| 2 | 0.25 | 21 | 3 | 21 |
| 3 | 0.25 | 21 | 7 | 20 |
| 4 | 0.75 | 21 | 0 | 20 |
| 5 | 0.75 | 21 | 1 | 21 |

Thus this was not a successful realization of a well-populated three-band
instrument. A nominal selector setting is not evidence of the selected
geometry. Every fit still had full coordinate rank and projected Weyl rank
ten because the remaining time bands have width; rank does not protect
against near-degeneracy or large statistical errors.

## Both observed spectra

The expected electric AND magnetic eigenvalue spectra in these plane-wave
populations are `(-A,0,A)`, in squared-tip-time units. The estimator is not
given that form or A.

| Case | Electric spectrum | Magnetic spectrum |
|---|---|---|
| 0 | (-0.0784, -0.0669, 0.1454) | (-0.6141, -0.0230, 0.6370) |
| 1 | (-0.2839, 0.0683, 0.2155) | (-0.5624, 0.0166, 0.5458) |
| 2 | (-0.3421, -0.0379, 0.3800) | (-1.4293, 0.1930, 1.2362) |
| 3 | (-0.2315, 0.0118, 0.2197) | (-0.9470, 0.3093, 0.6377) |
| 4 | (-0.8795, 0.1400, 0.7395) | (-1.1338, 0.0500, 1.0839) |
| 5 | (-0.8147, 0.1721, 0.6426) | (-0.7512, 0.0134, 0.7377) |

Case 5 has a close magnetic point estimate; the other cases rule out presenting
that one result as a consistently successful reconstruction. No accuracy,
coverage, or variance-improvement verdict is assigned to these twelve pilot
cases across the two designs. They are not paired comparison data.

All optimizers reported success, all inferred times stayed inside the tip
interval, and all data-generation ambiguity counters were zero. The six
generation-and-fit times ranged from about 2.2 to 5.2 seconds. These checks
do not repair the missing central landmarks.

## Input integrity and next action

The new archives, under `outputs/intrinsic_weyl_joint_pilot_v2`, have the
previous five order-only fields plus `selection_design="three"`, a frozen
algorithm choice. All six archive hashes and all generating source hashes
were verified against the reachable pre-run commit. Case 4 was replayed in
a separate estimator process with no generator import; all ten coefficients
matched to `1e-11`. Prior five-field archives also replay unchanged.

The next experiment must fix the observed selection failure, not simply raise
the reconstruction event count or relabel this pilot as successful. Two
specific issues need to be addressed before another frozen run:

1. Selection resolution: in the flat clock approximation, only 1024 selection
   events make an interval with expected count two have proper-time ratio
   `(2/1024)^(1/4) ~= 0.2102`. Many center-to-outer separations are on this
   scale or smaller. This is an expectation calculation, not a deterministic
   lower bound on intervals with two observed events.
2. The current candidate score is computed before restricting to the selected
   set. Its degree may therefore poorly predict connectivity to the actual
   retained outer landmarks. Central coverage must be measured after selection
   and enforced prospectively; insufficient coverage must be a recorded design
   failure, not silently accepted as a three-band population.

A denser, still independent selection stage and a score based on connections
to the selected outer landmarks are concrete candidates. Their population,
coverage conditions, cost, and failure handling must be fixed before new data.
Independent accuracy/coverage confirmation and genuinely varying-profile
validation also remain outstanding. The active goal remains open.
