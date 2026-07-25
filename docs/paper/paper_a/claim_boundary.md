# Paper A claim boundary

Complements manuscript Section 7. Every claim is a controlled validation in a
known model and tied to a cited experiment output. The clean-checkout evidence
boundary is recorded separately in `artifact_manifest.md`.

## Claimed (supported, as controlled validations)

- **Dimension** is estimated from causal-order statistics in flat
  Poisson-sprinkled Alexandrov intervals, with lower endpoint RMSE at N=2400
  than N=300 but non-monotonic finite-sample fluctuations in 3D/4D (exp10).
- **Timelike proper time** is recoverable from Alexandrov interval cardinality
  once an event density is supplied, with error consistent with finite-sampling
  noise (exp07, exp03, exp08).
- **Raw longest-chain length** is order-theoretic; its conversion to timelike
  proper time follows Brightwell-Gregory scaling only after supplied density
  and dimension-dependent normalization (exp09).
- **Radar time and unsigned distance** are recoverable from an observer chain
  with clock labels; error falls with tick density (exp11).
- **Signed coordinates and the Lorentz map** are recoverable with an added
  orientation reference (exp13).
- **Atlas transition maps** (Poincare) are approximately consistent across
  overlapping oriented, calibrated charts (exp14, exp15).
- **Volume under the conformal ambiguity** is recoverable with supplied measure
  information, implemented in exp19 as local weights; in the tested
  random-thinning protocol, reconstruction is stable after density rescaling
  (exp18-20, exp23).
- **Rindler wedge** is the reconstructible region for an accelerated observer; a
  horizon appears as a reconstruction-inaccessibility boundary (exp16, exp17).

## Negative results (bounding the ladder)

- Causal order alone does NOT fix conformal/absolute scale (needs a measure).
- A single observer gives only UNSIGNED distance — reflection degeneracy (needs
  an orientation reference) (exp12).
- Finite signal speed alone does NOT imply Lorentzian structure (finite-speed
  lattice counterexample) (exp05).
- The spacelike-distance proxy is exploratory and boundary-dependent, NOT a
  validated estimator (exp06).

## NOT claimed

- No claim that spacetime is reducible to, or emerges from, causal order.
- No claim that causal order alone yields absolute scale, the conformal factor,
  signed coordinates, or a unique observer atlas.
- No claim that finite signal speed implies relativity.
- Reconstructing a geometry that was put into the model (by sprinkling from
  Minkowski or by a supplied protocol) is NOT deriving geometry from order.
- The emergence/representability question (order not built from a geometry) is
  explicitly out of scope and is the subject of the companion Paper B.
- Results are controlled and mostly 1+1D (dimension checked to 4D); constants
  are convention-dependent and stated where used.
- Paper B's N=600 action-weighted samples are not certified equilibrium draws;
  no equilibrium transition or finite-size scaling is imported into Paper A.
- Fourteen of the 19 cited summary tables are not yet archived in the clean
  checkout; manuscript submission requires the complete evidence package.
