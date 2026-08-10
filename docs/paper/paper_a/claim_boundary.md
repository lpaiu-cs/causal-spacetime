# Paper A claim boundary

Complements manuscript Section 8. Every claim is a controlled validation in a
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
- **Capstone C1 (preregistered, 3+1D).** At the frozen operating point of the
  pure-Weyl plane-wave construction, the paired ensemble mean of the global
  relation-fraction change exceeds its frozen margin: 0.0502929
  [0.0501046, 0.0504812] vs epsilon_Delta = 3.579e-4
  (`docs/prereg/p14_prereg_results.json`).
- **Capstone C2 (preregistered replication, 3+1D).** A frozen single-poset
  classifier separates curved from flat ensembles under three frozen interval
  rules (s = 11.199 [11.035, 11.362]; AUC = 1.0 [0.999232, 1.0];
  BA = 1.0 [0.986, 1.014]), independently replicating the probe chain's
  confirmation (same file).
- **Type-D extension C1 (preregistered, paired).** On the frozen
  Schwarzschild coordinates and S1 domain, the same-point-set paired mean
  shift passes the frozen detection gate (identified CI95
  [-0.036211, -0.035953] vs threshold -0.0036) and quantitatively replicates
  the exploration block (Welch CI95 of the difference inside +-0.0012):
  stage verdict CONFIRMED, as a program-internal statement
  (`docs/prereg/p14_s4_results.json`).
- **Type-D extension C2 (preregistered, single-poset, separate stage).** On
  the same frozen domain with two independent arms, the per-reading global
  relation fraction discriminates flat from Schwarzschild ensembles above
  the independently declared 0.60 threshold: AUC 0.9734, DeLong CI95
  [0.9630, 0.9837] — outcome DETECTED, with incomplete separation; secondary
  out-of-sample BA 0.9033, CP-Bonferroni CI95 [0.8356, 0.9500]
  (`docs/prereg/p14_s5_results.json`).

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
- The submission gate — all 19 cited legacy summary tables committed and
  digest-locked plus the finalized Section 6 capstone evidence bundle (which
  itself carries no regeneration debt) — is met; `artifact_manifest.md` is
  the authority for the archived state.
- Type-D extension non-claims: the C1 and C2 results arise from SEPARATE
  preregistered stages on the same frozen domain — there is no joint primary
  verdict, and "confirmed" applies to C1 only (the C2 outcome is DETECTED);
  NOT complete separation (AUC ≈ 0.973, unlike the plane-wave C2, and no
  completeness gate was preregistered); the secondary BA verdict does not
  strengthen or combine with the primary; NOT mass-general (single frozen
  M = 1); NOT prediction-anchored — no diamond-volume oracle is used, and the
  oracle path is only partially derived
  (`docs/theory/schwarzschild_volume_oracle_note.md`).
- Capstone non-claims: NOT a general-Weyl discriminator (Petrov type N, exact
  volume form is unrepresentative); NOT Weyl-tensor recovery (no reconstructed
  quantity — a boundary experiment, not a ladder rung); NOT box- or
  density-independent (licensed at the frozen operating point only); NO verdict
  on the Schwarzschild path (S1 prices one solver path only; the diamond-volume
  oracle is open). Priority is claimed only within this program — no
  literature-wide first before a priority search.
