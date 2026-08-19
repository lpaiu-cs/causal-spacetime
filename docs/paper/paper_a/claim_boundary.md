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
- **Auxiliary O4b instrument audit (frozen, recovered).** At the single
  frozen Schwarzschild configuration, the sampler + S1 volume response +
  certified oracle stack is CONCORDANT under the frozen composite error
  budget of 3.25% with simultaneous coverage >= 95%: G1 interval
  [56.448806185841875, 56.9822829864225] vs the certified
  [56.212737, 57.348019] (identified discrepancy
  [-0.8992123405928396, 0.7695461186219887], band 1.7034113309135284);
  G2 leak upper 0.14195058753928652 within budget 0.14195094424279403,
  zero leaking points. Published as `run_kind =
  recovered_completed_campaign` — the frozen decision functions re-applied
  to the completed run's preserved statistics: no new seed, no solver call,
  no resampling, no gate change (`docs/prereg/p14_o4b_results.json`). An
  instrument statement only; see the non-claims below.
- **Prediction-anchored Poisson count, per rung (preregistered, 3+1D).**
  On each preregistered rung of the mass ladder independently — `mu = 2M/r_c`
  in {0.1333, 0.1867, 0.2400} on the shell `[10, 20]` with anchors `(12, 18)`
  fixed in ABSOLUTE coordinates, so no rung is an isometric copy — the
  certified-membership count of a Poisson sprinkling realized that rung's
  independently certified continuum volume within `tau = 2.5%`, under a rule
  frozen before any count data existed. The gate is the identified
  discrepancy of a conservative outer Garwood interval (lower limit at
  `K_certain`, upper at `K_certain + U_amb`, each tail at level 0.025)
  rescaled by the intensity, required to lie in `+- tau*V_ref`; three-way
  verdict, 96-bit throughout. CONCORDANT at every rung:
  `mu = 0.1333` K = 53,285, U = 0, D = [-0.854796, +0.676594] in +-1.419420
  (`docs/prereg/p14_o5_count.json`);
  `mu = 0.1867` K = 53,829, U = 1, D = [-0.675256, +1.069450] in +-1.616487
  (`docs/prereg/p14_s6_m14_count.json`);
  `mu = 0.2400` K = 53,513, U = 0, D = [-1.016709, +0.987400] in +-1.858507
  (`docs/prereg/p14_s6_m18_count.json`).

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
  M = 1 for THESE stages — the mass ladder is the separate Section 6.8
  stage); NOT prediction-anchored — no diamond-volume oracle is used in the
  S4/S5 verdicts themselves. The oracle is now fully certified and executed
  (`docs/theory/schwarzschild_volume_oracle_note.md`,
  `docs/prereg/p14_o3_volume.json`), and the auxiliary O4b audit above
  consumed it — but that audit does not promote S4/S5 to
  prediction-anchored (the S4/S5 verdicts themselves use no oracle). The
  Poisson-count stage is a SEPARATE preregistered stage, now executed
  across the mass ladder (see the claimed block above); it does not
  promote S4/S5.
- Capstone non-claims: NOT a general-Weyl discriminator (Petrov type N, exact
  volume form is unrepresentative); NOT Weyl-tensor recovery (no reconstructed
  quantity — a boundary experiment, not a ladder rung); NOT box- or
  density-independent (licensed at the frozen operating point only). The
  Schwarzschild path now carries its own preregistered verdicts and
  non-claims (the Type-D blocks above); the diamond-volume oracle is
  certified and instrument-audited there (auxiliary), while the
  prediction-anchored Poisson-count stage remains open. A literature
  priority search (August 2026) identified partial overlap in four
  clusters — see `literature_priority_oracle_arc.md` for the full
  assessment. The directed-rounding volume enclosure and the
  preregistered TOST gate against it have no direct prior; Poisson
  sprinkle verification on Schwarzschild has prior art in Homšak–Veroni
  (arXiv:2404.11670 §V.1).
- Auxiliary O4b audit non-claims: the CONCORDANT verdict is a statement that
  the instrument stack agrees with the certified continuum volume within the
  frozen margin at ONE frozen configuration. It is NOT a prediction-anchored
  promotion of S4/S5, NOT a C1/C2 joint verdict, NOT a Poisson causal-set
  count verification (sprinkle counts vs rho*V is a separate stage, executed
  separately across the mass ladder and not by this audit),
  NOT mass- or domain-generality (single frozen M = 1 and box), NOT a
  confirmation of causal-set theory, and NOT complete separation or general
  volume accuracy. Nothing in Sections 6-6.7 is upgraded by it.
- Prediction-anchored Poisson-count non-claims: the three rung verdicts are
  SEPARATE preregistered stages — there is NO joint or composite cross-rung
  verdict and no pooled statistic. NOT a statement at any compactness between
  or beyond the tested rungs: no interpolation, no extrapolation, no fitted
  trend in `mu` is claimed (the monotone growth of the certified volume with
  `mu` is geometry, not a result). NOT a confirmation of causal-set theory
  and NOT evidence that spacetime is discrete — a count agreeing with a
  volume says the instrument reproduces the continuum measure it was pointed
  at. NOT a promotion of the Section 6-6.7 verdicts or of the auxiliary O4b
  audit, none of which use it. NOT density- or domain-generality: one fixed
  diamond and one fixed intensity per rung, on a Schwarzschild exterior
  shell. The tolerance `tau = 2.5%` was chosen to be FEASIBLE at the frozen
  sizing, not to be sharp, so CONCORDANT means "inside that band", not
  "agreement to 2.5% is the instrument's accuracy".
