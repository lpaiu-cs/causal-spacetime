# Paper B claim boundary

A crisp checklist complementing manuscript Section 10. Every claim below is
backed by a frozen artifact under `docs/prereg/frozen/` (empirical claims)
or by a proved, CI-pinned statement in
`docs/theory/t1_parallax_identifiability.md` (theory claims, Section 8).

## Claimed (supported)

- The pipeline is a **discriminator** on 1+1D Minkowski causal sets at the
  stated scale: it passes the fixed gates on geometric order (PC-V1 H-SENS,
  20/20 confirmatory seeds) and blocks on matched geometry-free order (PC-V1
  H-SPEC, both control families 10/10).
- On geometric scenes the fit **recovers true spatial order** (calibration
  truth-order error mean 0.097; endpoint reproduction 19/20 in P1-B).
- Recovery **degrades monotonically** as geometry is diluted at fixed relation
  density (P1 H-MONO, 20/20 seeds with rho >= 0.85).
- The transition is **graded and identifiable**, not a cliff (P1 H-THRESH,
  median epsilon* = 0.31).
- **Embeddability alone over-reports geometry**: there is a false-pass window
  where profiles embed in 1D but no longer recover true space (P1 H-LAG,
  median crossing gap 0.19). A recovery/specificity check is therefore
  required.
- The degradation is **not explained by global relation density** (achieved
  density held at 0.566-0.584 across all epsilon; density is the only order
  statistic held fixed).
- The result is **robust to 2+1D**: the same pipeline passes and recovers true
  2D position, blocks on geometry-free order, and **selects the correct spatial
  dimension** (d=2: d=1 underfits, d=3 adds nothing) on 20/20 confirmatory
  seeds (P2-v2). The initial 2+1D run (P2) missed the sensitivity bar for a
  scene-generation reason and was remediated by preregistration, not retuning.
- **A geometry-free growth dynamics produces no reconstructable geometry**:
  transitive percolation blocks 100/100 across the preregistered density
  sweep — structurally at low density, at the geometry gate (null gap ~ 0)
  at higher density — while the same order-intrinsic protocol passes on
  sprinkled geometry (P3; separation d = 10.4).
- **The unrestricted ensemble has an exact crystalline low-action
  obstruction**: the complete bipartite crystal has smeared 2D BD action
  eps N (2 - eps N) exactly (test-locked identity), below the sprinkled value
  whenever the nonlocality scale resolves the system; exploratory sampling
  found no geometric window (this is a low-action state plus a sampling
  null, not a free-energy theorem). Corollary: Myrheim-Meyer dimension ~ 2
  for this crystal — a dimension estimator alone cannot exclude it.
- **The action-weighted restricted ensemble (2D orders) has a continuum
  phase** at beta below a hysteretic, first-order-like transition
  (P4: continuum 30/30 chains at beta <= 2 with dual-start agreement;
  crystal 10/10 at beta = 6 reaching the exact bipartite action; frozen pass
  criteria all met, sharp predictions 5/6 with the miss pre-flagged).
- **Post-burn-in samples of that continuum phase pass the full
  discriminator** (P5: 9/9 at beta = 2 and 9/9 at beta = 8, truth-order error
  ~ 0.13-0.14 vs 0.5 chance, gate statistics within the range of the
  uniform-ensemble calibration; the crystal control blocks 4/4 — structurally,
  at the chain-extraction stage, before any geometry gate). Emergence
  content: reconstructable geometry *survives* action weighting below beta_c
  and is *destroyed* above it.
- **Rejection happens at the geometry gates, not only at chain
  extraction**: eight preregistered chain-rich layered cells (N = 600,
  seeds 100-119) all met their frozen expectation — every seed reached
  the numerical gates (no structural blocks) and >= 16/20 per cell
  blocked there (149/160 total; P6a). The proposed local-shuffle family
  was retired by construction audit (coordinate remapping, not geometry
  destruction: original-coordinate gate 0/10 vs current-coordinate
  10/10 at 600 moves).
- **On identical data the instrument has the highest ROC AUC** (0.993 on
  362 labelled orders across P1/P3/P5/P6; height 0.967, MM 0.939,
  abundance 0.933), and MM dimension false-passes 25/27 of the P1
  false-pass window (P6b, all labels/metrics frozen before computation).
  NOT claimed: superiority over height as a ranking statistic on P1 —
  height is marginally more monotone there (median rho 0.994 vs 0.976).
- **The profile observable's identifiability is settled (1+1D)**: spatial
  order up to global reversal is decodable from the parallax dissimilarity
  alone (strict Robinson structure; already for two observers) — exactly
  for exact profiles, and on measured Model-D data for every comparison
  above the proved resolution/margin bounds (sub-resolution pairs may
  invert, and the instrument check records exactly those); spacings
  are provably NOT decodable from it (same-dissimilarity counterexample
  with affinely inequivalent targets); radial-distance error on the
  instrument's deterministic clock is bounded by half a tick spacing
  (position, via flanking differences, by one tick) with a 1/K RMSE
  law (measured slope -1.017); Poisson-clock order recovery concentrates
  at exponential rates (same-slice pairs never strictly invert, pathwise).
  Every proved statement is a deterministic CI regression -- Model-D
  claims checked against the instrument itself, Model-P claims by direct
  seeded simulation of the stated model (no frozen instrument realizes
  Model P).
- **The inverse-root density law is protocol-dependent**: it holds for
  Poisson-thinned clocks (lambda ~ rho, proved corollary; measured RMSE
  exponent -0.463) and fails to transfer to sprinkling-harvested clocks,
  whose rate couples at the discreteness scale (lambda ~ sqrt(rho),
  measured; error exponent -0.32). Any density-scaling statement must name
  its clock protocol.

## NOT claimed

- No claim of spacetime emergence, or that spacetime "is" information/order.
- No claim that causal order alone yields metric scale, signed coordinates,
  an affine structure, a manifold, or a unique observer atlas.
- No continuum-limit claim and no quantum dynamics. The dynamical orders are
  classical: sequential growth (P3) and classical statistical weighting
  exp(-beta S) on a restricted ensemble (P4/P5), a standard analytic
  continuation — not a quantum theory of gravity.
- No claim that geometry arises from nothing: the 2D-orders restriction is
  itself sprinkling-equivalent at beta = 0. The claims are survival below
  beta_c, destruction above it, and the two controlled nulls (dynamics
  alone; unrestricted action) showing each ingredient is necessary.
- "Pass" is gate-satisfiability by measured geometric order, NOT proof of a
  physical theory. "Block" is the expected null outcome on geometry-free
  order, NOT falsification of any theory.
- The absolute gate/threshold values are scale-dependent; only the contrast
  (pass vs block) and the monotone dose-response transfer. epsilon* is a
  property of this generator family and scale.
- Core results are 1+1D, single diamond, fixed observer layout, single
  dilution family; the emergence chain is one dynamics, one action family,
  one restricted ensemble at N <= 600. No curvature, quantum, or matter
  claim; beta_c and hysteresis strength are scale- and eps-dependent.
- Theory claims are 1+1D and instrument-relative. The Poisson-clock
  concentration statements describe an idealization no frozen instrument
  realizes; no inverse-root density law is claimed for harvested-chain
  clocks; the harvested protocol's tube selection reads embedded
  coordinates (coordinate-tube, not order-intrinsic), and its fluctuation
  class is an open question, not a result.
