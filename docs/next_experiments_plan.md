# Next experiments plan — post Paper B v0.5

Date: 2026-07-15. Status: **P6 COMPLETE; P7 FROZEN AT N=600
CHARACTERIZATION; T1 ACTIVE.** This remains a roadmap rather than a
preregistration. Frozen protocols and outcomes live in the
experiment-specific `docs/prereg/p*.md` files.

Execution update (2026-07-14):

- P6a constructed hard negatives completed on fresh seeds; all eight frozen
  cells remained chain-rich and blocked as predicted. The local-shuffle arm
  was audited as coordinate relabelling and was not promoted to a false hard
  negative.
- P6b same-data comparison completed. The instrument had ROC AUC 0.993 versus
  0.967 for height, 0.939 for Myrheim-Meyer, and 0.933 for interval abundance;
  height remained slightly more monotone on P1, so no blanket superiority is
  claimed.
- P6a/P7 N=600 characterization completed. The continuous score G was frozen
  from existing gate margins and validated on frozen data. After the
  review-driven Geyer IAT correction (logged in the P6a/P7 deviations logs),
  only beta 12 passes the frozen mixing screen; beta 14 is marginal (minimum
  ESS 19.7 against the frozen 20) with full phase agreement, and should be
  re-judged with longer chains rather than dropped. Beta 16--28 showed severe
  metastability and a broad start-dependent G hysteresis window.
- P7 N=900/1200 confirmatory FSS is not yet justified. A validated enhanced-
  sampling method with overlap/convergence evidence must precede scaling;
  ordinary local chains or a coarse replica ladder do not resolve the observed
  action barrier.

Execution update (2026-07-15):

- Wang-Landau / multicanonical sampling was built, validated against exact
  enumeration at N=6, and shown to tunnel the N=600-blocking action barrier
  at N <= 80 -- then eliminated at N=600 by a measured tunneling exponent,
  tau ~ N^5.88 with 90% CI [5.44, 6.38]: even the optimistic end costs
  ~278,000 core-hours for a single converged ln g(S). Details and data in
  `docs/p7_enhanced_sampling.md`. This closes the third quantified
  elimination for the N=600 barrier (after local Metropolis and the replica
  ladder arithmetic).
- **Decision: P7 is frozen at its N=600 characterization.** No further P7
  compute until either a qualitatively different sampler (nonlocal/cluster
  moves) is found, or the FSS axis is re-scoped to the demonstrably
  samplable range N <= ~120 under a fresh prereg freeze. Neither choice is
  made now.
- **T1 (parallax identifiability, theory) is the active track.**

Execution update (2026-07-16):

- T1 completed through v0.5: Theorems 1-2 and Lemmas 1-4 all [PROVED]
  and pinned in CI (PRs #5-#7); G1/G3 closed; G2 instrumented with the
  rho^{-1/2} law shown protocol-dependent (PR #8). Remaining T1 items
  are open questions (harvested-chain fluctuation class, order-only
  harvest design) and G4 (2+1D), not blockers.
- **The P7 FSS re-scope route is now closed** (fourth quantified
  elimination): reconnaissance showed the sampler-feasible window
  (N <= ~120) and the instrument-operable window (N >= ~500 for the
  protocol as frozen, greedy selector included — an optimal selector
  could lower that boundary but not below the ~156 envelope, under
  which even one 25-tick chain typically does not exist; N >= ~160
  marginal for the most permissive re-scoped spec) are disjoint. The compute budget was
  not the binding constraint. Details and dual-criterion verdict in
  `docs/p7_fss_rescope_recon.md`; probe and tracked results committed.
  P7 reopens only with a sampler reaching N >= 500 or a qualitatively
  new small-N geometry instrument (own calibration + freeze).

Context: the second external review judged the current Paper B arXiv-ready
after claim calibration (done, v0.5) and packaging, a strong candidate for
CQG / SciPost Physics Core as-is, and a defensible SciPost Physics submission
only with one "A-tier lever" completed. It also identified two validity gaps
that are cheap to close: the crystal control is structurally easy (blocks at
chain extraction, never exercising the geometry gates), and the discriminator
has not been compared head-to-head with existing manifoldlikeness diagnostics
on the same data. Both are recorded as future work in manuscript Section 10.

## P6 — hard negatives and diagnostics head-to-head

**Question.** Does the discriminator's geometry verdict (not just its
chain-extraction gate) reject chain-rich non-geometric order, and what does
it add over cheaper diagnostics on identical data?

### P6a: chain-rich hard negatives

Three candidate families, all within or adjacent to existing machinery
(`positive_control/two_orders.py`, the frozen order-intrinsic protocol):

1. **Near-critical 2D orders** — beta in {12, 16, 24} at N = 600
   (reconnaissance located beta_c in (8, 32)). Chain-rich but possibly
   partially layered. The verdict here is *not* known in advance: states
   below beta_c may legitimately pass. This arm is therefore a
   characterization arm (expected verdict conditional on the measured phase),
   and it doubles as the first beta_c refinement at N = 600 — reusable
   directly by P7.
2. **Constructed partially layered permutations** — k-block layered
   permutations interpolating continuum (k -> N) to bipartite crystal
   (k = 2), analogous to PC-V1's hand-built controls. For moderate k these
   are chain-rich. Prediction to test in Stage A and freeze: a k window
   exists where chain extraction succeeds but the geometry gates (held-out,
   null gap) block — the missing "blocks at the geometry gate, not
   structurally" negative.
3. **Locally shuffled continuum states** — frozen P5 continuum
   configurations with windowed transpositions that destroy local geometric
   consistency at approximately fixed global statistics: the in-ensemble
   analogue of P1's dilution.

Design: Stage A on families 2-3 (characterize which parameters give
chain-rich states; calibrate), freeze per-family expected verdicts, Stage B
on fresh seeds. Family 1 is the only MCMC arm (3 chains x 3 beta at N = 600).
Cost estimate: families 2-3 are constructions, hours; family 1 about one
desktop day with the incremental sampler (the P5 confirmatory stage at
N = 600 took roughly 5-7 h wall-clock with shard parallelism).

### P6b: same-data comparison with existing diagnostics

Analysis-only; no new sampling. All inputs regenerate deterministically from
recorded seeds (Section 11 of the manuscript).

- **Inputs:** P1 epsilon-sweep scenes, P3 percolation orders, P5
  configurations, plus the P6a negatives.
- **Diagnostics:** Myrheim-Meyer dimension, interval-abundance profile
  distance to the sprinkled reference, order height; optional (heavier,
  literature-cited): spacelike-distance constructions, homology-based
  manifoldlikeness.
- **Metrics:** ROC/AUC for geometric-vs-not (labels: P1 low-epsilon vs
  high-epsilon endpoints, P3, crystal), truth-rank correlation against
  epsilon on P1, and overlap of each diagnostic's false-pass region with the
  frozen H-LAG window.
- **Anti-cherry-picking rule:** the diagnostic list and the metric set are
  fixed *before* computing anything on confirmatory data (a lightweight
  freeze, since this is comparative analysis rather than a new confirmatory
  claim about the instrument).

Deliverable: one table + one figure; fulfills the Section 10 promise and the
review's "what does the new discriminator add" demand. Candidate Paper B
v0.6 section or supplement.

## P7 — geometry order parameter G(beta, N, eps) and finite-size scaling

The recommended A-tier lever. **Question:** does instrument-level geometry
loss coincide with the thermodynamic crystallization transition, or precede
it? Either answer is new physics: the reconstruction verdict becomes a
continuous observable of the 2D-orders phase diagram, which existing studies
characterize only by macroscopic observables.

- **G definition (choose and freeze in Stage A):** a per-configuration
  continuous score built from the already-frozen gate quantities — e.g. a
  clipped margin composite of (gate - heldout), (null_gap - gate), and
  (0.5 - truth)/0.5. Requirements to validate on *existing frozen data*
  before freezing: G ~ 1 on uniform-ensemble calibration, G ~ 0 on crystal
  and column-shuffle, monotone under P1-style dilution. Convention: G = 0
  when six 25-tick chains cannot be extracted (makes G total, including the
  structural-block regime).
- **Design axes:** N in {600, 900, 1200} at fixed eps N = 12
  (eps = 0.02, 0.0133, 0.01); per N, a reconnaissance pass brackets
  beta_c(N) (dual-start bisection), then a dense grid of ~8 beta spanning
  the bracket; dual starts x >= 3 seeds per point.
- **Chain-quality instrumentation (closes the review's ESS gap):** measured
  integrated autocorrelation time / effective sample size per chain
  (blocking analysis), dual-start convergence check per point, and a
  replica-exchange fallback if mixing collapses near beta_c.
- **Thermodynamic observables alongside:** <S>, n0 susceptibility, a Binder
  cumulant (n0 or height), action histograms (double-peak = first-order
  signature), hysteresis width vs N.
- **Confirmatory hypotheses (to freeze after Stage A):** H-G-MONO
  (G decreases in beta at each N), H-COINCIDE (relative position of
  beta_c^geo — where G crosses its frozen threshold — versus beta_c^thermo,
  stated as an ordered prediction), H-HYST (first-order character persists
  at N >= 900).
- **Analysis caution:** the P4 transition is first-order-like; the analysis
  plan must test first-order scaling (hysteresis persistence, double-peaked
  histograms, interface arguments) rather than assume continuous-transition
  power laws.
- **Cost (estimate):** incremental sampler is O(N^2)/move, so N = 1200
  chains cost ~4x the N = 600 chains; the full grid
  (3 N x ~10 beta x 2 starts x 3 seeds) is multi-week on one desktop.
  Staged budget: (i) N = 600 full grid first — its beta_c refinement is
  shared with P6a family 1; (ii) N = 900/1200 on a reduced grid around
  beta_c(N). Parallel shards as in P5.

## T1 — parallax identifiability and stability (theory track)

A-tier lever 2; no compute, 1-2 focused weeks; can run in parallel.

- **Goal:** state and prove conditions under which bracket-width echo
  profiles on observer chains determine the targets' spatial order up to
  reflection (and scale), plus a stability bound under Poisson sprinkling
  noise (expected error decreasing in sprinkling density).
- **Sketch:** in a 1+1D diamond the expected bracket width is a monotone
  function of |x_target - x_observer|; parallax centering removes the shared
  time scalar; monotonicity across >= 3 non-collinear observers pins the
  spatial order; concentration of Poisson counts gives an error bound
  ~ 1/sqrt(rho * area). Numerical verification harness on the existing
  generators (predicted vs measured error scaling).
- **Positioning:** complements the 2026 embedding-uniqueness theorem
  (uniqueness *given* an embedding) with a finite, order-intrinsic
  reconstruction criterion — the pairing the review called natural.

## Recommended sequence

1. **P6** (families 2-3 construction + P6b metric freeze and analysis;
   family 1 MCMC last): closes both review validity gaps cheaply ->
   Paper B v0.6.
2. **P7** Stage A at N = 600 (reusing P6a family 1) -> freeze -> staged
   Stage B across N: the SciPost Physics lever.
3. **T1** in parallel as bandwidth allows.
4. **Packaging** (user-gated, independent of experiments): LaTeX/PDF build,
   authors/affiliations, LICENSE, CITATION.cff, tagged release + Zenodo DOI,
   arXiv submission.

## Non-goals this cycle

Curvature recovery (review lever 4: high risk, premature before P7);
quantum/complex weights; any claim upgrades in Paper B without a
corresponding confirmatory pass; moving any frozen gate.
