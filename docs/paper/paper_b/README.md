# Paper B (draft)

Working title: **A validated discriminator for latent geometry in discrete
causal order: dilution response and survival in an action-weighted
ensemble.**

Status: DRAFT v0.6 (2026-07-17). Empirical prose grounded entirely in the
frozen preregistration artifacts under `docs/prereg/frozen/` (PC-V1, P1,
P2/P2-v2, P3, P4, P5, P6a/P6b); the Section 8 theory results are analysis-only and
grounded in the committed, CI-pinned theory artifacts (`docs/theory/`,
nothing frozen there by design). No number in the manuscript is from
memory; each traces to a committed registry, CSV, or tracked theory
artifact.

## Scope

This paper reports the validated-instrument line of the program:

1. **PC-V1** validates the response-profile representability pipeline as a
   *discriminator* for latent low-dimensional geometry in discrete causal
   order: it passes on measured Minkowski-sprinkled causal sets and blocks on
   matched geometry-free order.
2. **P1** makes geometry the manipulated variable and measures the
   dose-response: at fixed relation density, geometry recovery degrades
   monotonically as the order is diluted from Minkowski to geometry-free, with
   an identifiable graded transition and a "false-pass" window where
   embeddability outlasts true-geometry recovery.
3. **P2/P2-v2** extend the discriminator to 2+1D, where it selects the
   correct spatial dimension (d = 2) on 20/20 confirmatory seeds.
4. **The emergence chain (P3, P4, P5, plus an exact obstruction)** carries the
   validated instrument to orders not built by hand: a geometry-free growth
   dynamics blocks 100/100; action weighting of unrestricted orders meets an
   exact crystalline low-action obstruction (bipartite crystal,
   S_eps = eps N (2 - eps N)) with no geometric window found in exploratory
   sampling; the action-weighted restricted ensemble (2D orders) has a
   continuum phase with a hysteretic crystallization transition; and
   post-burn-in samples of that continuum phase pass the full discriminator
   (18/18) while the crystal control blocks structurally (4/4).
5. **Validity hardening (P6a/P6b, new in v0.6)** closes the two review
   validity gaps: constructed chain-rich layered negatives block at the
   geometry gates themselves (8/8 preregistered cells, 149/160 fresh
   seeds; the local-shuffle candidate was retired by construction audit
   as coordinate remapping), and a preregistered same-data head-to-head
   gives the instrument the highest ROC AUC with its full gate set
   (0.990 vs height 0.967, MM 0.939, abundance 0.933), with MM
   dimension false-passing 25/27 of the P1 false-pass window — while
   height stays slightly more monotone on P1, so no blanket superiority
   is claimed. Removing the margin's truth term drops it to 0.968, a
   tie with height: the aggregate ranking *does* depend on truth
   assistance, correcting an earlier claim that rested on a
   recomputation which had dropped a second, truth-independent gate.
   Manuscript Sections 7.5-7.6.
6. **The identifiability theory (T1, new in v0.6)** proves what the profile
   observable can and cannot identify: spatial order up to global reversal
   is decodable from the parallax dissimilarity alone (already for two
   observers; margin-qualified on measured data), spacings provably are not
   (explicit same-dissimilarity counterexample), radial-distance error on
   the instrument's deterministic clock is resolution-limited (delta/2
   pointwise; position via flanking differences, one tick), Poisson-clock
   order recovery concentrates at exponential rates, and the inverse-root
   density law for the error is protocol-dependent (holds for thinned
   clocks, fails for sprinkling-harvested clocks; the order-only
   anchored harvest is measured wandering-dominated at the KPZ exponent
   -1/6, closing the order-only design question). Every proved statement
   is a CI regression — Model-D claims against the instrument itself,
   Model-P claims by direct seeded simulation of the stated model;
   every fitted scaling exponent carries a residual-based interval and
   a split-half check. The last open question, the count-fluctuation
   class, is closed: measured from chain lengths rather than from the
   distance error it comes out protocol-dependent -- Tracy-Widom for
   the order-only chain, Poisson for the tube-confined one.
   Its dimensional reach is split and stated (Section 8.5): the
   identity, band, resolution law and concentration hold in any spatial
   dimension and are verified in 2+1D on the frozen P2/P2-v2 scene
   builder, and labeled identifiability there is multilateration from
   three non-collinear observers. The unlabeled result changes
   character with dimension: ordinal in 1+1D, metric in 2+1D, where
   four or more observers make the dissimilarity determine the whole
   scene up to congruence while three provably do not (characterization
   at stated strength -- infinitesimal rigidity, exact model).
   That threshold is not a 2+1D fact but a general-dimension law
   (Section 8.6): the observable rigidifies exactly from `R >= d + 2`
   observers, so physical 3+1 needs five, from 19 targets, and 1+1D is
   the lone exception because its profile surface is a polyline of zero
   curvature. Measured across spatial dimensions one through six
   against 20 predictions preregistered in four rounds and never
   amended afterwards, 20/20 -- including the two that carry
   information: a fibre slope of 2 where every prior measurement had
   slope 1, and the `R = d + 2` threshold formula carried out of sample
   to `d = 5` (predicted 41, measured 41). Most of that law is now
   *derived* rather than observed: the flex condition reduces exactly to
   infinitesimal rigidity of the centered profile cloud, which yields a
   closed form below the threshold and a proved lower bound above it,
   both with no free parameter and both reproducing every measured cell.
   The bound is attained with no slack in all 15 observer counts tried
   across six spatial dimensions. Attainment is proved too, by way of
   the profile surface's conical singularities -- one per observer,
   which no continuous ambient symmetry can permute, and which span
   because the un-squared observer distance matrix is nonsingular. That
   brackets the threshold from above, and the bracket closes exactly at
   `R = d + 2`: **physical 3+1 is pinned at five observers and 19
   targets, both theorems**, modulo one explicit rank hypothesis on
   observer configurations (verified 18/18 across `d = 2..7`).
   Sharpness for `R > d + 2` stays open, though rigidity there does not.
   `d = 1` is exempt in the proof for the same reason it was exempt in
   the measurements: its cone degenerates to two rays. Scope is
   unchanged by the extra dimensions, the derivation, or the proof.
   Section 8.7 adds the operating conditions, measured rather than
   assumed: at the pipeline's own configuration a half-tick readout
   bound costs about two of itself in position and falls linearly, so
   the metric result survives finite resolution -- but the Jacobian
   margin IS the error budget, crowding observers destroys it, and
   `R = d + 2`, the corner whose threshold is pinned, is the worst place
   to operate, worst of all in 3+1D. The theorem bounds how few
   observers are possible, not how many are wise. Descriptive, and
   narrow: only the readout error's *support* is proved, while its
   uniformity and its independence across target-observer pairs are
   assumptions -- the second load-bearing, since the parallax centering
   annihilates common-mode error outright and the measured cost runs
   from about two readout bounds down to zero across correlation
   structures. Nothing preregistered there and no gate consumes it.
   Source:
   `docs/theory/t1_parallax_identifiability.md` (v1.5),
   `docs/theory/t1_g4c_predictions{,_round2,_round3,_round4}.json`, and manuscript
   Section 8.

It is deliberately conservative: it does not claim continuum spacetime
emergence or quantum dynamics; the emergence claims are survival/destruction
of reconstructable geometry with each necessity backed by a controlled null.
See `manuscript.md` Section 10 (claim boundary).

## Files

- `manuscript.md` — the v0.6 draft (all sections incl. the Section 7 emergence
  chain and the Section 8 identifiability theory, a Conclusion, and Appendix
  tables; reproducibility in Section 12). Inline citations are pandoc keys
  matching the bibliography.
- `claim_boundary.md` — crisp claim/non-claim checklist (complements Section 10).
- `figures/make_figures.py` — regenerates every figure from committed
  artifacts: the frozen CSVs and summaries (Figs 1-5, 7; Fig 7 also reads
  `confound_data.csv`) and the two tracked theory tables (Fig 6). It
  asserts the gate-complete and order-only AUCs against the dated P6b
  correction table, re-derives the historical frozen proxy, and checks
  the zero band violations and the slope window at read time.
- `figures/compute_confound_data.py` — recomputes the Fig 7 confound data (raw
  vs parallax dissimilarity on the Stage C seeds); the parallax column
  reproduces the frozen Stage C registry exactly.
- `figures/*.png` — generated figures (committed for convenience). Fig 1
  discriminator separation, Fig 2 dose-response, Fig 3 2+1D dimension
  selection, Fig 4 P6b diagnostics head-to-head (ROC + H-LAG safety),
  Fig 5 emergence chain (P3 | P4 | P5), Fig 6 theory (quantization band,
  1/K law, density scaling by protocol), Fig 7 shared-scalar confound.
- `citations/references.bib` — verified bibliography (every entry confirmed
  against an authoritative source; no fabricated field).
- `citations/citation_verification_report.md` — per-entry source + corrections.

## Provenance

- PC-V1 frozen at commit `b77f588` (calibration `9162e8e`).
- P1 frozen at commit `a218d9a` (calibration `6b21bb7`).
- Confirmatory results: PC-V1 `891498f`, P1 `4c05cf2`.
