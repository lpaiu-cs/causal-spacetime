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
   gives the instrument the highest ROC AUC (0.993 vs height 0.967,
   MM 0.939, abundance 0.933; unchanged with the margin's truth term
   removed), with MM dimension false-passing 25/27 of
   the P1 false-pass window — while height stays slightly more monotone
   on P1, so no blanket superiority is claimed. Manuscript Sections
   7.5-7.6.
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
   at stated strength -- infinitesimal rigidity, exact model). Source:
   `docs/theory/t1_parallax_identifiability.md` (v1.0) and manuscript
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
  asserts the frozen registry AUCs, the order-only AUC, the zero band
  violations, and the slope window at read time.
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
