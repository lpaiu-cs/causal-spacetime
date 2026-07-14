# Paper B (draft)

Working title: **A validated discriminator for latent geometry in discrete
causal order, and its response to geometry dilution.**

Status: DRAFT v0.5 (2026-07-14). Prose grounded entirely in the frozen
preregistration artifacts under `docs/prereg/frozen/` (PC-V1, P1, P2/P2-v2,
P3, P4, P5). No number in the manuscript is from memory; each traces to a
committed registry or CSV.

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

It is deliberately conservative: it does not claim continuum spacetime
emergence or quantum dynamics; the emergence claims are survival/destruction
of reconstructable geometry with each necessity backed by a controlled null.
See `manuscript.md` Section 9 (claim boundary).

## Files

- `manuscript.md` — the v0.5 draft (all sections incl. the Section 7 emergence chain, a Conclusion, and Appendix tables;
  reproducibility in Section 11). Inline citations are pandoc keys matching the
  bibliography.
- `claim_boundary.md` — crisp claim/non-claim checklist (complements Section 9).
- `figures/make_figures.py` — regenerates every figure from the frozen CSVs
  (and, for Fig 3, from `confound_data.csv`).
- `figures/compute_confound_data.py` — recomputes the Fig 3 confound data (raw
  vs parallax dissimilarity on the Stage C seeds); the parallax column
  reproduces the frozen Stage C registry exactly.
- `figures/*.png` — generated figures (committed for convenience). Fig 1
  discriminator separation, Fig 2 dose-response, Fig 3 2+1D dimension
  selection, Fig 4 shared-scalar confound.
- `citations/references.bib` — verified bibliography (every entry confirmed
  against an authoritative source; no fabricated field).
- `citations/citation_verification_report.md` — per-entry source + corrections.

## Provenance

- PC-V1 frozen at commit `b77f588` (calibration `9162e8e`).
- P1 frozen at commit `a218d9a` (calibration `6b21bb7`).
- Confirmatory results: PC-V1 `891498f`, P1 `4c05cf2`.
