# Paper B (draft)

Working title: **A validated discriminator for latent geometry in discrete
causal order, and its response to geometry dilution.**

Status: DRAFT v0.1 (2026-07-09). Prose grounded entirely in the frozen
preregistration artifacts under `docs/prereg/frozen/` (PC-V1 and P1). No number
in the manuscript is from memory; each traces to a committed registry or CSV.

## Scope

This paper reports two controlled results on the `restart/m17-baseline` line:

1. **PC-V1** validates the response-profile representability pipeline as a
   *discriminator* for latent low-dimensional geometry in discrete causal
   order: it passes on measured Minkowski-sprinkled causal sets and blocks on
   matched geometry-free order.
2. **P1** makes geometry the manipulated variable and measures the
   dose-response: at fixed relation density, geometry recovery degrades
   monotonically as the order is diluted from Minkowski to geometry-free, with
   an identifiable graded transition and a "false-pass" window where
   embeddability outlasts true-geometry recovery.

It is deliberately conservative: it does not claim spacetime emergence, a
continuum limit, or dynamics. See `manuscript.md` Section 7 (claim boundary).

## Files

- `manuscript.md` — the v0.2 draft (all sections + Appendix A parameter table;
  reproducibility in Section 9). Inline citations are pandoc keys matching the
  bibliography.
- `claim_boundary.md` — crisp claim/non-claim checklist (complements Section 7).
- `figures/make_figures.py` — regenerates every figure from the frozen CSVs.
- `figures/*.png` — generated figures (committed for convenience).
- `citations/references.bib` — verified bibliography (every entry confirmed
  against an authoritative source; no fabricated field).
- `citations/citation_verification_report.md` — per-entry source + corrections.

## Provenance

- PC-V1 frozen at commit `b77f588` (calibration `9162e8e`).
- P1 frozen at commit `a218d9a` (calibration `6b21bb7`).
- Confirmatory results: PC-V1 `891498f`, P1 `4c05cf2`.
