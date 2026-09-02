# Paper A — journal-format draft (Classical and Quantum Gravity target)

A submission-format LaTeX rendition of `../manuscript.md`, typeset to the
IOP/CQG iopart 12pt preprint conventions. The Markdown manuscript remains
the source of truth — it is what the repository's contract tests
(`tests/test_paper_a_results_integration.py`,
`tests/test_paper_a_count_integration.py`) hold to the committed
artifacts. This directory renders it; it must never disagree with it.

## Build

```sh
./build.sh          # latexmk -xelatex twice: main.tex, then si.tex
```

Output: `main.pdf` (the article) and `si.pdf` (the supplementary-material
document). TeX Live 2025 suffices; the font is Latin Modern, the iopart
preprint face.

## Layout

- `main.tex` — article assembly only; all content in `sections/*.tex`.
- `si.tex` — supplementary-material assembly (sections S1-S2 = the
  manuscript's appendices A-B, same files under `sections/`); its
  references into the article resolve from `main.aux` via `xr-hyper`, so
  it must build second.
- `cqgmimic.sty` — every journal-specific typesetting decision (iopart
  look-alike: title block, run-in Abstract, `1.`/`1.1.` headings, 5-pica
  equation indent, `Figure 1.`/`Table 1.` captions, running heads,
  appendix `A1` numbering). Swapping to the real `iopart.cls` later means
  changing `\documentclass` and deleting this file; sections stay.
- `sections/front.tex` — title block; the abstract is the manuscript's
  abstract verbatim (unified 2026-08-22, under the CQG 300-word limit),
  rendered in LaTeX notation only.
- `sections/references.tex` — hand-built IOP Vancouver bibliography,
  transcribed from `../citations/references.bib` (note fields omitted),
  entry order pinned to first-citation order. Re-derive the order from
  `main.aux` if a `\cite` is added or moved.
- `figures/make_journal_figures.py` — regenerates all seven vector PDFs and
  matching 300 dpi PNGs from the
  same committed inputs the manuscript figures use (`../figures/data/*.csv`
  and `docs/prereg/*.json`); the two new result figures (capstone,
  Schwarzschild) draw the raw per-reading arrays stored inside the frozen
  artifacts, and the setup figure is an illustration (fixed seed, no
  measured quantity). Run from the repository root with the project venv.
- `check_numbers.py` — fidelity gate: every numeric token in the LaTeX
  must exist in `manuscript.md` and vice versa (normalized, string-exact).
  Run from anywhere: `python docs/paper/paper_a/latex/check_numbers.py`.

## Deliberate deviations from the manuscript (all register-level)

- British spelling with `-ize` endings, per IOP house style.
- Three method citations the manuscript text does not carry inline
  (`delong1988`, `garwood1936`, `clopperpearson1934`) are added at first
  mention — all from the verified bibliography.
- Two new result figures (plane-wave capstone; Schwarzschild S4/S5) and a
  genre-standard setup illustration, none previously in the manuscript.
- The two adjacent "provenance is in Appendix B" sentences at the end of
  section 6.7 are merged into one.
- A CQG-required Data availability statement is added (repository is
  public; no DOI claimed).
- The manuscript's appendices A-B are packaged as a separate
  supplementary-material document (`si.pdf`, sections S1-S2) per the
  submission plan; the manuscript keeps them as appendices, and body
  references render as "the supplementary material (section S1/S2)".
  Note CQG's own default is in-article appendices — reverting is two
  `\input` lines plus the appendix preamble (see git history).
