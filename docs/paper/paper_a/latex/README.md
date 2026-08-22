# Paper A — journal-format draft (Classical and Quantum Gravity target)

A submission-format LaTeX rendition of `../manuscript.md`, typeset to the
IOP/CQG iopart 12pt preprint conventions. The Markdown manuscript remains
the source of truth — it is what the repository's contract tests
(`tests/test_paper_a_results_integration.py`,
`tests/test_paper_a_count_integration.py`) hold to the committed
artifacts. This directory renders it; it must never disagree with it.

## Build

```sh
./build.sh          # latexmk -xelatex; needs XeLaTeX (Korean frozen sentences)
```

Output: `main.pdf` (25 pp). TeX Live 2025 suffices; the font is Latin
Modern, the iopart preprint face.

## Layout

- `main.tex` — assembly only; all content in `sections/*.tex`.
- `cqgmimic.sty` — every journal-specific typesetting decision (iopart
  look-alike: title block, run-in Abstract, `1.`/`1.1.` headings, 5-pica
  equation indent, `Figure 1.`/`Table 1.` captions, running heads,
  appendix `A1` numbering). Swapping to the real `iopart.cls` later means
  changing `\documentclass` and deleting this file; sections stay.
- `sections/front.tex` — title block; the abstract is a compressed
  rendition of the manuscript's (CQG limit 300 words), claiming nothing
  the manuscript does not.
- `sections/references.tex` — hand-built IOP Vancouver bibliography,
  transcribed from `../citations/references.bib` (note fields omitted),
  entry order pinned to first-citation order. Re-derive the order from
  `main.aux` if a `\cite` is added or moved.
- `figures/make_journal_figures.py` — regenerates all seven PDFs from the
  same committed inputs the manuscript figures use (`../figures/data/*.csv`
  and `docs/prereg/*.json`); the two new result figures (capstone,
  Schwarzschild) draw the raw per-reading arrays stored inside the frozen
  artifacts, and the setup figure is an illustration (fixed seed, no
  measured quantity). Run from the repository root with the project venv.
- `check_numbers.py` — fidelity gate: every numeric token in the LaTeX
  must exist in `manuscript.md` and vice versa (normalized, string-exact).
  Run from anywhere: `python docs/paper/paper_a/latex/check_numbers.py`.

## Deliberate deviations from the manuscript (all register-level)

- Abstract compressed to journal length; body claims untouched.
- British spelling with `-ize` endings, per IOP house style.
- Three method citations added at first mention (`delong1988`,
  `garwood1936`, `clopperpearson1934`) — all from the verified
  bibliography, previously cited only in the claim-boundary section.
- Two new result figures (plane-wave capstone; Schwarzschild S4/S5) and a
  genre-standard setup illustration, none previously in the manuscript.
- The two adjacent "provenance is in Appendix B" sentences at the end of
  section 6.7 are merged into one.
- A CQG-required Data availability statement is added (repository is
  public; no DOI claimed).
