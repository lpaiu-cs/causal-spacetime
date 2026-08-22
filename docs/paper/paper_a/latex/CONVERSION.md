# Conversion brief: manuscript.md -> journal LaTeX (CQG target)

Source of truth: `docs/paper/paper_a/manuscript.md`. That file is
contract-tested against the committed artifacts
(`tests/test_paper_a_results_integration.py`,
`tests/test_paper_a_count_integration.py`), so THIS conversion must be a
faithful rendition of it. The LaTeX lives in
`docs/paper/paper_a/latex/sections/*.tex`, included by `main.tex`
(class: article + `cqgmimic.sty`, XeLaTeX).

## Hard rules (violations are defects)

1. **Every number is copied verbatim.** Never round, recompute, or
   reformat digits. Allowed representation change only:
   `3.579e-4` -> `$3.579\times10^{-4}$` (digits unchanged).
2. **No claim may strengthen or weaken.** Grade words are frozen
   vocabulary: CONFIRMED, DETECTED, CONCORDANT, POSITIVE, REPLICATED,
   confirmed, detected, "program-internal statement", "incomplete
   separation", "controlled validation". Keep each where it stands.
3. **No Korean in the paper.** Frozen sentences appear as English
   renderings only; their Korean originals live in the results artifacts
   and are bound to those renderings by `FROZEN_RENDERINGS` in
   `tests/test_paper_a_results_integration.py`. Appendix B states this
   once. Two tests enforce it — every rendering must still match its
   original and reach its section, and no Hangul may reach the
   manuscript. Do not reintroduce quoted originals into either file.
4. **Citations:** `[@key1; @key2]` -> `\cite{key1,key2}`; `\cite{key}`
   stays. Every key must exist in `../citations/references.bib`. Add no
   new citations.
5. **No new content.** Polish is register-level only (see below). Do not
   reorder results, merge claims, or drop hedges/bounds ("at the tested
   settings", "in the tested protocol" etc. all stay).

## Structure and labels

Sections map 1:1 to the manuscript's; `\section`/`\subsection`/
`\subsubsection` numbering is automatic. Use sentence case for headings
(iopart style), e.g. "The reconstruction ladder".

Fixed labels (use exactly these; cross-reference with `\sref{}` /
`\Sref{}` (sentence start) / `\fref{}` / `\Fref{}` / `\tref{}` /
`\Tref{}` defined in cqgmimic.sty):

| label | object |
| --- | --- |
| `sec:intro` | §1 Introduction |
| `sec:ladder` | §2 The reconstruction ladder |
| `sec:methods` | §3 Methods |
| `sec:results` | §4 Results |
| `sec:r0`..`sec:r5` (4.1-4.6) and `sec:rindler` (4.7) | §4 subsections |
| `sec:negative` | §5 Negative results |
| `sec:capstone` | §6 Capstone |
| `sec:ceiling`, `sec:pureweyl`, `sec:design`, `sec:preregresult`, `sec:establishes`, `sec:notestablish` | §6.1-6.6 |
| `sec:schwarzschild` | §6.7 |
| `sec:count` | §6.8 |
| `sec:oracle` | §6.9 |
| `sec:discussion` §7, `sec:claims` §8, `sec:repro` §9 |
| `app:conventions` App A, `app:provenance` App B |
| `fig:ladder`, `fig:convergence`, `fig:measure`, `fig:capstone`, `fig:schwarzschild`, `fig:count` | figures 1-6 |
| `tab:rungs` (§2), `tab:capstone` (§6.4), `tab:s4` and `tab:s5` (§6.7), `tab:count` (§6.8), `tab:claims` (§8), `tab:conventions` (App A) | tables |

Manuscript cross-references like "Section 4.6" become `\sref{sec:r5}`
etc. "Sections 6--6.7" -> `sections~\ref{sec:capstone}.1--\ref{sec:schwarzschild}`
is wrong — write it as `sections~\ref{sec:capstone} and~\ref{sec:schwarzschild}`
or keep the range in prose ("the capstone and its extensions,
sections~\ref{sec:capstone}--\ref{sec:count}") matching the meaning.
"Appendix A/B" -> `\ref{app:conventions}` / `\ref{app:provenance}` via
`appendix~\ref{...}`.

Figure numbering intent (order of appearance): 1 ladder schematic (§2),
2 convergence panel (§4.4), 3 measure dependence (§4.6),
4 capstone separation (§6.4, NEW), 5 Schwarzschild extension (§6.7, NEW),
6 mass-ladder count (§6.8). Graphics files:
`figures/fig1_ladder.pdf` ... `figures/fig6_ladder_count.pdf`
(generated separately; reference them with
`\includegraphics[width=...]{figures/figN_....pdf}` even if the PDF does
not exist yet). Figure environments: `[t]` placement, `\centering`,
caption BELOW for figures. Table captions ABOVE (caption package handles
spacing; put `\caption` before the tabular). Use booktabs
(`\toprule`/`\midrule`/`\bottomrule`).

## Typography conventions

- Math mode for every symbolic quantity: $N = 300$, $\rho$, $\tau$,
  $\mu = 2M/r_c$, $\beta$, $\varepsilon_\Delta$ (eps_delta),
  $\varepsilon_{\mathrm{det}}$, $\varepsilon_{\mathrm{rep}}$, $K$, $U_{\mathrm{amb}}$,
  $D$, $B$, $V$, $C$, $L \sim \sqrt{2\rho}\,\tau$,
  $\tau_{\mathrm{est}} = \sqrt{2K/\rho}$, $\det g = -1$, $f(d)$,
  $r \in [10, 20]$, $wT = 1$, AUC stays roman text.
- Intervals verbatim: `[-0.036211, -0.035953]` ->
  $[-0.036211, -0.035953]$.
- Dimension notation: write `1+1D` as `$(1{+}1)$D`, `2+1D` as
  `$(2{+}1)$D`, `3+1D` as `$(3{+}1)$D`; standalone "4D" -> `4D` (text).
- `->` in prose -> `$\to$` or reword; `x` as multiplication sign
  ("140x") -> `$140\times$`.
- Backticked file paths and artifact names -> `\artifact{...}` with `_`
  escaped (`\_`). Experiment ids (exp07 etc.) -> `\expid{exp07}`.
- Names: Myrheim--Meyer, Brightwell--Gregory, Clopper--Pearson,
  He--Rideout, Homšak--Veroni (typeset `Hom\v{s}ak--Veroni`),
  Berthiere--Gibbons--Solodukhin, Alfyorov--Shnyukov: en-dash `--`
  between distinct persons. Ranges: `sections 4--5`, `PR \#38--\#48`.
- Em dash: use `---` closed up ("word---word") sparingly, mirroring the
  manuscript's dashes.
- Quotation marks: ``...'' (curly). Percent: `2.5\%`. `#` -> `\#`,
  `&` -> `\&`.
- Petrov "type N"/"type D" as in source. "Poincare" -> `Poincar\'e`.
  "Alexandrov", "Rindler", "Schwarzschild", "Minkowski" unchanged.
- ASCII source only, EXCEPT em/en dashes in copied prose and
  `Hom\v{s}ak` handled via macro. Unicode −, ≤, ± in copied English text
  must become math (`$-$`, `$\le$`, `$\pm$`).

## Register polish (allowed and wanted)

- Resolve markdown artifacts: bullet lists may stay `itemize` where the
  manuscript uses them structurally (methods, negative results), but
  narrative lists inside paragraphs become prose.
- Smooth telegraphic openers into full sentences where needed
  ("Contributions: (1) ..." -> "The contributions are: (i) ...").
- **British spelling with `-ize` endings** (IOP house style): behaviour,
  neighbouring, colour; but normalize/organize/recognize keep `-ize`. An
  earlier version of this brief said American spelling, which was wrong
  and had to be undone.
- Keep the author's voice: precise, spare, occasionally pointed. Do not
  flatten deliberate constructions ("three ways out and no fourth",
  "by certification rather than by taste") — they are intentional.
- Sentence-per-claim structure must survive: a sentence that carries a
  gate/verdict keeps its content in one sentence.

## What NOT to convert

- Do not modify `../manuscript.md`, `../claim_boundary.md`, or anything
  under `../figures/` or `docs/prereg/`.
- Do not renumber or re-grade anything.
- Do not invent keywords, PACS, or acknowledgements beyond what the
  manuscript carries (front matter is handled separately).
