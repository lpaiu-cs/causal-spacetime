# Paper A submission set (CQG target)

The pieces the journal receives, and how each is produced. The LaTeX
article (`../latex/main.pdf`) remains the typeset reference; the Word
files exist because the submission plan calls for editable manuscripts.

| File | What it is | Produced by |
| --- | --- | --- |
| `paper_a_article.docx` | The article, editable (native Word equations, 300 dpi figures, literal cross-reference and citation numbers) | `make_docx.py` |
| `paper_a_supplementary.docx` | Supplementary material S1-S2 (= manuscript appendices A-B) | `make_docx.py` |
| `cover_letter.md` / `.docx` | One-page editor-facing cover letter | hand-written; docx via `make_docx.py` |
| `submission_form_notes.md` | Internal article-type and referee notes; do not upload as the cover letter | hand-written |

`make_docx.py` flattens `../latex/sections/*.tex` into pandoc-friendly
TeX: cross-references and citation numbers are read from `main.aux` /
`si.aux` (never guessed, assert-guarded), captions get their printed
"Figure N."/"Table N." prefixes, custom macros (`\artifact`, `\expid`,
`\sref`...) are resolved, table column specs are simplified, and figure
PDFs are swapped for the 300 dpi PNGs in `../latex/figures/png300/`.
Run it AFTER a clean `../latex/build.sh` (it consumes the aux files):

```sh
cd docs/paper/paper_a/latex && ./build.sh && cd ../submission && python make_docx.py
```

Word-rendition caveats (all cosmetic, none numeric): commit hashes
inside math (e.g. `P' = 51875a2`) render in math italics; table column
widths are Word rather than LaTeX picas. The 11-column mass-ladder table is
split into two stacked Word tables (design/counts, then interval comparison)
for portrait-page legibility; its values are parsed from the LaTeX table and
assert-checked rather than retyped. Every number is still governed by
`../latex/check_numbers.py` at the LaTeX layer, which this conversion consumes
verbatim.

Submission itself (uploading, naming referees on the form) is performed
by the author.
