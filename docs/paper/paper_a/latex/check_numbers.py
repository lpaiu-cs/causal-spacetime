"""Number-fidelity check: the LaTeX rendition against manuscript.md.

The manuscript is the contract-tested consumer of the artifacts; the
LaTeX draft is a rendition of it. This check holds the rendition to two
promises:

1. NO INVENTED NUMBERS: every numeric token in the LaTeX body appears in
   manuscript.md (after normalization), or is explicitly whitelisted as
   layout/markup.
2. NO DROPPED RESULTS: every numeric token in manuscript.md's body
   appears somewhere in the LaTeX.
3. COUNTS MATCH: each token occurs the same number of times on both
   sides, up to the justified per-token differences in COUNT_DELTAS
   (each bound to its reason; a stale entry is itself a defect). This
   closes the review-found gap where deleting one occurrence of a
   repeated number, or mis-transcribing a result into a number that
   already exists elsewhere, still passed the membership checks.
   (A compensating swap that keeps every count equal remains invisible
   to any counter; position-level truth is held by the manuscript-side
   claim contracts and the FROZEN_RENDERINGS LaTeX pins.)

Numbers are compared as normalized strings (scientific notation unified,
digit-group commas removed), never parsed as floats, so 0.10 and 0.1
remain distinct -- transcription must be verbatim.

Usage: python docs/paper/paper_a/latex/check_numbers.py   (exit 1 on defect)
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
SECTIONS = HERE / "sections"
MANUSCRIPT = HERE.parent / "manuscript.md"

# Layout values that legitimately exist only in the LaTeX markup.
LATEX_MARKUP_PATTERNS = [
    r"\\includegraphics\[[^]]*\]\{[^}]*\}",   # widths and file names
    r"\\begin\{tabular\*?\}\{[^}]*\}",        # column specs
    r"\\begin\{tabularx\}\{[^}]*\}\{[^}]*\}",
    r"\\(?:hspace|vspace|rule)\*?\{[^}]*\}",
    r"\\setlength\{[^}]*\}\{[^}]*\}",       # lengths (e.g. tabcolsep 2.5pt)
    r"\\\\\[[0-9.]+(?:pt|em|ex|pc)\]",         # line-break spacing (\\[2pt])
    r"\\multicolumn\{\d+\}",                # span counts are layout
    r"\\label\{[^}]*\}",
    r"\\ref\{[^}]*\}",
    r"\\cite\{[^}]*\}",
    r"\\(?:begin|end)\{[^}]*\}",
    r"\\bibliography(?:style)?\{[^}]*\}",
    r"^%.*$",                                  # comments
    r"(?<!\\)%.*$",                            # trailing comments
    r"\\thebibliography\{[^}]*\}",
    r"p\{[0-9.]+(?:cm|in|pc|pt|em)\}",
    r"m\{[0-9.]+(?:cm|in|pc|pt|em)\}",
    r"[0-9.]+\\(?:textwidth|linewidth|columnwidth)",
]

# Tokens that appear in the LaTeX for legitimate presentational reasons
# and are not manuscript numbers. Keep this list SHORT and justified.
TEX_WHITELIST: dict[str, str] = {
    # (none yet)
}

# Per-token occurrence differences (LaTeX count minus manuscript count)
# that are DELIBERATE, each bound to its reason. Anything not listed must
# occur the same number of times on both sides; a listed token whose
# actual difference deviates from the declared one is a defect either
# way, so entries cannot go stale silently.
CAPTION_RESTATED = ("restated once in the caption of a LaTeX-only figure "
                    "(fig0 setup / fig4 capstone / fig5 Schwarzschild; "
                    "deliberate deviations listed in README.md)")
CONTRIB_ROMAN = ("the manuscript numbers its contribution list (1)-(6); "
                 "the LaTeX renders the labels as roman (i)-(vi) per the "
                 "conversion brief")
COUNT_DELTAS: dict[str, tuple[int, str]] = {
    "-0.036211": (+1, CAPTION_RESTATED),
    "-0.035953": (+1, CAPTION_RESTATED),
    "0.0036":    (+1, CAPTION_RESTATED),
    "0.9734":    (+1, CAPTION_RESTATED),
    "0.9630":    (+1, CAPTION_RESTATED),
    "0.9837":    (+1, CAPTION_RESTATED),
    "0.999232":  (+1, CAPTION_RESTATED),
    "0.0501046": (+1, CAPTION_RESTATED),
    "0.0502929": (+1, CAPTION_RESTATED),
    "0.0504812": (+1, CAPTION_RESTATED),
    "3.579e-4":  (+1, CAPTION_RESTATED),
    "3000":      (+1, CAPTION_RESTATED),
    "4800":      (+1, CAPTION_RESTATED),
    "140":       (+1, CAPTION_RESTATED),
    "95":        (+1, CAPTION_RESTATED),
    "1.0":       (+1, CAPTION_RESTATED),
    "1":         (+1, "the fig0 setup caption's two occurrences minus the "
                      "manuscript's (1) contribution label"),
    "2":         (-1, CONTRIB_ROMAN),
    "3":         (-1, CONTRIB_ROMAN),
    "4":         (-1, CONTRIB_ROMAN),
    "5":         (-1, CONTRIB_ROMAN),
    "6":         (-1, CONTRIB_ROMAN),
}


def _normalize_sci(text: str) -> str:
    # 3.579\times10^{-4} / 3.579 \times 10^{-4}  ->  3.579e-4
    text = re.sub(
        r"(\d+(?:\.\d+)?)\s*\\times\s*10\^\{?(-?\d+)\}?",
        lambda m: f"{m.group(1)}e{int(m.group(2))}",
        text,
    )
    # markdown scientific notation stays as-is (3.579e-4 -> 3.579e-4)
    text = re.sub(
        r"(\d+(?:\.\d+)?)[eE](-?\d+)",
        lambda m: f"{m.group(1)}e{int(m.group(2))}",
        text,
    )
    return text


NUM = re.compile(
    r"(?<![A-Za-z0-9_])"        # not inside an identifier (exp07, p14, S4)
    r"-?"
    r"\d{1,3}(?:,\d{3})+"       # digit-grouped integers (26,831,117)
    r"|(?<![A-Za-z0-9_.])-?\d+\.\d+(?:e-?\d+)?"  # decimals / sci
    r"|(?<![A-Za-z0-9_.])\d+e-?\d+"              # bare sci (1e-12)
    # integers: reject a digit or decimal continuation, but ACCEPT a
    # sentence-final period ("c = 1." must count its 1 -- the earlier
    # (?![0-9.]) lookahead silently skipped every integer that ended a
    # sentence, asymmetrically between prose and math)
    r"|(?<![A-Za-z0-9_.#])\d+(?!\d)(?!\.\d)"
)


def _tokens(text: str) -> Counter:
    text = _normalize_sci(text)
    out: Counter = Counter()
    for m in NUM.finditer(text):
        tok = m.group(0).replace(",", "")
        out[tok] += 1
    return out


def _strip_latex(text: str) -> str:
    # \texorpdfstring{tex}{pdf} prints only its first argument; counting
    # both would double every numeral in it.
    text = re.sub(
        r"\\texorpdfstring\{((?:[^{}]|\{[^{}]*\})*)\}\{(?:[^{}]|\{[^{}]*\})*\}",
        r"\1", text)
    for pat in LATEX_MARKUP_PATTERNS:
        text = re.sub(pat, " ", text, flags=re.MULTILINE)
    # en-dash markup between numbers (0.93--1.07) is a range, not a sign
    text = text.replace("--", "–")
    return text


def _strip_markdown(text: str) -> str:
    # image links and heading numerals are markdown artifacts
    text = re.sub(r"!\[[^]]*\]\([^)]*\)", " ", text)
    text = re.sub(r"^#{1,6} [0-9.]+ ", "# ", text, flags=re.MULTILINE)
    text = re.sub(r"^#{1,6} ", "# ", text, flags=re.MULTILINE)
    # +- is the markdown spelling of \pm: the following number is unsigned
    text = text.replace("+-", "±")
    # literal cross-references become \ref{} in the LaTeX (auto-numbered),
    # so their numerals are markdown-only; result numbers never appear in
    # this "Section(s)/Figure N" shape.
    text = re.sub(
        r"(?:Sections?|Figures?|Table)\s+[0-9]+(?:\.[0-9]+)*"
        r"(?:\s*(?:-|and|,|through)\s*[0-9]+(?:\.[0-9]+)*)*",
        " ", text)
    text = re.sub(r"Section [0-9]+(?:\.[0-9]+)*'s", " ", text)
    return text


def main() -> int:
    # references.tex is excluded: its numbers (years, volumes, pages) are
    # bibliographic facts whose source of truth is ../citations/references.bib,
    # not manuscript.md.
    tex_text = "\n".join(
        p.read_text(encoding="utf-8") for p in sorted(SECTIONS.glob("*.tex"))
        if p.name != "references.tex"
    )
    tex_tokens = _tokens(_strip_latex(tex_text))
    md_text = MANUSCRIPT.read_text(encoding="utf-8")
    md_tokens = _tokens(_strip_markdown(md_text))

    defects = 0

    # The frozen sentences are printed as English renderings only; their
    # Korean originals live in the results artifacts, bound to those
    # renderings by FROZEN_RENDERINGS in
    # tests/test_paper_a_results_integration.py (which guards the same
    # property on the manuscript side). Appendix B says so once.
    # scanned for Hangul beyond the section files: the bibliography and
    # the assembly file are excluded from the NUMBER comparison (their
    # numbers answer to references.bib), but no rendition file may print
    # Korean.
    hangul_text = tex_text + "\n".join(
        p.read_text(encoding="utf-8")
        for p in [SECTIONS / "references.tex", HERE / "main.tex",
                  HERE / "si.tex"]
        if p.exists())
    hangul = re.findall(r"[가-힣]+", hangul_text)
    if hangul:
        print("== Korean in the LaTeX rendition ==")
        for run in dict.fromkeys(hangul[:8]):
            print(f"  {run}")
        defects += len(set(hangul))

    invented = {t: c for t, c in tex_tokens.items()
                if t not in md_tokens and t not in TEX_WHITELIST}
    if invented:
        print("== numeric tokens in LaTeX not found in manuscript.md ==")
        stripped = _strip_latex(tex_text)
        norm = _normalize_sci(stripped)
        for tok in sorted(invented):
            lines = [ln.strip() for ln in norm.splitlines()
                     if tok in ln.replace(",", "")]
            print(f"  {tok}  (x{invented[tok]})")
            for ln in lines[:2]:
                print(f"      | {ln[:110]}")
        defects += len(invented)

    miscounted = {
        t: (tex_tokens[t], md_tokens[t])
        for t in set(tex_tokens) | set(md_tokens)
        if t in md_tokens and t in tex_tokens
        and tex_tokens[t] - md_tokens[t] != COUNT_DELTAS.get(t, (0,))[0]
    }
    stale = {t: d for t, (d, _r) in COUNT_DELTAS.items()
             if tex_tokens[t] - md_tokens[t] != d}
    if miscounted:
        print("== tokens whose occurrence counts disagree "
              "(tex, manuscript; allowed delta from COUNT_DELTAS) ==")
        for tok in sorted(miscounted):
            tex_n, md_n = miscounted[tok]
            allowed = COUNT_DELTAS.get(tok, (0, "not listed"))
            print(f"  {tok}  tex x{tex_n}, md x{md_n}, "
                  f"allowed delta {allowed[0]:+d} ({allowed[1][:60]})")
        defects += len(miscounted)
    if stale:
        # every stale entry is already in miscounted unless the token
        # vanished from both sides entirely; report those too
        gone = {t for t in stale if t not in miscounted}
        for tok in sorted(gone):
            print(f"== COUNT_DELTAS entry stale: {tok} no longer occurs ==")
            defects += 1

    dropped = {t: c for t, c in md_tokens.items() if t not in tex_tokens}
    if dropped:
        print("== numeric tokens in manuscript.md not found in LaTeX ==")
        norm = _normalize_sci(_strip_markdown(md_text))
        for tok in sorted(dropped):
            lines = [ln.strip() for ln in norm.splitlines()
                     if tok in ln.replace(",", "")]
            print(f"  {tok}  (x{dropped[tok]})")
            for ln in lines[:2]:
                print(f"      | {ln[:110]}")
        defects += len(dropped)

    if defects == 0:
        print(f"OK: {sum(tex_tokens.values())} numeric tokens in LaTeX, "
              f"{sum(md_tokens.values())} in manuscript; no inventions, "
              f"no droppages, counts reconciled "
              f"({len(COUNT_DELTAS)} justified deltas).")
        return 0
    print(f"\n{defects} distinct token defects.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
