"""Flatten the Paper A LaTeX into pandoc-friendly TeX and build the Word
submission set (article, supplementary material, cover letter).

Deterministic, assert-guarded: every substitution that must happen a known
number of times is counted. Cross-references and citation numbers are read
from main.aux / si.aux, never guessed.
"""
import re
import subprocess
import sys
from pathlib import Path

LATEX = Path(__file__).resolve().parents[1] / 'latex'
SUB = LATEX.parent / 'submission'
BUILD = SUB / 'build'
BUILD.mkdir(exist_ok=True)


def aux_labels(aux_path):
    text = aux_path.read_text(encoding='utf-8')
    out = {}
    for m in re.finditer(r'\\newlabel\{([^}]*)\}\{\{([^{}]*)\}', text):
        out[m.group(1)] = m.group(2)
    return out


def aux_citations(aux_path):
    text = aux_path.read_text(encoding='utf-8')
    order, seen = [], set()
    for m in re.finditer(r'\\citation\{([^}]*)\}', text):
        for k in m.group(1).split(','):
            k = k.strip()
            if k and k not in seen:
                seen.add(k)
                order.append(k)
    return {k: i + 1 for i, k in enumerate(order)}


LABELS = aux_labels(LATEX / 'main.aux')
LABELS_SI = aux_labels(LATEX / 'si.aux')
CITES = aux_citations(LATEX / 'main.aux')
assert len(CITES) == 31, len(CITES)


def clean_colspec(spec):
    spec = re.sub(r'@\{[^{}]*\}', '', spec)
    spec = re.sub(r'>\{[^{}]*\}', '', spec)
    spec = re.sub(r'[pm]\{[^{}]*\}', 'l', spec)
    spec = spec.replace('X', 'l')
    out = ''.join(c for c in spec if c in 'lcr')
    assert out, spec
    return out


def transform(text, labels):
    text = text.replace('\\allowbreak{}', '')
    text = text.replace('\\artifact{', '\\texttt{')
    text = text.replace('\\expid{', '\\texttt{')

    def label_num(key):
        assert key in labels, key
        return labels[key]

    for macro, word in [('sref', 'section'), ('Sref', 'Section'),
                        ('fref', 'figure'), ('Fref', 'Figure'),
                        ('tref', 'table'), ('Tref', 'Table')]:
        text = re.sub(r'\\' + macro + r'\{([^}]*)\}',
                      lambda m, w=word: f'{w}~{label_num(m.group(1))}', text)
    text = re.sub(r'\\ref\{([^}]*)\}', lambda m: label_num(m.group(1)), text)

    def cite_repl(m):
        keys = [k.strip() for k in m.group(2).split(',')]
        for k in keys:
            assert k in CITES, k
        sep = ' ' if m.group(1) else ''
        return sep + '[' + ', '.join(str(CITES[k]) for k in keys) + ']'
    text = re.sub(r'(~)?\\cite\{([^}]*)\}', cite_repl, text)

    # caption numbering: prepend the printed "Figure N."/"Table N." using the
    # label that follows the caption inside the same float.
    def float_repl(m):
        block = m.group(0)
        lab = re.search(r'\\label\{((?:fig|tab):[^}]*)\}', block)
        if lab:
            num = label_num(lab.group(1))
            kind = 'Figure' if lab.group(1).startswith('fig:') else 'Table'
            block = block.replace('\\caption{', f'\\caption{{{kind} {num}. ', 1)
        return block
    text = re.sub(r'\\begin\{(figure|table)\}.*?\\end\{\1\}',
                  float_repl, text, flags=re.DOTALL)
    text = re.sub(r'\\label\{[^}]*\}', '', text)

    SPEC = r'((?:[^{}]|\{[^{}]*\})*)'
    def mc_repl(m):
        keep = ''.join(c for c in re.sub(r'@\{[^{}]*\}', '', m.group(2))
                       if c in 'lcr')
        assert keep, m.group(2)
        return '\\multicolumn{' + m.group(1) + '}{' + keep + '}'
    text = re.sub(r'\\multicolumn\{(\d+)\}\{' + SPEC + r'\}', mc_repl, text)
    text = re.sub(r'\\begin\{tabularx\}\{\\textwidth\}\{' + SPEC + r'\}',
                  lambda m: ('\\begin{tabular}{'
                             + clean_colspec(m.group(1)) + '}'),
                  text)
    text = text.replace('\\end{tabularx}', '\\end{tabular}')
    text = re.sub(r'\\begin\{tabular\}\{' + SPEC + r'\}',
                  lambda m: ('\\begin{tabular}{'
                             + clean_colspec(m.group(1)) + '}'),
                  text)
    text = re.sub(r'\\includegraphics(\[[^]]*\])?\{figures/([^}]*)\.pdf\}',
                  r'\\includegraphics\1{figures/png300/\2.png}', text)
    text = re.sub(r'^\\shorttitle\{[^}]*\}\s*$', '', text, flags=re.MULTILINE)
    return text


# ---------------------------------------------------------------- article
front = (LATEX / 'sections' / 'front.tex').read_text(encoding='utf-8')
abstract = re.search(r'\\begin\{paperabstract\}\n(.*?)\\end\{paperabstract\}',
                     front, re.DOTALL).group(1).strip()
keywords = re.search(r'\\paperkeywords\{(.*?)\}', front, re.DOTALL).group(1)
keywords = ' '.join(keywords.split())

order = ['s01_intro', 's02_ladder', 's03_methods', 's04_results',
         's05_negative', 's06_capstone', 's067_schwarzschild', 's068_count',
         's069_oracle', 's07_discussion', 's08_claims', 's10_repro',
         'backmatter']
body = '\n\n'.join(
    (LATEX / 'sections' / f'{name}.tex').read_text(encoding='utf-8')
    for name in order)
refs = (LATEX / 'sections' / 'references.tex').read_text(encoding='utf-8')
refs = refs.replace('\\begin{thebibliography}{99}', '\\section*{References}')
refs = refs.replace('\\end{thebibliography}', '')
n_bib = 0


def bib_repl(m):
    global n_bib
    n_bib += 1
    key = m.group(1)
    assert CITES[key] == n_bib, (key, n_bib)
    return f'[{n_bib}]'


refs = re.sub(r'\\bibitem\{([^}]*)\}', bib_repl, refs)
assert n_bib == 31

TITLE = ('Spacetime quantities from causal order: an operational ladder '
         'with preregistered validations at finite density')
article = f"""\\title{{{TITLE}}}
\\author{{Juneyoung Kim \\\\ Independent researcher \\\\ lpaiu.cs@gmail.com}}
\\date{{}}
\\maketitle

\\section*{{Abstract}}
{abstract}

\\noindent\\textbf{{Keywords:}} {keywords}

\\noindent Submitted to: \\emph{{Classical and Quantum Gravity}}

{body}

{refs}
"""
article = transform(article, LABELS)
(BUILD / 'article_flat.tex').write_text(article, encoding='utf-8')

# ----------------------------------------------------------------- SI doc
si_body = '\n\n'.join(
    (LATEX / 'sections' / f'{name}.tex').read_text(encoding='utf-8')
    for name in ['appendix_a', 'appendix_b'])
# S-numbered headings become literal (pandoc numbering would restart at 1)
si_body = si_body.replace(
    '\\section{Conventions and normalizations}',
    '\\section*{S1. Conventions and normalizations}', 1)
si_body = si_body.replace(
    '\\section{Capstone execution provenance}',
    '\\section*{S2. Capstone execution provenance}', 1)
si_labels = dict(LABELS)
si_labels.update(LABELS_SI)
si = f"""\\title{{Supplementary material for “{TITLE}”}}
\\author{{Juneyoung Kim \\\\ Independent researcher \\\\ lpaiu.cs@gmail.com}}
\\date{{}}
\\maketitle

\\noindent This document carries the article's conventions table
(section~S1) and the execution-provenance record of the preregistered
stages (section~S2). Section numbers without an S prefix refer to the
article.

{si_body}
"""
si = transform(si, si_labels)
(BUILD / 'si_flat.tex').write_text(si, encoding='utf-8')

# ------------------------------------------------------------------ build
SUB.mkdir(exist_ok=True)


def pandoc(args):
    r = subprocess.run(['pandoc'] + args, cwd=LATEX, capture_output=True,
                       text=True)
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)
    assert r.returncode == 0, r.returncode


pandoc(['-f', 'latex', '-t', 'docx', '--number-sections',
        str(BUILD / 'article_flat.tex'),
        '-o', str(SUB / 'paper_a_article.docx')])
pandoc(['-f', 'latex', '-t', 'docx',
        str(BUILD / 'si_flat.tex'),
        '-o', str(SUB / 'paper_a_supplementary.docx')])
pandoc(['-f', 'markdown', '-t', 'docx',
        str(SUB / 'cover_letter.md'),
        '-o', str(SUB / 'cover_letter.docx')])
print('docx set written to', SUB)
