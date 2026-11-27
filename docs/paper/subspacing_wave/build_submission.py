"""Assemble the submission package for the subspacing gravitational-wave manuscript.

Usage: build_submission.py OUTDIR [TECTONIC]
Compiles manuscript.pdf inside an isolated source/ directory from note.tex and
committed Fig*.pdf files, packages manuscript_source.zip, verifies the ASCII-only
single-file TeX contract and abstract word count, writes abstract.txt, and builds
cover_letter.pdf via pandoc if available.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
TECTONIC = Path.home() / "AppData/Local/Programs/tectonic/tectonic.exe"
PLAIN = (
    (r"\$([^$]*)\$", r"\1"),
    (r"\\times", " x "),
    (r"\\le", "<="),
    (r"\\ge", ">="),
    (r"\\%", "%"),
    (r"--", "-"),
    (r"\\[a-zA-Z]+|[{}~]", " "),
    (r"\s+", " "),
)


def abstract_text(tex: str) -> str:
    text = tex.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]
    for pattern, replacement in PLAIN:
        text = re.sub(pattern, replacement, text)
    return text.strip()


def main(out: str, tectonic: str | Path = TECTONIC) -> None:
    out_path = Path(out).resolve()
    source = out_path / "source"
    source.mkdir(parents=True, exist_ok=True)
    tex = (HERE / "note.tex").read_text(encoding="utf-8")
    if "\\input{" in tex or "\\include{" in tex:
        raise ValueError("note.tex must be a single source file")
    if any(ord(c) > 127 for c in tex):
        raise ValueError("note.tex must be pure ASCII for pdfLaTeX compatibility")
    (source / "main.tex").write_text(tex, encoding="utf-8", newline="\n")

    figs = sorted(HERE.glob("Fig*.pdf"))
    if not figs:
        raise FileNotFoundError(
            "No compiled Fig*.pdf assets found in manuscript directory"
        )
    for fig in figs:
        shutil.copyfile(fig, source / fig.name)

    subprocess.run(
        [str(tectonic), "--outdir", str(out_path), "main.tex"], cwd=source, check=True
    )
    (out_path / "main.pdf").replace(out_path / "manuscript.pdf")
    with zipfile.ZipFile(
        out_path / "manuscript_source.zip", "w", zipfile.ZIP_DEFLATED
    ) as bundle:
        bundle.write(source / "main.tex", "main.tex")
        for fig in figs:
            bundle.write(source / fig.name, fig.name)

    abstract = abstract_text(tex)
    words = len(abstract.split())
    if not 150 <= words <= 250:
        raise ValueError(f"Abstract has {words} words; must stay within 150-250")
    (out_path / "abstract.txt").write_text(
        abstract + "\n", encoding="utf-8", newline="\n"
    )
    if (HERE / "cover_letter.md").exists():
        shutil.copyfile(HERE / "cover_letter.md", out_path / "cover_letter.md")
        if shutil.which("pandoc"):
            subprocess.run(
                [
                    "pandoc",
                    str(HERE / "cover_letter.md"),
                    "-o",
                    str(out_path / "cover_letter.pdf"),
                    f"--pdf-engine={tectonic}",
                    "-V",
                    "papersize=a4",
                    "-V",
                    "geometry:margin=25mm",
                    "-V",
                    "fontsize=11pt",
                    "-V",
                    "pagestyle=empty",
                ],
                check=True,
            )
    print(
        f"abstract words: {words}; figures: {len(figs)}; package written to {out_path}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:3])
