"""Assemble the journal submission package from note.tex, Fig1.pdf and cover_letter.md.

Usage: build_submission.py OUTDIR [TECTONIC]
Writes manuscript.pdf (compiled from the packaged source, so the archive is
known to build on its own), manuscript_source.zip, abstract.txt and, when
pandoc is installed, cover_letter.pdf. Nothing here changes the manuscript.
"""

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
    (r"\\%", "%"),
    (r"--", "-"),
    (r"\\[a-zA-Z]+|[{}~]", " "),
    (r"\s+", " "),
)


def abstract_text(tex):
    text = tex.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]
    for pattern, replacement in PLAIN:
        text = re.sub(pattern, replacement, text)
    return text.strip()


def main(out, tectonic=TECTONIC):
    out = Path(out).resolve()
    source = out / "source"
    source.mkdir(parents=True, exist_ok=True)
    tex = (HERE / "note.tex").read_text(encoding="utf-8")
    if any(ord(c) > 127 for c in tex):
        raise ValueError("The submission server uses pdfLaTeX: keep note.tex ASCII")
    (source / "main.tex").write_text(tex, encoding="utf-8", newline="\n")
    shutil.copyfile(HERE / "Fig1.pdf", source / "Fig1.pdf")
    subprocess.run(
        [str(tectonic), "--outdir", str(out), "main.tex"], cwd=source, check=True
    )
    (out / "main.pdf").replace(out / "manuscript.pdf")
    with zipfile.ZipFile(
        out / "manuscript_source.zip", "w", zipfile.ZIP_DEFLATED
    ) as bundle:
        for name in ("main.tex", "Fig1.pdf"):
            bundle.write(source / name, name)
    abstract = abstract_text(tex)
    words = len(abstract.split())
    if not 150 <= words <= 250:
        raise ValueError(f"Abstract has {words} words; the journal allows 150-250")
    (out / "abstract.txt").write_text(abstract + "\n", encoding="utf-8", newline="\n")
    shutil.copyfile(HERE / "cover_letter.md", out / "cover_letter.md")
    if shutil.which("pandoc"):
        subprocess.run(
            [
                "pandoc",
                str(HERE / "cover_letter.md"),
                "-o",
                str(out / "cover_letter.pdf"),
                f"--pdf-engine={tectonic}",
                *("-V", "papersize=a4", "-V", "geometry:margin=25mm"),
                *("-V", "fontsize=11pt", "-V", "pagestyle=empty"),
            ],
            check=True,
        )
    print(f"abstract words: {words}; package written to {out}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
