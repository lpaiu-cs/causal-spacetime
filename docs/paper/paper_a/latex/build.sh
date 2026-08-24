#!/bin/sh
# Build the CQG-style Paper A draft + its supplementary-material document.
# Run from this directory. main.tex must build first: si.tex resolves its
# cross-references into the article from main.aux (xr-hyper).
set -e
latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex
latexmk -xelatex -interaction=nonstopmode -halt-on-error si.tex
