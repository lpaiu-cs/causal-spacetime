#!/bin/sh
# Build the CQG-style Paper A draft. Run from this directory.
set -e
latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex
