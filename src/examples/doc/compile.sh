#!/usr/bin/env bash
# Compile the column report to PDF.
#
# Usage:  ./compile.sh            (from examples/doc/)
# Output: column_report.pdf
#
# Runs pdflatex three times so the table of contents and cross-references resolve.
# Figures are read from ./figures/ and ../output/column/ (set in the .tex \graphicspath).

set -euo pipefail
cd "$(dirname "$0")"

DOC=column_report

if ! command -v pdflatex >/dev/null 2>&1; then
  echo "error: pdflatex not found (install a TeX distribution, e.g. MacTeX/TeX Live)." >&2
  exit 1
fi

pdflatex -interaction=nonstopmode -halt-on-error "$DOC.tex"
pdflatex -interaction=nonstopmode -halt-on-error "$DOC.tex"
pdflatex -interaction=nonstopmode -halt-on-error "$DOC.tex"

# tidy auxiliary files
rm -f "$DOC".{aux,log,out,toc}

echo "built $DOC.pdf"
