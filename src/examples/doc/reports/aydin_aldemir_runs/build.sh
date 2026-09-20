#!/usr/bin/env bash
# Build the Aldemir-wall run sheets.
#
#   ./build.sh              regenerate every run sheet, then typeset
#   ./build.sh --no-gen     typeset only (source already generated)
#   ./generate.py --only <substr> ...   regenerate a subset by hand
#
# Layout, deliberately flat:
#   src/      LaTeX source — report.tex, preamble.tex, runs/*.tex, figures/*.pdf (generated)
#   build/    every temporary LaTeX product (aux, log, toc, out)
#   ./        report.pdf and the two scripts, nothing else
set -euo pipefail
cd "$(dirname "$0")"
command -v pdflatex >/dev/null 2>&1 || { echo "error: pdflatex not found" >&2; exit 1; }

[ "${1:-}" = "--no-gen" ] || (cd ../../../.. && uv run python \
    examples/doc/reports/aydin_aldemir_runs/generate.py)

SHARED="$PWD/../../shared"
mkdir -p build
# twice, so the table of contents and cross-references resolve
for _ in 1 2; do
  TEXINPUTS="$SHARED:$PWD/src:" pdflatex -interaction=nonstopmode \
      -output-directory build -jobname report src/report.tex >/dev/null 2>&1 || true
done
if [ -f build/report.pdf ]; then
  cp build/report.pdf report.pdf
  miss=$(grep -c "not found" build/report.log 2>/dev/null || true)
  pages=$(grep -oE "Output written .*\(([0-9]+) page" build/report.log | grep -oE "[0-9]+ page" | head -1)
  echo "built report.pdf  ($pages${miss:+, WARNING: $miss missing input(s)})"
else
  echo "FAILED — see build/report.log" >&2; exit 1
fi
