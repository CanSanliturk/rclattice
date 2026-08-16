#!/usr/bin/env bash
# Build an rclattice report to PDF.
#
#   ./build.sh                 build every report under reports/
#   ./build.sh kutay_wall      build one
#
# Each report lives in reports/<name>/report.tex and may keep its own figures/.
# The shared preamble is shared/rclattice-report.sty; PDFs land in out/<name>_report.pdf.
set -euo pipefail
cd "$(dirname "$0")"
command -v pdflatex >/dev/null 2>&1 || { echo "error: pdflatex not found" >&2; exit 1; }
mkdir -p out

build_one() {
  local name=$1 dir="reports/$1"
  [ -f "$dir/report.tex" ] || { echo "error: no $dir/report.tex" >&2; return 1; }
  local tmp; tmp=$(mktemp -d)
  # Run twice so cross-references resolve; keep aux files out of the tree. A missing figure is
  # reported as a warning rather than a hard stop, so one study whose analysis output has not been
  # regenerated does not block the other reports.
  local rc=0
  for _ in 1 2; do
    TEXINPUTS="$PWD/shared:$PWD/$dir:" pdflatex -interaction=nonstopmode \
        -output-directory "$tmp" "$dir/report.tex" >/dev/null 2>&1 || rc=$?
  done
  local missing
  missing=$(grep -c "not found" "$tmp/report.log" 2>/dev/null || true)
  if [ -f "$tmp/report.pdf" ]; then
    mv "$tmp/report.pdf" "out/${name}_report.pdf"
    if [ "${missing:-0}" -gt 0 ]; then
      echo "built out/${name}_report.pdf  (WARNING: ${missing} figure(s) not found — run the study first)"
    else
      echo "built out/${name}_report.pdf"
    fi
  else
    echo "FAILED $name (rc=$rc); see $tmp/report.log" >&2
    return 1
  fi
  rm -rf "$tmp"
}

if [ $# -eq 0 ]; then
  found=0; failed=0
  for d in reports/*/; do
    # a directory without a report.tex is a placeholder for a study not yet written up;
    # skip it rather than failing the whole build
    if [ -f "$d/report.tex" ]; then
      build_one "$(basename "$d")" || failed=1
      found=1
    fi
  done
  [ "$found" -eq 1 ] || echo "no reports found under reports/"
  [ "$failed" -eq 0 ] || exit 1
else
  build_one "$1"
fi
