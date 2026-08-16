# Reports

Model and analysis documentation, one directory per study.

```
doc/
  build.sh                 build one report, or all of them
  shared/                  preamble shared by every report
    rclattice-report.sty
  reports/
    column/                RC cantilever column study
      report.tex
      figures/             figures committed with the report
      make_model_figures.py
    kutay_wall/            Kutay's wall — SW-NC-FF (Sahinkaya et al. 2025)
      report.tex
    katrin_wall/           Katrin's wall — WSH3 (Dazio, Beyer & Bachmann 2009)
  out/                     built PDFs (not committed)
```

## Building

```bash
./build.sh                 # every report
./build.sh kutay_wall      # just one
```

Output: `out/<name>_report.pdf`. Auxiliary files are written to a temporary directory, so the
tree stays clean.

## Adding a report

1. `mkdir reports/<name>` and write `reports/<name>/report.tex`.
2. Start it with `\documentclass[11pt,a4paper]{article}` and `\usepackage{rclattice-report}` —
   `build.sh` puts `shared/` on `TEXINPUTS`, so the package resolves without a path.
3. Point `\graphicspath` at the analysis output, e.g.
   `\graphicspath{{../../../output/<study>/}{figures/}}` (relative to `reports/<name>/`).

The shared style provides `\unit{}` and three table environments — `spectable` (label / value /
note), `deftable` (label / definition) and `numtable{cols}` (small, centred, for numeric tables).
