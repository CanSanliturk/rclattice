"""The shared parametric-study harness (D103).

Built for the Aldemir wall (`examples/aydin_aldemir_wall/study/`, D83-D102) and lifted here so a
specimen CONFIGURES it instead of copying it. A specimen supplies a `StudySpec` — its parameter
registry, its `build(params)`, its node selectors, its measured references, its protocols and its
reporting text — and the harness supplies everything else: the CLI, the run-directory contract, the
per-run report, the metrics and rescoring, the master matrix, the advisor page and the LaTeX run
sheets.

What lives where:

* `registry.py`   the parameter registry — one `Param` per parameter, and everything (CLI, run
                  name, `params.json`, comparability) derived from it
* `spec.py`       the `StudySpec` contract a specimen fills in
* `runner.py`     `run.py`'s body: one run = one timestamped directory
* `metrics.py`    smoothed peak, plateau, drift capacity on a monotonic or a tip-envelope basis
* `protocols.py`  cyclic ladders and their cost
* `references.py` the digitized record and matched-displacement comparison
* `report.py`     the per-run `report.md`
* `master.py`     the matrix over every run, with `variant_note`
* `rescore.py`    re-score finished runs from their own stored series
* `page.py`       the advisor-facing HTML page and its JSON payload
* `runsheet.py`   the LaTeX run-sheet generator

Nothing here imports openseespy: the runners come from `rclattice.opensees`, which stays the only
module that does.
"""
