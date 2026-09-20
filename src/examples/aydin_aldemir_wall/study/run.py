"""Single entry point for every analysis in the study (PLAN.md §9; harness lifted in D103).

    run.py --analysis pushover --comp linear --tail solved --bond perfect --drift 0.003

Every parameter comes from `params.REGISTRY`, so the CLI, the run-directory name, `params.json` and
the master report never drift apart. One run = one timestamped directory containing `command.txt`,
`params.json` (EVERY parameter, defaults included, plus a schema version), a teed `console.log`,
`data.json`, `figures/` and a self-contained `report.md`. The body is `rclattice.study.runner`.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from study_spec import SPEC                       # noqa: E402
from rclattice.study import runner                # noqa: E402

if __name__ == "__main__":
    runner.main(SPEC)
