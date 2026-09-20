"""The matrix-aware master report (PLAN.md §10; harness lifted in D103) — `rclattice.study.master`."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from study_spec import SPEC                       # noqa: E402
from rclattice.study import master as shared      # noqa: E402

AXES = SPEC.axes
MEASURED = SPEC.measured_F_kN


def runs(root: Path) -> list[dict]:
    return shared.runs(SPEC, root)


def variant_note(params: dict) -> str:
    return shared.variant_note(SPEC, params)


def build(root: Path) -> str:
    return shared.build(SPEC, root)


if __name__ == "__main__":
    shared.main(SPEC)
