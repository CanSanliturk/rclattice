"""Build the JSON the published master-report page renders itself from (harness lifted in D103).

    python master.py && python report_payload.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from study_spec import SPEC                       # noqa: E402
from rclattice.study import page                  # noqa: E402

OUT = SPEC.page.out_payload


def build(root: Path) -> dict:
    return page.payload(SPEC, root)


if __name__ == "__main__":
    page.write_payload(SPEC)
