"""Render the published master-report page from the template (RW2, D105).

    python master.py && python report_page.py
    # then republish the emitted file to the same artifact URL
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from study_spec import SPEC                       # noqa: E402
from rclattice.study import page                  # noqa: E402

TEMPLATE = SPEC.page.template
OUT = SPEC.page.out_html

if __name__ == "__main__":
    page.render(SPEC)
