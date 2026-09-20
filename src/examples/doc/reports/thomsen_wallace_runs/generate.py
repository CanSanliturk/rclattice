#!/usr/bin/env python
"""Generate one LaTeX run sheet, and one figure, per Thomsen & Wallace RW2 pushover or cyclic run.

    python generate.py                     every pushover and cyclic run that produced data
    python generate.py --only 2026-09-10   just the runs whose directory name contains this
    python generate.py --list              show what would be generated, write nothing

The generator is `rclattice.study.runsheet` (lifted in D103); the RW2 vocabulary — axis
descriptions, title fragments, record labels — is `study_spec.SPEC.runsheet`.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]                              # .../src
sys.path.insert(0, str(ROOT / "examples" / "thomsen_wallace_wall" / "study"))

from study_spec import SPEC                         # noqa: E402
from rclattice.study import runsheet                # noqa: E402

if __name__ == "__main__":
    raise SystemExit(runsheet.main(SPEC, HERE))
