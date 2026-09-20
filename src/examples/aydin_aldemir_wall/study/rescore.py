"""Rescore finished runs from their own stored series (D92; harness lifted in D103).

    python rescore.py [--force] [--dry-run] [<dir> ...]
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from study_spec import SPEC                       # noqa: E402
from rclattice.study import rescore as shared     # noqa: E402

ROOT = SPEC.out_root

if __name__ == "__main__":
    raise SystemExit(shared.main(SPEC))
