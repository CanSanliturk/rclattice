"""The digitized Fig. 10(b) and Table 4 (PLAN.md §10; harness lifted in D103).

Kept as a compatibility surface: the record itself is described once in `study_spec.REFERENCES`
(a `rclattice.study.references.ReferenceSet`), and the functions below delegate to it.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from study_spec import REFERENCES                 # noqa: E402

DATA = REFERENCES.path
MEASURED_F_KN = REFERENCES.measured_F_kN
MEASURED_K_KNMM = REFERENCES.measured_K_kNmm
CLIP_MM = REFERENCES.clip_mm
AYDIN = REFERENCES.author_by_horizon

load = REFERENCES.load
envelope_at = REFERENCES.envelope_at
compare = REFERENCES.compare
comparison_points = REFERENCES.comparison_points
