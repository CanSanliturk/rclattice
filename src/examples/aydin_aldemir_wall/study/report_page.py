"""Render the published master-report page from the template and the current run directories.

The page is an advisor-facing deliverable on a PUBLICLY SHARED artifact, and a public artifact
cannot be granted the `db` capability — so it cannot fetch anything at load time and instead
carries its own data. That makes rebuilding it a step in the study's own workflow rather than a
one-off piece of hand-written HTML:

    python master.py && python report_page.py
    # then republish the emitted file to the same artifact URL

Everything numeric comes from `report_payload.build`, which reads the runs; the template holds only
prose and layout. Nothing here is hand-maintained, so a run that lands tomorrow appears without the
page being edited.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import references                                                          # noqa: E402
import report_payload                                                      # noqa: E402

TEMPLATE = HERE / "report_page.template.html"
OUT = HERE / "report_page.html"
H = 2250.0            # drift denominator, mm (PLAN §2)


def reference_data() -> dict:
    """The digitized experiment — fixed data, so it is baked in and never regenerated."""
    d = references.load()
    if d is None:
        return {"cloud": [], "env": []}
    cloud = np.asarray(d["cloud"], float)
    mx, my = np.asarray(d["mono_x"], float), np.asarray(d["mono_y"], float)
    order = np.argsort(mx)
    return {"cloud": [[round(float(a / H * 100), 4), round(float(b), 1)] for a, b in cloud[::5]],
            "env": [[round(float(mx[k] / H * 100), 4), round(float(my[k]), 1)] for k in order]}


def main() -> None:
    import os
    os.environ.setdefault("ALDEMIR_TW", "210.0")
    import specimen
    payload = report_payload.build(specimen.OUT / "study")
    html = (TEMPLATE.read_text(encoding="utf8")
            .replace("__PAYLOAD__", json.dumps(payload))
            .replace("__REFDATA__", json.dumps(reference_data())))
    OUT.write_text(html, encoding="utf8")
    print(f"{OUT}  {OUT.stat().st_size:,} bytes")
    print(f"  {payload['run_count']} runs, {len(payload['rows'])} rows, "
          f"{len(payload['series'])} curves, generated {payload['generated']}")


if __name__ == "__main__":
    main()
