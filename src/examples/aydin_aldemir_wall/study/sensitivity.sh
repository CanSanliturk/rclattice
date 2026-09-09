#!/bin/sh
# Drift capacity as a FUNCTION of the assumed rupture strain, plus the mechanism split.
# Every cell is monotonic to 3.0% drift; the existing eps0.05/res0 run at 2.0% is the fifth point.
#   a  crushing only          -- does the wall fail without any bar fracture?
#   b  fracture only          -- does it fail with the crushing floor left at D22's 0.2?
#   c  eps_su 0.025, both     -- welded mesh, low ductility
#   d  eps_su 0.075, both     -- hot-rolled deformed bar
PY=../../../.venv/bin/python
set -x
$PY run.py --analysis pushover --comp c02 --tail solved --drift 0.03 --concrete-residual 0.0                        --progress-every 40000 || true
$PY run.py --analysis pushover --comp c02 --tail solved --drift 0.03 --steel-rupture 0.05                           --progress-every 40000 || true
$PY run.py --analysis pushover --comp c02 --tail solved --drift 0.03 --steel-rupture 0.025 --concrete-residual 0.0  --progress-every 40000 || true
$PY run.py --analysis pushover --comp c02 --tail solved --drift 0.03 --steel-rupture 0.075 --concrete-residual 0.0  --progress-every 40000 || true
$PY master.py
