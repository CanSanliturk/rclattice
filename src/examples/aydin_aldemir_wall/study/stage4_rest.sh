#!/bin/sh
# Restart of the Stage 4 nobond cells the kill left undone. The four Stage 3 static cells and
# pushover lincomp-solved already completed and are NOT repeated.
PY=../../../.venv/bin/python
set -x
$PY run.py --analysis pushover --comp lincomp --tail paper  --drift 0.003 --progress-every 10000 || true
$PY run.py --analysis pushover --comp eppcomp --tail solved --drift 0.003 --progress-every 10000 || true
$PY run.py --analysis pushover --comp eppcomp --tail paper  --drift 0.003 --progress-every 10000 || true
$PY master.py
