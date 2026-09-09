#!/bin/sh
# 1. eppcomp/solved re-run PAST its target: its Stage 4 peak landed at 0.2997% of a 0.3000% target,
#    i.e. still ascending, so that peak is a lower bound and not a capacity (the D78 trap).
# 2. The bonded half, SOLVED TAIL ONLY. The nobond half measured `paper` as 5.26 x Gf — a mesh
#    artefact of transplanting his 20 mm fit — so bonded paper-tail cells would cost ~2.4 h to
#    vary a knob already known to be fictitious. Stage 4's stated job is to narrow before spending.
PY=../../../.venv/bin/python
set -x
$PY run.py --analysis pushover --comp eppcomp --tail solved --drift 0.005 --progress-every 20000 || true
$PY run.py --analysis pushover --comp lincomp --tail solved --bond bond-a06 --drift 0.003 --progress-every 20000 || true
$PY run.py --analysis pushover --comp eppcomp --tail solved --bond bond-a06 --drift 0.003 --progress-every 20000 || true
$PY master.py
