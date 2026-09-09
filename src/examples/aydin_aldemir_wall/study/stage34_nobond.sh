#!/bin/sh
# Stage 3 (static diagnostics) then Stage 4 (quasi-static, decisive) — the NOBOND half of the
# matrix: comp {lincomp, eppcomp} x tail {solved, paper}. The bonded half waits on PLAN §4.
PY=../../../.venv/bin/python
set -x
for comp in lincomp eppcomp; do
  for tail in solved paper; do
    $PY run.py --analysis static --comp $comp --tail $tail --drift 0.003 || true
  done
done
for comp in lincomp eppcomp; do
  for tail in solved paper; do
    $PY run.py --analysis pushover --comp $comp --tail $tail --drift 0.003 \
        --progress-every 10000 || true
  done
done
$PY master.py
