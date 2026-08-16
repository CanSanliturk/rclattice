#!/bin/sh
# Status + results of the wall cyclic analysis.  Usage: sh examples/wall/check.sh [logfile]
LOG=${1:-/tmp/cyc4.log}
if pgrep -f "cyclic.py --drift" >/dev/null 2>&1; then
    echo "STATUS: RUNNING"
    grep "step " "$LOG" | tail -1
    grep "step " "$LOG" | tail -1 | awk -F'[/ ]+' '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/ && $(i+1) ~ /^[0-9]+$/){d=$i;t=$(i+1);break}} END{if(t>0) printf "       %.1f%% complete\n", 100*d/t}'
else
    echo "STATUS: FINISHED"
fi
echo "convergence failures: $(grep -c 'failed to converge' "$LOG" 2>/dev/null)"
echo
grep -E "completed|drift reached|peak base|paper SW|residual|saved" "$LOG" 2>/dev/null
exit 0
