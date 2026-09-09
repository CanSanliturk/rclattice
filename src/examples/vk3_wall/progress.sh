#!/bin/sh
# Progress of the newest running/finished VK3 analysis.
#   sh examples/vk3_wall/progress.sh          # one line
#   sh examples/vk3_wall/progress.sh -w       # refresh every 10 s
#
# Looks in both places a console log can live: /tmp while a run is in flight (the shell wrapper
# moves it into the run directory only on exit), and the run directory once it has landed.
OUT="$(dirname "$0")/../output/vk3_wall"
show() {
  LOG=$(ls -t /tmp/vk3_*.log "$OUT"/runs/*/console.log 2>/dev/null | head -1)
  [ -z "$LOG" ] && { echo "no VK3 log found"; return; }
  grep -E '^ +step' "$LOG" | tail -1 | awk -v log="$LOG" '{
    split($2,a,"/"); pct=a[1]/a[2]*100
    el=$NF; gsub(/[^0-9]/,"",el)
    eta=(a[2]-a[1])*el/a[1]/60
    printf "%5.1f%%  step %s/%s  drift %s  shear %s kN  elapsed %ds  ETA ~%.0f min\n",
           pct, a[1], a[2], $4, $6, el, eta
  }'
  pgrep -f "vk3_wall/(preflight|cyclic|pushover|elastic)\.py" >/dev/null \
    && echo "        RUNNING   $(basename "$(dirname "$LOG")")" \
    || echo "        not running (log: $LOG)"
}
if [ "$1" = "-w" ]; then while :; do show; sleep 10; done; else show; fi
