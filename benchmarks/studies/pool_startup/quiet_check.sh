#!/usr/bin/env bash
# Is the remote benchmark host quiet? One line per check; exit 0 when it passes.
#   bash quiet_check.sh <ssh-alias>
# Pass: no julia process for any user, 1-min load < 3, no open remote job (a
# job.meta without END whose PID is alive), no other process above 50% CPU.
ssh "$1" 'set -u
julia=$(pgrep -c julia); load=$(cut -d" " -f1 /proc/loadavg)
open=0
for m in ~/spaceagora_remote/jobs/*/job.meta; do
  grep -q "^END=" "$m" && continue
  p=$(grep "^PID=" "$m" | cut -d= -f2)
  [ -n "$p" ] && kill -0 "$p" 2>/dev/null && open=$((open + 1))
done
heavy=$(ps -eo user:20,pcpu,comm --no-headers | awk "\$2 > 50 {printf \"%s:%s:%s \", \$1, \$3, \$2}")
pass=1
[ "$julia" -eq 0 ] || pass=0
awk "BEGIN{exit !($load < 3)}" || pass=0
[ "$open" -eq 0 ] || pass=0
[ -z "$heavy" ] || pass=0
echo "$(date -Iseconds) pass=$pass julia=$julia load1=$load open_jobs=$open heavy=[${heavy% }]"
[ "$pass" -eq 1 ]'
