#!/bin/bash
# Probe GC chain on one host: pilot -> frozen horizon rule -> untimed calls -> main. Usage: probe_gc_chain.sh <results-root> <host-stamp>   (run inside the repo root; julia on PATH)
set -e
ROOT=$1; HOST=$2; D=benchmarks/studies/periapsis_border_switchability
run() { systemd-run --user --scope -p MemoryMax=11G julia --project=. $D/probe_gc_run.jl "$@" --results-root="$ROOT" --host="$HOST"; }
uptime; ps -eo pcpu,nlwp,args | grep '[j]ulia' | cut -c1-120
run --phase=pilot
H=$(python3 -I - "$ROOT/probe_gc/summary_pilot_$HOST.csv" <<'PY'
import csv, math, sys
r = [x for x in csv.DictReader(open(sys.argv[1])) if x["run_id"] == "pilot_direct20_pert"][0]
W = float(r["run_wall_s"]) / 3.0          # wall seconds per day
print(30 if 30 * W <= 900 else max(10, math.floor(900 / W)))
PY
)
echo "HORIZON_DAYS=$H (frozen rule)"; echo "$H" > "$ROOT/probe_gc/horizon_days.txt"
run --phase=calls
run --phase=main --horizon-days="$H"
uptime
echo CHAIN_DONE
