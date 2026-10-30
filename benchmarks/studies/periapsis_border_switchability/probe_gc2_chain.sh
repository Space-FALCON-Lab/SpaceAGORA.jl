#!/bin/bash
# Probe GC-2 chain on one host. Usage: probe_gc2_chain.sh <results-root> <host-stamp> <smoke|full|r4x>   (run inside the repo root; julia on PATH)
set -e
ROOT=$1; HOST=$2; MODE=$3; D=benchmarks/studies/periapsis_border_switchability
run() { systemd-run --user --scope -p MemoryMax=11G julia --project=. $D/probe_gc2_run.jl "$@" --results-root="$ROOT" --host="$HOST"; }
echo "host check:"; uptime; ps -eo pcpu,nlwp,args | grep '[j]ulia' | cut -c1-100
if [ "$MODE" = r4x ]; then run --phase=states --horizon-days=15; run --phase=r4 --workers=3; echo CHAIN_DONE; exit 0; fi
if [ "$MODE" = explore ]; then run --phase=r1v; run --phase=r4b --workers=3; run --phase=r4v --workers=3; uptime; run --phase=timedv; uptime; echo CHAIN_DONE; exit 0; fi
if [ "$MODE" = smoke ]; then run --phase=smoke; echo CHAIN_DONE; exit 0; fi
run --phase=states
run --phase=r1
run --phase=r4 --workers=3
run --phase=ens --workers=3
find "$ROOT/probe_gc2/runs" -name knots.csv -exec gzip -f {} \;
uptime; ps -eo pcpu,nlwp,args | grep '[j]ulia' | cut -c1-100
run --phase=timed
uptime
echo CHAIN_DONE
