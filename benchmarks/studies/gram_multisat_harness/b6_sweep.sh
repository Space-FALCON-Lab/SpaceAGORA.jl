#!/bin/bash
# Exercise the B6 GRAM cases across B6's mode ladder at low/high thread budget.
# B6 phase definition: modes = serial, outer_threads, outer_process,
# outer_inner_adaptive, full_smart; thread_mode = :low_high.
SCRATCH="$1"; CAP="${2:-900}"
REPO=/home/space-falcon-1/Documents/SpaceAGORA.jl
cd "$REPO"
OUT="$SCRATCH/b6_sweep.csv"
echo "case,mode,threads,success,retcode,wall_s,terminal_s,pos_norm_m,vel_norm_mps,peak_rss_mb" > "$OUT"
for CASE in multi_16_gram_live montecarlo_mars_gram_live; do
for MODE in serial outer_threads outer_process outer_inner_adaptive full_smart; do
for THR in 1 24; do
  F="$SCRATCH/sw_${CASE}_${MODE}_${THR}"
  timeout $CAP julia --project=. -t $THR benchmarks/studies/parallelization_performance.jl full \
    --worker --case="$CASE" --mode="$MODE" --thread-count=$THR \
    --repeat=1 --outfile="$F.json" > "$F.log" 2>&1 &
  PID=$!
  PEAK=0
  while kill -0 $PID 2>/dev/null; do
    R=$(awk '/VmRSS/{print $2}' /proc/$PID/status 2>/dev/null)
    R=$(( ${R:-0} / 1024 )); [ $R -gt $PEAK ] && PEAK=$R
    sleep 3
  done
  wait $PID; RC=$?
  if [ -f "$F.json" ]; then
    python3 - "$F.json" "$CASE" "$MODE" "$THR" "$PEAK" >> "$OUT" <<'PY'
import csv,sys
p,case,mode,thr,peak = sys.argv[1:6]
try:
    r=list(csv.DictReader(open(p)))[0]
    print(f"{case},{mode},{thr},{r['success']},{r['retcode']},{float(r['wall_time_s']):.2f},"
          f"{r['terminal_time_s']},{r['final_primary_pos_norm_m']},{r['final_primary_vel_norm_mps']},{peak}")
except Exception as e:
    print(f"{case},{mode},{thr},PARSE_FAIL,{e},,,,,{peak}")
PY
  else
    echo "$CASE,$MODE,$THR,false,NO_OUTPUT_rc=$RC,,,,,$PEAK" >> "$OUT"
  fi
done; done; done
echo "SWEEP_DONE" >> "$OUT"
