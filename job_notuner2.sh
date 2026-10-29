#!/bin/bash
# Arms: a=predictive_notuner+SWEEP_NS_FLOOR=1, b=predictive_notuner, c=inner_only, d=predictive_nocalib_notuner.
# Fresh process per (round, case, arm); 1 warm-up + 10 repeats; fresh store each; launch order rotated by round.
export OPENBLAS_NUM_THREADS=1 GKSwstype=100 SPACEAGORA_PPB_NO_PLOTS=1
J=$HOME/.juliaup/bin/julia
R=$HOME/spaceagora_remote
cd $R/pol_notuner
export JULIA_DEPOT_PATH=$R/shared_depot
O=results_notuner2; mkdir -p $O
echo "=== host $(date)"; uptime; ps -eo pcpu,nlwp,args | grep [j]ulia
BUDGET=16
export JULIA_NUM_THREADS="" SPACEAGORA_PERF_PROCS="" SPACEAGORA_PPB_WORKLOAD=0
export SPACEAGORA_CORE_BUDGET=$BUDGET SPACEAGORA_PPC_PHYSICAL_CORES=$BUDGET SPACEAGORA_PPC_EQUAL_CORE_BUDGET=1
export SPACEAGORA_PPB_SHARED_PROCESS=1 SPACEAGORA_PPB_WARMUP_SPAN=0.1 SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=auto
export SPACEAGORA_CAMPAIGN_CORRECTIONS=off
ORD[1]="a b c d"; ORD[2]="b c d a"; ORD[3]="c d a b"
declare -A MODE=([a]=predictive_notuner [b]=predictive_notuner [c]=inner_only [d]=predictive_nocalib_notuner)
for i in 1 2 3; do for C in gravity_4096sat_l50_vacuum_5800s gravity_1024sat_l50_vacuum_24600s; do for arm in ${ORD[$i]}; do
  m=${MODE[$arm]}
  D=$O/${C}/${arm}_$i; rm -rf $D; mkdir -p $D
  unset SPACEAGORA_RHS_REGION_WIDTH SPACEAGORA_RHS_CALIBRATE SPACEAGORA_RHS_CALIBRATE_SEED_TOPK SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO SPACEAGORA_RHS_CALIBRATE_SEED_ACCEPT_WINNER SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR SPACEAGORA_RHS_CALIBRATE_PIN_MARGIN
  [ $arm = a ] && export SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR=1
  export SPACEAGORA_RHS_CALIBRATION_PATH=$R/pol_notuner/$D/store.toml SPACEAGORA_PPB_SNAPSHOT_DIR=$R/pol_notuner/$D/snap
  env | grep '^SPACEAGORA_RHS' | sort > $D/env_before_launch.txt
  echo "=== $C arm=$arm mode=$m round $i $(date)"; uptime
  systemd-run --user --scope -p MemoryMax=16G $J --project=examples benchmarks/studies/parallelization_performance.jl full \
    --outdir=$D/ppb --cases=$C --modes=$m --threads=$BUDGET --process-workers=$BUDGET --repeats=10 --warmup=1 --parity-cases=none > $D/ppb.log 2>&1
  echo "EXIT $?"
done; done; done
unset SPACEAGORA_RHS_CALIBRATE_SWEEP_NS_FLOOR
echo "=== DONE $(date)"
