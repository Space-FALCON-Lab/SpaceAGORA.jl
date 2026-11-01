#!/bin/bash
# Generic Julia job for Expanse. Submitted by `expanse.sh run`; resources come from
# sbatch flags, the command from $JOB_CMD. Runs inside the synced repo directory.
#SBATCH --job-name=spaceagora
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
echo "Job $SLURM_JOB_ID on $(hostname) | commit ${GIT_COMMIT:-unknown} | cpus ${SLURM_CPUS_PER_TASK:-1} | $(date)"

module purge
module load cpu gcc/10.2.0
export PATH="$HOME/julia/current/bin:$PATH"
export JULIA_DEPOT_PATH="$HOME/.julia"
# Memory and CPU are bounded by Slurm (cgroup); process-route workers read the CPU count.
export SPACEAGORA_PPB_PROCESS_WORKERS="${SLURM_CPUS_PER_TASK:-1}"

# Quiet-node gate (BENCH=1 only). The node is exclusive, but a previous tenant's processes can still be winding down.
# The gate passes on the INSTANTANEOUS reading: whole-node busy % over a short window, and no process of another user
# using CPU. It does not wait for the 1-minute load average, which lags minutes behind (it held jobs 55386812 and
# 55363544 for 146 s and 205 s while busy read 0.0-0.15% from the first sample); load1 is only logged. The gate polls
# until the node is quiet or the wait is used up, logs what it saw, and aborts the job (exit 75) rather than time
# anything on a busy node. The thresholds are ASSUMED values, not measured ones; override them in the environment of `expanse.sh run`.
QUIET_OTHER_PCPU_MAX="${QUIET_OTHER_PCPU_MAX:-5}"  # no process of another user (root excluded) above this %cpu
QUIET_BUSY_PCT_MAX="${QUIET_BUSY_PCT_MAX:-2.0}"    # whole-node CPU busy %, over a QUIET_SAMPLE_S window, must be below this
QUIET_SAMPLE_S="${QUIET_SAMPLE_S:-5}"              # length of the busy-% window
QUIET_WAIT_S="${QUIET_WAIT_S:-420}"                # give up after this many seconds
QUIET_POLL_S="${QUIET_POLL_S:-15}"                 # pause between checks

cpu_busy_pct() {   # whole-node busy % over QUIET_SAMPLE_S, from /proc/stat (instantaneous, unlike the load average)
  local a b
  a=($(head -1 /proc/stat)); sleep "$QUIET_SAMPLE_S"; b=($(head -1 /proc/stat))
  local ia=$(( a[4] + a[5] )) ib=$(( b[4] + b[5] )) ta=0 tb=0 i
  for i in 1 2 3 4 5 6 7 8; do ta=$(( ta + a[i] )); tb=$(( tb + b[i] )); done
  awk -v dt=$(( tb - ta )) -v di=$(( ib - ia )) 'BEGIN { v = (dt > 0) ? 100 * (dt - di) / dt : 0; printf "%.2f", v }'
}
log_node_state() {
  echo "== load: $(cat /proc/loadavg)"
  echo "== busiest processes (user pid %cpu %mem etime command):"
  ps -eo user:16,pid,pcpu,pmem,etime,args --sort=-pcpu | head -12 | cut -c1-160
}
other_users_procs() {   # processes of other users (root excluded) above QUIET_OTHER_PCPU_MAX %cpu; empty when none
  ps -eo user:20,pcpu,pid,args --sort=-pcpu | awk -v me="$(id -un)" -v m="$QUIET_OTHER_PCPU_MAX" 'NR>1 && $2+0>m && $1!=me && $1!="root"' | head -5 | cut -c1-160 || true
}
quiet_gate() {
  local t0=$SECONDS load1 busy others
  while :; do
    load1=$(cut -d' ' -f1 /proc/loadavg); busy=$(cpu_busy_pct); others=$(other_users_procs)
    if [[ -z "$others" ]] && awk -v b="$busy" -v bm="$QUIET_BUSY_PCT_MAX" 'BEGIN { exit !(b < bm) }'; then
      echo "== quiet gate passed after $((SECONDS - t0)) s: busy=${busy}% (< $QUIET_BUSY_PCT_MAX%), no other-user process above ${QUIET_OTHER_PCPU_MAX}% cpu; load1=$load1 (logged, not gated)"
      return 0
    fi
    if (( SECONDS - t0 >= QUIET_WAIT_S )); then
      echo "== NODE NOT QUIET after $((SECONDS - t0)) s: busy=${busy}% (limit $QUIET_BUSY_PCT_MAX%) load1=$load1 (not gated)"
      [[ -z "$others" ]] || printf '== other-user processes:\n%s\n' "$others"
      log_node_state
      echo "== aborting without running the command; resubmit, or raise QUIET_WAIT_S / the thresholds if the node is known good" >&2
      return 1
    fi
    echo "== waiting for a quiet node: busy=${busy}% other_users=$([[ -z "$others" ]] && echo none || echo present) load1=$load1 ($((SECONDS - t0)) s of $QUIET_WAIT_S s)"
    sleep "$QUIET_POLL_S"
  done
}

if [[ "${BENCH:-0}" == 1 ]]; then
  echo "== bench node: $(hostname) | $(lscpu | sed -n 's/^Model name: *//p') | $(nproc) cpus"
  echo "== load before: $(cat /proc/loadavg)"
  log_node_state
  quiet_gate || exit 75
fi
bash -c "$JOB_CMD"
if [[ "${BENCH:-0}" == 1 ]]; then echo "== load after: $(cat /proc/loadavg)"; fi
