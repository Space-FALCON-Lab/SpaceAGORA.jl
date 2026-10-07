#!/bin/sh
# Run one case end to end: Basilisk on sf3 (accuracy run), SpaceAGORA locally, then compare.
#   sh run_case.sh <caseA|caseB|caseC> <duration_s> <local-julia-env>
# Preconditions: ssh alias sf3 is quiet; harness synced to ~/xval_mujoco/harness there (rsync, see README.md).
set -e
case=$1; T=$2; ENV=$3; here=$(cd "$(dirname "$0")" && pwd)
tag=${case}_T$T
# Basilisk replaces MuJoCo's warning handler with its logger, so warnings appear on stdout, not in MUJOCO_LOG.TXT.
bsklog=$(mktemp)
ssh sf3 "cd ~/xval_mujoco/harness && BSK_DIST=\$HOME/xval_mujoco/bsk_mj/dist3 BSK_PY=\$HOME/Documents/basilisk_2_12/.venv/bin/python sh run_bsk.sh cases/$case.toml ~/xval_mujoco/results/sf3/$tag $T" 2>&1 | tee "$bsklog"
if grep -qiE "MuJoCo internal warning|unstable|Traceback|WARNING" "$bsklog"; then echo "ABORT: Basilisk run reported a warning or instability (see above)"; exit 1; fi
mkdir -p "$here/results/sf3" "$here/results/$(hostname)"
rsync -a sf3:xval_mujoco/results/sf3/$tag "$here/results/sf3/"
# The SpaceAGORA driver aborts a level if MUJOCO_LOG.TXT grows in its working directory, if mjData.time diverges from
# the scene time (MuJoCo reset the state), or if a recorded value is non-finite.
flock /tmp/claude-1000/spaceagora-agents-julia.lock systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
  julia --project="$ENV" "$here/sagora_driver.jl" "$here/cases/$case.toml" "$here/results/$(hostname)/$tag" $T
python3 "$here/compare.py" "$here/cases/$case.toml" "$here/results/$(hostname)/$tag" "$here/results/sf3/$tag" "$here/results/summary_$tag"
