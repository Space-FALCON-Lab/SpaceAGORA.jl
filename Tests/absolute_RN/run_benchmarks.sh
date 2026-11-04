#!/usr/bin/env bash
set -euo pipefail

test_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project="$(cd -- "$test_dir/../.." && pwd)"
julia_bin="${JULIA_BIN:-$(command -v julia)}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_PYTHONCALL_EXE="$project/9_SpaceAGORA.jl-main_full integration/.CondaPkg/.pixi/envs/default/bin/python"
test -x "$JULIA_PYTHONCALL_EXE"

# Reuse the installed Python runtime; do not provision a second environment
# merely because the frozen source snapshot lives at a different path.
while read -r kind label altitude inclination helpers policy orbits power; do
  for path in "$test_dir/PilotOutput/$label" "$test_dir/results/$label.toml" \
    "$test_dir/results/$label.log" "$test_dir/results/$label.time" "$test_dir/results/$label.verify.log"; do
    if [[ -e "$path" ]]; then
      echo "Refusing to overwrite benchmark evidence: $path" >&2
      exit 1
    fi
  done
  echo "Benchmark: $label"
  /usr/bin/time -v -o "$test_dir/results/$label.time" \
    "$julia_bin" --startup-file=no --threads=1 --gcthreads=1 \
    "$test_dir/pilot.jl" "$kind" "$label" "$altitude" "$inclination" \
    "$helpers" "$policy" "$orbits" "$power" benchmark \
    >"$test_dir/results/$label.log" 2>&1
  "$julia_bin" --startup-file=no --threads=1 --gcthreads=1 \
    "$test_dir/validate_pilots.jl" verify "$label" \
    >"$test_dir/results/$label.verify.log" 2>&1
  echo "Verified benchmark: $label"
done <<'CASES'
baseline baseline_sma800 1000 1 800 gve_sma 10 10000
new new_sma800 1000 1 800 gve_sma 10 10000
baseline baseline_inc800 1000 1 800 gve_inc 10 10000
new new_inc800 1000 1 800 gve_inc 10 10000
baseline baseline_zero 1000 1 100 gve_sma 0.1 0
CASES
