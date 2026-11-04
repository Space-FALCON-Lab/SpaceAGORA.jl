#!/usr/bin/env bash
# Build the paper harness's precompile workload image: create the workload
# environment, then precompile SpaceAGORAPaperWorkload (which runs the workload)
# in a memory-capped scope, then confirm the image is current for the way the
# harness launches its workers.
#
# Usage: build_workload.sh [setup|precompile|check|all]      (default: all)
# Env:   SPACEAGORA_PPB_WORKLOAD_ENV         workload environment
#                                            (default: <repo>/output/paper_workload/env)
#        SPACEAGORA_PPB_WORKLOAD_MEMORY_MAX  MemoryMax of the precompile scope (default: 16G)
#        SPACEAGORA_PPB_WORKLOAD_SCOPE       0 runs without a systemd scope (no memory cap)
#        SPACEAGORA_PPB_WORKLOAD_LOG         log directory (default: dirname of the env)
#
# The precompile step also (re)builds SpaceAGORA's own image for the repository
# project when that is stale, which is the image the harness itself loads, so
# running this once before a benchmark run pays both compilations up front.
set -uo pipefail
STAGE="${1:-all}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
ENVDIR="${SPACEAGORA_PPB_WORKLOAD_ENV:-$REPO/output/paper_workload/env}"
LOGDIR="${SPACEAGORA_PPB_WORKLOAD_LOG:-$(dirname "$ENVDIR")}"
MEMMAX="${SPACEAGORA_PPB_WORKLOAD_MEMORY_MAX:-16G}"
JULIA_BIN="$(command -v julia || echo "$HOME/.juliaup/bin/julia")"
PKG=SpaceAGORAPaperWorkload
LOAD_PATH_STACK="@:${ENVDIR}:@v#.#:@stdlib"
mkdir -p "$LOGDIR"

# Run a command inside a memory-capped transient scope when systemd is there to
# provide one, and under /usr/bin/time -v for wall time and peak RSS.
capped() {
  local log="$1"; shift
  local timer=()
  [[ -x /usr/bin/time ]] && timer=(/usr/bin/time -v)
  if [[ "${SPACEAGORA_PPB_WORKLOAD_SCOPE:-1}" != "0" ]] && command -v systemd-run >/dev/null 2>&1 \
      && systemd-run --user --scope -q -- true >/dev/null 2>&1; then
    systemd-run --user --scope -q -p MemoryMax="$MEMMAX" -p MemorySwapMax=0 -- \
      "${timer[@]}" "$@" >> "$log" 2>&1
  else
    echo "build_workload: no systemd user scope available; running without a memory cap" | tee -a "$log"
    "${timer[@]}" "$@" >> "$log" 2>&1
  fi
}

stage_setup() {
  echo "== setup $ENVDIR" | tee -a "$LOGDIR/setup.log"
  capped "$LOGDIR/setup.log" "$JULIA_BIN" --startup-file=no "$HERE/setup_workload_env.jl" "$ENVDIR"
}

stage_precompile() {
  echo "== precompile $PKG (log: $LOGDIR/precompile.log)" | tee -a "$LOGDIR/precompile.log"
  (cd "$REPO" && capped "$LOGDIR/precompile.log" \
    env JULIA_LOAD_PATH="$LOAD_PATH_STACK" JULIA_NUM_PRECOMPILE_TASKS=1 \
        OPENBLAS_NUM_THREADS=1 GKSwstype=100 \
    "$JULIA_BIN" --project="$REPO" --startup-file=no \
      -e "t = time(); using $PKG; println(\"precompile_workload_loaded_s=\", round(time() - t; digits=1))")
  local rc=$?
  grep -E 'precompile_workload_loaded_s|Elapsed \(wall clock\)|Maximum resident set size|failed' "$LOGDIR/precompile.log" | tail -n 5
  return $rc
}

# The same probe the harness controller runs before it launches workers.
stage_check() {
  JULIA_LOAD_PATH="$LOAD_PATH_STACK" "$JULIA_BIN" --project="$REPO" --startup-file=no -e "
    id = Base.identify_package(\"$PKG\")
    ok = id !== nothing && Base.isprecompiled(id)
    println(ok ? \"$PKG image is current for $ENVDIR\" : \"$PKG image is NOT current for $ENVDIR\")
    exit(ok ? 0 : 3)"
}

case "$STAGE" in
  setup) stage_setup ;;
  precompile) stage_precompile ;;
  check) stage_check ;;
  all) stage_setup && stage_precompile && stage_check ;;
  *) echo "unknown stage $STAGE" >&2; exit 2 ;;
esac
