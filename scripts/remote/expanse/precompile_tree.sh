#!/usr/bin/env bash
# Untimed: precompile this source tree on the `shared` partition so the exclusive benchmark node does not pay for it.
# Run once per source tree, after the sync and before the --bench job (do not change the tree in between):
#   EXPANSE_ACCOUNT=... ./expanse.sh run --cpus 8 --mem 24G --time 00:45:00 -- bash scripts/remote/expanse/precompile_tree.sh
# Both lines, as the benchmark script's own prelude (. is the children's env, examples the controller's).
set -eu
cd "$(dirname "$0")/../../.." || exit 1
julia --project=. -e 'using SpaceAGORA'
julia --project=examples -e 'using SpaceAGORA'
