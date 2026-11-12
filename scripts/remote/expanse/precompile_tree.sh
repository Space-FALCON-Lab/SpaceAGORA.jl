#!/usr/bin/env bash
# Untimed: precompile this source tree on the `shared` partition so the exclusive benchmark node does not pay for it.
# Run once per source tree, after the sync and before the --bench job (do not change the tree in between):
#   EXPANSE_ACCOUNT=... ./expanse.sh run --cpus 8 --mem 24G --time 00:45:00 -- bash scripts/remote/expanse/precompile_tree.sh
# Both lines, as the benchmark script's own prelude (. is the children's env, examples the controller's).
set -eu
cd "$(dirname "$0")/../../.." || exit 1
# examples/Manifest.toml is gitignored: seed it from the root manifest and instantiate it when it is missing.
[ -f examples/Manifest.toml ] || julia --project=. examples/setup.jl
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using SpaceAGORA'
julia --project=examples -e 'using SpaceAGORA'
# GRAM extension (SpaceAGORAGRAMSuiteExt): GRAMSuite is a weak dependency, so it cannot be imported from the root project
# directly; load it the way the harness does (cases.jl ppc_ensure_gramsuite_loaded!: vendored environment first on LOAD_PATH).
julia --project=. -e 'using SpaceAGORA; SpaceAGORA.RuntimeServices._prepare_gramsuite_load_path!(joinpath(pwd(), "data", "GRAMSuite.jl"); prepend=true); @eval import GRAMSuite; ext = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing; println("GRAMSuite loaded, extension attached: ", ext); ext || exit(4)'
