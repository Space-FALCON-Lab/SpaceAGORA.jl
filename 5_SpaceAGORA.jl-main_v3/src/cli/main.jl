const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"))

exit(SpaceAGORA.run_cli(copy(ARGS)))
