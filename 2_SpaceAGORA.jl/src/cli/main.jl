const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(@__DIR__, "..", "SpaceAGORA.jl"))

exit(SpaceAGORA.run_cli(copy(ARGS)))
