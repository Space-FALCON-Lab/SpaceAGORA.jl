repo_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(@__DIR__, "..", "..", "src", "SpaceAGORA.jl"))

SpaceAGORA.SpaceAGORACLI.setup_open_assets(; repo_root=repo_root)
