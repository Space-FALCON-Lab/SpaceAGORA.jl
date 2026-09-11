repo_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(repo_root, "src", "SpaceAGORA.jl"))

SpaceAGORA.SpaceAGORACLI.setup_open_assets(; repo_root=repo_root)
