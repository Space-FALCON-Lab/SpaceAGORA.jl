repo_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(repo_root, "src", "SpaceAGORA.jl"))

entries = SpaceAGORA.SpaceAGORACLI.load_asset_manifest(; repo_root=repo_root)
SpaceAGORA.SpaceAGORACLI.render_asset_manifest(entries)
