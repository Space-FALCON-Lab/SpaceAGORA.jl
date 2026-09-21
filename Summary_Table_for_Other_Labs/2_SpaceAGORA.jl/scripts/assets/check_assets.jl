repo_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(repo_root, "src", "SpaceAGORA.jl"))

report = SpaceAGORA.check_assets(; repo_root=repo_root)
SpaceAGORA.render_asset_report(report)
