repo_root = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(@__DIR__, "..", "..", "src", "SpaceAGORA.jl"))

report = SpaceAGORA.check_assets(; repo_root=repo_root)
SpaceAGORA.render_asset_report(report)
