const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEMO_OUT_ROOT = get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos"))
using SpaceAGORA, Arrow, DataFrames, Statistics
dir = joinpath(DEMO_OUT_ROOT, "odyssey_ed")
println(readdir(dir))
df = DataFrame(Arrow.Table(joinpath(dir, "simulation_results.feather")))
lp = filter(n -> startswith(n, "sc1_link_pose"), names(df))
println("rows=", nrow(df), " link_pose cols=", length(lp))
for k in 1:14
    c = df[!, "sc1_link_pose_$k"]
    println("  col $k: min=", round(minimum(c); digits=4), " max=", round(maximum(c); digits=4), " changes=", count(!=(0.0), diff(c)))
end
alpha = df[!, "sc1_panel_alpha_rad_1"]
chg = findall(!=(0.0), diff(alpha))
println("panel alpha changes: n=", length(chg), " first t=", isempty(chg) ? NaN : df.time[chg[1]], " last t=", isempty(chg) ? NaN : df.time[chg[end]], " tEnd=", df.time[end])
println("alt min t=", df.time[argmin(df.sc1_altitude)], " alt min=", minimum(df.sc1_altitude))
scene = read_visualization_scene(joinpath(dir, "simulation_results_scene.json"))
println("planet=", scene.planet.name, " links=", length(scene.spacecraft[1].links), " html=", filesize(joinpath(dir, "simulation_results_viewer.html")))
