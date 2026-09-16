# HyPR/PSO rendezvous path planning around the ISS: the Earth_RPO_CubeSat_MPC
# demo with NASA's ISS model as the station geometry, tracked by LQ-MPC.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEMO_OUT_ROOT = get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos"))
include(joinpath(REPO_ROOT, "examples", "Earth_RPO_CubeSat_MPC.jl"))
using Random, StaticArrays, JSON
T = DEMO_OUT_ROOT
outdir = joinpath(T, "iss_hypr_out")
iss = joinpath(REPO_ROOT, "data", "models", "iss_nasa_3d_resources_b.glb")
points = sample_model_pointcloud(iss; n_points=12000, rng=MersenneTwister(7), scale=2.4, rotation_deg=(-90, 0, -90))
println("station cloud extents: ", [(round(minimum(points[k, :]); digits=1), round(maximum(points[k, :]); digits=1)) for k in 1:3])
t_plan = @elapsed demo = build_rpo_cubesat_mpc_demo(;
    mission_time=900.0, results_directory=outdir, seed=741,
    start_rtn=SVector{3, Float64}(-170.0, -100.0, 60.0), goal_rtn=SVector{3, Float64}(70.0, 20.0, 0.0),
    station_points=points, station_keepout_radius_m=3.0, station_name="iss", station_dims_m=(73.0, 109.0, 20.0), station_mass_kg=420_000.0,
    safe_distance_m=2.0, cost_ref_distance_m=150.0, search_margin_m=80.0, sample_ds_m=0.5,
    pso_n_particles=80, pso_n_iters=25, data_rate_s=2.0, verbose=false)
path = demo.plan_result.path
ref = Matrix{Float64}(demo.initial_plan.r_ref_rtn)   # dense retimed reference, 3 x N
println("plan: waypoints ", size(path), " reference ", size(ref), " cost=", round(demo.plan_result.cost; digits=3), " t_ref_end=", demo.initial_plan.t_ref_s[end], " mission=", demo.args.mission_configuration.mission_time, " planned in ", round(t_plan; digits=1), " s")
if !isfile(joinpath(outdir, "simulation_results_scene.json"))
    @time run_simulation(demo.args; visualization=true)
else
    println("results present; skipping the simulation")
end
html = export_visualization(joinpath(outdir, "simulation_results"); max_frames=3000, trail_orbits=Inf,
    title="AGORA RPO · HyPR path to the ISS", models=Dict(201 => iss), model_scale=2.4, model_rotation_deg=Dict(201 => (-90, 0, -90)),
    paths=[(name="HyPR/PSO reference (RTN)", points_m=ref, frame=:rtn, target=2, color="#7fe0ff", dashed=true),
           (name="PSO waypoints", points_m=Matrix{Float64}(path), frame=:rtn, target=2, color="#ffd166", dashed=false)])
println("html: ", html, " ", filesize(html))
