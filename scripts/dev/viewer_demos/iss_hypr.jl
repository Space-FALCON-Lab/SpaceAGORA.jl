# HyPR/PSO rendezvous path planning around the ISS: the Earth_RPO_CubeSat_MPC
# example with NASA's ISS display model sampled as the station point cloud,
# tracked by LQ-MPC and exported as a standalone viewer page with the planned
# path overlaid on the model.
#
#   julia --project=. scripts/dev/viewer_demos/iss_hypr.jl
#
# Environment:
#   SPACEAGORA_VIEWER_DEMO_OUT   output root (default: output/viewer_demos)
#   SPACEAGORA_DEMO_SMOKE=1      bounded run: a short hop beside the station,
#                                a small swarm and coarse sampling; it
#                                exercises the pipeline, not the approach
#
# Results go to <root>/iss_hypr_<inputs digest>/. A finished run is reused
# only when the provenance sidecar written beside it lists exactly the inputs
# of the current invocation and the plan it flew is on disk; the page is then
# rebuilt from that recorded plan, never from a new one, because the PSO is
# reproducible only on a single Julia thread. Otherwise the run is fresh.
# Nothing is deleted.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "Earth_RPO_CubeSat_MPC.jl"))
using Dates
using JSON
using Random
using SHA
using StaticArrays

const DEMO_OUT_ROOT = get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos"))
const DEMO_SMOKE = get(ENV, "SPACEAGORA_DEMO_SMOKE", "0") == "1"
const ISS_MODEL = joinpath(REPO_ROOT, "data", "models", "iss_nasa_3d_resources_b.glb")
const ISS_SCALE = 2.4                           # meters per model unit; truss near 109 m
const ISS_ROTATION_DEG = (-90.0, 0.0, -90.0)    # truss along body y, modules along body x
const STATION_ID = 201                          # the example's station spacecraft id
const STATION_INDEX = 2                         # its position in the run's spacecraft list

"""Every input that shapes the run, in one record; its digest names the output directory."""
function iss_hypr_inputs(; smoke::Bool=DEMO_SMOKE)
    isfile(ISS_MODEL) || error("ISS display model missing at $(ISS_MODEL)")
    common = (
        model=relpath(ISS_MODEL, REPO_ROOT),
        model_sha256=bytes2hex(open(sha256, ISS_MODEL)),
        model_scale=ISS_SCALE,
        model_rotation_deg=ISS_ROTATION_DEG,
        cloud_seed=7,
        station_keepout_radius_m=3.0,
        station_name="iss",
        station_dims_m=(73.0, 109.0, 20.0),
        station_mass_kg=420_000.0,
        safe_distance_m=2.0,
        cost_ref_distance_m=150.0,
        search_margin_m=80.0,
        seed=741,
        smoke=smoke,
    )
    bounded = smoke ?
        (n_points=2000, start_rtn=(-90.0, -40.0, 20.0), goal_rtn=(-75.0, -30.0, 15.0),
         mission_time=120.0, sample_ds_m=1.0, pso_n_particles=12, pso_n_iters=3, data_rate_s=5.0) :
        (n_points=12000, start_rtn=(-170.0, -100.0, 60.0), goal_rtn=(70.0, 20.0, 0.0),
         mission_time=900.0, sample_ds_m=0.5, pso_n_particles=80, pso_n_iters=25, data_rate_s=2.0)
    return merge(common, bounded)
end

_json_roundtrip(x) = JSON.parse(JSON.json(x))

iss_hypr_outdir(inputs) = joinpath(DEMO_OUT_ROOT, "iss_hypr_" * bytes2hex(sha256(JSON.json(inputs)))[1:10])

"""Sample the ISS model into the station point cloud and build the example scenario for it."""
function iss_hypr_build(inputs, outdir::AbstractString)
    points = sample_model_pointcloud(
        ISS_MODEL;
        n_points=inputs.n_points,
        rng=MersenneTwister(inputs.cloud_seed),
        scale=inputs.model_scale,
        rotation_deg=inputs.model_rotation_deg,
    )
    extents = [(minimum(points[k, :]), maximum(points[k, :])) for k in 1:3]
    demo = build_rpo_cubesat_mpc_demo(;
        mission_time=inputs.mission_time,
        results_directory=outdir,
        seed=inputs.seed,
        start_rtn=SVector{3, Float64}(inputs.start_rtn),
        goal_rtn=SVector{3, Float64}(inputs.goal_rtn),
        station_points=points,
        station_keepout_radius_m=inputs.station_keepout_radius_m,
        station_name=inputs.station_name,
        station_dims_m=inputs.station_dims_m,
        station_mass_kg=inputs.station_mass_kg,
        safe_distance_m=inputs.safe_distance_m,
        cost_ref_distance_m=inputs.cost_ref_distance_m,
        search_margin_m=inputs.search_margin_m,
        sample_ds_m=inputs.sample_ds_m,
        pso_n_particles=inputs.pso_n_particles,
        pso_n_iters=inputs.pso_n_iters,
        data_rate_s=inputs.data_rate_s,
        verbose=false,
    )
    return demo, extents
end

"""True when a finished run with exactly these inputs, its scene sidecar and its flown plan already sit in `outdir`."""
function iss_hypr_matching_run(sidecar::AbstractString, inputs, prefix::AbstractString, planfile::AbstractString)
    isfile(sidecar) && isfile(planfile) && isfile(prefix * "_scene.json") || return false
    recorded = JSON.parsefile(sidecar)
    return get(recorded, "inputs", nothing) == _json_roundtrip(inputs)
end

"""The flown plan as written by a fresh run: waypoints, retimed reference and cost."""
function iss_hypr_plan_record(demo)
    return Dict(
        "path_rtn" => Matrix{Float64}(demo.plan_result.path),
        "t_ref_s" => Vector{Float64}(demo.initial_plan.t_ref_s),
        "r_ref_rtn" => Matrix{Float64}(demo.initial_plan.r_ref_rtn),
        "v_ref_rtn" => Matrix{Float64}(demo.initial_plan.v_ref_rtn),
        "cost" => demo.plan_result.cost,
    )
end

"""Read a 3 x N matrix back from the JSON column list."""
_matrix3(columns) = Matrix{Float64}(reduce(hcat, [Float64.(c) for c in columns]))

_file_record(path) = isfile(path) ? (path=relpath(path, REPO_ROOT), bytes=filesize(path), sha256=bytes2hex(open(sha256, path))) : nothing

function main()
    inputs = iss_hypr_inputs()
    outdir = iss_hypr_outdir(inputs)
    prefix = joinpath(outdir, "simulation_results")
    sidecar = joinpath(outdir, "iss_hypr_provenance.json")
    planfile = joinpath(outdir, "iss_hypr_plan.json")
    mkpath(outdir)

    reused = iss_hypr_matching_run(sidecar, inputs, prefix, planfile)
    if reused
        recorded = JSON.parsefile(sidecar)
        saved = JSON.parsefile(planfile)
        path = _matrix3(saved["path_rtn"])
        ref = _matrix3(saved["r_ref_rtn"])
        station = recorded["station"]
        plan = recorded["plan"]
        extents = recorded["cloud_extents_m"]
        println("results with matching provenance present in ", outdir, "; reusing the recorded plan and simulation")
    else
        t_plan = @elapsed (demo, extents) = iss_hypr_build(inputs, outdir)
        path = Matrix{Float64}(demo.plan_result.path)
        ref = Matrix{Float64}(demo.initial_plan.r_ref_rtn)   # dense retimed reference, 3 x N
        clearance = SM.rpo_path_clearance_stats(ref, demo.geometry; safe_distance_m=inputs.safe_distance_m)
        plan = _json_roundtrip((
            waypoints=size(path, 2),
            reference_samples=size(ref, 2),
            cost=demo.plan_result.cost,
            t_ref_end_s=demo.initial_plan.t_ref_s[end],
            simulation_time_s=demo.args.mission_configuration.mission_time,
            min_clearance_m=clearance.min_clearance,
            violation_count=clearance.violation_count,
            planning_seconds=t_plan,
        ))
        station = _json_roundtrip(demo.station)
        println("station cloud extents (m): ", [(round(lo; digits=1), round(hi; digits=1)) for (lo, hi) in extents])
        println("plan: ", plan["waypoints"], " waypoints, ", plan["reference_samples"], " reference samples, cost=", round(plan["cost"]; digits=3),
            ", reference ends at ", round(plan["t_ref_end_s"]; digits=1), " s, simulated ", round(plan["simulation_time_s"]; digits=1),
            " s, minimum clearance ", round(plan["min_clearance_m"]; digits=2), " m, violations ", plan["violation_count"])
        open(planfile, "w") do io
            JSON.print(io, iss_hypr_plan_record(demo))
        end
        println("simulating into ", outdir)
        # The engine simulates a copy of the configuration whose LQ-MPC
        # controller shares the native OSQP workspace of `demo.control`; keep
        # the builder result alive for the whole run so that workspace is not
        # finalized underneath the copy.
        GC.@preserve demo run_simulation(demo.args; visualization=true)
    end

    html = export_visualization(
        prefix;
        max_frames=3000,
        trail_s=plan["simulation_time_s"],
        title="AGORA RPO · HyPR path to the ISS",
        models=Dict(STATION_ID => ISS_MODEL),
        model_scale=inputs.model_scale,
        model_rotation_deg=Dict(STATION_ID => inputs.model_rotation_deg),
        paths=[
            (name="HyPR/PSO reference (RTN)", points_m=ref, frame=:rtn, target=STATION_INDEX, color="#7fe0ff", dashed=true),
            (name="PSO waypoints", points_m=path, frame=:rtn, target=STATION_INDEX, color="#ffd166", dashed=false),
        ],
    )

    provenance = Dict(
        "inputs" => _json_roundtrip(inputs),
        "station" => station,
        "cloud_extents_m" => extents,
        "plan" => plan,
        "simulation_reused" => reused,
        "outputs" => Dict(
            "plan" => _file_record(planfile),
            "scene" => _file_record(prefix * "_scene.json"),
            "results_csv" => _file_record(prefix * ".csv"),
            "viewer_html" => _file_record(html),
        ),
        "software" => Dict("julia" => string(VERSION), "julia_threads" => Threads.nthreads(), "spaceagora" => string(pkgversion(SpaceAGORA))),
        "written_at" => Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") * "Z",
    )
    open(sidecar, "w") do io
        JSON.print(io, provenance, 2)
    end
    println("html: ", html, " (", filesize(html), " bytes)")
    println("provenance: ", sidecar)
    return provenance
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
