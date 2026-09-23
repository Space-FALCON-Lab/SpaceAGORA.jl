# HyPR relocation around the ISS: the Earth_RPO_CubeSat_MPC example with
# NASA's ISS display model in a flight-like attitude held in LVLH (truss along
# N, pressurized modules along T, smallest extent along R), planned by HyPR in
# its `:manuscript` mode, retimed with an acceleration limit, tracked by LQ-MPC
# and exported as a standalone viewer page with the plan overlaid on the model.
#
#   julia --project=. scripts/dev/viewer_demos/iss_hypr.jl
#
# Environment:
#   SPACEAGORA_VIEWER_DEMO_OUT       output root (default: output/viewer_demos)
#   SPACEAGORA_DEMO_SMOKE=1          bounded run: a short hop beside the
#                                    station, a small swarm and a coarse
#                                    cloud; it exercises the pipeline, not
#                                    the relocation
#   SPACEAGORA_ISS_HYPR_PLAN_ONLY=1  plan and write the plan record without
#                                    simulating, for comparing plans across
#                                    thread counts
#
# The scenario relocates the chaser between two V-bar hold points, from 100 m
# behind the station's aft end to 30 m ahead of its forward end, at rest at
# both; the straight segment between them crosses the station, so the planner
# has to route around it. The display model also holds seven disc-shaped
# meshes, about 37 m beyond the arrays and apart from the station, that are
# not ISS hardware. The run writes a copy of the model without them and uses
# that copy for the point cloud, the hold points and the viewer.
#
# Results go to <root>/iss_hypr_<inputs digest>/. A finished run is reused
# only when the provenance sidecar written beside it lists exactly the inputs
# of the current invocation and the plan, scene, Feather results and control
# log match their recorded hashes; the page is then rebuilt from that recorded
# plan, never from a new one. With fixed iteration budgets and no wall-clock
# limits, seeded planning is reproducible across thread counts. Otherwise the
# run is fresh. Nothing is deleted.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "Earth_RPO_CubeSat_MPC.jl"))
using Dates
using JSON
using Random
using SHA
using StaticArrays

const DEMO_OUT_ROOT = get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos"))
const DEMO_SMOKE = get(ENV, "SPACEAGORA_DEMO_SMOKE", "0") == "1"
const PLAN_ONLY = get(ENV, "SPACEAGORA_ISS_HYPR_PLAN_ONLY", "0") == "1"
const ISS_MODEL = joinpath(REPO_ROOT, "data", "models", "iss_nasa_3d_resources_b.glb")
const ISS_SCALE = 2.4                           # meters per model unit; truss near 109 m
# The model's own axes run along the station's height (x), the pressurized
# modules (y) and the truss (z), so the identity holds the truss along N and
# the modules along T, with the smallest extent along R.
const ISS_ROTATION_DEG = (0.0, 0.0, 0.0)
# Disc-shaped meshes about 37 m beyond the arrays, apart from the station.
const ISS_EXCLUDED_NODES = ("bendedtru1", "bendedtru2", "bendedtru3", "bendedtrus", "pCylinder1", "pCylinder2", "pCylinder8")
const STATION_ID = 201                          # the example's station spacecraft id
const STATION_INDEX = 2                         # its position in the run's spacecraft list

"""
Write a copy of a GLB without the named nodes: they are dropped from their
parents and from the scene roots, and the binary chunk is kept as is, so
every other node draws and samples exactly as in the source.
"""
function glb_without_nodes(src::AbstractString, dst::AbstractString, names)
    bytes = read(src)
    String(bytes[1:4]) == "glTF" || error("$(src) is not a binary glTF file")
    json_length = Int(reinterpret(UInt32, bytes[13:16])[1])
    String(bytes[17:20]) == "JSON" || error("$(src): the first chunk is not JSON")
    doc = JSON.parse(String(bytes[21:(20 + json_length)]))
    rest = bytes[(21 + json_length):end]
    drop = Set(i - 1 for (i, node) in enumerate(doc["nodes"]) if get(node, "name", "") in names)
    length(drop) == length(names) || error("$(src): expected the nodes $(names), found $(length(drop)) of them")
    for node in doc["nodes"]
        haskey(node, "children") && (node["children"] = [c for c in node["children"] if !(c in drop)])
    end
    for scene in doc["scenes"]
        scene["nodes"] = [c for c in scene["nodes"] if !(c in drop)]
    end
    chunk = Vector{UInt8}(JSON.json(doc))
    while length(chunk) % 4 != 0
        push!(chunk, 0x20)
    end
    mkpath(dirname(dst))
    open(dst, "w") do io
        write(io, b"glTF", UInt32(2), UInt32(12 + 8 + length(chunk) + length(rest)))
        write(io, UInt32(length(chunk)), b"JSON", chunk, rest)
    end
    return dst
end

"""The station model the run uses: the display model without the excluded nodes, written under `root`."""
function iss_station_model(root::AbstractString=joinpath(DEMO_OUT_ROOT, "models"))
    source = bytes2hex(open(sha256, ISS_MODEL))
    return glb_without_nodes(ISS_MODEL, joinpath(root, "iss_station_model_$(source[1:10]).glb"), ISS_EXCLUDED_NODES)
end

"""Half extents and extents of the station model (scaled, rotated and centred as the viewer draws it), per body axis."""
function iss_model_extents(model_path::AbstractString)
    tris = load_model_triangles(model_path; scale=ISS_SCALE, rotation_deg=ISS_ROTATION_DEG)
    lo = [minimum(tris[k, :]) for k in 1:3]
    hi = [maximum(tris[k, :]) for k in 1:3]
    return (half_extent_m=Tuple(0.5 .* (hi .- lo)), extent_m=Tuple(hi .- lo))
end

"""Every input that shapes the run, in one record; its digest names the output directory."""
function iss_hypr_inputs(; smoke::Bool=DEMO_SMOKE)
    isfile(ISS_MODEL) || error("ISS display model missing at $(ISS_MODEL)")
    station_model = iss_station_model()
    extents = iss_model_extents(station_model)
    t_half = extents.half_extent_m[2]
    common = (
        model=relpath(ISS_MODEL, REPO_ROOT),
        model_sha256=bytes2hex(open(sha256, ISS_MODEL)),
        model_excluded_nodes=ISS_EXCLUDED_NODES,
        station_model_sha256=bytes2hex(open(sha256, station_model)),
        model_scale=ISS_SCALE,
        model_rotation_deg=ISS_ROTATION_DEG,
        station_half_extent_m=extents.half_extent_m,
        implementation="hypr-manuscript-accel-limited-v1",
        hypr_mode=:manuscript,
        cloud_seed=7,
        station_keepout_radius_m=3.0,
        station_name="iss",
        station_dims_m=(20.0, 73.0, 109.0),
        station_mass_kg=420_000.0,
        safe_distance_m=2.0,
        # Eq. 5 and 6 with the weights and sharpness of Table I. The tolerance
        # covers the cloud's sampling of the surface and the path's sampling.
        w_obs=1.0e6,
        w_fuel=1.0,
        obstacle_sigmoid_k=1.0e5,
        obstacle_sigmoid_tol_m=1.5,
        # Retiming of Table I, with the acceleration-limited passes and no global cap.
        retime_a_max_mps2=0.00625,
        retime_reaction_time_s=0.25,
        retime_speed_scale=0.5,
        reference_max_speed_mps=nothing,
        retime_accel_limit=true,
        retime_dt_s=0.1,
        mpc_terminal_weight_scale=20.0,
        search_margin_m=80.0,
        seed=741,
        smoke=smoke,
    )
    bounded = if smoke
        (
            n_points=20_000,
            start_rtn=(-30.0, -(t_half + 40.0), 0.0),
            goal_rtn=(-30.0, -(t_half + 20.0), 5.0),
            adaptive_sampling_max_ds_m=4.0,
            adaptive_sampling_far_clearance_m=20.0,
            rrt=(iters=200, step_size_m=5.0, goal_sample_rate=0.05, collision_sample_ds_m=0.5, shortcut_iters=20, box_margin_m=20.0),
            ranges=(waypoints=(3, 4), particles=(8, 12), iters=(2, 3)),
            data_rate_s=1.0,
            post_reference_hold_s=10.0,
        )
    else
        (
            n_points=100_000,
            start_rtn=(0.0, -(t_half + 100.0), 0.0),
            goal_rtn=(0.0, t_half + 30.0, 0.0),
            adaptive_sampling_max_ds_m=4.0,
            adaptive_sampling_far_clearance_m=20.0,
            rrt=(iters=1000, step_size_m=5.0, goal_sample_rate=0.05, collision_sample_ds_m=0.5, shortcut_iters=80, box_margin_m=30.0),
            ranges=(waypoints=(3, 8), particles=(60, 160), iters=(10, 60)),
            data_rate_s=0.1,
            post_reference_hold_s=60.0,
        )
    end
    return merge(common, bounded)
end

"""HyPR settings of the run: the `:manuscript` mode with the Table I ranges, scaled to a 100 m structure."""
function iss_hypr_configurator(inputs)
    return RPOPSOConfigurator(
        hypr_mode=:manuscript,
        swarm=RPOPSOSwarmSettings(
            n_waypoints=inputs.ranges.waypoints[1],
            n_particles=inputs.ranges.particles[1],
            n_iters=inputs.ranges.iters[1],
            curve_type=:bezier,
            search_margin_m=inputs.search_margin_m,
            sample_ds_m=inputs.safe_distance_m,
        ),
        objective=RPOPSOObjectiveSettings(
            w_obs=inputs.w_obs,
            w_fuel=inputs.w_fuel,
            obstacle_sigmoid_k=inputs.obstacle_sigmoid_k,
            obstacle_sigmoid_tol_m=inputs.obstacle_sigmoid_tol_m,
            fuel_proxy_dt_s=inputs.retime_dt_s,
        ),
        adaptive=RPOPSOAdaptiveSettings(
            enabled=true,
            n_waypoints_min=inputs.ranges.waypoints[1],
            n_waypoints_max=inputs.ranges.waypoints[2],
            n_particles_min=inputs.ranges.particles[1],
            n_particles_max=inputs.ranges.particles[2],
            n_iters_min=inputs.ranges.iters[1],
            n_iters_max=inputs.ranges.iters[2],
            w_inertia_min=0.4,
            w_inertia_max=0.75,
            c1_min=1.2,
            c1_max=1.8,
            c2_min=1.2,
            c2_max=2.2,
        ),
        adaptive_sampling=SM.GuidanceHooks.RPOAdaptiveSamplingSettings(
            max_ds_m=inputs.adaptive_sampling_max_ds_m,
            far_clearance_m=inputs.adaptive_sampling_far_clearance_m,
        ),
        cull=RPOPSOCullSettings(fraction_max=0.25),
        schedule=RPOPSOScheduleSettings(transition_fraction=0.5),
        stagnation=RPOPSOStagnationSettings(stagnation_learning_threshold=10),
        early_stopping=RPOPSOEarlyStoppingSettings(enabled=false),
        reexplore=RPOPSOReexploreSettings(trigger_iter=10),
        rrt_warmstart=RPOPSORRTConnectWarmstartSettings(
            enabled=true,
            n_iters=inputs.rrt.iters,
            step_size_m=inputs.rrt.step_size_m,
            goal_sample_rate=inputs.rrt.goal_sample_rate,
            collision_sample_ds_m=inputs.rrt.collision_sample_ds_m,
            shortcut_iters=inputs.rrt.shortcut_iters,
            box_margin_m=inputs.rrt.box_margin_m,
        ),
        retiming=RPOPSORetimingSettings(
            dt_s=inputs.retime_dt_s,
            reaction_time_s=inputs.retime_reaction_time_s,
            a_max_mps2=inputs.retime_a_max_mps2,
            speed_scale=inputs.retime_speed_scale,
            accel_limit_enable=inputs.retime_accel_limit,
            initial_speed_mps=0.0,
        ),
    )
end

"""Replace non-finite numbers (which JSON cannot hold) by `nothing` and symbols by strings, recursively."""
_json_safe(x::AbstractFloat) = isfinite(x) ? x : nothing
_json_safe(x::Union{Tuple, AbstractVector}) = map(_json_safe, collect(x))
_json_safe(x::AbstractMatrix) = [_json_safe(collect(c)) for c in eachcol(x)]
_json_safe(x::NamedTuple) = NamedTuple{keys(x)}(map(_json_safe, values(x)))
_json_safe(x::AbstractDict) = Dict(string(k) => _json_safe(v) for (k, v) in x)
_json_safe(x::Symbol) = String(x)
_json_safe(x) = x

_json_roundtrip(x) = JSON.parse(JSON.json(_json_safe(x)))

iss_hypr_outdir(inputs) = joinpath(DEMO_OUT_ROOT, "iss_hypr_" * bytes2hex(sha256(JSON.json(_json_safe(inputs))))[1:10])

"""Sample the station model into the point cloud and build the example scenario for it."""
function iss_hypr_build(inputs, outdir::AbstractString)
    station_model = iss_station_model()
    points = sample_model_pointcloud(
        station_model;
        n_points=inputs.n_points,
        rng=MersenneTwister(inputs.cloud_seed),
        scale=inputs.model_scale,
        rotation_deg=inputs.model_rotation_deg,
    )
    extents = [(minimum(points[k, :]), maximum(points[k, :])) for k in 1:3]
    # The simulated time is the reference plus the hold at the goal; the
    # :manuscript objective does not use the example's `tf_s = mission_time`.
    demo = build_rpo_cubesat_mpc_demo(;
        mission_time=1.0,
        results_directory=outdir,
        seed=inputs.seed,
        start_rtn=SVector{3, Float64}(inputs.start_rtn),
        goal_rtn=SVector{3, Float64}(inputs.goal_rtn),
        pso_configurator=iss_hypr_configurator(inputs),
        station_points=points,
        station_keepout_radius_m=inputs.station_keepout_radius_m,
        station_name=inputs.station_name,
        station_dims_m=inputs.station_dims_m,
        station_mass_kg=inputs.station_mass_kg,
        safe_distance_m=inputs.safe_distance_m,
        search_margin_m=inputs.search_margin_m,
        retime_accel_mps2=inputs.retime_a_max_mps2,
        mpc_terminal_weight_scale=inputs.mpc_terminal_weight_scale,
        post_reference_hold_s=inputs.post_reference_hold_s,
        record_control_commands=true,
        data_rate_s=inputs.data_rate_s,
        verbose=false,
    )
    return demo, extents, station_model
end

"""
Point-cloud sampling of the station surface: distances from dense uniform
surface samples of the model to the nearest cloud point. The largest is an
estimate of the cloud's covering radius (a lower bound; the true value can
exceed it by about the dense samples' spacing).
"""
function iss_cloud_spacing(station_model::AbstractString, geometry; n_dense::Integer=400_000, seed::Integer=99)
    dense = sample_model_pointcloud(station_model; n_points=n_dense, rng=MersenneTwister(seed), scale=ISS_SCALE, rotation_deg=ISS_ROTATION_DEG)
    d = sort!([sqrt(nearest_station_distance_sq(SVector{3, Float64}(dense[:, j]), geometry.station)) for j in 1:size(dense, 2)])
    tris = load_model_triangles(station_model; scale=ISS_SCALE, rotation_deg=ISS_ROTATION_DEG)
    area = 0.0
    for i in 1:(size(tris, 2) ÷ 3)
        a = SVector{3, Float64}(tris[:, 3i - 2])
        area += 0.5 * norm(cross(SVector{3, Float64}(tris[:, 3i - 1]) - a, SVector{3, Float64}(tris[:, 3i]) - a))
    end
    quantile_of(p) = d[clamp(ceil(Int, p * length(d)), 1, length(d))]
    n_cloud = size(geometry.station.points_body, 2)
    return (
        cloud_points=n_cloud,
        dense_samples=n_dense,
        dense_seed=seed,
        surface_area_m2=area,
        mean_spacing_m=sqrt(area / n_cloud),
        dense_spacing_m=sqrt(area / n_dense),
        median_distance_m=quantile_of(0.5),
        p99_distance_m=quantile_of(0.99),
        covering_radius_estimate_m=d[end],
    )
end

"""Clearance along the straight segment between the hold points, to show that it crosses the station."""
function iss_straight_segment_clearance(inputs, geometry)
    straight = hcat(collect(inputs.start_rtn), collect(inputs.goal_rtn))
    samples = SM.GuidanceHooks.rpo_sample_path_polyline(straight, 0.25)
    distances = [sqrt(nearest_station_distance_sq(SVector{3, Float64}(samples[:, j]), geometry.station)) for j in 1:size(samples, 2)]
    stats = rpo_path_clearance_stats(samples, geometry; safe_distance_m=inputs.safe_distance_m)
    return (
        min_distance_to_cloud_m=minimum(distances),
        min_clearance_m=stats.min_clearance,
        samples_below_safe_distance=stats.violation_count,
        samples=size(samples, 2),
        closest_rtn=samples[:, argmin(distances)],
    )
end

"""Which limit sets the reference speed at each interior retiming sample: clearance, curvature, the global cap or the acceleration passes."""
function iss_reference_limit_fractions(profile, cfg)
    n = length(profile.v)
    counts = Dict("clearance" => 0, "curvature" => 0, "acceleration" => 0, "global_cap" => 0)
    for j in 2:(n - 1)
        if profile.v[j] < profile.v_point[j] * (1.0 - 1.0e-9)
            counts["acceleration"] += 1
            continue
        end
        κ = profile.curvature[j]
        v_curve = κ > 0.0 ? cfg.retime_speed_scale * sqrt(cfg.retime_a_max_mps2 / κ) : Inf
        v_cap = cfg.retime_speed_scale * cfg.retime_max_speed_mps
        if isfinite(v_cap) && profile.v_point[j] >= v_cap * (1.0 - 1.0e-9)
            counts["global_cap"] += 1
        elseif profile.v_point[j] >= v_curve * (1.0 - 1.0e-9)
            counts["curvature"] += 1
        else
            counts["clearance"] += 1
        end
    end
    interior = max(n - 2, 1)
    return Dict(k => v / interior for (k, v) in counts)
end

"""The flown plan and the diagnostics of every planning stage (warm start, exploration score, PSO, refinement, retiming)."""
function iss_hypr_plan_record(demo, inputs)
    plan = demo.plan_result
    cfg = plan.config
    ref = SM.GuidanceHooks.rpo_retimed_reference(plan.path, demo.geometry, cfg; safe_distance_m=inputs.safe_distance_m)
    ref.r_rtn == demo.initial_plan.r_ref_rtn || error("the recorded reference differs from the flown one")
    demand = SM.GuidanceHooks.rpo_reference_accel_demand(ref, demo.orbit.mean_motion_radps)
    a_axis = demo.chaser.max_axis_accel_mps2
    per_axis = [maximum(abs, demand.u_rtn[:, k]) for k in axes(demand.u_rtn, 2)]
    magnitude = [norm(demand.u_rtn[:, k]) for k in axes(demand.u_rtn, 2)]
    ka = argmax(per_axis)
    km = argmax(magnitude)
    warm = plan.warmstart
    profile = ref.profile
    return Dict(
        "path_rtn" => Matrix{Float64}(plan.path),
        "cost" => plan.cost,
        "components" => plan.components,
        "refinement_improved" => plan.refinement_improved,
        "pso_best_before_refinement_rtn" => Matrix{Float64}(plan.pso_path),
        "initial_control_points_rtn" => Matrix{Float64}(plan.initial_control_points),
        "pso_bounds_rtn" => plan.initial_bounds,
        "cost_history" => plan.cost_history,
        "counts" => Dict(
            "particles" => cfg.n_particles,
            "iterations" => cfg.n_iters,
            "iterations_run" => length(plan.cost_history),
            "internal_control_points" => cfg.n_waypoints,
            "internal_control_points_from_eta" => plan.adaptive.n_waypoints,
        ),
        "exploration" => plan.adaptive,
        "warm_start" => Dict(
            "attempted" => warm.attempted,
            "path_found" => warm.path_found,
            "iterations" => warm.iterations,
            "path_rtn" => warm.path,
            "raw_path_rtn" => warm.raw_path,
            "path_length_m" => warm.path_length_m,
            "raw_path_length_m" => warm.raw_path_length_m,
            "clearance_threshold_m" => warm.clearance_threshold_m,
            "cost" => warm.cost,
        ),
        "coefficients" => Dict("w_inertia" => cfg.w_inertia, "c1" => cfg.c1, "c2" => cfg.c2,
            "direction" => "equations of Sec. III.A: w and c1 increase with eta, c2 decreases"),
        "t_ref_s" => demo.initial_plan.t_ref_s,
        "r_ref_rtn" => demo.initial_plan.r_ref_rtn,
        "v_ref_rtn" => demo.initial_plan.v_ref_rtn,
        "reference" => Dict(
            "speed_mps" => ref.speed_mps,
            "arc_length_m" => ref.s_m,
            "tangential_accel_mps2" => ref.accel_tangential_mps2,
            "duration_s" => profile.duration_s,
            "length_m" => profile.length_m,
            "start_speed_mps" => ref.speed_mps[1],
            "end_speed_mps" => ref.speed_mps[end],
            "max_speed_mps" => maximum(ref.speed_mps),
            "max_tangential_accel_mps2" => maximum(abs, profile.a_seg),
            "global_speed_cap_mps" => cfg.retime_max_speed_mps,
            "limit_fractions" => iss_reference_limit_fractions(profile, cfg),
            "fallback_samples" => profile.fallback_count,
            "samples" => Dict(
                "s_m" => profile.s, "t_s" => profile.t, "speed_mps" => profile.v, "pointwise_speed_mps" => profile.v_point,
                "curvature_1pm" => profile.curvature, "clearance_m" => profile.clearance,
            ),
            "accel_demand" => Dict(
                "per_axis_capability_mps2" => a_axis,
                "max_per_axis_rtn_mps2" => per_axis[ka],
                "max_per_axis_rtn_ratio" => per_axis[ka] / a_axis,
                "max_per_axis_rtn_at_s" => ref.t_s[ka],
                "max_magnitude_mps2" => magnitude[km],
                "max_magnitude_ratio" => magnitude[km] / a_axis,
                "max_magnitude_at_s" => ref.t_s[km],
                "max_tangential_mps2" => maximum(demand.tangential_mps2),
                "max_centripetal_mps2" => maximum(demand.centripetal_mps2),
                "max_hcw_mps2" => maximum(demand.hcw_mps2),
                "u_rtn_mps2" => demand.u_rtn,
            ),
        ),
        "fuel_proxy" => Dict(
            "J_fuel_kg" => plan.components.J_fuel,
            "delta_v_eq_mps" => plan.components.delta_v_eq_mps,
            "dt_s" => SM.GuidanceHooks.rpo_fuel_proxy_dt_s(cfg),
            "mean_motion_radps" => cfg.mean_motion_radps,
            "mass_kg" => cfg.mass_kg,
            "isp_s" => cfg.isp_s,
        ),
    )
end

"""Write the controller's per-update record as CSV: time, RTN state and command, QP status, forces, mass and attitude."""
function iss_write_control_log(path::AbstractString, log)
    open(path, "w") do io
        println(io, "t_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,ax_cmd_mps2,ay_cmd_mps2,az_cmd_mps2,qp_status,",
            "fbx_des_n,fby_des_n,fbz_des_n,f1_n,f2_n,f3_n,f4_n,f5_n,f6_n,mass_kg,q1,q2,q3,q4")
        for k in eachindex(log.t_s)
            head = vcat([log.t_s[k]], collect(log.x_rel_rtn[k]), collect(log.accel_cmd_rtn[k]))
            tail = vcat(collect(log.force_body_desired_n[k]), collect(log.thruster_forces_n[k]), [log.mass_kg[k]], collect(log.q_chaser[k]))
            println(io, join(head, ","), ",", log.qp_status[k], ",", join(tail, ","))
        end
    end
    return path
end

"""
True when a finished run with exactly these inputs has intact recorded results.
The provenance must include the station, plan and cloud records used by reuse,
and matching SHA-256 hashes for the flown plan, scene, Feather results and the
control log. Missing, empty, altered or older unrecorded bundles cause a fresh run.
"""
function iss_hypr_matching_run(sidecar::AbstractString, inputs, prefix::AbstractString, planfile::AbstractString)
    isfile(sidecar) && filesize(sidecar) > 0 || return false
    recorded = try
        JSON.parsefile(sidecar)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    recorded isa AbstractDict || return false
    all(haskey(recorded, key) for key in ("inputs", "station", "plan", "cloud_extents_m", "outputs")) || return false
    recorded["inputs"] == _json_roundtrip(inputs) || return false
    outputs = recorded["outputs"]
    outputs isa AbstractDict || return false
    for (key, path) in (("plan", planfile), ("scene", prefix * "_scene.json"),
                        ("results_feather", prefix * ".feather"),
                        ("control_log", joinpath(dirname(planfile), "iss_hypr_control_log.csv")))
        isfile(path) && filesize(path) > 0 || return false
        metadata = get(outputs, key, nothing)
        metadata isa AbstractDict || return false
        get(metadata, "bytes", nothing) == filesize(path) || return false
        get(metadata, "sha256", nothing) == bytes2hex(open(sha256, path)) || return false
    end
    return true
end

"""Read a 3 x N matrix back from the JSON column list."""
_matrix3(columns) = Matrix{Float64}(reduce(hcat, [Float64.(c) for c in columns]))

_file_record(path) = isfile(path) ? (path=relpath(path, REPO_ROOT), bytes=filesize(path), sha256=bytes2hex(open(sha256, path))) : nothing

"""Commit of the checkout and whether its tracked files differ from it."""
function iss_git_state()
    commit = try
        readchomp(`git -C $(REPO_ROOT) rev-parse HEAD`)
    catch
        "unknown"
    end
    changed = try
        !isempty(readchomp(`git -C $(REPO_ROOT) status --porcelain --untracked-files=no`))
    catch
        nothing
    end
    return (commit=commit, tracked_changes=changed)
end

"""Every field of the flattened HyPR settings the plan used."""
iss_config_record(cfg) = Dict(String(name) => _json_safe(getfield(cfg, name)) for name in fieldnames(typeof(cfg)))

function main()
    inputs = iss_hypr_inputs()
    outdir = iss_hypr_outdir(inputs)
    prefix = joinpath(outdir, "simulation_results")
    sidecar = joinpath(outdir, "iss_hypr_provenance.json")
    planfile = joinpath(outdir, "iss_hypr_plan.json")
    logfile = joinpath(outdir, "iss_hypr_control_log.csv")
    mkpath(outdir)

    if PLAN_ONLY
        t_plan = @elapsed (demo, _, _) = iss_hypr_build(inputs, outdir)
        record = iss_hypr_plan_record(demo, inputs)
        open(io -> JSON.print(io, _json_safe(record)), planfile, "w")
        println("plan only: ", planfile, " sha256 ", bytes2hex(open(sha256, planfile)), " cost ", demo.plan_result.cost,
            " (", round(t_plan; digits=1), " s, ", Threads.nthreads(), " threads)")
        return record
    end

    reused = iss_hypr_matching_run(sidecar, inputs, prefix, planfile)
    station_model = iss_station_model()
    if reused
        recorded = JSON.parsefile(sidecar)
        saved = JSON.parsefile(planfile)
        path = _matrix3(saved["path_rtn"])
        ref = _matrix3(saved["r_ref_rtn"])
        warm = saved["warm_start"]["path_rtn"]
        warm_path = isempty(warm) ? zeros(3, 0) : _matrix3(warm)
        station = recorded["station"]
        plan = recorded["plan"]
        extents = recorded["cloud_extents_m"]
        println("results with matching provenance present in ", outdir, "; reusing the recorded plan and simulation")
    else
        t_plan = @elapsed (demo, extents, _) = iss_hypr_build(inputs, outdir)
        path = Matrix{Float64}(demo.plan_result.path)
        ref = Matrix{Float64}(demo.initial_plan.r_ref_rtn)   # dense retimed reference, 3 x N
        warm_path = Matrix{Float64}(demo.plan_result.warmstart.path)
        record = iss_hypr_plan_record(demo, inputs)
        open(io -> JSON.print(io, _json_safe(record)), planfile, "w")
        clearance = rpo_path_clearance_stats(ref, demo.geometry; safe_distance_m=inputs.safe_distance_m)
        plan = _json_roundtrip((
            control_points=size(path, 2),
            reference_samples=size(ref, 2),
            cost=demo.plan_result.cost,
            t_ref_end_s=demo.initial_plan.t_ref_s[end],
            simulation_time_s=demo.args.mission_configuration.mission_time,
            min_clearance_m=clearance.min_clearance,
            violation_count=clearance.violation_count,
            eta=demo.plan_result.adaptive.eta,
            warm_start_iterations=demo.plan_result.warmstart.iterations,
            reference_max_speed_mps=record["reference"]["max_speed_mps"],
            accel_demand_max_magnitude_ratio=record["reference"]["accel_demand"]["max_magnitude_ratio"],
            planning_seconds=t_plan,
        ))
        station = _json_roundtrip(demo.station)
        station["cloud_spacing"] = _json_roundtrip(iss_cloud_spacing(station_model, demo.geometry))
        station["straight_segment"] = _json_roundtrip(iss_straight_segment_clearance(inputs, demo.geometry))
        station["orbit"] = _json_roundtrip(demo.orbit)
        station["chaser"] = _json_roundtrip(demo.chaser)
        station["mpc"] = _json_roundtrip(demo.mpc)
        station["config"] = iss_config_record(demo.plan_result.config)
        println("station cloud extents (m): ", [(round(lo; digits=1), round(hi; digits=1)) for (lo, hi) in extents])
        println("plan: eta=", round(plan["eta"]; digits=3), ", ", plan["control_points"], " control points, cost=", plan["cost"],
            ", reference ends at ", round(plan["t_ref_end_s"]; digits=1), " s, minimum clearance ", round(plan["min_clearance_m"]; digits=2), " m")
        println("simulating into ", outdir)
        run_simulation(demo.args; visualization=true, isolate_state=false)
        iss_write_control_log(logfile, demo.control.command_log)
    end
    cp(station_model, joinpath(outdir, "iss_station_model.glb"); force=true)

    paths = Any[
        (name="HyPR reference (RTN)", points_m=ref, frame=:rtn, target=STATION_INDEX, color="#7fe0ff", dashed=true),
        (name="HyPR control points", points_m=path, frame=:rtn, target=STATION_INDEX, color="#ffd166", dashed=false),
    ]
    size(warm_path, 2) >= 2 &&
        push!(paths, (name="RRT-Connect warm start", points_m=warm_path, frame=:rtn, target=STATION_INDEX, color="#c792ea", dashed=true))
    html = export_visualization(
        prefix;
        max_frames=3000,
        trail_s=plan["simulation_time_s"],
        title="AGORA RPO · HyPR relocation around the ISS",
        models=Dict(STATION_ID => station_model),
        model_scale=inputs.model_scale,
        model_rotation_deg=Dict(STATION_ID => inputs.model_rotation_deg),
        paths=paths,
    )

    provenance = Dict(
        "inputs" => _json_roundtrip(inputs),
        "station" => station,
        "cloud_extents_m" => extents,
        "plan" => plan,
        "simulation_reused" => reused,
        "models" => Dict(
            "planet" => "Earth",
            "orbit" => "circular, equatorial, radius Rp_e + 420 km",
            "gravity" => "InverseSquaredGravityModel (point mass)",
            "atmosphere" => "NoAtmosphereModel",
            "ephemerides" => "SimpleEphemeridesModel",
            "spice_kernels" => relpath(SPICE_PATH, REPO_ROOT),
            "station_attitude" => "body frame = RTN: identity at t = 0 and spin n about the orbit normal",
        ),
        "outputs" => Dict(
            "plan" => _file_record(planfile),
            "scene" => _file_record(prefix * "_scene.json"),
            "results_csv" => _file_record(prefix * ".csv"),
            "results_feather" => _file_record(prefix * ".feather"),
            "control_log" => _file_record(logfile),
            "station_model" => _file_record(joinpath(outdir, "iss_station_model.glb")),
            "viewer_html" => _file_record(html),
        ),
        "software" => Dict("julia" => string(VERSION), "julia_threads" => Threads.nthreads(),
            "spaceagora" => string(pkgversion(SpaceAGORA)), "git" => _json_roundtrip(iss_git_state())),
        "seeds" => Dict("planner" => inputs.seed, "cloud" => inputs.cloud_seed,
            "particles" => "one MersenneTwister per particle, drawn from the planner stream after the warm start"),
        "written_at" => Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") * "Z",
    )
    open(sidecar, "w") do io
        JSON.print(io, _json_safe(provenance), 2)
    end
    println("html: ", html, " (", filesize(html), " bytes)")
    println("provenance: ", sidecar)
    return provenance
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
