ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

include(joinpath(@__DIR__, "Earth_RPO_CubeSat_MPC_Batch.jl"))

import Plots

function _rpo_hypr_video_help()
    return """
    Usage:
      julia --project=. examples/Earth_RPO_CubeSat_MPC_HYPR_Path_Video.jl [N_CASES] [options]

    Options:
      --cases, -n N             Number of seeded HYPR cases to plan and animate.
      --output PATH             Output MP4 path.
      --seed N                  Batch seed.
      --particles N             HYPR/PSO particle count.
      --iters N                 HYPR/PSO iteration count.
      --station-points N        Station point-cloud samples used for planning.
      --display-faces N         Station CAD mesh triangles shown in the video; 0 uses the full mesh.
      --frames-per-case N       Animation frames used to reveal each case.
      --hold-frames N           Extra frames after each completed case.
      --fps N                   Output MP4 frames per second.
      --safe-distance X         Planner safe distance in meters.
      --width N                 Video width in pixels.
      --height N                Video height in pixels.
      --zoom X                  Tighten view padding without cropping; larger values fit closer.
      --allow-crop-zoom         Let --zoom crop the fitted mesh/path bounds.
      --camera-azimuth X        3D camera azimuth angle.
      --camera-elevation X      3D camera elevation angle.

    Example:
      julia --project=. examples/Earth_RPO_CubeSat_MPC_HYPR_Path_Video.jl 10 --output output/rpo_hypr_path_video/hypr_10_cases.mp4
    """
end

function _rpo_hypr_video_arg_value(args, i)
    token = args[i]
    eq = findfirst(==('='), token)
    if eq !== nothing
        return token[eq + 1:end], i
    end
    i < length(args) || error("Missing value after $(token).")
    return args[i + 1], i + 1
end

function _rpo_hypr_video_parse_args(args)
    opts = Dict{Symbol, Any}(
        :n_cases => _rpo_batch_smoke_mode() ? 1 : _env_int("SPACEAGORA_RPO_HYPR_VIDEO_N", 10),
        :seed => _env_int("SPACEAGORA_RPO_HYPR_VIDEO_SEED", _env_int("SPACEAGORA_RPO_BATCH_SEED", 740)),
        :output => get(ENV, "SPACEAGORA_RPO_HYPR_VIDEO_OUTPUT", joinpath(REPO_ROOT, "output", "rpo_hypr_path_video", "hypr_paths.mp4")),
        :particles => _rpo_batch_smoke_mode() ? 8 : _env_int("SPACEAGORA_RPO_HYPR_VIDEO_PSO_PARTICLES", 80),
        :iters => _rpo_batch_smoke_mode() ? 2 : _env_int("SPACEAGORA_RPO_HYPR_VIDEO_PSO_ITERS", 30),
        :station_points => _rpo_batch_smoke_mode() ? 800 : _env_int("SPACEAGORA_RPO_HYPR_VIDEO_STATION_POINTS", 6000),
        :display_faces => _rpo_batch_smoke_mode() ? 200 : _env_int("SPACEAGORA_RPO_HYPR_VIDEO_DISPLAY_FACES", 0),
        :frames_per_case => _rpo_batch_smoke_mode() ? 8 : _env_int("SPACEAGORA_RPO_HYPR_VIDEO_FRAMES_PER_CASE", 36),
        :hold_frames => _rpo_batch_smoke_mode() ? 2 : _env_int("SPACEAGORA_RPO_HYPR_VIDEO_HOLD_FRAMES", 8),
        :fps => _env_int("SPACEAGORA_RPO_HYPR_VIDEO_FPS", 12),
        :safe_distance_m => _env_float("SPACEAGORA_RPO_HYPR_VIDEO_SAFE_DISTANCE", 0.1),
        :width_px => _env_int("SPACEAGORA_RPO_HYPR_VIDEO_WIDTH", 900),
        :height_px => _env_int("SPACEAGORA_RPO_HYPR_VIDEO_HEIGHT", 720),
        :zoom => _env_float("SPACEAGORA_RPO_HYPR_VIDEO_ZOOM", 1.0),
        :allow_crop_zoom => _env_bool("SPACEAGORA_RPO_HYPR_VIDEO_ALLOW_CROP_ZOOM", false),
        :camera_azimuth => _env_float("SPACEAGORA_RPO_HYPR_VIDEO_CAMERA_AZIMUTH", 42.0),
        :camera_elevation => _env_float("SPACEAGORA_RPO_HYPR_VIDEO_CAMERA_ELEVATION", 24.0),
    )

    i = 1
    while i <= length(args)
        token = strip(args[i])
        if token in ("--help", "-h")
            println(_rpo_hypr_video_help())
            return nothing
        elseif !startswith(token, "-")
            opts[:n_cases] = parse(Int, token)
        elseif token in ("--cases", "--n-cases", "-n") || startswith(token, "--cases=") || startswith(token, "--n-cases=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:n_cases] = parse(Int, value)
        elseif token == "--output" || startswith(token, "--output=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:output] = value
        elseif token == "--seed" || startswith(token, "--seed=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:seed] = parse(Int, value)
        elseif token == "--particles" || startswith(token, "--particles=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:particles] = parse(Int, value)
        elseif token == "--iters" || startswith(token, "--iters=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:iters] = parse(Int, value)
        elseif token == "--station-points" || startswith(token, "--station-points=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:station_points] = parse(Int, value)
        elseif token in ("--display-faces", "--mesh-faces") || startswith(token, "--display-faces=") || startswith(token, "--mesh-faces=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:display_faces] = parse(Int, value)
        elseif token == "--frames-per-case" || startswith(token, "--frames-per-case=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:frames_per_case] = parse(Int, value)
        elseif token == "--hold-frames" || startswith(token, "--hold-frames=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:hold_frames] = parse(Int, value)
        elseif token == "--fps" || startswith(token, "--fps=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:fps] = parse(Int, value)
        elseif token == "--safe-distance" || startswith(token, "--safe-distance=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:safe_distance_m] = parse(Float64, value)
        elseif token == "--width" || startswith(token, "--width=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:width_px] = parse(Int, value)
        elseif token == "--height" || startswith(token, "--height=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:height_px] = parse(Int, value)
        elseif token == "--zoom" || startswith(token, "--zoom=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:zoom] = parse(Float64, value)
        elseif token == "--allow-crop-zoom"
            opts[:allow_crop_zoom] = true
        elseif token == "--no-allow-crop-zoom"
            opts[:allow_crop_zoom] = false
        elseif token == "--camera-azimuth" || startswith(token, "--camera-azimuth=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:camera_azimuth] = parse(Float64, value)
        elseif token == "--camera-elevation" || startswith(token, "--camera-elevation=")
            value, i = _rpo_hypr_video_arg_value(args, i)
            opts[:camera_elevation] = parse(Float64, value)
        else
            error("Unknown option $(token). Run with --help for usage.")
        end
        i += 1
    end

    return opts
end

function _rpo_hypr_video_config(; n_particles::Integer, n_iters::Integer, safe_distance_m::Real)
    cfg = rpo_740_mpc_final_pso_config(
        safe_distance_m=safe_distance_m;
        n_particles=Int(n_particles),
        n_iters=Int(n_iters),
        adaptive_enable=!_rpo_batch_smoke_mode(),
        adaptive_n_particles_min=Int(n_particles),
        adaptive_n_particles_max=Int(n_particles),
        adaptive_n_iters_min=Int(n_iters),
        adaptive_n_iters_max=Int(n_iters),
    )
    if _rpo_batch_smoke_mode()
        cfg = rpo_pso_config(
            cfg;
            cull_enable=false,
            schedule_enable=false,
            reexplore_enable=false,
            refinement_enable=false,
            sample_ds_m=0.5,
            retime_dt_s=0.5,
            retime_a_max_mps2=0.05,
            retime_max_steps=5_000,
        )
    end
    return cfg
end

function _rpo_hypr_video_geometry(seed::Integer, n_station_points::Integer)
    station_points = SpaceAGORA.load_rpo_station_cad_pointcloud(
        :gateway;
        n_points=n_station_points,
        rng=MersenneTwister(seed),
    )
    return RPOReferenceGeometry(
        RPOStationGeometry(station_points; keepout_radius_m=0.25, name="gateway_core");
        chaser=RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.3)),
    )
end

function _rpo_hypr_video_plan_cases(cases, geometry, pso_cfg, safe_distance_m::Real)
    planned = NamedTuple[]
    for case in cases
        label = _rpo_batch_case_label(case.case_id)
        print("Planning HYPR video case $(label) ... ")
        flush(stdout)
        runtime_s = @elapsed plan = rpo_pso_plan_path(
            case.start_rtn,
            case.goal_rtn,
            geometry,
            pso_cfg;
            safe_distance_m=safe_distance_m,
            rng=MersenneTwister(case.seed),
        )
        t_ref, r_ref, _ = rpo_reference_from_path(
            plan.path,
            geometry,
            plan.config;
            safe_distance_m=safe_distance_m,
        )
        println("done in $(round(runtime_s; digits=2)) s")
        push!(planned, (
            case_id=case.case_id,
            label=label,
            seed=case.seed,
            start_rtn=SVector{3, Float64}(case.start_rtn),
            goal_rtn=SVector{3, Float64}(case.goal_rtn),
            plan=plan,
            t_ref=collect(Float64.(t_ref)),
            r_ref=Matrix{Float64}(r_ref),
            runtime_s=runtime_s,
        ))
    end
    return planned
end

function _rpo_hypr_video_mesh(triangles; max_faces::Integer=0)
    tris = Matrix{Float64}(triangles)
    ntri = size(tris, 2) ÷ 3
    ntri > 0 || return (x=Float64[], y=Float64[], z=Float64[], i=Int[], j=Int[], k=Int[])
    tri_idxs = if max_faces <= 0 || max_faces >= ntri
        1:ntri
    else
        stride = max(1, Int(ceil(ntri / Int(max_faces))))
        1:stride:ntri
    end
    nfaces = length(tri_idxs)
    vertices = zeros(Float64, 3, 3 * nfaces)
    i = zeros(Int, nfaces)
    j = zeros(Int, nfaces)
    k = zeros(Int, nfaces)
    for (face_idx, tri_idx) in enumerate(tri_idxs)
        base = 3 * (tri_idx - 1)
        dst = 3 * (face_idx - 1)
        vertices[:, dst + 1] .= tris[:, base + 1]
        vertices[:, dst + 2] .= tris[:, base + 2]
        vertices[:, dst + 3] .= tris[:, base + 3]
        i[face_idx] = dst
        j[face_idx] = dst + 1
        k[face_idx] = dst + 2
    end
    return (x=vertices[1, :], y=vertices[2, :], z=vertices[3, :], i=i, j=j, k=k)
end

function _rpo_hypr_video_bounds(station_triangles, planned; zoom::Real=1.0, allow_crop_zoom::Bool=false)
    mins = vec(minimum(station_triangles; dims=2))
    maxs = vec(maximum(station_triangles; dims=2))
    for item in planned
        path = item.r_ref
        mins .= min.(mins, vec(minimum(path; dims=2)), collect(item.start_rtn), collect(item.goal_rtn))
        maxs .= max.(maxs, vec(maximum(path; dims=2)), collect(item.start_rtn), collect(item.goal_rtn))
    end
    center = 0.5 .* (mins .+ maxs)
    span = maximum(maxs .- mins)
    pad = max(0.12 * span, 0.5)
    zoom_factor = Float64(zoom)
    zoom_factor > 0.0 || throw(ArgumentError("zoom must be positive."))
    half = allow_crop_zoom ?
        (0.5 * span + pad) / zoom_factor :
        0.5 * span + pad / zoom_factor
    return (
        x=(center[1] - half, center[1] + half),
        y=(center[2] - half, center[2] + half),
        z=(center[3] - half, center[3] + half),
    )
end

function _rpo_hypr_video_case_color(case_idx::Integer)
    colors = (
        :dodgerblue3,
        :orangered3,
        :seagreen3,
        :purple3,
        :goldenrod3,
        :deepskyblue4,
        :tomato2,
        :darkolivegreen3,
        :magenta3,
        :sienna3,
    )
    return colors[mod1(Int(case_idx), length(colors))]
end

function _rpo_hypr_video_frame_plot(
    station_mesh,
    planned,
    case_idx::Integer,
    frame_idx::Integer,
    frames_per_case::Integer,
    bounds;
    camera=(42.0, 24.0),
    plot_size=(900, 720),
)
    item = planned[case_idx]
    n_path = size(item.r_ref, 2)
    frac = frames_per_case <= 1 ? 1.0 : (frame_idx - 1) / (frames_per_case - 1)
    prefix = clamp(Int(round(1 + frac * (n_path - 1))), 1, n_path)
    title_time = isempty(item.t_ref) ? 0.0 : item.t_ref[prefix]
    title_total = isempty(item.t_ref) ? 0.0 : item.t_ref[end]

    p = Plots.mesh3d(
        station_mesh.x,
        station_mesh.y,
        station_mesh.z;
        connections=(station_mesh.i, station_mesh.j, station_mesh.k),
        fillcolor=:gray70,
        fillalpha=0.42,
        linecolor=:gray45,
        linewidth=0.15,
        label=false,
        xlims=bounds.x,
        ylims=bounds.y,
        zlims=bounds.z,
        xlabel="radial (m)",
        ylabel="along-track (m)",
        zlabel="cross-track (m)",
        title="RPO HYPR path sequence: case $(item.label)  t=$(round(title_time; digits=1))/$(round(title_total; digits=1)) s",
        legend=false,
        camera=camera,
        size=plot_size,
        dpi=120,
    )

    for prior_idx in 1:(case_idx - 1)
        prior = planned[prior_idx]
        color = _rpo_hypr_video_case_color(prior_idx)
        Plots.plot!(
            p,
            prior.r_ref[1, :],
            prior.r_ref[2, :],
            prior.r_ref[3, :];
            linecolor=color,
            linewidth=2,
            linealpha=0.28,
            label=false,
        )
    end

    color = _rpo_hypr_video_case_color(case_idx)
    Plots.plot!(
        p,
        item.r_ref[1, 1:prefix],
        item.r_ref[2, 1:prefix],
        item.r_ref[3, 1:prefix];
        linecolor=color,
        linewidth=4,
        label=false,
    )
    Plots.scatter!(
        p,
        [item.start_rtn[1]],
        [item.start_rtn[2]],
        [item.start_rtn[3]];
        marker=:circle,
        markersize=5,
        markercolor=:green3,
        markerstrokecolor=:white,
        label=false,
    )
    Plots.scatter!(
        p,
        [item.goal_rtn[1]],
        [item.goal_rtn[2]],
        [item.goal_rtn[3]];
        marker=:diamond,
        markersize=6,
        markercolor=:red3,
        markerstrokecolor=:white,
        label=false,
    )
    Plots.scatter!(
        p,
        [item.r_ref[1, prefix]],
        [item.r_ref[2, prefix]],
        [item.r_ref[3, prefix]];
        marker=:circle,
        markersize=7,
        markercolor=color,
        markerstrokecolor=:white,
        label=false,
    )
    return p
end

function _rpo_hypr_video_write_summary(output_path::AbstractString, planned)
    rows = [
        (
            case_id=item.case_id,
            seed=item.seed,
            start_x=item.start_rtn[1],
            start_y=item.start_rtn[2],
            start_z=item.start_rtn[3],
            goal_x=item.goal_rtn[1],
            goal_y=item.goal_rtn[2],
            goal_z=item.goal_rtn[3],
            planner_runtime_s=item.runtime_s,
            planned_duration_s=isempty(item.t_ref) ? 0.0 : item.t_ref[end],
            path_points=size(item.r_ref, 2),
            planner_cost=item.plan.cost,
            min_clearance_m=item.plan.components.min_clearance,
            violation_count=item.plan.components.violation_count,
        )
        for item in planned
    ]
    csv_path = joinpath(dirname(output_path), "hypr_path_video_summary.csv")
    CSV.write(csv_path, DataFrame(rows))
    return csv_path
end

function run_rpo_hypr_path_video(;
    n_cases::Integer=10,
    seed::Integer=740,
    output::AbstractString=joinpath(REPO_ROOT, "output", "rpo_hypr_path_video", "hypr_paths.mp4"),
    pso_n_particles::Integer=80,
    pso_n_iters::Integer=30,
    n_station_points::Integer=6000,
    display_faces::Integer=0,
    frames_per_case::Integer=36,
    hold_frames::Integer=8,
    fps::Integer=12,
    safe_distance_m::Real=0.1,
    width_px::Integer=900,
    height_px::Integer=720,
    zoom::Real=1.0,
    allow_crop_zoom::Bool=false,
    camera_azimuth::Real=42.0,
    camera_elevation::Real=24.0,
)
    n_cases > 0 || throw(ArgumentError("n_cases must be positive."))
    frames_per_case > 0 || throw(ArgumentError("frames_per_case must be positive."))
    fps > 0 || throw(ArgumentError("fps must be positive."))

    output_path = abspath(output)
    mkpath(dirname(output_path))

    geometry = _rpo_hypr_video_geometry(seed, n_station_points)
    cases = generate_rpo_seeded_batch_cases(
        n_cases=n_cases,
        seed=seed,
        geometry_seed=seed,
        n_station_points=n_station_points,
    )
    pso_cfg = _rpo_hypr_video_config(
        n_particles=pso_n_particles,
        n_iters=pso_n_iters,
        safe_distance_m=safe_distance_m,
    )
    planned = _rpo_hypr_video_plan_cases(cases, geometry, pso_cfg, Float64(safe_distance_m))

    station_triangles = SpaceAGORA.load_rpo_station_cad_triangles(:gateway)
    station_mesh = _rpo_hypr_video_mesh(station_triangles; max_faces=display_faces)
    bounds = _rpo_hypr_video_bounds(station_triangles, planned; zoom=zoom, allow_crop_zoom=allow_crop_zoom)
    camera = (Float64(camera_azimuth), Float64(camera_elevation))
    plot_size = (Int(width_px), Int(height_px))

    Plots.gr()
    anim = Plots.Animation()
    for case_idx in eachindex(planned)
        for frame_idx in 1:frames_per_case
            plt = _rpo_hypr_video_frame_plot(
                station_mesh,
                planned,
                case_idx,
                frame_idx,
                frames_per_case,
                bounds;
                camera=camera,
                plot_size=plot_size,
            )
            Plots.frame(anim, plt)
        end
        for _ in 1:max(0, hold_frames)
            plt = _rpo_hypr_video_frame_plot(
                station_mesh,
                planned,
                case_idx,
                frames_per_case,
                frames_per_case,
                bounds;
                camera=camera,
                plot_size=plot_size,
            )
            Plots.frame(anim, plt)
        end
    end

    Plots.mp4(anim, output_path; fps=fps)
    summary_path = _rpo_hypr_video_write_summary(output_path, planned)
    println("HYPR path MP4: ", output_path)
    println("Summary CSV:   ", abspath(summary_path))
    return (
        video_path=output_path,
        summary_path=summary_path,
        cases=cases,
        planned=planned,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = _rpo_hypr_video_parse_args(ARGS)
    if opts !== nothing
        run_rpo_hypr_path_video(
            n_cases=opts[:n_cases],
            seed=opts[:seed],
            output=opts[:output],
            pso_n_particles=opts[:particles],
            pso_n_iters=opts[:iters],
            n_station_points=opts[:station_points],
            display_faces=opts[:display_faces],
            frames_per_case=opts[:frames_per_case],
            hold_frames=opts[:hold_frames],
            fps=opts[:fps],
            safe_distance_m=opts[:safe_distance_m],
            width_px=opts[:width_px],
            height_px=opts[:height_px],
            zoom=opts[:zoom],
            allow_crop_zoom=opts[:allow_crop_zoom],
            camera_azimuth=opts[:camera_azimuth],
            camera_elevation=opts[:camera_elevation],
        )
    end
end
