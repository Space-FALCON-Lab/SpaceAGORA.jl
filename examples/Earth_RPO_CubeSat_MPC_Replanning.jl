include(joinpath(@__DIR__, "Earth_RPO_CubeSat_MPC_Batch.jl"))

using Printf

const RPO_HYPR = SM.GuidanceHooks

function _rpo_replanning_env_int(name::AbstractString, default::Integer)
    token = strip(get(ENV, name, ""))
    isempty(token) && return Int(default)
    return parse(Int, token)
end

function _rpo_replanning_env_float(name::AbstractString, default::Real)
    token = strip(get(ENV, name, ""))
    isempty(token) && return Float64(default)
    return parse(Float64, token)
end

function _rpo_replanning_smoke_mode()
    return get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
end

function _rpo_replanning_progress_bar(frac::Real; width::Integer=30)
    f = clamp(Float64(frac), 0.0, 1.0)
    filled = Int(floor(f * width))
    return "[" * repeat("=", filled) * repeat(" ", width - filled) * "]"
end

function _rpo_replanning_progress_line!(label::AbstractString, completed::Integer, total::Integer; detail::AbstractString="")
    frac = total <= 0 ? 1.0 : completed / total
    percent = 100.0 * clamp(frac, 0.0, 1.0)
    suffix = isempty(detail) ? "" : "  $(detail)"
    @printf("\r%s %s %5.1f%%  %d/%d%s", label, _rpo_replanning_progress_bar(frac), percent, completed, total, suffix)
    flush(stdout)
    completed >= total && println()
    return nothing
end

function rpo_replanning_scenario_specs()
    return [
        (id=:baseline_static_map, label="Baseline static map", expected_action=:none),
        (id=:far_obstacle, label="Small obstacle-cost change far from path", expected_action=:none),
        (id=:near_path_feasible, label="Obstacle appears near path but outside keepout", expected_action=:retime),
        (id=:partial_block, label="Obstacle partially blocks path", expected_action=:replan),
        (id=:behind_vehicle, label="Obstacle appears behind vehicle", expected_action=:none),
        (id=:directly_ahead, label="Obstacle appears directly ahead", expected_action=:replan),
        (id=:goal_region_change, label="Desired end-position change", expected_action=:replan),
        (id=:tracking_disturbance_only, label="Tracking disturbance only", expected_action=:retime),
        (id=:sensor_noise_false_positive, label="Sensor noise false positive", expected_action=:none),
    ]
end

function _rpo_replanning_path_samples(demo)
    cfg = demo.plan_result.config
    return RPO_HYPR.rpo_sample_path(demo.initial_plan.path_rtn, cfg.sample_ds_m; curve_type=cfg.curve_type)
end

function _rpo_replanning_sample_index(samples, alpha)
    return clamp(Int(round(1 + clamp(Float64(alpha), 0.0, 1.0) * (size(samples, 2) - 1))), 1, size(samples, 2))
end

function _rpo_replanning_sample_normal(samples, idx)
    i0 = max(idx - 2, 1)
    i1 = min(idx + 2, size(samples, 2))
    tangent = samples[:, i1] - samples[:, i0]
    norm(tangent) < 1.0e-9 && (tangent = [1.0, 0.0, 0.0])
    axis = abs(dot(tangent ./ norm(tangent), [0.0, 0.0, 1.0])) < 0.9 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    normal = cross(tangent, axis)
    return normal ./ max(norm(normal), 1.0e-9)
end

function _rpo_sphere_path_clearance(samples, center, radius, geometry)
    buffer = geometry.station.keepout_radius_m + maximum(geometry.chaser.half_extents_body)
    min_clearance = Inf
    for j in 1:size(samples, 2)
        min_clearance = min(min_clearance, norm(samples[:, j] .- center) - radius - buffer)
    end
    return min_clearance
end

function _rpo_sphere_center_for_clearance(samples, point, normal, radius, geometry, target_clearance; center_on_path=false)
    center_on_path && return point
    buffer = geometry.station.keepout_radius_m + maximum(geometry.chaser.half_extents_body)
    hi = radius + buffer + max(Float64(target_clearance), 0.0) + 1.0
    while _rpo_sphere_path_clearance(samples, point .+ hi .* normal, radius, geometry) < target_clearance
        hi *= 2.0
        hi > 100.0 && break
    end
    lo = 0.0
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        center = point .+ mid .* normal
        if _rpo_sphere_path_clearance(samples, center, radius, geometry) < target_clearance
            lo = mid
        else
            hi = mid
        end
    end
    return point .+ hi .* normal
end

function _rpo_replanning_sphere_at_path(demo; alpha, appear_alpha, clearance, radius, label, disappear_alpha=Inf, normal_sign=1.0, center_on_path=false)
    samples = _rpo_replanning_path_samples(demo)
    idx = _rpo_replanning_sample_index(samples, alpha)
    point = samples[:, idx]
    normal = normal_sign .* _rpo_replanning_sample_normal(samples, idx)
    center = _rpo_sphere_center_for_clearance(samples, point, normal, radius, demo.geometry, clearance; center_on_path=center_on_path)
    tf_ref = max(demo.initial_plan.t_ref_s[end], demo.pso_config.retime_dt_s)
    return RPO_HYPR.RPOReplanningSphere(
        center,
        radius;
        appear_time_s=Float64(appear_alpha) * tf_ref,
        disappear_time_s=isfinite(disappear_alpha) ? Float64(disappear_alpha) * tf_ref : Inf,
        label=label,
    )
end

function _rpo_replanning_endpoint_buffer(demo)
    return demo.geometry.station.keepout_radius_m + maximum(demo.geometry.chaser.half_extents_body)
end

function _rpo_replanning_sphere_endpoint_clearance(demo, center, radius, point)
    return norm(SVector{3, Float64}(center) - SVector{3, Float64}(point)) - Float64(radius) - _rpo_replanning_endpoint_buffer(demo)
end

function _rpo_replanning_endpoint_safe_center(demo, center, radius, safe_distance_m; clearance_margin_m=1.0e-6)
    adjusted = SVector{3, Float64}(center)
    endpoints = (SVector{3, Float64}(demo.initial_relative_state_rtn[1:3]), SVector{3, Float64}(demo.goal_rtn))
    target = Float64(safe_distance_m) + Float64(clearance_margin_m)
    for _ in 1:24
        clearances = [_rpo_replanning_sphere_endpoint_clearance(demo, adjusted, radius, endpoint) for endpoint in endpoints]
        min_clearance, idx = findmin(clearances)
        min_clearance >= target && return adjusted
        endpoint = endpoints[idx]
        direction = adjusted - endpoint
        if norm(direction) < 1.0e-9
            direction = _rpo_replanning_sample_normal(_rpo_replanning_path_samples(demo), 1)
        end
        adjusted += (target - min_clearance) .* (direction ./ max(norm(direction), 1.0e-9))
    end
    return adjusted
end

function _rpo_replanning_endpoint_safe_sphere(demo, sphere::RPO_HYPR.RPOReplanningSphere, safe_distance_m)
    center = _rpo_replanning_endpoint_safe_center(demo, sphere.center_rtn, sphere.radius_m, safe_distance_m)
    return RPO_HYPR.RPOReplanningSphere(
        center,
        sphere.radius_m;
        appear_time_s=sphere.appear_time_s,
        disappear_time_s=sphere.disappear_time_s,
        velocity_rtn_mps=sphere.velocity_rtn_mps,
        label=sphere.label,
    )
end

function _rpo_replanning_directly_ahead_sphere(demo; current_alpha, appear_alpha, radius, label, safe_distance_m, clearance_margin_m=0.15)
    samples = _rpo_replanning_path_samples(demo)
    current_ref_idx = _rpo_replanning_sample_index(demo.initial_plan.r_ref_rtn, current_alpha)
    current = demo.initial_plan.r_ref_rtn[:, current_ref_idx]
    start_idx = _rpo_replanning_sample_index(samples, current_alpha)
    buffer = _rpo_replanning_endpoint_buffer(demo)
    min_center_distance = Float64(radius) + buffer + Float64(safe_distance_m) + Float64(clearance_margin_m)
    idx = size(samples, 2)
    for j in start_idx:size(samples, 2)
        if norm(samples[:, j] .- current) >= min_center_distance
            idx = j
            break
        end
    end
    tf_ref = max(demo.initial_plan.t_ref_s[end], demo.pso_config.retime_dt_s)
    return _rpo_replanning_endpoint_safe_sphere(demo, RPO_HYPR.RPOReplanningSphere(
        samples[:, idx],
        radius;
        appear_time_s=Float64(appear_alpha) * tf_ref,
        label=label,
    ), safe_distance_m)
end

function _rpo_replanning_changed_goal(
    demo,
    case_seed::Integer;
    endpoint_min_clearance_m::Real=0.26,
    endpoint_max_clearance_m::Real=1.0,
    min_separation_m::Real=1.5,
    surrounded_max_distance_m=2.0,
    max_sampling_tries::Integer=4000,
)
    rng = MersenneTwister(Int(case_seed) + 10_000)
    station_triangles = SpaceAGORA.load_rpo_station_cad_triangles(:gateway)
    _, tri_cdf = _rpo_triangle_areas_and_cdf(station_triangles)
    station_points = demo.geometry.station.points_body
    start = SVector{3, Float64}(demo.initial_relative_state_rtn[1:3])
    old_goal = SVector{3, Float64}(demo.goal_rtn)
    for _ in 1:max_sampling_tries
        candidate = _sample_rpo_near_surface_endpoint(
            rng,
            demo.geometry,
            station_triangles,
            tri_cdf;
            min_clearance_m=endpoint_min_clearance_m,
            max_clearance_m=endpoint_max_clearance_m,
            max_tries=max_sampling_tries,
        )
        if norm(candidate - start) >= min_separation_m &&
            norm(candidate - old_goal) >= min_separation_m &&
            !_rpo_is_surrounded_endpoint(
                candidate,
                station_points;
                max_distance=surrounded_max_distance_m,
            )
            return candidate
        end
    end
    error("Failed to sample a changed goal with batch endpoint constraints after $(max_sampling_tries) tries.")
end

function _rpo_replanning_target_goal(demo, config)
    config.desired_goal_rtn === nothing && return SVector{3, Float64}(demo.goal_rtn)
    return config.desired_goal_rtn
end

function rpo_replanning_case_config(demo, spec; safe_distance_m=0.1, desired_goal_rtn=nothing)
    radius = 0.22
    retime_clearance = 2.0 * safe_distance_m
    near_feasible_clearance = safe_distance_m + 0.80 * (retime_clearance - safe_distance_m)
    blocking_clearance = 0.20 * safe_distance_m
    spheres = RPO_HYPR.RPOReplanningSphere[]
    tracking_error_retime_m = Inf
    desired_goal = spec.id == :goal_region_change ? desired_goal_rtn : nothing

    if spec.id == :far_obstacle
        push!(spheres, _rpo_replanning_sphere_at_path(demo; alpha=0.40, appear_alpha=0.20, clearance=retime_clearance + 5.0, radius=radius, label="far"))
    elseif spec.id == :near_path_feasible
        push!(spheres, _rpo_replanning_sphere_at_path(demo; alpha=0.45, appear_alpha=0.20, clearance=near_feasible_clearance, radius=radius, label="near_feasible"))
    elseif spec.id == :partial_block
        push!(spheres, _rpo_replanning_sphere_at_path(demo; alpha=0.45, appear_alpha=0.20, clearance=blocking_clearance, radius=radius, label="partial_block"))
    elseif spec.id == :behind_vehicle
        push!(spheres, _rpo_replanning_sphere_at_path(demo; alpha=0.12, appear_alpha=0.55, clearance=0.0, radius=radius, label="behind"))
    elseif spec.id == :directly_ahead
        push!(spheres, _rpo_replanning_directly_ahead_sphere(
            demo;
            current_alpha=0.25,
            appear_alpha=0.12,
            radius=radius,
            label="direct_ahead",
            safe_distance_m=safe_distance_m,
        ))
    elseif spec.id == :goal_region_change
        desired_goal === nothing && error("goal_region_change requires desired_goal_rtn.")
    elseif spec.id == :tracking_disturbance_only
        tracking_error_retime_m = 0.75
    elseif spec.id == :sensor_noise_false_positive
        tf_ref = max(demo.initial_plan.t_ref_s[end], demo.pso_config.retime_dt_s)
        push!(spheres, _rpo_replanning_sphere_at_path(
            demo;
            alpha=0.42,
            appear_alpha=0.20,
            disappear_alpha=0.20 + 0.5 * demo.pso_config.retime_dt_s / tf_ref,
            clearance=-0.05,
            radius=radius,
            label="false_positive",
        ))
    elseif spec.id == :moving_obstacle_crossing
        samples = _rpo_replanning_path_samples(demo)
        idx = _rpo_replanning_sample_index(samples, 0.50)
        point = samples[:, idx]
        normal = _rpo_replanning_sample_normal(samples, idx)
        buffer = demo.geometry.station.keepout_radius_m + maximum(demo.geometry.chaser.half_extents_body)
        offset = radius + buffer + max(1.5, retime_clearance)
        tf_ref = max(demo.initial_plan.t_ref_s[end], demo.pso_config.retime_dt_s)
        appear_time = 0.28 * tf_ref
        disappear_time = 0.62 * tf_ref
        duration = max(disappear_time - appear_time, demo.pso_config.retime_dt_s)
        velocity = (-2.0 * offset - 0.02) / duration .* normal
        push!(spheres, RPO_HYPR.RPOReplanningSphere(point .+ offset .* normal, radius; appear_time_s=appear_time, disappear_time_s=disappear_time, velocity_rtn_mps=velocity, label="moving_crossing"))
    end

    spheres = [_rpo_replanning_endpoint_safe_sphere(demo, sphere, safe_distance_m) for sphere in spheres]

    return RPO_HYPR.RPOReplanningConfig(
        enabled=true,
        spheres=spheres,
        desired_goal_rtn=desired_goal,
        goal_change_tolerance_m=0.05,
        safe_distance_m=safe_distance_m,
        retime_clearance_m=retime_clearance,
        hysteresis_samples=2,
        remaining_sample_ds_m=max(demo.pso_config.sample_ds_m, 0.05),
        tracking_error_retime_m=tracking_error_retime_m,
    )
end

function _rpo_replanning_eval_state(demo, spec)
    tf_ref = max(demo.initial_plan.t_ref_s[end], demo.pso_config.retime_dt_s)
    alpha = spec.id == :behind_vehicle ? 0.58 :
        spec.id == :goal_region_change ? 0.46 :
        spec.id == :moving_obstacle_crossing ? 0.45 :
        spec.id == :sensor_noise_false_positive ? 0.21 :
        spec.id == :tracking_disturbance_only ? 0.30 :
        0.25
    idx = _rpo_replanning_sample_index(demo.initial_plan.r_ref_rtn, alpha)
    current = copy(demo.initial_plan.r_ref_rtn[:, idx])
    if spec.id == :tracking_disturbance_only
        current .+= 1.25 .* _rpo_replanning_sample_normal(demo.initial_plan.r_ref_rtn, idx)
    end
    return current, alpha * tf_ref
end

function _rpo_replanning_pso_cost_callback(label::AbstractString)
    return function (iter, cost, components)
        if components === nothing
            @printf("PSO cost [%s] iter=%03d total=%.9g\n", label, iter, cost)
        else
            @printf(
                "PSO cost [%s] iter=%03d total=%.9g len=%.6g obs=%.6g fuel=%.6g min_clearance=%.6g violations=%d\n",
                label,
                iter,
                cost,
                components.J_len_norm,
                components.J_obs,
                components.J_fuel_norm,
                components.min_clearance,
                components.violation_count,
            )
        end
        flush(stdout)
    end
end

function _rpo_replanning_changed_path(demo, config, decision, current, t_eval, case_seed; pso_iteration_callback=nothing)
    if decision.action == :retime
        plan = RPO_HYPR.rpo_retime_existing_plan(
            demo.initial_plan,
            current,
            decision.geometry,
            demo.pso_config,
            config.safe_distance_m,
            t_eval,
        )
        return plan.r_ref_rtn
    elseif decision.action == :replan
        replan_cfg = RPO_HYPR.rpo_pso_config(demo.pso_config; rrt_warmstart_enable=true)
        target_goal = _rpo_replanning_target_goal(demo, config)
        result = RPO_HYPR.rpo_pso_plan_path(
            SVector{3, Float64}(current),
            target_goal,
            decision.geometry,
            replan_cfg;
            safe_distance_m=config.safe_distance_m,
            rng=MersenneTwister(case_seed),
            iteration_callback=pso_iteration_callback,
        )
        return RPO_HYPR.rpo_sample_path(
            result.path,
            result.config.sample_ds_m;
            curve_type=result.config.curve_type,
        )
    end
    return zeros(3, 0)
end

function _rpo_replanning_result_path(demo, decision, current, changed_path)
    size(changed_path, 2) > 0 && return changed_path
    return RPO_HYPR.rpo_remaining_reference_path(
        demo.initial_plan,
        current;
        sample_ds_m=max(demo.pso_config.sample_ds_m, 0.05),
    )
end

function _rpo_replanning_result_metrics(demo, config, decision, current, changed_path; goal_tolerance_m::Real=0.25)
    path = _rpo_replanning_result_path(demo, decision, current, changed_path)
    stats = RPO_HYPR.rpo_path_clearance_stats(path, decision.geometry; safe_distance_m=config.safe_distance_m)
    target_goal = _rpo_replanning_target_goal(demo, config)
    goal_error = size(path, 2) == 0 ? Inf : norm(SVector{3, Float64}(path[:, end]) - target_goal)
    safety_ok = stats.min_clearance + 1.0e-9 >= config.safe_distance_m
    goal_reached = goal_error <= Float64(goal_tolerance_m)
    return (
        path=path,
        safety_ok=safety_ok,
        goal_reached=goal_reached,
        passed=safety_ok && goal_reached,
        min_result_clearance_m=stats.min_clearance,
        result_keepout_violations=stats.violation_count,
        result_goal_error_m=goal_error,
    )
end

function _rpo_replanning_active_obstacle_endpoint_clearances(demo, config, t_eval)
    active = RPO_HYPR.rpo_active_replanning_spheres(config, t_eval)
    isempty(active) && return (start=Inf, goal=Inf)
    start = SVector{3, Float64}(demo.initial_relative_state_rtn[1:3])
    goal = _rpo_replanning_target_goal(demo, config)
    start_clearance = minimum(_rpo_replanning_sphere_endpoint_clearance(demo, sphere.center_rtn, sphere.radius_m, start) for sphere in active)
    goal_clearance = minimum(_rpo_replanning_sphere_endpoint_clearance(demo, sphere.center_rtn, sphere.radius_m, goal) for sphere in active)
    return (start=start_clearance, goal=goal_clearance)
end

function _rpo_replanning_sphere_trace(sphere, t_eval; name_prefix="obstacle")
    center = RPO_HYPR.rpo_replanning_sphere_center(sphere, t_eval)
    radius = sphere.radius_m
    u = range(0.0, 2π; length=24)
    v = range(0.0, π; length=12)
    x = [center[1] + radius * cos(uu) * sin(vv) for vv in v, uu in u]
    y = [center[2] + radius * sin(uu) * sin(vv) for vv in v, uu in u]
    z = [center[3] + radius * cos(vv) for vv in v, uu in u]
    return surface(
        x=x,
        y=y,
        z=z,
        opacity=0.42,
        colorscale=[[0.0, "rgb(205,70,65)"], [1.0, "rgb(245,145,80)"]],
        showscale=false,
        name="$(name_prefix): $(sphere.label)",
        hoverinfo="name",
    )
end

function _rpo_replanning_case_plot(demo, config, decision, current, t_eval, changed_path, spec, case_id)
    baseline = _rpo_replanning_path_samples(demo)
    original_goal = SVector{3, Float64}(demo.goal_rtn)
    target_goal = _rpo_replanning_target_goal(demo, config)
    traces = [
        _station_mesh_trace(),
        scatter3d(
            x=baseline[1, :],
            y=baseline[2, :],
            z=baseline[3, :],
            mode="lines",
            line=attr(width=6, color="rgb(95,105,115)", dash="dash"),
            name="baseline PSO path",
        ),
        scatter3d(
            x=demo.initial_plan.r_ref_rtn[1, :],
            y=demo.initial_plan.r_ref_rtn[2, :],
            z=demo.initial_plan.r_ref_rtn[3, :],
            mode="lines",
            line=attr(width=4, color="rgb(40,105,180)", dash="dot"),
            name="baseline retimed reference",
        ),
        scatter3d(
            x=[current[1]],
            y=[current[2]],
            z=[current[3]],
            mode="markers",
            marker=attr(size=8, color="rgb(25,25,25)", symbol="x"),
            name="retime/replan evaluation point",
        ),
        scatter3d(
            x=[demo.initial_relative_state_rtn[1]],
            y=[demo.initial_relative_state_rtn[2]],
            z=[demo.initial_relative_state_rtn[3]],
            mode="markers",
            marker=attr(size=7, color="rgb(45,150,90)", symbol="circle"),
            name="start",
        ),
        scatter3d(
            x=[target_goal[1]],
            y=[target_goal[2]],
            z=[target_goal[3]],
            mode="markers",
            marker=attr(size=8, color="rgb(190,60,55)", symbol="diamond"),
            name=config.desired_goal_rtn === nothing ? "goal" : "updated goal",
        ),
    ]

    if config.desired_goal_rtn !== nothing && norm(target_goal - original_goal) > config.goal_change_tolerance_m
        push!(traces, scatter3d(
            x=[original_goal[1]],
            y=[original_goal[2]],
            z=[original_goal[3]],
            mode="markers",
            marker=attr(size=7, color="rgb(120,120,120)", symbol="diamond-open"),
            name="original goal",
        ))
    end

    for sphere in RPO_HYPR.rpo_active_replanning_spheres(config, t_eval)
        push!(traces, _rpo_replanning_sphere_trace(sphere, t_eval))
    end

    if size(changed_path, 2) > 0
        push!(traces, scatter3d(
            x=changed_path[1, :],
            y=changed_path[2, :],
            z=changed_path[3, :],
            mode="lines",
            line=attr(width=7, color=decision.action == :replan ? "rgb(215,80,40)" : "rgb(35,150,120)"),
            name=decision.action == :replan ? "replanned path" : "retimed remaining path",
        ))
    end

    title = "Case $(lpad(case_id, 3, '0')) | $(spec.label) | action=$(decision.action), clearance=$(round(decision.min_clearance, digits=4)) m"
    return Plot(
        traces,
        Layout(
            title=title,
            scene=attr(
                aspectmode="data",
                xaxis_title="radial (m)",
                yaxis_title="along-track (m)",
                zaxis_title="cross-track (m)",
            ),
            legend=attr(orientation="h"),
        ),
    )
end

function _rpo_replanning_save_case_plot(results_directory, demo, config, decision, current, t_eval, changed_path, spec, case_id)
    plot_dir = joinpath(results_directory, "replanning_path_plots", "case_$(lpad(case_id, 3, '0'))")
    plot = _rpo_replanning_case_plot(demo, config, decision, current, t_eval, changed_path, spec, case_id)
    return _save_plot(plot, plot_dir, "replanning_$(spec.id).html")
end

function run_rpo_cubesat_mpc_replanning_cases(;
    n_cases::Integer=_rpo_replanning_env_int("SPACEAGORA_RPO_REPLANNING_N", 10),
    results_directory=joinpath(REPO_ROOT, "output", "rpo_replanning_cases"),
    seed::Integer=_rpo_replanning_env_int("SPACEAGORA_RPO_BATCH_SEED", 740),
    mission_time::Real=_rpo_replanning_env_float("SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME", 180.0),
    pso_n_particles::Integer=_rpo_replanning_smoke_mode() ? 10 : _rpo_replanning_env_int("SPACEAGORA_RPO_BATCH_PSO_PARTICLES", 200),
    pso_n_iters::Integer=_rpo_replanning_smoke_mode() ? 4 : _rpo_replanning_env_int("SPACEAGORA_RPO_BATCH_PSO_ITERS", 55),
    pso_config=nothing,
    pso_configurator=nothing,
    n_station_points::Integer=_rpo_replanning_smoke_mode() ? 800 : _rpo_replanning_env_int("SPACEAGORA_RPO_BATCH_STATION_POINTS", 10000),
    save_path_plots::Bool=true,
    show_progress::Bool=true,
    print_pso_costs::Bool=false,
    pso_iteration_runtime_limit_s::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_REPLANNING_PSO_ITERATION_LIMIT", 10.0),
    goal_tolerance_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_REPLANNING_GOAL_TOL", 0.25),
)
    n_cases >= 1 || throw(ArgumentError("n_cases must be >= 1."))
    rows = NamedTuple[]
    demos = NamedTuple[]
    specs = rpo_replanning_scenario_specs()
    cases = generate_rpo_seeded_batch_cases(
        n_cases=n_cases,
        seed=seed,
        geometry_seed=seed,
        n_station_points=n_station_points,
    )
    total_rows = Int(n_cases) * length(specs)
    completed_rows = 0
    show_progress && _rpo_replanning_progress_line!("RPO replanning cases", completed_rows, total_rows)
    for case in cases
        case_id = case.case_id
        case_seed = case.seed
        case_dir = joinpath(results_directory, "case_$(lpad(case_id, 3, '0'))")
        pso_initial_callback = print_pso_costs ? _rpo_replanning_pso_cost_callback("case $(lpad(case_id, 3, '0')) initial") : nothing
        show_progress && print_pso_costs && println()
        demo = build_rpo_cubesat_mpc_demo(
            mission_time=mission_time,
            results_directory=case_dir,
            seed=case_seed,
            station_geometry_seed=seed,
            start_rtn=case.start_rtn,
            goal_rtn=case.goal_rtn,
            pso_n_particles=pso_n_particles,
            pso_n_iters=pso_n_iters,
            pso_config=pso_config,
            pso_configurator=pso_configurator,
            n_station_points=n_station_points,
            data_rate_s=10.0,
            pso_iteration_runtime_limit_s=pso_iteration_runtime_limit_s,
            pso_iteration_callback=pso_initial_callback,
            verbose=false,
        )
        push!(demos, demo)
        for spec in specs
            if show_progress
                detail = "case $(lpad(case_id, 3, '0')) $(spec.id)"
                _rpo_replanning_progress_line!("RPO replanning cases", completed_rows, total_rows; detail=detail)
            end
            desired_goal = spec.id == :goal_region_change ?
                _rpo_replanning_changed_goal(demo, case_seed) :
                nothing
            config = rpo_replanning_case_config(demo, spec; desired_goal_rtn=desired_goal)
            target_goal = _rpo_replanning_target_goal(demo, config)
            current, t_eval = _rpo_replanning_eval_state(demo, spec)
            decision_time_s = @elapsed decision = RPO_HYPR.rpo_replanning_decision(demo.initial_plan, current, demo.geometry, config, t_eval)
            tracking_error = RPO_HYPR.rpo_reference_tracking_error(demo.initial_plan, current)
            changed_path = zeros(3, 0)
            changed_path_time_s = 0.0
            plot_time_s = 0.0
            plot_path = ""
            pso_replan_callback = print_pso_costs && decision.action == :replan ?
                _rpo_replanning_pso_cost_callback("case $(lpad(case_id, 3, '0')) $(spec.id) replan") :
                nothing
            if decision.action in (:retime, :replan)
                show_progress && pso_replan_callback !== nothing && println()
                changed_path_time_s = @elapsed changed_path = _rpo_replanning_changed_path(
                    demo,
                    config,
                    decision,
                    current,
                    t_eval,
                    case_seed,
                    pso_iteration_callback=pso_replan_callback,
                )
            end
            result_metrics = _rpo_replanning_result_metrics(
                demo,
                config,
                decision,
                current,
                changed_path;
                goal_tolerance_m=goal_tolerance_m,
            )
            endpoint_clearance = _rpo_replanning_active_obstacle_endpoint_clearances(demo, config, t_eval)
            if save_path_plots
                plot_time_s = @elapsed plot_path = _rpo_replanning_save_case_plot(
                    results_directory,
                    demo,
                    config,
                    decision,
                    current,
                    t_eval,
                    changed_path,
                    spec,
                    case_id,
                )
            end
            push!(rows, (
                case_id=case_id,
                seed=case_seed,
                start_rtn_1_m=case.start_rtn[1],
                start_rtn_2_m=case.start_rtn[2],
                start_rtn_3_m=case.start_rtn[3],
                goal_rtn_1_m=case.goal_rtn[1],
                goal_rtn_2_m=case.goal_rtn[2],
                goal_rtn_3_m=case.goal_rtn[3],
                desired_goal_rtn_1_m=target_goal[1],
                desired_goal_rtn_2_m=target_goal[2],
                desired_goal_rtn_3_m=target_goal[3],
                goal_change_distance_m=norm(target_goal - SVector{3, Float64}(case.goal_rtn)),
                scenario_id=string(spec.id),
                scenario_label=spec.label,
                expected_action=string(spec.expected_action),
                action=string(decision.action),
                reason=string(decision.reason),
                measured_trigger_distance_m=decision.min_clearance,
                replan_trigger_distance_m=config.safe_distance_m,
                retime_trigger_distance_m=config.retime_clearance_m,
                distance_to_replan_threshold_m=decision.min_clearance - config.safe_distance_m,
                distance_to_retime_threshold_m=decision.min_clearance - config.retime_clearance_m,
                tracking_error_m=tracking_error,
                tracking_error_retime_m=config.tracking_error_retime_m,
                obstacle_start_clearance_m=endpoint_clearance.start,
                obstacle_goal_clearance_m=endpoint_clearance.goal,
                result_min_clearance_m=result_metrics.min_result_clearance_m,
                result_keepout_violations=result_metrics.result_keepout_violations,
                result_goal_error_m=result_metrics.result_goal_error_m,
                safety_ok=result_metrics.safety_ok,
                goal_reached=result_metrics.goal_reached,
                active_spheres=length(decision.spheres),
                decision_runtime_s=decision_time_s,
                changed_path_runtime_s=changed_path_time_s,
                plot_runtime_s=plot_time_s,
                plot_path=plot_path,
                passed=result_metrics.passed,
            ))
            completed_rows += 1
            show_progress && _rpo_replanning_progress_line!("RPO replanning cases", completed_rows, total_rows)
        end
    end
    mkpath(results_directory)
    df = DataFrame(rows)
    CSV.write(joinpath(results_directory, "rpo_replanning_cases.csv"), df)
    return (cases=cases, demos=demos, results=df, csv_path=joinpath(results_directory, "rpo_replanning_cases.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_rpo_cubesat_mpc_replanning_cases()
    println("RPO replanning cases complete. Results: $(abspath(result.csv_path))")
    show(result.results; allrows=true, allcols=true)
    println()
end
