include(joinpath(@__DIR__, "Earth_RPO_CubeSat_MPC_Batch.jl"))

using Printf

const RPO_HYPR = SM.GuidanceHooks
const RPO_REPLANNING_DEFAULT_PLANNERS = [:hypr, :rrt_connect, :rrt_star]

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

function _rpo_replanning_env_bool(name::AbstractString, default::Bool)
    token = lowercase(strip(get(ENV, name, "")))
    isempty(token) && return default
    token in ("1", "true", "yes", "on") && return true
    token in ("0", "false", "no", "off") && return false
    error("Environment variable $name must be boolean-like, got '$token'.")
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
        (id=:vehicle_displacement_retime, label="Vehicle displacement requires retiming", expected_action=:retime),
        (id=:vehicle_displacement_replan, label="Vehicle displacement requires replanning", expected_action=:replan),
        (id=:sensor_noise_false_positive, label="Sensor noise false positive skip", expected_action=:none),
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

function _rpo_replanning_safe_displaced_current(
    demo,
    idx::Integer;
    desired_offset_m::Real,
    min_tracking_error_m::Real,
    safe_distance_m::Real,
    clearance_margin_m::Real=1.0e-6,
)
    ref = SVector{3, Float64}(demo.initial_plan.r_ref_rtn[:, idx])
    i0 = max(Int(idx) - 2, 1)
    i1 = min(Int(idx) + 2, size(demo.initial_plan.r_ref_rtn, 2))
    tangent = SVector{3, Float64}(demo.initial_plan.r_ref_rtn[:, i1] - demo.initial_plan.r_ref_rtn[:, i0])
    tangent_norm = norm(tangent)
    tangent_hat = tangent_norm > 1.0e-9 ? tangent / tangent_norm : SVector{3, Float64}(1.0, 0.0, 0.0)
    normal = SVector{3, Float64}(_rpo_replanning_sample_normal(demo.initial_plan.r_ref_rtn, idx))
    clearance = RPO_HYPR.rpo_clearance_to_station(ref, demo.geometry)
    away = ref - SVector{3, Float64}(clearance.nearest_point)
    away = norm(away) > 1.0e-9 ? away / norm(away) : normal
    outward_normal = away - dot(away, tangent_hat) * tangent_hat
    outward_normal = norm(outward_normal) > 1.0e-9 ? outward_normal / norm(outward_normal) : away
    directions = SVector{3, Float64}[]
    function add_direction!(direction)
        direction_norm = norm(direction)
        direction_norm > 1.0e-9 || return nothing
        candidate = SVector{3, Float64}(direction ./ direction_norm)
        any(existing -> abs(dot(existing, candidate)) > 0.985, directions) && return nothing
        push!(directions, candidate)
        return nothing
    end
    add_direction!(outward_normal)
    add_direction!(-outward_normal)
    add_direction!(normal)
    add_direction!(-normal)
    add_direction!(away)
    for sample_idx in unique(clamp.((idx - 12, idx - 6, idx - 3, idx + 3, idx + 6, idx + 12), 1, size(demo.initial_plan.r_ref_rtn, 2)))
        sample_ref = SVector{3, Float64}(demo.initial_plan.r_ref_rtn[:, sample_idx])
        sample_i0 = max(Int(sample_idx) - 2, 1)
        sample_i1 = min(Int(sample_idx) + 2, size(demo.initial_plan.r_ref_rtn, 2))
        sample_tangent = SVector{3, Float64}(demo.initial_plan.r_ref_rtn[:, sample_i1] - demo.initial_plan.r_ref_rtn[:, sample_i0])
        sample_tangent_norm = norm(sample_tangent)
        sample_tangent_hat = sample_tangent_norm > 1.0e-9 ? sample_tangent / sample_tangent_norm : tangent_hat
        sample_clearance = RPO_HYPR.rpo_clearance_to_station(sample_ref, demo.geometry)
        sample_away = sample_ref - SVector{3, Float64}(sample_clearance.nearest_point)
        sample_away = norm(sample_away) > 1.0e-9 ? sample_away / norm(sample_away) : normal
        sample_outward = sample_away - dot(sample_away, sample_tangent_hat) * sample_tangent_hat
        add_direction!(sample_outward)
    end
    target_clearance = Float64(safe_distance_m) + Float64(clearance_margin_m)
    best_current = nothing
    best_tracking_error = -Inf
    for offset in range(Float64(desired_offset_m), 0.0; length=24)
        for direction in directions
            current = ref + offset * direction
            RPO_HYPR.rpo_clearance_distance_to_station(current, demo.geometry) >= target_clearance || continue
            tracking_error = RPO_HYPR.rpo_reference_tracking_error(demo.initial_plan, current)
            tracking_error > Float64(min_tracking_error_m) && return current
            if tracking_error > best_tracking_error
                best_current = current
                best_tracking_error = tracking_error
            end
        end
    end
    return best_current === nothing ? ref : best_current
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

function rpo_replanning_case_config(
    demo,
    spec;
    safe_distance_m=0.1,
    desired_goal_rtn=nothing,
    tracking_error_retime_m_override=nothing,
    tracking_error_replan_m_override=nothing,
)
    radius = 0.22
    retime_clearance = safe_distance_m + 0.10
    near_feasible_clearance = safe_distance_m + 0.80 * (retime_clearance - safe_distance_m)
    blocking_clearance = 0.20 * safe_distance_m
    spheres = RPO_HYPR.RPOReplanningSphere[]
    tracking_error_retime_m = Inf
    tracking_error_replan_m = Inf
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
    elseif spec.id == :vehicle_displacement_retime
        tracking_error_retime_m = 0.75
    elseif spec.id == :vehicle_displacement_replan
        tracking_error_retime_m = 0.75
        tracking_error_replan_m = 1.10
    elseif spec.id == :sensor_noise_false_positive
        push!(spheres, _rpo_replanning_sphere_at_path(
            demo;
            alpha=0.42,
            appear_alpha=0.18,
            disappear_alpha=0.20,
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

    tracking_error_retime_m_override !== nothing && (tracking_error_retime_m = Float64(tracking_error_retime_m_override))
    tracking_error_replan_m_override !== nothing && (tracking_error_replan_m = Float64(tracking_error_replan_m_override))

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
        tracking_error_replan_m=tracking_error_replan_m,
    )
end

function _rpo_replanning_eval_state(demo, spec; safe_distance_m::Real=0.1)
    tf_ref = max(demo.initial_plan.t_ref_s[end], demo.pso_config.retime_dt_s)
    alpha = spec.id == :behind_vehicle ? 0.58 :
        spec.id == :goal_region_change ? 0.46 :
        spec.id == :moving_obstacle_crossing ? 0.45 :
        spec.id == :sensor_noise_false_positive ? 0.21 :
        spec.id in (:vehicle_displacement_retime, :vehicle_displacement_replan) ? 0.30 :
        0.25
    idx = _rpo_replanning_sample_index(demo.initial_plan.r_ref_rtn, alpha)
    current = copy(demo.initial_plan.r_ref_rtn[:, idx])
    if spec.id == :vehicle_displacement_retime
        current = _rpo_replanning_safe_displaced_current(
            demo,
            idx;
            desired_offset_m=0.95,
            min_tracking_error_m=0.75,
            safe_distance_m=safe_distance_m,
        )
    elseif spec.id == :vehicle_displacement_replan
        current = _rpo_replanning_safe_displaced_current(
            demo,
            idx;
            desired_offset_m=1.35,
            min_tracking_error_m=1.10,
            safe_distance_m=safe_distance_m,
        )
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

function _rpo_replanning_planner_symbols(planners)
    return [RPO_HYPR.normalize_rpo_comparison_planner_type(planner) for planner in planners]
end

function _rpo_replanning_sample_planned_path(plan)
    return RPO_HYPR.rpo_sample_path(
        plan.path,
        plan.config.sample_ds_m;
        curve_type=plan.config.curve_type,
    )
end

function _rpo_replanning_changed_path(
    demo,
    config,
    decision,
    current,
    t_eval,
    case_seed;
    planner::Symbol=:hypr,
    planner_index::Integer=1,
    pso_iteration_callback=nothing,
    rrt_connect_settings::RPO_HYPR.RPORRTConnectSettings=RPO_HYPR.RPORRTConnectSettings(),
    rrt_star_settings::RPO_HYPR.RPORRTStarSettings=RPO_HYPR.RPORRTStarSettings(),
    planner_runtime_limit_s::Real=Inf,
)
    if decision.action == :retime
        plan = RPO_HYPR.rpo_retime_existing_plan(
            demo.initial_plan,
            current,
            decision.geometry,
            demo.pso_config,
            config.safe_distance_m,
            t_eval,
        )
        return (
            path=plan.r_ref_rtn,
            planner_executed=false,
            planner_path_found=true,
            planner_iteration_count=0,
            planner_cost=plan.cost,
        )
    elseif decision.action == :replan
        target_goal = _rpo_replanning_target_goal(demo, config)
        rng = MersenneTwister(Int(case_seed) + 1009 * Int(planner_index))
        if planner == :hypr
            replan_cfg = RPO_HYPR.rpo_pso_config(demo.pso_config; rrt_warmstart_enable=true)
            result = RPO_HYPR.rpo_pso_plan_path(
                SVector{3, Float64}(current),
                target_goal,
                decision.geometry,
                replan_cfg;
                safe_distance_m=config.safe_distance_m,
                rng=rng,
                iteration_callback=pso_iteration_callback,
            )
            return (
                path=_rpo_replanning_sample_planned_path(result),
                planner_executed=true,
                planner_path_found=true,
                planner_iteration_count=length(result.cost_history),
                planner_cost=result.cost,
            )
        elseif planner == :rrt_connect
            replan_cfg = RPO_HYPR.rpo_pso_config(demo.pso_config; curve_type=:polyline)
            result = RPO_HYPR.rpo_rrt_connect_plan_path(
                SVector{3, Float64}(current),
                target_goal,
                decision.geometry,
                replan_cfg;
                safe_distance_m=config.safe_distance_m,
                settings=rrt_connect_settings,
                max_runtime_s=planner_runtime_limit_s,
                rng=rng,
            )
            return (
                path=_rpo_replanning_sample_planned_path(result),
                planner_executed=true,
                planner_path_found=result.path_found,
                planner_iteration_count=result.iterations,
                planner_cost=result.cost,
            )
        elseif planner == :rrt_star
            replan_cfg = RPO_HYPR.rpo_pso_config(demo.pso_config; curve_type=:polyline)
            result = RPO_HYPR.rpo_rrt_star_plan_path(
                SVector{3, Float64}(current),
                target_goal,
                decision.geometry,
                replan_cfg;
                safe_distance_m=config.safe_distance_m,
                settings=rrt_star_settings,
                max_runtime_s=planner_runtime_limit_s,
                rng=rng,
            )
            return (
                path=_rpo_replanning_sample_planned_path(result),
                planner_executed=true,
                planner_path_found=result.path_found,
                planner_iteration_count=result.iterations,
                planner_cost=result.cost,
            )
        end
        throw(ArgumentError("Unsupported replanning planner $(planner)."))
    end
    return (
        path=zeros(3, 0),
        planner_executed=false,
        planner_path_found=true,
        planner_iteration_count=0,
        planner_cost=NaN,
    )
end

function _rpo_replanning_tracking_metrics(demo, config, decision, target_goal, path; tracking_settings, goal_tolerance_m::Real)
    size(path, 2) > 0 || return (
        success=false,
        fuel_used=NaN,
        fuel_used_pct=NaN,
        planned_travel_duration=NaN,
        actual_travel_duration=NaN,
    )
    tracking = RPO_HYPR.RPOLQMPCTrackingSettings(
        dt_s=tracking_settings.dt_s,
        mean_motion_radps=tracking_settings.mean_motion_radps,
        horizon=tracking_settings.horizon,
        mass_kg=tracking_settings.mass_kg,
        propellant_mass_kg=tracking_settings.propellant_mass_kg,
        isp_s=tracking_settings.isp_s,
        g0_mps2=tracking_settings.g0_mps2,
        u_max_mps2=tracking_settings.u_max_mps2,
        q_pos=tracking_settings.q_pos,
        q_vel=tracking_settings.q_vel,
        r_accel=tracking_settings.r_accel,
        qf_pos=tracking_settings.qf_pos,
        qf_vel=tracking_settings.qf_vel,
        settle_time_s=tracking_settings.settle_time_s,
        final_position_tol_m=Float64(goal_tolerance_m),
    )
    return RPO_HYPR.rpo_track_retimed_path_lqmpc(
        path,
        target_goal,
        decision.geometry,
        demo.pso_config,
        tracking;
        safe_distance_m=config.safe_distance_m,
    )
end

function _rpo_replanning_baseline_for_planner(
    demo,
    planner::Symbol,
    case_seed,
    planner_index::Integer;
    safe_distance_m::Real=0.1,
    rrt_connect_settings::RPO_HYPR.RPORRTConnectSettings=RPO_HYPR.RPORRTConnectSettings(),
    rrt_star_settings::RPO_HYPR.RPORRTStarSettings=RPO_HYPR.RPORRTStarSettings(),
    planner_runtime_limit_s::Real=Inf,
)
    start = SVector{3, Float64}(demo.initial_relative_state_rtn[1:3])
    goal = SVector{3, Float64}(demo.goal_rtn)
    rng = MersenneTwister(Int(case_seed) + 50_000 + 1009 * Int(planner_index))
    if planner == :hypr
        base_cfg = RPO_HYPR.rpo_pso_config(demo.pso_config)
        runtime_s = @elapsed result = RPO_HYPR.rpo_pso_plan_path(
            start,
            goal,
            demo.geometry,
            base_cfg;
            safe_distance_m=safe_distance_m,
            rng=rng,
        )
        path_found = true
        iteration_count = length(result.cost_history)
    elseif planner == :rrt_connect
        base_cfg = RPO_HYPR.rpo_pso_config(demo.pso_config; curve_type=:polyline)
        runtime_s = @elapsed result = RPO_HYPR.rpo_rrt_connect_plan_path(
            start,
            goal,
            demo.geometry,
            base_cfg;
            safe_distance_m=safe_distance_m,
            settings=rrt_connect_settings,
            max_runtime_s=planner_runtime_limit_s,
            rng=rng,
        )
        path_found = result.path_found
        iteration_count = result.iterations
    elseif planner == :rrt_star
        base_cfg = RPO_HYPR.rpo_pso_config(demo.pso_config; curve_type=:polyline)
        runtime_s = @elapsed result = RPO_HYPR.rpo_rrt_star_plan_path(
            start,
            goal,
            demo.geometry,
            base_cfg;
            safe_distance_m=safe_distance_m,
            settings=rrt_star_settings,
            max_runtime_s=planner_runtime_limit_s,
            rng=rng,
        )
        path_found = result.path_found
        iteration_count = result.iterations
    else
        throw(ArgumentError("Unsupported replanning baseline planner $(planner)."))
    end

    plan = RPO_HYPR.rpo_plan_from_path(
        result.path,
        demo.geometry,
        result.config,
        safe_distance_m,
        0.0;
        cost=result.cost,
        diagnostics=(
            planner=planner,
            path_found=path_found,
            iterations=iteration_count,
        ),
    )
    return (
        demo=merge(demo, (initial_plan=plan, plan_result=result, pso_config=result.config)),
        runtime_s=runtime_s,
        runtime_limit_s=Float64(planner_runtime_limit_s),
        path_found=path_found,
        iteration_count=iteration_count,
        cost=result.cost,
    )
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
    collision_free = stats.violation_count == 0
    goal_reached = goal_error <= Float64(goal_tolerance_m)
    return (
        path=path,
        safety_ok=collision_free,
        goal_reached=goal_reached,
        passed=collision_free && goal_reached,
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

function _rpo_replanning_case_plot(demo, config, decision, current, t_eval, changed_path, spec, case_id, planner)
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

    title = "Case $(lpad(case_id, 3, '0')) | $(spec.label) | $(RPO_HYPR.rpo_comparison_planner_label(planner)) | action=$(decision.action), clearance=$(round(decision.min_clearance, digits=4)) m"
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

function _rpo_replanning_save_case_plot(results_directory, demo, config, decision, current, t_eval, changed_path, spec, case_id, planner)
    plot_dir = joinpath(results_directory, "replanning_path_plots", "case_$(lpad(case_id, 3, '0'))")
    plot = _rpo_replanning_case_plot(demo, config, decision, current, t_eval, changed_path, spec, case_id, planner)
    return _save_plot(plot, plot_dir, "replanning_$(spec.id)_$(planner).html")
end

function _rpo_replanning_terminal_log_path(results_directory)
    return joinpath(results_directory, "rpo_replanning_terminal_output.log")
end

function _rpo_replanning_failed(result)
    :passed in propertynames(result.results) || return false
    return any(pass -> pass != true, result.results.passed)
end

function _rpo_replanning_drain_terminal_pipe(pipe::Pipe, terminal_io::IO, log_io::IO, log_lock::ReentrantLock)
    while !eof(pipe)
        bytes = readavailable(pipe)
        isempty(bytes) && continue
        write(terminal_io, bytes)
        flush(terminal_io)
        lock(log_lock)
        try
            write(log_io, bytes)
            flush(log_io)
        finally
            unlock(log_lock)
        end
    end
end

function _rpo_replanning_capture_terminal_output(f, log_io::IO)
    original_stdout = stdout
    original_stderr = stderr
    log_lock = ReentrantLock()
    stdout_pipe = redirect_stdout()
    stderr_pipe = redirect_stderr()
    stdout_task = @async _rpo_replanning_drain_terminal_pipe(stdout_pipe, original_stdout, log_io, log_lock)
    stderr_task = @async _rpo_replanning_drain_terminal_pipe(stderr_pipe, original_stderr, log_io, log_lock)
    try
        return f()
    finally
        flush(stdout)
        flush(stderr)
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        close(stdout_pipe)
        close(stderr_pipe)
        wait(stdout_task)
        wait(stderr_task)
        flush(log_io)
    end
end

function run_rpo_cubesat_mpc_replanning_cases(;
    save_failure_terminal_output::Bool=true,
    failure_terminal_output_path=nothing,
    keep_success_terminal_output::Bool=false,
    kwargs...,
)
    if !save_failure_terminal_output
        return _run_rpo_cubesat_mpc_replanning_cases(; kwargs...)
    end

    results_directory = get(kwargs, :results_directory, joinpath(REPO_ROOT, "output", "rpo_replanning_cases"))
    log_path = something(failure_terminal_output_path, _rpo_replanning_terminal_log_path(results_directory))
    mkpath(dirname(log_path))

    result = open(log_path, "w") do log_io
        println(log_io, "RPO replanning terminal output")
        println(log_io, "results_directory=$(abspath(results_directory))")
        println(log_io)
        flush(log_io)
        try
            _rpo_replanning_capture_terminal_output(log_io) do
                _run_rpo_cubesat_mpc_replanning_cases(; kwargs...)
            end
        catch err
            bt = catch_backtrace()
            println(log_io)
            println(log_io, "RPO replanning crashed before completion.")
            Base.showerror(log_io, err, bt)
            println(log_io)
            flush(log_io)
            println(stderr, "RPO replanning terminal output saved to: $(abspath(log_path))")
            rethrow()
        end
    end

    if _rpo_replanning_failed(result)
        open(log_path, "a") do log_io
            println(log_io)
            println(log_io, "RPO replanning completed with failed rows.")
            println(log_io, "csv_path=$(abspath(result.csv_path))")
            show(log_io, result.results; allrows=true, allcols=true)
            println(log_io)
        end
        println(stderr, "RPO replanning terminal output saved to: $(abspath(log_path))")
    elseif !keep_success_terminal_output && isfile(log_path)
        rm(log_path)
    end
    return result
end

function _run_rpo_cubesat_mpc_replanning_cases(;
    n_cases::Integer=_rpo_replanning_smoke_mode() ? 1 : _rpo_replanning_env_int("SPACEAGORA_RPO_REPLANNING_N", _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_N", 10)),
    results_directory=joinpath(REPO_ROOT, "output", "rpo_replanning_cases"),
    seed::Integer=_rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_SEED", _rpo_replanning_env_int("SPACEAGORA_RPO_BATCH_SEED", 740)),
    mission_time::Real=_rpo_replanning_env_float("SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME", 180.0),
    pso_n_particles::Integer=_rpo_replanning_smoke_mode() ? 8 : _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_PSO_PARTICLES", 100),
    pso_n_iters::Integer=_rpo_replanning_smoke_mode() ? 2 : _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_PSO_ITERS", 60),
    pso_n_waypoints::Integer=_rpo_replanning_smoke_mode() ? 1 : _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_PSO_WAYPOINTS", 5),
    pso_config=nothing,
    pso_configurator=nothing,
    n_station_points::Integer=_rpo_replanning_smoke_mode() ? 800 : _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_STATION_POINTS", 10000),
    save_path_plots::Bool=true,
    show_progress::Bool=true,
    print_pso_costs::Bool=false,
    pso_iteration_runtime_limit_s::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_REPLANNING_PSO_ITERATION_LIMIT", 10.0),
    safe_distance_m::Real=RPO_HYPR.RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M,
    endpoint_clearance_margin_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_ENDPOINT_CLEARANCE_MARGIN", 0.05),
    endpoint_max_clearance_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_ENDPOINT_MAX_CLEARANCE", 1.0),
    min_separation_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MIN_SEPARATION", 1.5),
    surrounded_max_distance_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_SURROUNDED_MAX_DISTANCE", 2.0),
    max_sampling_tries::Integer=_rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_MAX_SAMPLING_TRIES", 4000),
    goal_tolerance_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_REPLANNING_GOAL_TOL", _rpo_replanning_smoke_mode() ? 0.75 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_FINAL_TOL", 0.25)),
    planners=RPO_REPLANNING_DEFAULT_PLANNERS,
    rrt_connect_iters::Integer=_rpo_replanning_smoke_mode() ? 25 : _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_RRT_CONNECT_ITERS", 1000),
    rrt_connect_step_size_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_RRT_CONNECT_STEP_SIZE", 0.75),
    rrt_star_iters::Integer=_rpo_replanning_smoke_mode() ? 25 : _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_RRT_STAR_ITERS", rrt_connect_iters),
    rrt_star_step_size_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_RRT_STAR_STEP_SIZE", rrt_connect_step_size_m),
    rrt_star_neighbor_radius_m::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_RRT_STAR_NEIGHBOR_RADIUS", 2.0),
    planner_runtime_limit_s::Real=_rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_RUNTIME_LIMIT_S", 30.0),
    match_hypr_runtime::Bool=_rpo_replanning_env_bool("SPACEAGORA_RPO_COMPARISON_MATCH_HYPR_RUNTIME", true),
    mpc_horizon::Integer=_rpo_replanning_smoke_mode() ? 4 : _rpo_replanning_env_int("SPACEAGORA_RPO_COMPARISON_MPC_HORIZON", 80),
    mpc_u_max_mps2::Real=_rpo_replanning_smoke_mode() ? 0.05 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_U_MAX", 0.01875),
    mpc_q_pos::Real=_rpo_replanning_smoke_mode() ? 10.0 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_Q_POS", 30.0),
    mpc_q_vel::Real=_rpo_replanning_smoke_mode() ? 1.0 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_Q_VEL", 4.0),
    mpc_r_accel::Real=_rpo_replanning_smoke_mode() ? 0.1 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_R_ACCEL", 0.025),
    mpc_qf_pos::Real=_rpo_replanning_smoke_mode() ? 50.0 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_QF_POS", 300.0),
    mpc_qf_vel::Real=_rpo_replanning_smoke_mode() ? 5.0 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_QF_VEL", 25.0),
    mpc_settle_time_s::Real=_rpo_replanning_smoke_mode() ? 2.0 : _rpo_replanning_env_float("SPACEAGORA_RPO_COMPARISON_MPC_SETTLE_TIME", 30.0),
)
    n_cases >= 1 || throw(ArgumentError("n_cases must be >= 1."))
    planner_symbols = _rpo_replanning_planner_symbols(planners)
    if match_hypr_runtime && (:hypr in planner_symbols)
        planner_symbols = vcat([:hypr], [planner for planner in planner_symbols if planner != :hypr])
    end
    comparison_pso_cfg = if pso_config === nothing && pso_configurator === nothing
        cfg = RPO_HYPR.rpo_740_mpc_final_pso_config(
            safe_distance_m=Float64(safe_distance_m);
            n_particles=Int(pso_n_particles),
            n_iters=Int(pso_n_iters),
            n_waypoints=Int(pso_n_waypoints),
            adaptive_enable=!_rpo_replanning_smoke_mode(),
            refinement_enable=true,
        )
        if _rpo_replanning_smoke_mode()
            cfg = RPO_HYPR.rpo_pso_config(
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
        cfg
    else
        pso_config
    end
    tracking_settings = RPO_HYPR.RPOLQMPCTrackingSettings(
        dt_s=(comparison_pso_cfg === nothing ? 0.1 : comparison_pso_cfg.retime_dt_s),
        horizon=Int(mpc_horizon),
        settle_time_s=Float64(mpc_settle_time_s),
        final_position_tol_m=Float64(goal_tolerance_m),
        u_max_mps2=SVector{3, Float64}(Float64(mpc_u_max_mps2), Float64(mpc_u_max_mps2), Float64(mpc_u_max_mps2)),
        q_pos=Float64(mpc_q_pos),
        q_vel=Float64(mpc_q_vel),
        r_accel=Float64(mpc_r_accel),
        qf_pos=Float64(mpc_qf_pos),
        qf_vel=Float64(mpc_qf_vel),
    )
    rrt_connect_settings = RPO_HYPR.RPORRTConnectSettings(
        n_iters=Int(rrt_connect_iters),
        step_size_m=Float64(rrt_connect_step_size_m),
        collision_sample_ds_m=max((comparison_pso_cfg === nothing ? 0.05 : comparison_pso_cfg.sample_ds_m), 0.10),
        shortcut_iters=_rpo_replanning_smoke_mode() ? 4 : 80,
    )
    rrt_star_settings = RPO_HYPR.RPORRTStarSettings(
        n_iters=Int(rrt_star_iters),
        step_size_m=Float64(rrt_star_step_size_m),
        collision_sample_ds_m=max((comparison_pso_cfg === nothing ? 0.05 : comparison_pso_cfg.sample_ds_m), 0.10),
        neighbor_radius_m=Float64(rrt_star_neighbor_radius_m),
        shortcut_iters=_rpo_replanning_smoke_mode() ? 4 : 80,
    )
    rows = NamedTuple[]
    demos = NamedTuple[]
    specs = rpo_replanning_scenario_specs()
    endpoint_min_clearance_m = Float64(safe_distance_m) + max(Float64(endpoint_clearance_margin_m), 0.0)
    cases = generate_rpo_seeded_batch_cases(
        n_cases=n_cases,
        seed=seed,
        geometry_seed=seed,
        n_station_points=n_station_points,
        endpoint_min_clearance_m=endpoint_min_clearance_m,
        endpoint_max_clearance_m=endpoint_max_clearance_m,
        min_separation_m=min_separation_m,
        surrounded_max_distance_m=surrounded_max_distance_m,
        max_sampling_tries=max_sampling_tries,
    )
    total_rows = Int(n_cases) * length(specs) * length(planner_symbols)
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
            pso_config=comparison_pso_cfg,
            pso_configurator=pso_configurator,
            n_station_points=n_station_points,
            data_rate_s=10.0,
            pso_iteration_runtime_limit_s=pso_iteration_runtime_limit_s,
            pso_iteration_callback=pso_initial_callback,
            verbose=false,
        )
        push!(demos, demo)
        baseline_by_planner = Dict{Symbol, NamedTuple}()
        hypr_baseline_runtime_s = NaN
        for (planner_index, planner) in enumerate(planner_symbols)
            baseline_runtime_limit_s = planner == :hypr ?
                Inf :
                (match_hypr_runtime && isfinite(hypr_baseline_runtime_s) ? hypr_baseline_runtime_s : planner_runtime_limit_s)
            baseline_by_planner[planner] = _rpo_replanning_baseline_for_planner(
                demo,
                planner,
                case_seed,
                planner_index;
                safe_distance_m=safe_distance_m,
                rrt_connect_settings=rrt_connect_settings,
                rrt_star_settings=rrt_star_settings,
                planner_runtime_limit_s=baseline_runtime_limit_s,
            )
            planner == :hypr && (hypr_baseline_runtime_s = baseline_by_planner[planner].runtime_s)
        end
        for spec in specs
            if show_progress
                detail = "case $(lpad(case_id, 3, '0')) $(spec.id)"
                _rpo_replanning_progress_line!("RPO replanning cases", completed_rows, total_rows; detail=detail)
            end
            hypr_replan_runtime_s = NaN
            for (planner_index, planner) in enumerate(planner_symbols)
                planner_baseline = baseline_by_planner[planner]
                planner_demo = planner_baseline.demo
                if show_progress
                    detail = "case $(lpad(case_id, 3, '0')) $(spec.id) $(planner)"
                    _rpo_replanning_progress_line!("RPO replanning cases", completed_rows, total_rows; detail=detail)
                end
                desired_goal = spec.id == :goal_region_change ?
                    _rpo_replanning_changed_goal(
                        planner_demo,
                        case_seed;
                        endpoint_min_clearance_m=endpoint_min_clearance_m,
                        endpoint_max_clearance_m=endpoint_max_clearance_m,
                        min_separation_m=min_separation_m,
                        surrounded_max_distance_m=surrounded_max_distance_m,
                        max_sampling_tries=max_sampling_tries,
                    ) :
                    nothing
                config = rpo_replanning_case_config(planner_demo, spec; safe_distance_m=safe_distance_m, desired_goal_rtn=desired_goal)
                current, t_eval = _rpo_replanning_eval_state(planner_demo, spec; safe_distance_m=config.safe_distance_m)
                tracking_error = RPO_HYPR.rpo_reference_tracking_error(planner_demo.initial_plan, current)
                tracking_error_retime_m_override = nothing
                tracking_error_replan_m_override = nothing
                if spec.id == :vehicle_displacement_retime && tracking_error <= config.tracking_error_retime_m
                    tracking_error_retime_m_override = max(0.0, 0.95 * tracking_error)
                elseif spec.id == :vehicle_displacement_replan && tracking_error <= config.tracking_error_replan_m
                    tracking_error_replan_m_override = max(0.0, 0.95 * tracking_error)
                end
                if tracking_error_retime_m_override !== nothing || tracking_error_replan_m_override !== nothing
                    config = rpo_replanning_case_config(
                        planner_demo,
                        spec;
                        safe_distance_m=safe_distance_m,
                        desired_goal_rtn=desired_goal,
                        tracking_error_retime_m_override=tracking_error_retime_m_override,
                        tracking_error_replan_m_override=tracking_error_replan_m_override,
                    )
                end
                target_goal = _rpo_replanning_target_goal(planner_demo, config)
                decision_time_s = @elapsed decision = RPO_HYPR.rpo_replanning_decision(
                    planner_demo.initial_plan,
                    current,
                    planner_demo.geometry,
                    config,
                    t_eval,
                )
                endpoint_clearance = _rpo_replanning_active_obstacle_endpoint_clearances(planner_demo, config, t_eval)
                pso_replan_callback = print_pso_costs && decision.action == :replan && planner == :hypr ?
                    _rpo_replanning_pso_cost_callback("case $(lpad(case_id, 3, '0')) $(spec.id) $(planner) replan") :
                    nothing
                show_progress && pso_replan_callback !== nothing && println()
                update = (
                    path=zeros(3, 0),
                    planner_executed=false,
                    planner_path_found=true,
                    planner_iteration_count=0,
                    planner_cost=NaN,
                )
                changed_path_time_s = 0.0
                effective_planner_runtime_limit_s = planner_runtime_limit_s
                if decision.action == :replan && match_hypr_runtime && planner != :hypr && isfinite(hypr_replan_runtime_s)
                    effective_planner_runtime_limit_s = hypr_replan_runtime_s
                end
                if decision.action in (:retime, :replan)
                    changed_path_time_s = @elapsed update = _rpo_replanning_changed_path(
                        planner_demo,
                        config,
                        decision,
                        current,
                        t_eval,
                        case_seed;
                        planner=planner,
                        planner_index=planner_index,
                        pso_iteration_callback=pso_replan_callback,
                        rrt_connect_settings=rrt_connect_settings,
                        rrt_star_settings=rrt_star_settings,
                        planner_runtime_limit_s=effective_planner_runtime_limit_s,
                    )
                end
                if decision.action == :replan && planner == :hypr
                    hypr_replan_runtime_s = changed_path_time_s
                end
                changed_path = update.path
                result_metrics = _rpo_replanning_result_metrics(
                    planner_demo,
                    config,
                    decision,
                    current,
                    changed_path;
                    goal_tolerance_m=goal_tolerance_m,
                )
                tracking_metrics = _rpo_replanning_tracking_metrics(
                    planner_demo,
                    config,
                    decision,
                    target_goal,
                    result_metrics.path;
                    tracking_settings=tracking_settings,
                    goal_tolerance_m=goal_tolerance_m,
                )
                plot_time_s = 0.0
                plot_path = ""
                if save_path_plots
                    plot_time_s = @elapsed plot_path = _rpo_replanning_save_case_plot(
                        results_directory,
                        planner_demo,
                        config,
                        decision,
                        current,
                        t_eval,
                        changed_path,
                        spec,
                        case_id,
                        planner,
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
                    planner=string(planner),
                    planner_label=RPO_HYPR.rpo_comparison_planner_label(planner),
                    planner_executed=update.planner_executed,
                    planner_path_found=update.planner_path_found,
                    planner_iteration_count=update.planner_iteration_count,
                    planner_cost=update.planner_cost,
                    baseline_path_found=planner_baseline.path_found,
                    baseline_iteration_count=planner_baseline.iteration_count,
                    baseline_planner_runtime_s=planner_baseline.runtime_s,
                    baseline_planner_runtime_limit_s=planner_baseline.runtime_limit_s,
                    baseline_planner_cost=planner_baseline.cost,
                    baseline_path_duration_s=planner_demo.initial_plan.t_ref_s[end],
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
                    planner_runtime_s=changed_path_time_s,
                    planner_runtime_limit_s=effective_planner_runtime_limit_s,
                    match_hypr_runtime=match_hypr_runtime,
                    fuel_used_kg=tracking_metrics.fuel_used,
                    fuel_used_pct=tracking_metrics.fuel_used_pct,
                    updated_path_duration_s=tracking_metrics.planned_travel_duration,
                    actual_path_duration_s=tracking_metrics.actual_travel_duration,
                    success_rate=result_metrics.passed ? 1.0 : 0.0,
                    plot_runtime_s=plot_time_s,
                    plot_path=plot_path,
                    passed=result_metrics.passed,
                ))
                completed_rows += 1
                show_progress && _rpo_replanning_progress_line!("RPO replanning cases", completed_rows, total_rows)
            end
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
