function rpo_estimate_geometry_complexity(start_rtn, goal_rtn, geometry; sample_ds_m::Real=0.25, safe_distance_m::Real=0.0)
    path = hcat(SVector{3, Float64}(start_rtn), SVector{3, Float64}(goal_rtn))
    samples = rpo_sample_path_polyline(path, sample_ds_m)
    stats = rpo_path_clearance_stats(samples, geometry; safe_distance_m=safe_distance_m)
    clearance_term = stats.min_clearance <= 0.0 ? 1.0 : 1.0 / (1.0 + stats.min_clearance)
    return clamp(0.7 * stats.violation_fraction + 0.3 * clearance_term, 0.0, 1.0)
end

function rpo_probe_geometry_metrics(start_rtn, goal_rtn, geometry; sample_ds_m::Real=0.25, safe_distance_m::Real=0.0)
    straight = hcat(SVector{3, Float64}(start_rtn), SVector{3, Float64}(goal_rtn))
    samples = rpo_sample_path_polyline(straight, sample_ds_m)
    stats = rpo_path_clearance_stats(samples, geometry; safe_distance_m=safe_distance_m)
    return (
        detour_ratio=1.0,
        min_clearance=stats.min_clearance,
        violation_fraction=stats.violation_fraction,
        success=stats.violation_count == 0,
    )
end

function rpo_adaptive_pso_config(base::RPOPSOConfig, start_rtn, goal_rtn, geometry; safe_distance_m::Real=0.0)
    dist = norm(SVector{3, Float64}(goal_rtn) - SVector{3, Float64}(start_rtn))
    if !base.adaptive_enable
        return validate_rpo_pso_config(base), (distance_m=dist, complexity=0.0, explore=0.0, enabled=false)
    end

    complexity = rpo_estimate_geometry_complexity(
        start_rtn,
        goal_rtn,
        geometry;
        sample_ds_m=base.sample_ds_m,
        safe_distance_m=safe_distance_m,
    )
    dist_norm = clamp(dist / max(base.cost_ref_distance_m, 1.0e-6), 0.0, 1.0)
    weight_sum = max(base.adaptive_complexity_weight + base.adaptive_distance_weight, 1.0e-9)
    explore = clamp(
        (base.adaptive_complexity_weight * complexity + base.adaptive_distance_weight * dist_norm) / weight_sum,
        0.0,
        1.0,
    )
    n_waypoints_max = base.adaptive_allow_downscale ?
        max(base.adaptive_n_waypoints_max, base.adaptive_n_waypoints_min) :
        max(base.adaptive_n_waypoints_max, base.n_waypoints)
    n_waypoints_min = base.adaptive_allow_downscale ?
        min(base.adaptive_n_waypoints_min, n_waypoints_max) :
        min(max(base.adaptive_n_waypoints_min, base.n_waypoints), n_waypoints_max)
    n_particles_max = base.adaptive_n_particles_max > 0 ? base.adaptive_n_particles_max : max(20, 3 * base.n_particles)
    n_iters_max = base.adaptive_n_iters_max > 0 ? base.adaptive_n_iters_max : max(5, 3 * base.n_iters)
    n_particles_min = min(base.adaptive_n_particles_min, n_particles_max)
    n_iters_min = min(base.adaptive_n_iters_min, n_iters_max)
    effort_scale = base.adaptive_effort_min_fraction +
        (base.adaptive_effort_max_fraction - base.adaptive_effort_min_fraction) * explore
    cfg = rpo_pso_config(
        base;
        n_waypoints=clamp(Int(round(base.n_waypoints + base.adaptive_waypoint_gain * complexity)), n_waypoints_min, n_waypoints_max),
        n_particles=clamp(Int(round(base.n_particles * effort_scale)), n_particles_min, n_particles_max),
        n_iters=clamp(Int(round(base.n_iters * effort_scale)), n_iters_min, n_iters_max),
        w_len=clamp(base.w_len * (1.25 - 0.5 * complexity), base.adaptive_w_len_min, base.adaptive_w_len_max),
        w_obs=clamp(base.w_obs, base.adaptive_w_obs_min, base.adaptive_w_obs_max),
        w_inertia=clamp(0.45 + 0.3 * explore, base.adaptive_w_inertia_min, base.adaptive_w_inertia_max),
        c1=clamp(1.2 + 0.4 * (1.0 - explore), base.adaptive_c1_min, base.adaptive_c1_max),
        c2=clamp(1.2 + 0.8 * explore, base.adaptive_c2_min, base.adaptive_c2_max),
        spread_scale=clamp(base.spread_scale * (0.75 + explore), base.adaptive_spread_scale_min, base.adaptive_spread_scale_max),
    )
    return validate_rpo_pso_config(cfg), (distance_m=dist, complexity=complexity, explore=explore, enabled=true)
end
