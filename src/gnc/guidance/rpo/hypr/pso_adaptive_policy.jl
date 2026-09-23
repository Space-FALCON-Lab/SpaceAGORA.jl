"""Estimate obstacle complexity along the straight-line start-to-goal corridor."""
function rpo_estimate_geometry_complexity(start_rtn, goal_rtn, geometry; sample_ds_m::Real=0.25, safe_distance_m::Real=0.0)
    path = hcat(SVector{3, Float64}(start_rtn), SVector{3, Float64}(goal_rtn))
    samples = rpo_sample_path_polyline(path, sample_ds_m)
    stats = rpo_path_clearance_stats(samples, geometry; safe_distance_m=safe_distance_m)
    clearance_term = stats.min_clearance <= 0.0 ? 1.0 : 1.0 / (1.0 + stats.min_clearance)
    buffer_fraction = stats.min_clearance < Float64(safe_distance_m) ? 1.0 : 0.0
    return clamp(0.7 * buffer_fraction + 0.3 * clearance_term, 0.0, 1.0)
end

"""Measure straight-line path length, clearance, and violation counts for adaptive planning."""
function rpo_probe_geometry_metrics(start_rtn, goal_rtn, geometry; sample_ds_m::Real=0.25, safe_distance_m::Real=0.0)
    straight = hcat(SVector{3, Float64}(start_rtn), SVector{3, Float64}(goal_rtn))
    samples = rpo_sample_path_polyline(straight, sample_ds_m)
    stats = rpo_path_clearance_stats(samples, geometry; safe_distance_m=safe_distance_m)
    return (
        detour_ratio=1.0,
        min_clearance=stats.min_clearance,
        violation_fraction=stats.violation_fraction,
        success=stats.min_clearance + 1.0e-9 >= Float64(safe_distance_m),
    )
end

"""Adjust an RPO PSO config based on quick geometry-probe metrics."""
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
    raw_n_particles_max = base.adaptive_n_particles_max > 0 ? base.adaptive_n_particles_max : max(20, 3 * base.n_particles)
    raw_n_iters_max = base.adaptive_n_iters_max > 0 ? base.adaptive_n_iters_max : max(5, 3 * base.n_iters)
    n_particles_max = base.adaptive_allow_downscale ? raw_n_particles_max : max(raw_n_particles_max, base.n_particles)
    n_iters_max = base.adaptive_allow_downscale ? raw_n_iters_max : max(raw_n_iters_max, base.n_iters)
    n_particles_min = base.adaptive_allow_downscale ?
        min(base.adaptive_n_particles_min, n_particles_max) :
        min(max(base.adaptive_n_particles_min, base.n_particles), n_particles_max)
    n_iters_min = base.adaptive_allow_downscale ?
        min(base.adaptive_n_iters_min, n_iters_max) :
        min(max(base.adaptive_n_iters_min, base.n_iters), n_iters_max)
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

"""
    rpo_manuscript_exploration_score(path_length_m, direct_distance_m, iterations; detour_eps_m=1e-6, effort_scale=100.0)

Exploration score of Sec. III.A of the HyPR manuscript:
D = clip[0,1](L_RRT / max(||r_f - r_0||, ε_d) - 1), S(N) = 1 - exp(-N / N_s)
and η = (2D + S(N)) / 3.
"""
function rpo_manuscript_exploration_score(
    path_length_m::Real,
    direct_distance_m::Real,
    iterations::Integer;
    detour_eps_m::Real=1.0e-6,
    effort_scale::Real=100.0,
)
    D = clamp(Float64(path_length_m) / max(Float64(direct_distance_m), Float64(detour_eps_m)) - 1.0, 0.0, 1.0)
    S = 1.0 - exp(-Float64(iterations) / Float64(effort_scale))
    return (eta=(2.0 * D + S) / 3.0, detour_score=D, search_effort_score=S)
end

"""
    rpo_manuscript_adaptive_pso_config(base, start_rtn, goal_rtn, warmstart)

Map the exploration score η of the RRT-Connect warm start to the initial PSO
coefficients and the particle, iteration and control-point counts (Sec.
III.A). The equations are used as written: w0 = (1-η) w_min + η w_max,
c1,0 = (1-η) c1,min + η c1,max and c2,0 = (1-η) c2,max + η c2,min, so the
inertia and cognitive coefficients grow with η and the social coefficient
shrinks. (The manuscript's prose states the opposite direction for c1 and
c2.) Counts interpolate linearly between the `adaptive_*_min` and
`adaptive_*_max` settings and round to the nearest integer. A warm start that
is disabled or finds no path counts as the largest detour, D = 1. With
`adaptive_enable` false the counts and coefficients of `base` are kept and η
is only reported.
"""
function rpo_manuscript_adaptive_pso_config(base::RPOPSOConfig, start_rtn, goal_rtn, warmstart)
    direct = norm(SVector{3, Float64}(goal_rtn) - SVector{3, Float64}(start_rtn))
    iterations = Int(warmstart.iterations)
    found = Bool(warmstart.path_found)
    length_m = found ? Float64(warmstart.path_length_m) : NaN
    score = if found
        rpo_manuscript_exploration_score(
            length_m,
            direct,
            iterations;
            detour_eps_m=base.adaptive_detour_eps_m,
            effort_scale=base.adaptive_search_effort_scale,
        )
    else
        S = 1.0 - exp(-Float64(iterations) / base.adaptive_search_effort_scale)
        (eta=(2.0 + S) / 3.0, detour_score=1.0, search_effort_score=S)
    end
    η = score.eta
    lerp(lo, hi) = (1.0 - η) * Float64(lo) + η * Float64(hi)
    cfg = if base.adaptive_enable
        rpo_pso_config(
            base;
            w_inertia=lerp(base.adaptive_w_inertia_min, base.adaptive_w_inertia_max),
            c1=lerp(base.adaptive_c1_min, base.adaptive_c1_max),
            c2=lerp(base.adaptive_c2_max, base.adaptive_c2_min),
            n_waypoints=round(Int, lerp(base.adaptive_n_waypoints_min, base.adaptive_n_waypoints_max)),
            n_particles=round(Int, lerp(base.adaptive_n_particles_min, base.adaptive_n_particles_max)),
            n_iters=round(Int, lerp(base.adaptive_n_iters_min, base.adaptive_n_iters_max)),
        )
    else
        validate_rpo_pso_config(base)
    end
    return cfg, (
        enabled=base.adaptive_enable,
        mode=:manuscript,
        eta=η,
        detour_score=score.detour_score,
        search_effort_score=score.search_effort_score,
        rrt_path_found=found,
        rrt_iterations=iterations,
        rrt_path_length_m=length_m,
        direct_distance_m=direct,
        w_inertia=cfg.w_inertia,
        c1=cfg.c1,
        c2=cfg.c2,
        n_particles=cfg.n_particles,
        n_iters=cfg.n_iters,
        n_waypoints=cfg.n_waypoints,
        coefficient_direction=:equations_c1_up_c2_down,
    )
end
