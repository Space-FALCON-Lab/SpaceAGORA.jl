"""Compute path-length and distance references used to normalize RPO cost terms."""
function rpo_path_cost_normalization_refs(points, cfg::RPOPSOConfig)
    pts = Matrix{Float64}(points)
    dx = pts[1, end] - pts[1, 1]
    dy = pts[2, end] - pts[2, 1]
    dz = pts[3, end] - pts[3, 1]
    straight_len = sqrt(dx * dx + dy * dy + dz * dz)
    len_ref = cfg.cost_ref_distance_m > 0.0 ? cfg.cost_ref_distance_m : straight_len
    len_ref = max(len_ref, cfg.sample_ds_m, 1.0e-6)
    v_ref = len_ref / max(cfg.tf_s, 1.0e-6)
    fuel_ref = cfg.mass_kg * v_ref / max(cfg.isp_s * cfg.g0_mps2, 1.0e-9)
    return (straight_len=straight_len, len_ref=len_ref, fuel_ref=max(fuel_ref, 1.0e-12))
end

"""Approximate fuel demand from sampled path increments."""
function rpo_fuel_proxy_from_samples(samples, cfg::RPOPSOConfig)
    pts = Matrix{Float64}(samples)
    size(pts, 2) < 3 && return 0.0
    dt = max(cfg.tf_s / max(size(pts, 2) - 1, 1), 1.0e-6)
    fuel = 0.0
    @inbounds for j in 1:(size(pts, 2) - 2)
        ax = (pts[1, j + 2] - 2.0 * pts[1, j + 1] + pts[1, j]) / (dt * dt)
        ay = (pts[2, j + 2] - 2.0 * pts[2, j + 1] + pts[2, j]) / (dt * dt)
        az = (pts[3, j + 2] - 2.0 * pts[3, j + 1] + pts[3, j]) / (dt * dt)
        fuel += cfg.mass_kg * sqrt(ax * ax + ay * ay + az * az) / max(cfg.isp_s * cfg.g0_mps2, 1.0e-9) * dt
    end
    return fuel
end

"""Summarize continuous capsule clearance and penalties along a sampled polyline."""
function rpo_clearance_stats_from_samples(
    samples,
    geometry,
    safe_distance_m::Real;
    cost_cutoff::Real=Inf,
    w_obs::Real=0.0,
    obstacle_sigmoid_k::Real=1.0e5,
    obstacle_sigmoid_tol_m::Real=0.0,
)
    min_clearance = Inf
    violation_count = 0
    obstacle_score = 0.0
    safe = Float64(safe_distance_m)
    threshold = safe - Float64(obstacle_sigmoid_tol_m)
    k = Float64(obstacle_sigmoid_k)
    @inbounds for j in 1:size(samples, 2)
        p = SVector{3, Float64}(samples[1, j], samples[2, j], samples[3, j])
        clearance = rpo_capsule_clearance_to_station(
            view(samples, :, max(1, j - 1)), p, geometry,
        )
        min_clearance = min(min_clearance, clearance)
        beta = rpo_obstacle_sigmoid_penalty(clearance, threshold, k)
        obstacle_score += beta
        if clearance < 0.0
            violation_count += 1
        end
        if w_obs > 0.0 && isfinite(cost_cutoff) && Float64(w_obs) * obstacle_score > Float64(cost_cutoff)
            return (
                min_clearance=min_clearance,
                violation_count=violation_count,
                violation_fraction=violation_count / max(size(samples, 2), 1),
                obstacle_score=obstacle_score,
                cutoff_exceeded=true,
            )
        end
    end
    return (
        min_clearance=min_clearance,
        violation_count=violation_count,
        violation_fraction=violation_count / max(size(samples, 2), 1),
        obstacle_score=obstacle_score,
        cutoff_exceeded=false,
    )
end

"""Evaluate the smooth obstacle penalty for a clearance value."""
@inline function rpo_obstacle_sigmoid_penalty(clearance::Real, threshold::Real, k::Real)
    x = Float64(k) * (Float64(clearance) - Float64(threshold))
    if x >= 0.0
        y = exp(-x)
        return y / (1.0 + y)
    else
        return 1.0 / (1.0 + exp(x))
    end
end

"""Compute normalized RPO objective components from caller-supplied path samples."""
function rpo_normalized_path_cost_components_from_samples(
    points,
    samples,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    cost_cutoff::Real=Inf,
    w_len::Real=cfg.w_len,
    w_obs::Real=cfg.w_obs,
    w_fuel::Real=cfg.w_fuel,
    clearance_stats=nothing,
    compute_fuel_proxy::Bool=true,
)
    local_w_len = Float64(w_len)
    local_w_obs = Float64(w_obs)
    local_w_fuel = Float64(w_fuel)
    stats = clearance_stats === nothing ? rpo_clearance_stats_from_samples(
        samples,
        geometry,
        safe_distance_m;
        cost_cutoff=cost_cutoff,
        w_obs=local_w_obs,
        obstacle_sigmoid_k=cfg.obstacle_sigmoid_k,
        obstacle_sigmoid_tol_m=cfg.obstacle_sigmoid_tol_m,
    ) : clearance_stats
    J_obs = stats.obstacle_score
    violation_count = stats.min_clearance + cfg.clearance_feasibility_tol_m <
        Float64(safe_distance_m) ? max(1, stats.violation_count) :
        stats.violation_count
    if stats.cutoff_exceeded
        return (
            total=Inf,
            J_len=0.0,
            J_len_norm=0.0,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=stats.min_clearance,
            violation_count=violation_count,
            len_ref=0.0,
            fuel_ref=0.0,
        )
    end
    refs = rpo_path_cost_normalization_refs(points, cfg)
    J_len = rpo_path_length(samples)
    J_len_norm = J_len / refs.len_ref
    partial_cost = local_w_obs * J_obs + local_w_len * J_len_norm^2
    if isfinite(cost_cutoff) && partial_cost > Float64(cost_cutoff)
        return (
            total=Inf,
            J_len=J_len,
            J_len_norm=J_len_norm,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=stats.min_clearance,
            violation_count=violation_count,
            len_ref=refs.len_ref,
            fuel_ref=refs.fuel_ref,
        )
    end
    J_fuel = compute_fuel_proxy ? rpo_fuel_proxy_from_samples(samples, cfg) : 0.0
    J_fuel_norm = J_fuel / refs.fuel_ref
    return (
        total=local_w_len * J_len_norm^2 + local_w_obs * J_obs + local_w_fuel * J_fuel_norm^2,
        J_len=J_len,
        J_len_norm=J_len_norm,
        J_obs=J_obs,
        J_fuel=J_fuel,
        J_fuel_norm=J_fuel_norm,
        min_clearance=stats.min_clearance,
        violation_count=violation_count,
        len_ref=refs.len_ref,
        fuel_ref=refs.fuel_ref,
    )
end


"""Compute normalized RPO objective components for a candidate path."""
function rpo_normalized_path_cost_components(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    cost_cutoff::Real=Inf,
)
    samples = rpo_sample_path(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_hypr_sampling_density_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    return rpo_normalized_path_cost_components_from_samples(
        points,
        samples,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        cost_cutoff=cost_cutoff,
    )
end


"""Sample, safety-check, and prepare one candidate for retimed evaluation."""
function rpo_prepare_retimed_candidate(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    retime_dt_s::Real=cfg.retime_dt_s,
    w_len::Real=cfg.w_len,
    w_obs::Real=cfg.w_obs,
    w_fuel::Real=cfg.w_fuel,
    compute_fuel_proxy::Bool=true,
    retime_mean_motion_radps=nothing,
    retime_command_limit_mps2=nothing,
)
    samples = rpo_sample_path(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=cfg.sample_ds_m,
        curve_type=cfg.curve_type,
    )
    samples = rpo_remove_near_duplicate_samples(samples; warn_removed=false)

    safe = Float64(safe_distance_m)
    threshold = safe - cfg.obstacle_sigmoid_tol_m
    obstacle_score = 0.0
    min_clearance = Inf
    violation_count = 0
    geometry_distance = zeros(size(samples, 2))
    clearance_offset = geometry.station.keepout_radius_m +
        norm(geometry.chaser.half_extents_body)
    @inbounds for j in axes(samples, 2)
        clearance = rpo_clearance_distance_to_station(
            SVector{3, Float64}(samples[1, j], samples[2, j], samples[3, j]),
            geometry,
        )
        geometry_distance[j] = clearance + clearance_offset
        min_clearance = min(min_clearance, clearance)
        obstacle_score += rpo_obstacle_sigmoid_penalty(
            clearance, threshold, cfg.obstacle_sigmoid_k,
        )
        violation_count += clearance < 0.0 ? 1 : 0
    end
    stats = (
        min_clearance=min_clearance,
        violation_count=violation_count,
        violation_fraction=violation_count / max(size(samples, 2), 1),
        obstacle_score=obstacle_score,
        cutoff_exceeded=false,
    )
    components = rpo_normalized_path_cost_components_from_samples(
        points,
        samples,
        geometry,
        cfg;
        safe_distance_m=safe,
        w_len=w_len,
        w_obs=w_obs,
        w_fuel=w_fuel,
        clearance_stats=stats,
        compute_fuel_proxy=compute_fuel_proxy,
    )
    required_violation_count = min_clearance + cfg.clearance_feasibility_tol_m < safe ?
        max(1, violation_count) : violation_count
    components = merge(components, (violation_count=required_violation_count,))
    retiming_samples = samples
    retiming_geometry_distance = geometry_distance
    if required_violation_count == 0 && isfinite(components.total) &&
       retime_mean_motion_radps !== nothing &&
       retime_command_limit_mps2 !== nothing && cfg.curve_type == :bezier
        maximum_speed = isfinite(cfg.retime_max_speed_mps) ?
            cfg.retime_max_speed_mps : sqrt(
                cfg.retime_a_max_mps2 * max(cfg.sample_ds_m, 1.0e-9),
            )
        kinematic_spacing = min(
            cfg.sample_ds_m,
            max(2.0 * Float64(retime_dt_s) * maximum_speed, 1.0e-3),
        )
        if kinematic_spacing + 1.0e-12 < cfg.sample_ds_m
            retiming_samples = rpo_sample_path(
                points, kinematic_spacing; curve_type=:bezier,
            )
            retiming_geometry_distance = [
                rpo_clearance_distance_to_station(view(retiming_samples, :, index), geometry) + clearance_offset
                for index in axes(retiming_samples, 2)
            ]
        end
    end
    # Certify the actual polyline used by retiming, not just its vertices.
    # Retiming may resample a Bezier curve, so the coarse checks cannot certify it.
    if required_violation_count == 0 && isfinite(components.total)
        capsule_min = rpo_clearance_distance_to_station(view(retiming_samples, :, 1), geometry)
        capsule_violations = 0
        capsule_score = 0.0
        for index in 1:(size(retiming_samples, 2) - 1)
            clearance = rpo_capsule_clearance_to_station(
                view(retiming_samples, :, index), view(retiming_samples, :, index + 1), geometry,
            )
            capsule_min = min(capsule_min, clearance)
            capsule_violations += clearance + cfg.clearance_feasibility_tol_m < safe ? 1 : 0
            capsule_score += rpo_obstacle_sigmoid_penalty(clearance, threshold, cfg.obstacle_sigmoid_k)
        end
        required_violation_count = capsule_violations
        components = merge(components, (
            min_clearance=capsule_min,
            violation_count=capsule_violations,
            J_obs=max(components.J_obs, capsule_score),
            total=components.total + Float64(w_obs) * (max(components.J_obs, capsule_score) - components.J_obs),
        ))
    end
    profile = required_violation_count == 0 && isfinite(components.total) ?
        rpo_retiming_profile_from_samples(
            retiming_samples,
            geometry,
            cfg;
            geometry_distance=retiming_geometry_distance,
            samples_are_clean=true,
            retime_dt_s=retime_dt_s,
            retime_mean_motion_radps=retime_mean_motion_radps,
            retime_command_limit_mps2=retime_command_limit_mps2,
        ) : nothing
    # Commands interpolate successive time-reference positions; check those
    # chords too, since a time step can straddle a polyline corner.
    if profile !== nothing
        reference, _, _ = rpo_retime_path_from_profile(profile)
        minimum_clearance = components.min_clearance
        violations = 0
        score = 0.0
        for index in 1:(size(reference, 2) - 1)
            first_segment = searchsortedlast(profile.s_samples, profile.s_ref[index])
            last_segment = searchsortedlast(profile.s_samples, profile.s_ref[index + 1])
            # Chords within one certified segment need no second query.
            first_segment == last_segment && continue
            clearance = rpo_capsule_clearance_to_station(
                view(reference, :, index), view(reference, :, index + 1), geometry,
            )
            minimum_clearance = min(minimum_clearance, clearance)
            violations += clearance + cfg.clearance_feasibility_tol_m < safe ? 1 : 0
            score += rpo_obstacle_sigmoid_penalty(clearance, threshold, cfg.obstacle_sigmoid_k)
        end
        components = merge(components, (
            min_clearance=minimum_clearance,
            violation_count=violations,
            J_obs=max(components.J_obs, score),
            total=components.total + Float64(w_obs) * (max(components.J_obs, score) - components.J_obs),
        ))
        violations > 0 && (profile = nothing)
    end
    return (components=components, profile=profile)
end

"""Return the scalar RPO path objective, optionally cutting off expensive candidates early."""
function rpo_path_cost(points, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0, cost_cutoff::Real=Inf)
    return rpo_normalized_path_cost_components(points, geometry, cfg; safe_distance_m=safe_distance_m, cost_cutoff=cost_cutoff).total
end

"""Evaluate a path with the standard HYPR objective or a caller-supplied component evaluator."""
function rpo_path_objective_components(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    cost_cutoff::Real=Inf,
    objective_evaluator=nothing,
)
    objective_evaluator === nothing && return rpo_normalized_path_cost_components(
        points,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        cost_cutoff=cost_cutoff,
    )
    components = objective_evaluator(
        points,
        geometry,
        cfg,
        Float64(safe_distance_m),
        Float64(cost_cutoff),
    )
    for field in (:total, :J_obs, :violation_count)
        hasproperty(components, field) || throw(ArgumentError(
            "custom RPO objective components must include .$field",
        ))
    end
    return components
end
