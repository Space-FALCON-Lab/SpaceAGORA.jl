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

"""Summarize obstacle clearance, violations, and penalties along sampled RPO path points."""
function rpo_clearance_stats_from_samples(
    samples,
    geometry,
    safe_distance_m::Real;
    cost_cutoff::Real=Inf,
    w_obs::Real=0.0,
    obstacle_sigmoid_k::Real=1.0e6,
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
        clearance = rpo_clearance_distance_to_station(p, geometry)
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
    stats = rpo_clearance_stats_from_samples(
        samples,
        geometry,
        safe_distance_m;
        cost_cutoff=cost_cutoff,
        w_obs=cfg.w_obs,
        obstacle_sigmoid_k=cfg.obstacle_sigmoid_k,
        obstacle_sigmoid_tol_m=cfg.obstacle_sigmoid_tol_m,
    )
    J_obs = stats.obstacle_score
    if stats.cutoff_exceeded
        return (
            total=Inf,
            J_len=0.0,
            J_len_norm=0.0,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=stats.min_clearance,
            violation_count=stats.violation_count,
            len_ref=0.0,
            fuel_ref=0.0,
        )
    end
    refs = rpo_path_cost_normalization_refs(points, cfg)
    J_len = rpo_path_length(samples)
    J_len_norm = J_len / refs.len_ref
    partial_cost = cfg.w_obs * J_obs + cfg.w_len * J_len_norm^2
    if isfinite(cost_cutoff) && partial_cost > Float64(cost_cutoff)
        return (
            total=Inf,
            J_len=J_len,
            J_len_norm=J_len_norm,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=stats.min_clearance,
            violation_count=stats.violation_count,
            len_ref=refs.len_ref,
            fuel_ref=refs.fuel_ref,
        )
    end
    J_fuel = rpo_fuel_proxy_from_samples(samples, cfg)
    J_fuel_norm = J_fuel / refs.fuel_ref
    return (
        total=cfg.w_len * J_len_norm^2 + cfg.w_obs * J_obs + cfg.w_fuel * J_fuel_norm^2,
        J_len=J_len,
        J_len_norm=J_len_norm,
        J_obs=J_obs,
        J_fuel=J_fuel,
        J_fuel_norm=J_fuel_norm,
        min_clearance=stats.min_clearance,
        violation_count=stats.violation_count,
        len_ref=refs.len_ref,
        fuel_ref=refs.fuel_ref,
    )
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
