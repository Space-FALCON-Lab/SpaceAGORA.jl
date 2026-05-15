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

function rpo_clearance_stats_from_samples(
    samples,
    geometry,
    safe_distance_m::Real;
    cost_cutoff::Real=Inf,
    w_obs::Real=0.0,
)
    min_clearance = Inf
    violation_count = 0
    safe = Float64(safe_distance_m)
    @inbounds for j in 1:size(samples, 2)
        p = SVector{3, Float64}(samples[1, j], samples[2, j], samples[3, j])
        clearance = rpo_clearance_distance_to_station(p, geometry)
        min_clearance = min(min_clearance, clearance)
        if clearance < safe
            violation_count += 1
            if w_obs > 0.0 && isfinite(cost_cutoff) && Float64(w_obs) * violation_count > Float64(cost_cutoff)
                return (
                    min_clearance=min_clearance,
                    violation_count=violation_count,
                    violation_fraction=violation_count / max(size(samples, 2), 1),
                    cutoff_exceeded=true,
                )
            end
        end
    end
    return (
        min_clearance=min_clearance,
        violation_count=violation_count,
        violation_fraction=violation_count / max(size(samples, 2), 1),
        cutoff_exceeded=false,
    )
end

function rpo_normalized_path_cost_components(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    cost_cutoff::Real=Inf,
)
    samples = rpo_sample_path(points, cfg.sample_ds_m; curve_type=cfg.curve_type)
    stats = rpo_clearance_stats_from_samples(
        samples,
        geometry,
        safe_distance_m;
        cost_cutoff=cost_cutoff,
        w_obs=cfg.w_obs,
    )
    J_obs = stats.violation_count
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

function rpo_path_cost(points, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0, cost_cutoff::Real=Inf)
    return rpo_normalized_path_cost_components(points, geometry, cfg; safe_distance_m=safe_distance_m, cost_cutoff=cost_cutoff).total
end
