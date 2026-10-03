
"""Fixed step at which `:manuscript` mode evaluates each candidate's retimed reference."""
rpo_fuel_proxy_dt_s(cfg::RPOPSOConfig) = cfg.fuel_proxy_dt_s > 0.0 ? cfg.fuel_proxy_dt_s : cfg.retime_dt_s

"""
Objective components of the `:manuscript` mode for one candidate: Eq. 5,
J = w_obs J_obs + w_fuel J_fuel, with J_obs from Eq. 6 (sigmoid centred at
d_safe + τ_tol) over the adaptive samples and J_fuel the HCW fuel proxy of
Sec. III.B on the candidate's retimed reference at `rpo_fuel_proxy_dt_s`.
The candidate is retimed by the same retimer as the final reference; with
`retime_accel_limit_enable` the samples used for J_obs are reused for it.
There is no length term: `J_len` and the normalized values are reported only.
A violation is a sample below the Eq. 6 threshold; `keepout_violation_count`
counts samples inside the keep-out surface itself.
"""
function rpo_manuscript_path_cost_components(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    cost_cutoff::Real=Inf,
)
    samples, params, clearances = rpo_sample_path_with_params(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_hypr_sampling_density_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    threshold = rpo_obstacle_sigmoid_threshold(safe_distance_m, cfg.obstacle_sigmoid_tol_m, :manuscript)
    k = Float64(cfg.obstacle_sigmoid_k)
    w_obs = Float64(cfg.w_obs)
    cutoff = Float64(cost_cutoff)
    min_clearance = Inf
    violation_count = 0
    keepout_violation_count = 0
    J_obs = 0.0
    cut = false
    @inbounds for j in 1:size(samples, 2)
        c = clearances[j]
        if isnan(c)
            c = rpo_clearance_distance_to_station(SVector{3, Float64}(samples[1, j], samples[2, j], samples[3, j]), geometry)
            clearances[j] = c
        end
        min_clearance = min(min_clearance, c)
        J_obs += rpo_obstacle_sigmoid_penalty(c, threshold, k)
        c < threshold && (violation_count += 1)
        c < 0.0 && (keepout_violation_count += 1)
        if w_obs > 0.0 && isfinite(cutoff) && w_obs * J_obs > cutoff
            cut = true
            break
        end
    end
    refs = rpo_path_cost_normalization_refs(points, cfg)
    J_len = rpo_path_length(samples)
    if cut
        return (
            total=Inf,
            J_len=J_len,
            J_len_norm=J_len / refs.len_ref,
            J_obs=J_obs,
            J_fuel=0.0,
            J_fuel_norm=0.0,
            min_clearance=min_clearance,
            violation_count=violation_count,
            len_ref=refs.len_ref,
            fuel_ref=refs.fuel_ref,
            keepout_violation_count=keepout_violation_count,
            delta_v_eq_mps=0.0,
            reference_duration_s=0.0,
        )
    end
    dt = rpo_fuel_proxy_dt_s(cfg)
    fuel = if cfg.retime_accel_limit_enable
        profile = rpo_retime_profile(
            RPORetimeCurve(points, cfg.curve_type),
            samples,
            params,
            clearances,
            geometry,
            cfg;
            safe_distance_m=safe_distance_m,
            warn=false,
        )
        rpo_profile_hcw_fuel_proxy(profile, dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
    else
        step_cfg = dt == cfg.retime_dt_s ? cfg : rpo_pso_config(cfg; retime_dt_s=dt)
        r_ref, _, _ = Logging.with_logger(Logging.NullLogger()) do
            rpo_retime_path(points, geometry, step_cfg; safe_distance_m=safe_distance_m)
        end
        proxy = rpo_hcw_fuel_proxy(r_ref, dt, cfg.mean_motion_radps, cfg.mass_kg, cfg.isp_s, cfg.g0_mps2)
        (J_fuel=proxy.J_fuel, delta_v_eq_mps=proxy.delta_v_eq_mps, steps=size(r_ref, 2) - 1, duration_s=(size(r_ref, 2) - 1) * dt)
    end
    return (
        total=w_obs * J_obs + Float64(cfg.w_fuel) * fuel.J_fuel,
        J_len=J_len,
        J_len_norm=J_len / refs.len_ref,
        J_obs=J_obs,
        J_fuel=fuel.J_fuel,
        J_fuel_norm=fuel.J_fuel / refs.fuel_ref,
        min_clearance=min_clearance,
        violation_count=violation_count,
        len_ref=refs.len_ref,
        fuel_ref=refs.fuel_ref,
        keepout_violation_count=keepout_violation_count,
        delta_v_eq_mps=fuel.delta_v_eq_mps,
        reference_duration_s=fuel.duration_s,
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
    if cfg.hypr_mode === :manuscript
        return rpo_manuscript_path_cost_components(
            points,
            geometry,
            cfg;
            safe_distance_m=safe_distance_m,
            cost_cutoff=cost_cutoff,
        )
    end
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
