

"""
Distance available for the stopping-distance limit of the retimer.

In `:manuscript` mode this is Sec. III.E's d_avail = max(0, c - d_safe), where c
is the clearance to the keep-out surface. The legacy retimer uses the raw
distance to the nearest station point.
"""
@inline function rpo_retime_available_distance(cfg::RPOPSOConfig, clearance::Real, distance::Real, safe_distance_m::Real)
    cfg.hypr_mode === :manuscript && return max(0.0, Float64(clearance) - Float64(safe_distance_m))
    return max(0.0, Float64(distance))
end

"""
Pointwise retiming speed limit of Sec. III.E at one sample: the clearance-limited
speed from the stopping-distance inequality and the curvature-limited speed,
scaled by the safety factor k_v. In `:manuscript` mode the global cap enters
the minimum before scaling, v = k_v min(v_clear, v_curv, v_max); the legacy
retimer applies the cap after scaling.
"""
@inline function rpo_retime_pointwise_speed(cfg::RPOPSOConfig, d_avail::Real, curvature::Real)
    amax = Float64(cfg.retime_a_max_mps2)
    reaction_time = Float64(cfg.retime_reaction_time_s)
    d = Float64(d_avail)
    κ = Float64(curvature)
    v_clear = if d <= 0.0 || amax <= 0.0
        0.0
    else
        -amax * reaction_time + sqrt((amax * reaction_time)^2 + 2.0 * amax * d)
    end
    v_clear = max(0.0, v_clear)
    v_curve = κ <= 0.0 || amax <= 0.0 ? Inf : sqrt(amax / κ)
    if cfg.hypr_mode === :manuscript
        return Float64(cfg.retime_speed_scale) * min(v_clear, v_curve, Float64(cfg.retime_max_speed_mps))
    end
    v = Float64(cfg.retime_speed_scale) * min(v_clear, v_curve)
    max_speed = Float64(cfg.retime_max_speed_mps)
    isfinite(max_speed) && (v = min(v, max_speed))
    return v
end

"""Collision-sampling spacing the retimer uses: the planner's density in `:manuscript` mode, `sample_ds_m` otherwise."""
function rpo_retime_sampling_ds_m(cfg::RPOPSOConfig, safe_distance_m::Real)
    cfg.hypr_mode === :manuscript && return rpo_hypr_sampling_density_m(cfg, safe_distance_m)
    return cfg.sample_ds_m
end

"""Retiming path samples into position, velocity, acceleration, and time references."""
function rpo_retime_path(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
)
    if cfg.retime_accel_limit_enable
        ref = rpo_retimed_reference(
            points,
            geometry,
            cfg;
            safe_distance_m=safe_distance_m,
            fallback_speed_mps=fallback_speed_mps,
            duplicate_tol_m=duplicate_tol_m,
        )
        return ref.r_rtn, ref.s_m, ref.speed_mps
    end
    raw_samples = rpo_sample_path(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_retime_sampling_ds_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    return rpo_retime_samples(
        raw_samples, geometry;
        max_speed_mps=cfg.retime_max_speed_mps,
        min_speed_mps=cfg.retime_min_speed_mps,
        dt_s=cfg.retime_dt_s,
        max_steps=cfg.retime_max_steps,
        available_distance=(clearance, distance, safe) -> rpo_retime_available_distance(cfg, clearance, distance, safe),
        pointwise_speed=(distance, curvature) -> rpo_retime_pointwise_speed(cfg, distance, curvature),
        safe_distance_m=safe_distance_m,
        fallback_speed_mps=fallback_speed_mps,
        duplicate_tol_m=duplicate_tol_m,
    )
end

"""Build the shared profile with the existing HYPR-configured speed policy."""
function rpo_retime_profile(
    curve::RPORetimeCurve,
    samples,
    params,
    clearances,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
    warn::Bool=true,
)
    return rpo_retime_profile(
        curve, samples, params, clearances, geometry;
        max_speed_mps=cfg.retime_max_speed_mps,
        min_speed_mps=cfg.retime_min_speed_mps,
        initial_speed_mps=cfg.retime_initial_speed_mps,
        a_max_mps2=cfg.retime_a_max_mps2,
        available_distance=(clearance, distance, safe) -> rpo_retime_available_distance(cfg, clearance, distance, safe),
        pointwise_speed=(distance, curvature) -> rpo_retime_pointwise_speed(cfg, distance, curvature),
        safe_distance_m=safe_distance_m,
        fallback_speed_mps=fallback_speed_mps,
        duplicate_tol_m=duplicate_tol_m,
        warn=warn,
    )
end

"""
    rpo_retimed_reference(points, geometry, cfg; safe_distance_m=0.0)

Acceleration-limited reference for a control polygon: sample the path as the
planner does, build `rpo_retime_profile`, and sample it at `retime_dt_s`.
"""
function rpo_retimed_reference(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
    warn::Bool=true,
)
    samples, params, clearances = rpo_sample_path_with_params(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_retime_sampling_ds_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    profile = rpo_retime_profile(
        RPORetimeCurve(points, cfg.curve_type),
        samples,
        params,
        clearances,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        fallback_speed_mps=fallback_speed_mps,
        duplicate_tol_m=duplicate_tol_m,
        warn=warn,
    )
    return rpo_retimed_reference_from_profile(profile, cfg.retime_dt_s)
end
