# Existing HYPR signatures preserve caller defaults and forward only sampling inputs.
@inline function _rpo_sampling_settings(cfg::RPOPSOConfig)
    return RPOAdaptiveSamplingSettings(
        cfg.adaptive_sampling_enable,
        cfg.adaptive_sampling_max_ds_m,
        cfg.adaptive_sampling_far_clearance_m,
        cfg.adaptive_sampling_power,
        cfg.adaptive_sampling_safe_distance_fraction,
        cfg.adaptive_sampling_obstacle_guard_fraction,
    )
end

"""Apply the existing HYPR sampling defaults through the shared sampling kernel."""
function rpo_adaptive_sampling_min_ds_m(
    base_ds::Real,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
)
    return rpo_adaptive_sampling_min_ds_m(base_ds, geometry, _rpo_sampling_settings(cfg); safe_distance_m=safe_distance_m)
end

"""Apply the existing HYPR sampling defaults through the shared sampling kernel."""
function rpo_sample_path_polyline_adaptive(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    base_ds_m::Real=cfg.sample_ds_m,
)
    return rpo_sample_path_polyline_adaptive(
        points, geometry, _rpo_sampling_settings(cfg);
        safe_distance_m=safe_distance_m, base_ds_m=base_ds_m,
    )
end

"""Apply the existing HYPR sampling defaults through the shared sampling kernel."""
function rpo_sample_path_bezier_adaptive_with_params(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    base_ds_m::Real=cfg.sample_ds_m,
)
    return rpo_sample_path_bezier_adaptive_with_params(
        points, geometry, _rpo_sampling_settings(cfg);
        safe_distance_m=safe_distance_m, base_ds_m=base_ds_m,
    )
end

"""Apply the existing HYPR sampling defaults through the shared sampling kernel."""
function rpo_sample_path_bezier_adaptive(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    base_ds_m::Real=cfg.sample_ds_m,
)
    return rpo_sample_path_bezier_adaptive(
        points, geometry, _rpo_sampling_settings(cfg);
        safe_distance_m=safe_distance_m, base_ds_m=base_ds_m,
    )
end

"""Apply the existing HYPR sampling defaults through the shared sampling kernel."""
function rpo_sample_path_with_params(
    points,
    cfg::RPOPSOConfig,
    geometry;
    safe_distance_m::Real=cfg.safe_distance_m,
    base_ds_m::Real=cfg.sample_ds_m,
    curve_type::Symbol=cfg.curve_type,
)
    return rpo_sample_path_with_params(
        points, _rpo_sampling_settings(cfg), geometry;
        safe_distance_m=safe_distance_m, base_ds_m=base_ds_m, curve_type=curve_type,
    )
end

"""Apply the existing HYPR sampling defaults through the shared sampling kernel."""
function rpo_sample_path(
    points,
    cfg::RPOPSOConfig,
    geometry;
    safe_distance_m::Real=cfg.safe_distance_m,
    base_ds_m::Real=cfg.sample_ds_m,
    curve_type::Symbol=cfg.curve_type,
)
    return rpo_sample_path(
        points, _rpo_sampling_settings(cfg), geometry;
        safe_distance_m=safe_distance_m, base_ds_m=base_ds_m, curve_type=curve_type,
    )
end
