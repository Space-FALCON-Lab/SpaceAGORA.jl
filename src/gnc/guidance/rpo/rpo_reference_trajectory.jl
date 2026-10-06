"""Materialize a rest-to-rest reference with interior forward-difference velocities."""
function rpo_reference_from_profile(profile::RPORetimingProfile)
    r_ref, _, _ = rpo_retime_path_from_profile(profile)
    n = size(r_ref, 2)
    v_ref = zeros(3, n)
    if n > 1
        @inbounds for j in 2:(n - 1)
            v_ref[:, j] .= (r_ref[:, j + 1] - r_ref[:, j]) / profile.dt_s
        end
    end
    t_ref = collect(0.0:profile.dt_s:(profile.dt_s * (n - 1)))
    return t_ref, r_ref, v_ref
end


"""Convert a geometric RPO path into a retimed reference trajectory and plan object."""
function rpo_reference_from_path(path_rtn, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0)
    raw_samples = rpo_sample_path(
        path_rtn,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=cfg.sample_ds_m,
        curve_type=cfg.curve_type,
    )
    profile = rpo_retiming_profile_from_samples(raw_samples, geometry, cfg)
    return rpo_reference_from_profile(profile)
end
