"""
Convert a geometric RPO path into a retimed reference trajectory and plan object.

With `retime_accel_limit_enable` the reference comes from the
acceleration-limited profile (`rpo_retimed_reference`): velocities are the
profile speed along the curve tangent, the first sample moves at
`retime_initial_speed_mps` and the last is the goal at rest. Otherwise
velocities are forward differences and the last one repeats the previous.
"""
function rpo_reference_from_path(path_rtn, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0)
    if cfg.retime_accel_limit_enable
        ref = rpo_retimed_reference(path_rtn, geometry, cfg; safe_distance_m=safe_distance_m)
        return ref.t_s, ref.r_rtn, ref.v_rtn
    end
    r_ref, _, _ = rpo_retime_path(path_rtn, geometry, cfg; safe_distance_m=safe_distance_m)
    n = size(r_ref, 2)
    v_ref = zeros(3, n)
    if n > 1
        @inbounds for j in 1:(n - 1)
            v_ref[:, j] .= (r_ref[:, j + 1] - r_ref[:, j]) / cfg.retime_dt_s
        end
        v_ref[:, end] .= v_ref[:, end - 1]
    end
    t_ref = collect(0.0:cfg.retime_dt_s:(cfg.retime_dt_s * (n - 1)))
    return t_ref, r_ref, v_ref
end
