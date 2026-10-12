"""Convert a geometric RPO path into a retimed reference trajectory and plan object."""
function rpo_reference_from_path(path_rtn, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0)
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

