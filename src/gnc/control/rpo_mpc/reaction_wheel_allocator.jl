"""Return body torque commanded through the RPO reaction wheel allocator."""
function rpo_reaction_wheel_torque_command(model::RPOMPCControlModel, u, sat_idx::Int)
    if model.attitude_kp == 0.0 && model.rate_kd == 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    q = SVector{4, Float64}(u.sc[sat_idx].q)
    ω = SVector{3, Float64}(u.sc[sat_idx].ω)
    q_vec = q[4] < 0.0 ? -q[1:3] : q[1:3]
    τ = -model.attitude_kp * q_vec - model.rate_kd * ω
    limit = model.max_rw_torque_nm
    if isfinite(limit)
        τ = clamp.(τ, -limit, limit)
    end
    return SVector{3, Float64}(τ)
end
