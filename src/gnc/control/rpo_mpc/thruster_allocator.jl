function rpo_allocate_six_axis_thrusters(force_body, thrusters::SixAxisThrusterModel)
    f = zeros(6)
    desired = SVector{3, Float64}(force_body)
    @inbounds for axis in 1:3
        pos_idx = 2 * axis - 1
        neg_idx = 2 * axis
        if desired[axis] >= 0.0
            f[pos_idx] = min(desired[axis], thrusters.max_thrust_n[pos_idx])
        else
            f[neg_idx] = min(-desired[axis], thrusters.max_thrust_n[neg_idx])
        end
    end
    return SVector{6, Float64}(f)
end

function rpo_thruster_wrench_body(thruster_forces, thrusters::SixAxisThrusterModel)
    F = SVector{3, Float64}(0.0, 0.0, 0.0)
    τ = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for j in 1:6
        dir = SVector{3, Float64}(thrusters.directions_body[:, j])
        loc = SVector{3, Float64}(thrusters.locations_body[:, j])
        fj = Float64(thruster_forces[j])
        F += fj * dir
        τ += cross(loc, fj * dir)
    end
    return F, τ
end

