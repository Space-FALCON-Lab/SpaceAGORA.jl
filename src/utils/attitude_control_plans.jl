include("../physical_models/Attitude_control_models.jl")
include("quaternion_utils.jl")
using Rotations, LinearAlgebra
using ControlSystemsBase

# Function to compute the rotation quaternion between two vectors
function rotation_between(v1::SVector{3, Float64}, v2::SVector{3, Float64})
    v1 = normalize(v1)
    v2 = normalize(v2)
    dot_prod = dot(v1, v2)

    if isapprox(dot_prod, 1.0; atol=1e-8)
        return SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)  # identity rotation
    elseif isapprox(dot_prod, -1.0; atol=1e-8)
        # Choose orthogonal axis
        axis = normalize(cross(v1, SVector(1.0, 0.0, 0.0)))
        if norm(axis) < 1e-8
            axis = normalize(cross(v1, SVector(0.0, 1.0, 0.0)))
        end
        return SVector{4, Float64}(axis[1], axis[2], axis[3], 0.0)
    else
        c = cross(v1, v2)
        s = sqrt((1 + dot_prod) * 2)
        q = SVector{4, Float64}(c[1] / s, c[2] / s, c[3] / s, s / 2)
        return q
    end
end

function constant_α_β(m, t0, b, bodies, root_index, args, vel_pp_rw)
"""
    Generate a constant attitude control plan with fixed angles α and β set to 0.
    
    # Arguments
    - `m`: Mission object containing spacecraft and other parameters.
    - `t0`: Initial time for the control plan.
    - `args`: Simulation parameters.
    
    # Returns
    - A tuple containing the time vector, attitude vector, and angular velocity vector.
"""
    # α = 0.0  # Angle of attack in radians
    # β = 0.0  # Angle of sideslip in radians

    # R = config.rotate_to_inertial(m.body, b, root_index)
    # body_frame_velocity = R' * m.planet.L_PI' * vel_pp_rw
    # current_α = atan(body_frame_velocity[2], body_frame_velocity[1])
    # current_β = atan(body_frame_velocity[2], norm(body_frame_velocity[1:2:3]))
    # ω = b.ω
    # Calculate the body frame x axis in the inertial frame
    x_axis_inertial = rot(m.body.roots[root_index].q)' * SVector(1.0, 0.0, 0.0)
    # Calculate the wind-relative velocity in the inertial frame
    wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
    # Calculate the orientation quaternion from the inertial x-axis to the wind-relative velocity
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
    error_quat = error_quaternion(m.body.roots[root_index].q, orientation_quat)
    δq = error_quat[1:3]
    δω = b.ω
    R = config.rotate_to_inertial(m.body, b, root_index)
    inertia_tensor = R * config.get_inertia_tensor(m.body, b) * R'
    kp = 0.5*100
    kd = 1.0*100
    b.rw_τ = R' * (cross(b.ω, inertia_tensor * b.ω) - kp*δq - kd*b.ω)

    return R * b.rw_τ
end

function lqr_constant_α_β(m, t0::Float64, b::config.Link, bodies::Vector{config.Link}, root_index::Int, args, vel_pp_rw::SVector{3, Float64}, aerobraking_phase::Int)
    """
    Generate a constant attitude control plan using LQR with fixed angles α and β set to 0.
    
    # Arguments
    - `m`: Mission object containing spacecraft and other parameters.
    - `t0`: Initial time for the control plan.
    - `args`: Simulation parameters.
    
    # Returns
    - A tuple containing the time vector, attitude vector, and angular velocity vector.
    """
    # x_axis_inertial = rot(m.body.roots[root_index].q)' * SVector(1.0, 0.0, 0.0)
    Rot = config.rotate_to_inertial(m.body, b, root_index)
    was_lqr = false
    was_pid = false
    # Calculate the wind-relative velocity in the inertial frame
    wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
    # Calculate the orientation quaternion from the inertial x-axis to the wind-relative velocity
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
    error_quat = error_quaternion(m.body.roots[root_index].q, orientation_quat)
    G = (q) -> [q[4] -q[3] q[2];
                q[3] q[4] -q[1];
                -q[2] q[1] q[4]]
    
    n = b.gyro
    J = Rot * config.get_inertia_tensor(m.body, root_index) * Rot'
    state = SVector{6+n, Float64}([error_quat[1:3]; b.ω; b.rw])

    A = ([zeros(3, 3) 0.5*I(3) zeros(3, n);
        zeros(3+n, 6+n)]) - 1e-6*I #SMatrix{length(state), length(state)}

    J_rw_inertial = zeros(3, n)
    for i in 1:n
        J_rw_inertial[:, i] = Rot * b.J_rw[:, i]
    end
    B = SMatrix{length(state), n, Float64}([zeros(3, n); inv(J)*(J_rw_inertial); I])# + 1e-6*ones(length(state), n)

    Q = Diagonal(vcat(1e2*ones(3), ones(3), 1e-6*ones(n)))

    R = 1e-1*I(n)

    K = lqr(A, B, Q, R)#ControlSystemsBase.Discrete, 
    # Return the wheel momentum derivatives
    return -K * state
end