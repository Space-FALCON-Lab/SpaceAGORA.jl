include(joinpath(@__DIR__, "..", "..", "gnc", "control", "attitude_control_models.jl"))
include(joinpath(@__DIR__, "..", "..", "core", "numerics", "quaternion_utils.jl"))
using Rotations, LinearAlgebra
using ControlSystemsBase
using MatrixEquations
using CSV
using .SimulationModel: Link

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

function constant_α_β(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int)
"""
    Generate a constant attitude control plan with fixed angles α and β set to 0.
    
    # Arguments
    - `m`: Mission object containing spacecraft and other parameters.
    - `t0`: Initial time for the control plan.
    - `args`: Simulation parameters.
    
    # Returns
    - A tuple containing the time vector, attitude vector, and angular velocity vector.
"""
    wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
    error_quat = error_quaternion(SVector{4, Float64}(m.body.roots[root_index].q), orientation_quat)
    δq = error_quat[1:3]
    δω = b.ω
    R = rotate_to_inertial(m.body, b, root_index)
    inertia_tensor = R * get_inertia_tensor(m.body, b) * R'
    kp = 0.5*100
    kd = 1.0*100
    b.rw_τ = R' * (cross(b.ω, inertia_tensor * b.ω) - kp*δq - kd*b.ω)
    ω_wheel_derivatives = MVector{m.body.n_reaction_wheels, Float64}(zeros(m.body.n_reaction_wheels))
    for i in 1:m.body.n_reaction_wheels
        ω_wheel_derivatives[i] = ((1/b.J_rw[:,i]) * R' * b.rw_τ)[1]
    end
    return ω_wheel_derivatives
end

function lqr_constant_α_β(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int)
    """
    Generate a constant attitude control plan using LQR with fixed angles α and β set to 0.
    
    # Arguments
    - `m`: Mission object containing spacecraft and other parameters.
    - `t0`: Initial time for the control plan.
    - `args`: Simulation parameters.
    
    # Returns
    - A tuple containing the time vector, attitude vector, and angular velocity vector.
    """ 
        Rot = rotate_to_inertial(m.body, b, root_index)
        wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
        orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
        error_quat = error_quaternion(SVector{4, Float64}(m.body.roots[root_index].q), orientation_quat)
        
        n = b.gyro
        state = SVector{6+n, Float64}([error_quat[1:3]; b.ω; b.rw])

        J = Rot * get_inertia_tensor(m.body, root_index) * Rot'
        A = SMatrix{length(state), length(state)}(Float64[zeros(3, 3) 0.5*I(3) zeros(3, n);
                                                    zeros(3+n, 6+n)]) - 1e-3*I 

        J_rw_inertial = MMatrix{3, n, Float64}(zeros(3, n))
        @inbounds for i in 1:n
            J_rw_inertial[:, i] .= Rot * b.J_rw[:, i]
        end
        B = SMatrix{length(state), n, Float64}([zeros(Float64, 3, n); J\J_rw_inertial; I(n)])# + 1e-6*ones(length(state), n)

        Q = Diagonal(SVector{6+n, Float64}([1.0e5*ones(3); 1.0*ones(3); 1.0e-6*ones(n)]))

        R = SMatrix{n, n, Float64}(0.1*I(n))
        # Recompute from scratch when the cached Riccati solution is missing or the
        # spacecraft is still far from the target attitude; otherwise iterate from
        # the previous solution for a cheaper update.
        if cnf.P == zeros(3, 3) || acosd(error_quat[4])*2 > 0.5 # if the error is larger than 0.5 degrees, recalculate using full LQR method
            K = SMatrix{n, length(state), Float64}(lqr(A, B, Q, R))
            cnf.P = pinv(B')*R*K
            P = SMatrix{length(state), length(state), Float64}(cnf.P)
        else
            P = solve_care_newton(A, B, Q, R; P0=cnf.P, max_iter=100, tol=1e-8)
            cnf.P .= P
        end
        return -R \ B' * P * state
end

function lqr_constant_α_β_σ(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int)
    """
    Generate a constant attitude control plan using LQR with fixed angles α and β set to 0.
    
    # Arguments
    - `m`: Mission object containing spacecraft and other parameters.
    - `t0`: Initial time for the control plan.
    - `args`: Simulation parameters.
    
    # Returns
    - A tuple containing the time vector, attitude vector, and angular velocity vector.
    """
    Rot = rotate_to_inertial(m.body, b, root_index)

    wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
    error_quat = error_quaternion(SVector{4, Float64}(m.body.roots[root_index].q), orientation_quat)

    # Fold the lift-vector bank error into the same quaternion state used by the
    # nominal attitude LQR so one controller regulates both objectives.
    h_ii_hat = m.planet.L_PI' * h_pp_hat
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), h_ii_hat)
    sc_orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), Rot*SVector{3, Float64}([0.0, 1.0, 0.0]))
    sc_orientation_error_quat = error_quaternion(sc_orientation_quat, orientation_quat)
    error_quat_2 = error_quaternion(error_quat, sc_orientation_error_quat)

    error_quat = SMatrix{4, 4, Float64}([Ψ(error_quat_2) error_quat_2])*error_quat
    n = b.gyro
    state = SVector{6+n, Float64}([error_quat[1:3]; b.ω; b.rw])

    J = Rot * get_inertia_tensor(m.body, root_index) * Rot'
    A = SMatrix{length(state), length(state)}(Float64[zeros(3, 3) 0.5*I(3) zeros(3, n);
                                                zeros(3+n, 6+n)]) - 1e-3*I 

    J_rw_inertial = MMatrix{3, n, Float64}(zeros(3, n))
    @inbounds for i in 1:n
        J_rw_inertial[:, i] .= Rot * b.J_rw[:, i]
    end
    B = SMatrix{length(state), n, Float64}([zeros(Float64, 3, n); J\J_rw_inertial; I(n)])# + 1e-6*ones(length(state), n)

    Q = Diagonal(SVector{6+n, Float64}([1.0e5*ones(3); 1.0*ones(3); 1.0e-6*ones(n)]))

    R = SMatrix{n, n, Float64}(0.1*I(n))
    if cnf.P == zeros(3, 3)
        K = SMatrix{n, length(state), Float64}(lqr(A, B, Q, R))
        cnf.P = pinv(B')*R*K
        P = SMatrix{length(state), length(state), Float64}(cnf.P)
    else
        P = solve_care_newton(A, B, Q, R; P0=cnf.P, max_iter=100, tol=1e-6)
        cnf.P .= P
    end
    return -R \ B' * P * state
end

function solve_care_newton(A::AbstractMatrix, B::AbstractMatrix, Q::AbstractMatrix, R::AbstractMatrix;
                          P0::Union{AbstractMatrix, Nothing}=nothing,
                          max_iter::Int=100, tol::Float64=1e-8,
                          verbose::Bool=false)
    @fastmath begin
        n = size(A, 1)

        BRinvBt = SMatrix{size(B, 1), size(B, 1), Float64}(B * (R \ B')) # B * inv(R) * B'

        P_k = isnothing(P0) ? MMatrix{n, n, Float64}(zeros(n, n)) : MMatrix{n, n, Float64}(P0)
        if !issymmetric(P_k)
            P_k .= Symmetric(P_k + P_k') / 2 # Ensure symmetry
        end

        converged = false
        residual_norm = Inf
        F_Pk = MMatrix{n, n, Float64}(zeros(n, n))
        Ak = MMatrix{n, n, Float64}(zeros(n, n))
        Delta_P = MMatrix{n, n, Float64}(zeros(n, n))
        @inbounds for _ in 1:max_iter
            F_Pk .= A'P_k + P_k * A - P_k * BRinvBt * P_k + Q
            residual_norm = norm(F_Pk)

            if residual_norm < tol
                converged = true
                break
            end

            Ak .= A - BRinvBt * P_k

            # Newton-Kleinman step for the continuous-time algebraic Riccati equation.
            Delta_P .= MatrixEquations.lyapc(Ak', F_Pk)
            if !issymmetric(Delta_P)
                Delta_P .= Symmetric(Delta_P + Delta_P') / 2 # Ensure symmetry for ΔP
            end

            P_k .+= Delta_P
            if !issymmetric(P_k)
                P_k .= Symmetric(P_k + P_k') / 2 # Re-symmetrize P_k to combat numerical errors
            end
        end
    end
    return SMatrix{n, n, Float64}(P_k)
end

function constant_thruster!(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    for thruster in b.thrusters
        thruster.thrust = 1.0
    end
end

function rw_mrp_feedback_control!(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    """
    MRP feedback control for comparison with Basilisk
    """
    q = m.body.links[root_index].q
    ω = m.body.links[root_index].ω
    current_MRP = -q[1:3]/(1+q[4])
    if norm(current_MRP) > 1
        current_MRP .= -current_MRP/norm(current_MRP)^2 # If the rotation is larger than 180 degrees, switch to shadow set
    end
    
    # Determine the error rotation
    R = rotate_to_inertial(m.body, b, root_index)
    J = get_inertia_tensor(m.body, b)
    ω_body = R' * ω
    P = 30.0
    K = 3.5
    Lr = -K*current_MRP - P*ω_body + cross(ω_body, J*ω_body) # Total control torque, body frame
    b.ω_wheel_derivatives .= -pinv(b.J_rw)*Lr
end

function rw_polyfit_control!(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    rw1 = [-5.60895312e-36,  2.67405595e-32, -5.75833882e-29,  7.39858321e-26,
       -6.30995548e-23,  3.75879226e-20, -1.60171291e-17,  4.91432034e-15,
       -1.07698378e-12,  1.64628922e-10, -1.67735734e-08,  1.04916101e-06,
       -3.43892609e-05,  4.42626398e-04, -8.56999390e-03,  3.83306056e-01]
    rw2 = [ 1.12500939e-35, -5.19204031e-32,  1.07916760e-28, -1.33381746e-25,
        1.08997484e-22, -6.19243499e-20,  2.50278492e-17, -7.23567055e-15,
        1.48276059e-12, -2.10179799e-10,  1.97229165e-08, -1.13827583e-06,
        3.48996218e-05, -2.05866906e-04, -2.64800295e-02,  4.02809724e-01]
    rw3 = [ 8.71762704e-36, -3.92648983e-32,  7.94120288e-29, -9.51632719e-26,
        7.50680130e-23, -4.09401651e-20,  1.57679295e-17, -4.29964414e-15,
        8.18216619e-13, -1.04948256e-10,  8.50458803e-09, -3.93126674e-07,
        1.09634047e-05, -5.65654872e-04,  2.88211721e-02, -1.41382754e-01]
    function polyfit(coeffs, t)
        exponent = length(coeffs) - 1
        output = 0
        for i in eachindex(coeffs)
            output += coeffs[i] * t^exponent
            exponent -= 1
        end
        return output
    end
    b.ω_wheel_derivatives .= t < 375 ? [polyfit(rw1, t), polyfit(rw2, t), polyfit(rw3, t)] : [0.0, 0.0, 0.0]
end

function basilisk_rw_read_csv!(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    data = CSV.File("/workspaces/SpaceAGORA.jl/basilisk_rw_torque.csv", delim=',', header=true) |> DataFrame
    data = Matrix(data)
    times = data[:, 1]
    rw_values = data[:, 2:end]

    idx = findmin(abs.(times .- t))[2]
    if t < times[idx]
        idx -= 1
    end

    # Set the wheel momentum derivatives
    b.ω_wheel_derivatives .= rw_values[idx, :]
end

function basilisk_thruster_read_csv!(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    data = CSV.File("/workspaces/SpaceAGORA.jl/basilisk_thruster_force.csv", delim=',', header=true) |> DataFrame
    data = Matrix(data)
    times = data[:, 1]
    thruster_values = data[:, 2:end]
    idx = findmin(abs.(times .- t))[2]
    for (i, thruster) in enumerate(b.thrusters)
        thruster.thrust = thruster_values[idx, i]
        thruster.stop_firing_time = Inf
    end
end

function basilisk_thruster_torque_read_csv!(m, b::Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    data = CSV.File("/workspaces/SpaceAGORA.jl/basilisk_thruster_torque.csv", delim=',', header=true) |> DataFrame
    data = Matrix(data)
    times = data[:, 1]
    thruster_torques = data[:, 2:end]
    idx = findmin(abs.(times .- t))[2]
    if t < times[idx]
        idx -= 1
    end
    update_thrusters!(b, thruster_torques[idx, :], t)
end
