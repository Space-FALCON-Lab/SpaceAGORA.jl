include("../physical_models/Attitude_control_models.jl")
include("quaternion_utils.jl")
using Rotations, LinearAlgebra
using ControlSystemsBase
using MatrixEquations

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

function constant_α_β(m, b::config.Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int)
"""
    Generate a constant attitude control plan with fixed angles α and β set to 0.
    
    # Arguments
    - `m`: Mission object containing spacecraft and other parameters.
    - `t0`: Initial time for the control plan.
    - `args`: Simulation parameters.
    
    # Returns
    - A tuple containing the time vector, attitude vector, and angular velocity vector.
"""
    # Rot = config.rotate_to_inertial(m.body, b, root_index)

    # Calculate the wind-relative velocity in the inertial frame
    wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
    # Calculate the orientation quaternion from the inertial x-axis to the wind-relative velocity
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
    error_quat = error_quaternion(SVector{4, Float64}(m.body.roots[root_index].q), orientation_quat)
    δq = error_quat[1:3]
    δω = b.ω
    R = config.rotate_to_inertial(m.body, b, root_index)
    inertia_tensor = R * config.get_inertia_tensor(m.body, b) * R'
    kp = 0.5*100
    kd = 1.0*100
    b.rw_τ = R' * (cross(b.ω, inertia_tensor * b.ω) - kp*δq - kd*b.ω)
    ω_wheel_derivatives = MVector{m.body.n_reaction_wheels, Float64}(zeros(m.body.n_reaction_wheels))
    # println("b.rw_τ: $(size(b.rw_τ))")
    # println("b.J_rw: $(size((1/b.J_rw[:, 1])))")
    # println("r: $(size(R))")
    for i in 1:m.body.n_reaction_wheels
        ω_wheel_derivatives[i] = ((1/b.J_rw[:,i]) * R' * b.rw_τ)[1]
    end
    return ω_wheel_derivatives
end

function lqr_constant_α_β(m, b::config.Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int)
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

        # Calculate the wind-relative velocity in the inertial frame
        wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
        # Calculate the orientation quaternion from the inertial x-axis to the wind-relative velocity
        orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
        error_quat = error_quaternion(SVector{4, Float64}(m.body.roots[root_index].q), orientation_quat)
        # println("error_quat: $(error_quat)")
        # G = (q) -> [q[4] -q[3] q[2];
        # q[3] q[4] -q[1];
        # -q[2] q[1] q[4]]
        
        n = b.gyro
        state = SVector{6+n, Float64}([error_quat[1:3]; b.ω; b.rw])

        J = Rot * config.get_inertia_tensor(m.body, root_index) * Rot'
        A = SMatrix{length(state), length(state)}(Float64[zeros(3, 3) 0.5*I(3) zeros(3, n);
                                                    zeros(3+n, 6+n)]) - 1e-3*I 

        J_rw_inertial = MMatrix{3, n, Float64}(zeros(3, n))
        @inbounds for i in 1:n
            J_rw_inertial[:, i] .= Rot * b.J_rw[:, i]
        end
        B = SMatrix{length(state), n, Float64}([zeros(Float64, 3, n); J\J_rw_inertial; I(n)])# + 1e-6*ones(length(state), n)

        Q = Diagonal(SVector{6+n, Float64}([1.0e5*ones(3); 1.0*ones(3); 1.0e-6*ones(n)]))

        R = SMatrix{n, n, Float64}(0.1*I(n))
        if config.cnf.P == zeros(3, 3) || acosd(error_quat[4])*2 > 0.5 # if the error is larger than 0.5 degrees, recalculate using full LQR method
        #     # Use LQR to compute the gain matrix K
            K = SMatrix{n, length(state), Float64}(lqr(A, B, Q, R))
            # println(typeof(K))
            config.cnf.P = pinv(B')*R*K
            P = SMatrix{length(state), length(state), Float64}(config.cnf.P)
        # elseif acosd(error_quat[4])*2 < 0.1 # 2 degrees
        #     # Use LQR to compute the gain matrix K
        #     P = config.cnf.P
        else
            P = solve_care_newton(A, B, Q, R; P0=config.cnf.P, max_iter=100, tol=1e-8)
            config.cnf.P .= P
            # config.cnf.K .= R \ B' * config.cnf.P
        end
        # K = lqr(A, B, Q, R)
        # # Return the wheel momentum derivatives
        return -R \ B' * P * state
        # return -K * state
end

function lqr_constant_α_β_σ(m, b::config.Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int)
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

    # Calculate the wind-relative velocity in the inertial frame
    wind_relative_velocity = m.planet.L_PI' * vel_pp_rw
    # Calculate the orientation quaternion from the inertial x-axis to the wind-relative velocity
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), wind_relative_velocity)
    error_quat = error_quaternion(SVector{4, Float64}(m.body.roots[root_index].q), orientation_quat)

    # Error quaternion for bank angle
    h_ii_hat = m.planet.L_PI' * h_pp_hat
    orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), h_ii_hat)
    sc_orientation_quat = rotation_between(SVector{3, Float64}([1.0, 0.0, 0.0]), Rot*SVector{3, Float64}([0.0, 1.0, 0.0]))
    sc_orientation_error_quat = error_quaternion(sc_orientation_quat, orientation_quat)
    error_quat_2 = error_quaternion(error_quat, sc_orientation_error_quat)

    error_quat = SMatrix{4, 4, Float64}([Ψ(error_quat_2) error_quat_2])*error_quat
    n = b.gyro
    state = SVector{6+n, Float64}([error_quat[1:3]; b.ω; b.rw])

    J = Rot * config.get_inertia_tensor(m.body, root_index) * Rot'
    A = SMatrix{length(state), length(state)}(Float64[zeros(3, 3) 0.5*I(3) zeros(3, n);
                                                zeros(3+n, 6+n)]) - 1e-3*I 

    J_rw_inertial = MMatrix{3, n, Float64}(zeros(3, n))
    @inbounds for i in 1:n
        J_rw_inertial[:, i] .= Rot * b.J_rw[:, i]
    end
    B = SMatrix{length(state), n, Float64}([zeros(Float64, 3, n); J\J_rw_inertial; I(n)])# + 1e-6*ones(length(state), n)

    Q = Diagonal(SVector{6+n, Float64}([1.0e5*ones(3); 1.0*ones(3); 1.0e-6*ones(n)]))

    R = SMatrix{n, n, Float64}(0.1*I(n))
    if config.cnf.P == zeros(3, 3)
    #     # Use LQR to compute the gain matrix K
        K = SMatrix{n, length(state), Float64}(lqr(A, B, Q, R))
        # println(typeof(K))
        config.cnf.P = pinv(B')*R*K
        P = SMatrix{length(state), length(state), Float64}(config.cnf.P)
    # elseif acosd(error_quat[4])*2 < 0.1 # 2 degrees
    #     # Use LQR to compute the gain matrix K
    #     P = config.cnf.P
    else
        P = solve_care_newton(A, B, Q, R; P0=config.cnf.P, max_iter=100, tol=1e-6)
        config.cnf.P .= P
        # config.cnf.K .= R \ B' * config.cnf.P
    end
    # K = lqr(A, B, Q, R) #ControlSystemsBase.Discrete,
    # # Return the wheel momentum derivatives
    return -R \ B' * P * state
end

function solve_care_newton(A::AbstractMatrix, B::AbstractMatrix, Q::AbstractMatrix, R::AbstractMatrix;
                          P0::Union{AbstractMatrix, Nothing}=nothing,
                          max_iter::Int=100, tol::Float64=1e-8,
                          verbose::Bool=false)
    @fastmath begin
        n = size(A, 1)

        # Pre-compute common terms
        BRinvBt = SMatrix{size(B, 1), size(B, 1), Float64}(B * (R \ B')) # B * inv(R) * B'

        # Initial guess for P
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
            # 1. Compute the residual F(P_k) = A'P_k + P_k A - P_k BRinvBt P_k + Q
            F_Pk .= A'P_k + P_k * A - P_k * BRinvBt * P_k + Q
            residual_norm = norm(F_Pk)

            # 2. Check for convergence
            if residual_norm < tol
                converged = true
                break
            end

            # 3. Compute A_k = A - BRinvBt P_k
            Ak .= A - BRinvBt * P_k

            # 4. Solve the Lyapunov equation: Ak' ΔP + ΔP Ak = -F(P_k) for ΔP
            # MatrixEquations.lyap solves AX + XB = C
            # Here: A -> Ak', X -> ΔP, B -> Ak, C -> -F_Pk
            Delta_P .= MatrixEquations.lyapc(Ak', F_Pk)
            if !issymmetric(Delta_P)
                Delta_P .= Symmetric(Delta_P + Delta_P') / 2 # Ensure symmetry for ΔP
            end


            # 5. Update P_k = P_k + ΔP
            P_k .+= Delta_P
            if !issymmetric(P_k)
                P_k .= Symmetric(P_k + P_k') / 2 # Re-symmetrize P_k to combat numerical errors
            end
        end
    end
    return SMatrix{n, n, Float64}(P_k)
end
