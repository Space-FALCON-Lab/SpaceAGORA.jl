# Attitude Determination Subsystem -- Multiplicative Extended Kalman Filter
# DRM 1 | Requirements: ADS-FLT-001 through ADS-FLT-003
# State: [d_theta (3, attitude error), d_b (3, gyro bias error)] = 6 states

using LinearAlgebra

mutable struct ADSMEKFState
    q_est::Vector{Float64}
    gyro_bias::Vector{Float64}
    P::Matrix{Float64}
end

function ads_mekf_init(;
    q0::AbstractVector = [1.0,0.0,0.0,0.0],
    bias0::AbstractVector = zeros(3),
    P0::AbstractMatrix = 1e-4*Matrix{Float64}(I,6,6))
    ADSMEKFState(normalize(q0), collect(Float64,bias0), collect(Float64,P0))
end

function ads_mekf_predict!(state::ADSMEKFState, omega_meas::AbstractVector, dt::Float64, Q::AbstractMatrix)
    omega = omega_meas .- state.gyro_bias
    theta = norm(omega)*dt
    if theta > 1e-12
        axis = omega./norm(omega)
        dq = [cos(theta/2); axis.*sin(theta/2)]
        state.q_est = normalize(quat_multiply(state.q_est, dq))
    end
    Om = skew(omega)
    F = [Matrix{Float64}(I,3,3)-Om*dt  -dt*Matrix{Float64}(I,3,3);
         zeros(3,3)                     Matrix{Float64}(I,3,3)]
    state.P = F*state.P*F' + Q
end

function ads_mekf_update_star_tracker!(state::ADSMEKFState, q_meas::AbstractVector, R_st::AbstractMatrix)
    q_inv = [state.q_est[1]; -state.q_est[2:4]]
    dq = quat_multiply(q_inv, q_meas)
    sign_w = dq[1] >= 0 ? 1.0 : -1.0
    z = 2.0*sign_w*dq[2:4]
    H = [Matrix{Float64}(I,3,3)  zeros(3,3)]
    _kf_update!(state, z, H, R_st)
end

function ads_mekf_update_vector!(state::ADSMEKFState, b_meas::AbstractVector, b_ref_I::AbstractVector, R_vec::AbstractMatrix)
    b_pred = quat_rotate(state.q_est, normalize(b_ref_I))
    z = b_meas .- b_pred
    H = [-skew(b_pred)  zeros(3,3)]
    _kf_update!(state, z, H, R_vec)
end

function _kf_update!(state::ADSMEKFState, z::AbstractVector, H::AbstractMatrix, R::AbstractMatrix)
    S = H*state.P*H' + R
    K = state.P*H'*(S\Matrix{Float64}(I,size(S)...))
    dx = K*z
    dphi = dx[1:3]
    angle = norm(dphi)
    if angle > 1e-12
        dq = quat_from_axis_angle(dphi./angle, angle)
        state.q_est = normalize(quat_multiply(state.q_est, dq))
    end
    state.gyro_bias .+= dx[4:6]
    I6 = Matrix{Float64}(I,6,6)
    IKH = I6 - K*H
    state.P = IKH*state.P*IKH' + K*R*K'
end

function compute_nees(state::ADSMEKFState, q_true::AbstractVector)
    q_inv = [state.q_est[1]; -state.q_est[2:4]]
    dq = quat_multiply(q_inv, q_true)
    sign_w = dq[1] >= 0 ? 1.0 : -1.0
    err = 2.0*sign_w*dq[2:4]
    P_att = state.P[1:3,1:3]
    err'*(P_att\err)
end

attitude_error_deg(s::ADSMEKFState, q_true::AbstractVector) = quat_angle_deg(s.q_est, q_true)

function skew(v::AbstractVector)
    [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
end
