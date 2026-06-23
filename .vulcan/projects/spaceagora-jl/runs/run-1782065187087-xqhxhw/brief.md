# Implement MEKF estimator (ads_mekf.jl)

Create src/gnc/navigation/ads_mekf.jl implementing a Multiplicative Extended Kalman Filter.

Requirements: ADS-FLT-001, ADS-FLT-002, ADS-FLT-003, ADS-SYS-002, ADS-SYS-005, ADS-SYS-006.

State: ADSMEKFState struct with fields: q_est (Quaternion), gyro_bias (SVector{3}), P (SMatrix{6,6}).
The 6-state error vector is [delta_phi (3, attitude error via MRP), delta_b (3, gyro bias error)].

Functions required:

ads_mekf_init(q0, bias0, P0) -> ADSMEKFState
  Initialize filter state.

ads_mekf_predict!(state::ADSMEKFState, omega_meas::SVector{3}, dt::Float64, Q::SMatrix{6,6})
  Gyro-based propagation:
  omega_corrected = omega_meas - state.gyro_bias
  q_new = q_est * delta_q(omega_corrected * dt)   [use axis-angle formula]
  F = linearized kinematics Jacobian (6x6)
  P_new = F * P * F' + Q

ads_mekf_update_star_tracker!(state::ADSMEKFState, q_meas::Quaternion, R_st::SMatrix{3,3})
  H = [I3 | 0_{3x3}]  (attitude block only)
  innovation = quat_error(q_meas, q_est) as 3-vec MRP
  Standard KF update: K = P*H'*(H*P*H'+R)^-1; delta_x = K*inn; P = (I-K*H)*P
  Apply correction: q_est = q_est * exp(delta_x[1:3]); gyro_bias += delta_x[4:6]

ads_mekf_update_vector!(state::ADSMEKFState, b_meas::SVector{3}, b_ref_I::SVector{3}, R_vec::SMatrix{3,3})
  For magnetometer and sun sensor measurements (unit vectors).
  b_ref_B = R(q_est) * b_ref_I   [rotate reference to body frame]
  H = [skew(b_ref_B) | 0_{3x3}]  [sensitivity matrix]
  innovation = b_meas - b_ref_B
  Standard KF update as above.

ads_mekf_reset_error_state!(state::ADSMEKFState)
  After update: delta_phi is reset to zero (q_est already updated). Reset P accordingly.