# Attitude Determination and Control System — Controller
# DRM 2 | Requirements: ADCS-CTL-001 through ADCS-INT-002
#
# Quaternion feedback attitude control + B-dot momentum desaturation.
# Uses existing ReactionWheelAssembly and Magnet hardware.
# Integrates via calcControlEffect! / calcControlForceTorque pattern.

using LinearAlgebra

# ── Configuration ─────────────────────────────────────────────────────────────

"""Configuration for the ADCS closed-loop controller (DRM 2)."""
@kwdef struct ADCSConfig
    kp_fine::Float64      = 0.05     # Proportional gain, Nm/rad (FINE_TRACK mode)
    kd_fine::Float64      = 0.5      # Derivative gain, Nm·s/rad
    kp_coarse::Float64    = 0.02     # Proportional gain (COARSE_ACQ mode)
    kd_coarse::Float64    = 0.3
    kd_detumble::Float64  = 0.4      # Rate-only damping gain, Nm·s/rad
    k_desaturation::Float64  = 5e5   # Dipole gain, Am²·s/(Nms·T)
    desaturation_threshold::Float64 = 0.70  # Activate at 70% saturation
    desaturation_exit::Float64      = 0.30  # Exit at 30% saturation
    max_dipole_Am2::Float64         = 30.0  # Max magnetorquer dipole, Am²
    max_body_torque_nm::Float64     = 0.3   # Nm
    detumble_omega_threshold_rads::Float64 = 5.0 * pi / 180
    coarse_to_fine_error_deg::Float64      = 5.0
    update_rate_hz::Float64  = 1.0
    sat_idx::Int             = 1
end

# ── ADCS operational modes ────────────────────────────────────────────────────

@enum ADCSMode::Int8 ADCS_SAFE=1 ADCS_DETUMBLE=2 ADCS_COARSE_POINT=3 ADCS_FINE_POINT=4 ADCS_DESATURATING=5

# ── State and held actuation ──────────────────────────────────────────────────

mutable struct ADCSHeldActuation
    torque_body::Vector{Float64}   # 3-vector, body frame
end
ADCSHeldActuation() = ADCSHeldActuation(zeros(3))

mutable struct ADCSControllerState
    mode::ADCSMode
    q_ref::Vector{Float64}         # Target attitude quaternion [x,y,z,w]
end
ADCSControllerState() = ADCSControllerState(ADCS_SAFE, [0.0, 0.0, 0.0, 1.0])

# ── Main effector struct ──────────────────────────────────────────────────────

mutable struct ADCSController
    config::ADCSConfig
    state::ADCSControllerState
    held::ADCSHeldActuation
end
ADCSController(; config::ADCSConfig = ADCSConfig()) =
    ADCSController(config, ADCSControllerState(), ADCSHeldActuation())

# ── Quaternion utilities ──────────────────────────────────────────────────────

"""Quaternion product p⊗q (scalar-last convention [x,y,z,w])."""
function quat_multiply_adcs(p::AbstractVector, q::AbstractVector)
    px, py, pz, pw = p[1], p[2], p[3], p[4]
    qx, qy, qz, qw = q[1], q[2], q[3], q[4]
    return [pw*qx + px*qw + py*qz - pz*qy,
            pw*qy - px*qz + py*qw + pz*qx,
            pw*qz + px*qy - py*qx + pz*qw,
            pw*qw - px*qx - py*qy - pz*qz]
end

"""Quaternion inverse for unit quaternion."""
quat_inv_adcs(q::AbstractVector) = [-q[1], -q[2], -q[3], q[4]]

"""Rotation matrix from quaternion (scalar-last [x,y,z,w])."""
function quat_to_rotmat(q::AbstractVector)
    qx, qy, qz, qw = q[1], q[2], q[3], q[4]
    [1-2*(qy^2+qz^2)  2*(qx*qy+qw*qz)  2*(qx*qz-qw*qy);
     2*(qx*qy-qw*qz)  1-2*(qx^2+qz^2)  2*(qy*qz+qw*qx);
     2*(qx*qz+qw*qy)  2*(qy*qz-qw*qx)  1-2*(qx^2+qy^2)]
end

"""Attitude error angle in degrees."""
function quat_error_deg_adcs(q_est::AbstractVector, q_ref::AbstractVector)
    dq = quat_multiply_adcs(quat_inv_adcs(q_ref), q_est)
    return 2.0 * acos(clamp(abs(dq[4]), 0.0, 1.0)) * 180.0 / pi
end

# ── Mode management ───────────────────────────────────────────────────────────

function update_adcs_mode!(ctrl::ADCSController, q_est, omega_est,
                            ads_valid::Bool, ads_fine::Bool,
                            rw_sat_fraction::Float64, t::Float64)
    cfg = ctrl.config
    current = ctrl.state.mode
    if !ads_valid
        ctrl.state.mode = ADCS_SAFE; return
    end
    omega_mag = norm(omega_est)
    if current == ADCS_SAFE
        ctrl.state.mode = omega_mag > cfg.detumble_omega_threshold_rads ? ADCS_DETUMBLE : ADCS_COARSE_POINT
    elseif current == ADCS_DETUMBLE
        omega_mag < cfg.detumble_omega_threshold_rads / 5.0 && (ctrl.state.mode = ADCS_COARSE_POINT)
    elseif current == ADCS_COARSE_POINT
        if ads_fine && quat_error_deg_adcs(q_est, ctrl.state.q_ref) < cfg.coarse_to_fine_error_deg
            ctrl.state.mode = ADCS_FINE_POINT
        end
    elseif current == ADCS_FINE_POINT
        if rw_sat_fraction > cfg.desaturation_threshold
            ctrl.state.mode = ADCS_DESATURATING
        elseif !ads_fine
            ctrl.state.mode = ADCS_COARSE_POINT
        end
    elseif current == ADCS_DESATURATING
        if rw_sat_fraction < cfg.desaturation_exit
            ctrl.state.mode = ADCS_FINE_POINT
        elseif !ads_valid
            ctrl.state.mode = ADCS_SAFE
        end
    end
end

# ── Quaternion feedback torque ────────────────────────────────────────────────

function quaternion_feedback_torque(ctrl::ADCSController,
                                     q_est::AbstractVector,
                                     omega_est::AbstractVector) :: Vector{Float64}
    cfg = ctrl.config
    mode = ctrl.state.mode
    mode == ADCS_SAFE && return zeros(3)
    if mode == ADCS_DETUMBLE
        return clamp.(-cfg.kd_detumble .* collect(omega_est), -cfg.max_body_torque_nm, cfg.max_body_torque_nm)
    end
    # Quaternion error: q_err = q_ref^{-1} ⊗ q_est
    q_err = quat_multiply_adcs(quat_inv_adcs(ctrl.state.q_ref), collect(q_est))
    # Short-path: flip if scalar negative
    q_err_vec = q_err[4] >= 0 ? q_err[1:3] : -q_err[1:3]
    kp = (mode == ADCS_FINE_POINT || mode == ADCS_DESATURATING) ? cfg.kp_fine : cfg.kp_coarse
    kd = (mode == ADCS_FINE_POINT || mode == ADCS_DESATURATING) ? cfg.kd_fine : cfg.kd_coarse
    tau = -kp .* q_err_vec - kd .* collect(omega_est)
    return clamp.(tau, -cfg.max_body_torque_nm, cfg.max_body_torque_nm)
end

# ── Reaction wheel allocation ─────────────────────────���───────────────────────

"""
Allocate body torque to reaction wheels via pseudo-inverse of J_rw.
Returns net body torque. Mutates rw.h_dot_wheels and rw.tau_body_net.
"""
function allocate_to_reaction_wheels!(rw, tau_body::AbstractVector) :: Vector{Float64}
    J = Matrix(rw.J_rw)
    J_pinv = J' * inv(J * J')
    h_dot = J_pinv * tau_body
    h_dot_clamped = clamp.(h_dot, -rw.max_wheel_torque, rw.max_wheel_torque)
    rw.h_dot_wheels .= h_dot_clamped
    tau_net = J * h_dot_clamped
    rw.tau_body_net .= tau_net
    return tau_net
end

# ── RW saturation fraction ────────────────────────────────────────────────────

function rw_saturation_fraction(rw) :: Float64
    h_cap = rw.max_wheel_h * rw.n_wheels
    h_stored = norm(Matrix(rw.J_rw) * Vector(rw.h_wheels))
    return h_stored / max(h_cap, 1e-10)
end

# ── Momentum desaturation ─────────────────────────────────────────────────────

"""
B-dot / cross-product desaturation.
m_cmd = -k_mag * cross(B_body, h_rw_total), clamped to max_dipole.
Sets Magnet.m to command the magnetorquer.
"""
function compute_desaturation_dipole!(ctrl::ADCSController, magnet, rw,
                                       q_est::AbstractVector, r_ecef::AbstractVector)
    cfg = ctrl.config
    B_I_nT = get_magnetic_field_dipole(r_ecef)       # nT, inertial
    B_I_T  = collect(B_I_nT) .* 1e-9                 # → Tesla
    R_ItoB = quat_to_rotmat(collect(q_est))
    B_body = R_ItoB * B_I_T
    h_total = Matrix(rw.J_rw) * Vector(rw.h_wheels)
    m_cmd = -cfg.k_desaturation .* cross(B_body, h_total)
    m_norm = norm(m_cmd)
    m_norm > cfg.max_dipole_Am2 && (m_cmd .*= cfg.max_dipole_Am2 / m_norm)
    magnet.m .= m_cmd
    return m_cmd
end

# ── calcControlEffect! integration ───────────────────────────────────────────

function calcControlEffect!(ctrl::ADCSController, u, p, t::Float64, sat_idx::Int)
    sat_idx == ctrl.config.sat_idx || return nothing
    q_est    = collect(u.sc[sat_idx].q)
    omega_est = collect(u.sc[sat_idx].omega)
    ads_valid = true; ads_fine = true       # TODO: read from ADS output when wired
    link   = u.sc[sat_idx].links[1]
    rw     = link.rw_assembly
    magnet = length(link.magnets) > 0 ? link.magnets[1] : nothing
    sat_frac = rw_saturation_fraction(rw)
    update_adcs_mode!(ctrl, q_est, omega_est, ads_valid, ads_fine, sat_frac, t)
    tau_cmd     = quaternion_feedback_torque(ctrl, q_est, omega_est)
    tau_body_net = allocate_to_reaction_wheels!(rw, tau_cmd)
    if (ctrl.state.mode == ADCS_DESATURATING || sat_frac > ctrl.config.desaturation_threshold) && magnet !== nothing
        compute_desaturation_dipole!(ctrl, magnet, rw, q_est, collect(u.sc[sat_idx].pos))
    elseif magnet !== nothing
        magnet.m .= zeros(3)
    end
    ctrl.held = ADCSHeldActuation(tau_body_net)
    return nothing
end

function calcControlForceTorque(ctrl::ADCSController, u, p, i::Int, t::Float64)
    i == ctrl.config.sat_idx || return zeros(3), zeros(3)
    return zeros(3), ctrl.held.torque_body
end

calcControlMassFlowRate(ctrl::ADCSController, u, p, i::Int, t::Float64) :: Float64 = 0.0

# ── Registration ─────────────────────────────────────────────────────────────

function register_adcs_controller!(sim_config; config::ADCSConfig = ADCSConfig())
    ctrl = ADCSController(config=config)
    # push!(sim_config.control_effectors, ctrl)  # adapt to actual SimConfig API
    return ctrl
end
