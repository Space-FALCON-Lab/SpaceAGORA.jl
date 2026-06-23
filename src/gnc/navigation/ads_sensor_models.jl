# Attitude Determination Subsystem — Sensor Measurement Models
# DRM 1 | Requirements: ADS-SNS-001 through ADS-SNS-004
# Implements measurement models for: StarTracker, Gyroscope, Magnetometer, SunSensor

# Quaternion convention: [w, x, y, z] (scalar first)

# ── Star Tracker Model ────────────────────────────────────────────────────────

struct StarTrackerModel
    noise_sigma_rad::Float64      # per-axis 1σ noise (rad)
    fov_half_angle_rad::Float64   # field-of-view half-angle (rad)
    sun_exclusion_rad::Float64    # minimum sun angle (rad)
end

function StarTrackerModel(;
    noise_arcsec::Float64 = 5.0,
    fov_deg::Float64 = 15.0,
    sun_exclusion_deg::Float64 = 30.0)
    StarTrackerModel(
        noise_arcsec * π / (180.0 * 3600.0),
        fov_deg * π / 180.0,
        sun_exclusion_deg * π / 180.0)
end

"""
    measure_star_tracker(model, q_true, sun_dir_body, eclipse_flag)

Return (q_meas, is_valid). ADS-SNS-001.
"""
function measure_star_tracker(model::StarTrackerModel,
                               q_true::AbstractVector,
                               sun_dir_body::AbstractVector,
                               eclipse_flag::Bool)
    eclipse_flag && return (copy(q_true), false)

    # Sun exclusion check
    s = normalize(sun_dir_body)
    boresight = [0.0, 0.0, 1.0]
    cos_angle = clamp(dot(s, boresight), -1.0, 1.0)
    abs(acos(cos_angle)) < model.sun_exclusion_rad && return (copy(q_true), false)

    # Add small rotation noise
    noise_axis = normalize(randn(3))
    noise_angle = model.noise_sigma_rad * randn()
    half = noise_angle / 2
    δq = [cos(half); noise_axis .* sin(half)]
    q_meas = quat_multiply(q_true, δq)
    return (normalize(q_meas), true)
end

# ── Gyroscope Model ───────────────────────────────────────────────────────────

mutable struct GyroModel
    arw_sigma::Float64        # angle random walk σ per root-Hz (rad/s/√Hz)
    bias_sigma::Float64       # bias random walk σ per step (rad/s)
    bias::Vector{Float64}     # current bias (truth, rad/s)
end

function GyroModel(;
    arw_deg_sqrthr::Float64 = 0.1,
    bias_deg_hr::Float64 = 1.0)
    arw = arw_deg_sqrthr * π / 180.0 / sqrt(3600.0)
    bias_rw = bias_deg_hr * π / 180.0 / 3600.0
    GyroModel(arw, bias_rw, zeros(3))
end

"""
    measure_gyro!(model, omega_true, dt) -> (omega_meas, is_valid)

ADS-SNS-002.
"""
function measure_gyro!(model::GyroModel, omega_true::AbstractVector, dt::Float64)
    arw_noise = (model.arw_sigma / sqrt(dt)) .* randn(3)
    model.bias .+= model.bias_sigma .* sqrt(dt) .* randn(3)
    return (omega_true .+ model.bias .+ arw_noise, true)
end

# ── Magnetometer Model ────────────────────────────────────────────────────────

struct MagnetometerModel
    noise_nT::Float64    # per-axis 1σ noise (nT)
end
MagnetometerModel(; noise_nT::Float64 = 50.0) = MagnetometerModel(noise_nT)

"""
    measure_magnetometer(model, q_true, B_I) -> (b_meas, is_valid)

B_I in inertial frame (nT). Returns normalized body-frame vector. ADS-SNS-003.
"""
function measure_magnetometer(model::MagnetometerModel,
                               q_true::AbstractVector,
                               B_I::AbstractVector)
    B_body = quat_rotate(q_true, B_I)
    B_noisy = B_body .+ model.noise_nT .* randn(3)
    n = norm(B_noisy)
    n < 1.0 && return (B_noisy, false)
    return (B_noisy ./ n, true)
end

# ── Sun Sensor Model ──────────────────────────────────────────────────────────

struct SunSensorModel
    noise_sigma_rad::Float64
    fov_half_angle_rad::Float64
end
function SunSensorModel(; noise_deg::Float64 = 1.0, fov_deg::Float64 = 70.0)
    SunSensorModel(noise_deg * π / 180.0, fov_deg * π / 180.0)
end

"""
    measure_sun_sensor(model, q_true, sun_I, eclipse_flag) -> (s_meas, is_valid)

ADS-SNS-004.
"""
function measure_sun_sensor(model::SunSensorModel,
                             q_true::AbstractVector,
                             sun_I::AbstractVector,
                             eclipse_flag::Bool)
    eclipse_flag && return ([0.0, 0.0, 1.0], false)

    sun_body = quat_rotate(q_true, normalize(sun_I))
    # FOV check against boresight [0, 0, 1]
    cos_angle = clamp(sun_body[3], -1.0, 1.0)
    acos(cos_angle) > model.fov_half_angle_rad && return (sun_body, false)

    noise = model.noise_sigma_rad .* randn(3)
    s_meas = normalize(sun_body .+ noise)
    return (s_meas, true)
end

# ── Quaternion utilities ──────────────────────────────────────────────────────

"""Hamilton product of two quaternions [w,x,y,z]."""
function quat_multiply(p::AbstractVector, q::AbstractVector)
    pw, px, py, pz = p[1], p[2], p[3], p[4]
    qw, qx, qy, qz = q[1], q[2], q[3], q[4]
    [pw*qw - px*qx - py*qy - pz*qz,
     pw*qx + px*qw + py*qz - pz*qy,
     pw*qy - px*qz + py*qw + pz*qx,
     pw*qz + px*qy - py*qx + pz*qw]
end

"""Rotate vector v by quaternion q (active rotation)."""
function quat_rotate(q::AbstractVector, v::AbstractVector)
    qv = [0.0; v]
    qc = [q[1]; -q[2]; -q[3]; -q[4]]
    r = quat_multiply(quat_multiply(q, qv), qc)
    r[2:4]
end

"""Angle in degrees between two quaternions."""
function quat_angle_deg(q1::AbstractVector, q2::AbstractVector)
    q1c = [q1[1]; -q1[2]; -q1[3]; -q1[4]]
    dq = quat_multiply(q1c, q2)
    2.0 * acos(clamp(abs(dq[1]), 0.0, 1.0)) * 180.0 / π
end

"""Quaternion from axis-angle."""
function quat_from_axis_angle(axis::AbstractVector, angle::Float64)
    n = norm(axis)
    n < 1e-12 && return [1.0, 0.0, 0.0, 0.0]
    a = axis ./ n
    [cos(angle/2); a .* sin(angle/2)]
end

