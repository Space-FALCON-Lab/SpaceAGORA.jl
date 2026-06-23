# Attitude Determination Subsystem — SpaceAGORA Navigation Effector
# DRM 1 | Requirements: ADS-SYS-001, ADS-INT-001, ADS-INT-002, ADS-INT-003
#
# Integration point: register as a calcNavigationEffect! effector.
# All sensor physics reuse existing SpaceAGORA models (no new environment code).

include("ads_sensor_models.jl")
include("ads_mekf.jl")
include("ads_mode_manager.jl")

# ── Configuration ─────────────────────────────────────────────────────────────

struct ADSConfig
    update_rate_hz::Float64
    star_tracker_noise_arcsec::Float64
    gyro_arw_deg_sqrthr::Float64
    gyro_bias_instability_deg_hr::Float64
    mag_noise_nT::Float64
    sun_sensor_noise_deg::Float64
    nees_threshold::Float64
    star_tracker_fov_deg::Float64
    sun_sensor_fov_deg::Float64
end

ADSConfig() = ADSConfig(1.0, 5.0, 0.1, 1.0, 50.0, 1.0, 9.49, 15.0, 70.0)

# ── Effector ──────────────────────────────────────────────────────────────────

mutable struct ADSNavigationEffector
    config::ADSConfig
    star_tracker::StarTrackerModel
    gyro::GyroModel
    magnetometer::MagnetometerModel
    sun_sensor::SunSensorModel
    mekf::ADSMEKFState
    mode_mgr::ADSModeManager
    last_t::Float64
    initialized::Bool
end

function ADSNavigationEffector(config::ADSConfig = ADSConfig())
    P0 = Diagonal([1e-3, 1e-3, 1e-3, 1e-5, 1e-5, 1e-5]) * 1.0
    ADSNavigationEffector(
        config,
        StarTrackerModel(noise_arcsec = config.star_tracker_noise_arcsec,
                         fov_deg      = config.star_tracker_fov_deg),
        GyroModel(arw_deg_sqrthr = config.gyro_arw_deg_sqrthr,
                  bias_deg_hr   = config.gyro_bias_instability_deg_hr),
        MagnetometerModel(noise_nT = config.mag_noise_nT),
        SunSensorModel(noise_deg = config.sun_sensor_noise_deg),
        ads_mekf_init(P0 = Matrix(P0)),
        ads_mode_init(nees_threshold = config.nees_threshold),
        0.0,
        false)
end

# Noise matrices (built once per config)
function _Q(config::ADSConfig, dt::Float64)
    arw_var  = (config.gyro_arw_deg_sqrthr * π/180 / sqrt(3600.0))^2 / dt
    bias_var = (config.gyro_bias_instability_deg_hr * π/180 / 3600.0)^2 * dt
    Diagonal([arw_var, arw_var, arw_var, bias_var, bias_var, bias_var]) * 1.0
end

_R_st(config::ADSConfig) =
    (config.star_tracker_noise_arcsec * π/180/3600)^2 * Matrix{Float64}(I, 3, 3)

_R_mag(config::ADSConfig) =
    (config.mag_noise_nT / 50000e-9 * 1e-9 / 50000e-9)^2 * Matrix{Float64}(I, 3, 3)

_R_sun(config::ADSConfig) =
    (config.sun_sensor_noise_deg * π/180)^2 * Matrix{Float64}(I, 3, 3)

# ── calcNavigationEffect! hook ────────────────────────────────────────────────

"""
    calcNavigationEffect!(effector::ADSNavigationEffector, model, u, p, t, sat_idx)

SpaceAGORA navigation effector hook. ADS-INT-002.

NOTE: Field names (model.body.roots[sat_idx].q, etc.) must match actual SpaceAGORA
model struct. Adapt if SpaceAGORA changes field naming conventions.
"""
function calcNavigationEffect!(effector::ADSNavigationEffector,
                                model, u, p, t::Float64, sat_idx::Int)
    dt = effector.initialized ?
        max(t - effector.last_t, 1e-6) :
        1.0 / effector.config.update_rate_hz

    # ── Truth state (ADS-INT-003) ─────────────────────────────────────────────
    # Adapt these fields to the actual SpaceAGORA model struct:
    q_true    = collect(Float64, model.body.roots[sat_idx].q)
    ω_true    = collect(Float64, model.body.roots[sat_idx].omega)

    # ── Reference vectors (reuse existing SpaceAGORA functions) ─────────────
    # These calls must match the actual SpaceAGORA API:
    sun_I_raw = get_sun_direction_inertial(model, t)        # from src/environment/ephemerides/
    B_I_raw   = get_magnetic_field_dipole(model, t, sat_idx) # from src/dynamics/coupled/perturbations.jl
    sun_I     = normalize(collect(Float64, sun_I_raw))
    B_I       = collect(Float64, B_I_raw)  # keep in nT for magnetometer noise model

    # Eclipse flag — use existing SpaceAGORA eclipse detection
    eclipse = is_in_eclipse(model, t, sat_idx)

    # ── Sensor measurements ───────────────────────────────────────────────────
    sun_body = quat_rotate(q_true, sun_I)
    (q_st,  st_valid)  = measure_star_tracker(effector.star_tracker,  q_true, sun_body, eclipse)
    (ω_meas, gyro_valid) = measure_gyro!(effector.gyro, ω_true, dt)
    (b_meas, mag_valid)  = measure_magnetometer(effector.magnetometer, q_true, B_I)
    (s_meas, sun_valid)  = measure_sun_sensor(effector.sun_sensor, q_true, sun_I, eclipse)

    # ── Mode update ───────────────────────────────────────────────────────────
    health = Dict(:star_tracker => st_valid, :gyro => gyro_valid,
                  :magnetometer => mag_valid, :sun_sensor => sun_valid)
    nees = effector.initialized ? compute_nees(effector.mekf, q_true) : 0.0
    err  = effector.initialized ? attitude_error_deg(effector.mekf, q_true) : 180.0
    mode = ads_mode_update!(effector.mode_mgr, health, eclipse, nees, err)
    active = ads_mode_get_active_sensors(effector.mode_mgr)

    # ── MEKF predict ──────────────────────────────────────────────────────────
    ads_mekf_predict!(effector.mekf, ω_meas, dt, Matrix(_Q(effector.config, dt)))

    # ── MEKF update (active sensors only) ────────────────────────────────────
    :star_tracker ∈ active && st_valid &&
        ads_mekf_update_star_tracker!(effector.mekf, q_st, _R_st(effector.config))
    :magnetometer ∈ active && mag_valid &&
        ads_mekf_update_vector!(effector.mekf, b_meas, normalize(B_I), _R_mag(effector.config))
    :sun_sensor ∈ active && sun_valid &&
        ads_mekf_update_vector!(effector.mekf, s_meas, sun_I, _R_sun(effector.config))

    # ── Output ────────────────────────────────────────────────────────────────
    valid = (mode != SAFE) && gyro_valid
    ω_est = ω_meas .- effector.mekf.gyro_bias

    effector.last_t = t
    effector.initialized = true

    # Return tuple; caller writes to NavigationEffect output struct
    return (effector.mekf.q_est, ω_est, valid, Int(mode))
end

"""
    register_ads_navigation!(sim_config; config=ADSConfig()) -> ADSNavigationEffector

Register ADS effector with SpaceAGORA SimulationConfiguration. ADS-INT-002.
"""
function register_ads_navigation!(sim_config; config::ADSConfig = ADSConfig())
    effector = ADSNavigationEffector(config)
    # Adapt to actual SpaceAGORA registration API:
    # push!(sim_config.navigation_effectors, (effector, 1.0 / config.update_rate_hz))
    return effector
end

