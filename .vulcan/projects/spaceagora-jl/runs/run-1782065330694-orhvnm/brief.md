# Implement ADS navigation effector (ads_navigation_effector.jl)

Create src/gnc/navigation/ads_navigation_effector.jl — the SpaceAGORA integration entry point.

Requirements: ADS-SYS-001, ADS-INT-001, ADS-INT-002, ADS-INT-003.

This file integrates sensor models + MEKF + mode manager and registers with SpaceAGORA.

struct ADSConfig
  update_rate_hz::Float64       # default 1.0
  star_tracker_noise_arcsec::Float64  # default 5.0
  gyro_arw_deg_sqrthr::Float64       # default 0.1
  gyro_bias_instability_deg_hr::Float64  # default 1.0
  mag_noise_nT::Float64              # default 50.0
  sun_sensor_noise_deg::Float64      # default 1.0
  nees_threshold::Float64            # default 9.49
  star_tracker_fov_deg::Float64      # default 15.0
  sun_sensor_fov_deg::Float64        # default 70.0
end
ADSConfig() = ADSConfig(1.0, 5.0, 0.1, 1.0, 50.0, 1.0, 9.49, 15.0, 70.0)

mutable struct ADSNavigationEffector
  config::ADSConfig
  sensor_suite   # tuple or struct of all four sensor models
  mekf_state::ADSMEKFState
  mode_manager::ADSModeManager
  last_update_time::Float64
  initialized::Bool
end

function calcNavigationEffect!(effector::ADSNavigationEffector, model, u, p, t, sat_idx)
  dt = t - effector.last_update_time
  q_true = model.body.roots[sat_idx].q
  omega_true = model.body.roots[sat_idx].omega  # or equivalent field name

  # Get reference vectors
  sun_dir_I = get_sun_direction(model, t)       # from existing ephemerides
  B_field_I = get_magnetic_field_dipole(model, t, sat_idx)

  # Measure all sensors
  (q_st, st_valid) = measure(effector.sensor_suite.star_tracker, q_true, omega_true, t)
  (omega_meas, gyro_valid) = measure(effector.sensor_suite.gyro, omega_true, t)
  (b_meas, mag_valid) = measure(effector.sensor_suite.magnetometer, q_true, B_field_I, t)
  (s_meas, sun_valid) = measure(effector.sensor_suite.sun_sensor, q_true, sun_dir_I, t)

  # Update mode
  eclipse_flag = is_eclipse(model, t, sat_idx)  # use existing SpaceAGORA eclipse check
  sensor_health = Dict(:star_tracker=>st_valid, :gyro=>gyro_valid, :magnetometer=>mag_valid, :sun_sensor=>sun_valid)
  nees = compute_nees(effector.mekf_state)
  mode = ads_mode_update!(effector.mode_manager, sensor_health, eclipse_flag, nees, estimate_error_deg(effector.mekf_state))

  # MEKF predict
  Q = build_process_noise(effector.config, dt)
  ads_mekf_predict!(effector.mekf_state, omega_meas, dt, Q)

  # MEKF update (only for active sensors in current mode)
  active = ads_mode_get_active_sensors(effector.mode_manager)
  if :star_tracker in active && st_valid
    R_st = build_st_noise(effector.config)
    ads_mekf_update_star_tracker!(effector.mekf_state, q_st, R_st)
  end
  if :magnetometer in active && mag_valid
    R_mag = build_mag_noise(effector.config)
    ads_mekf_update_vector!(effector.mekf_state, b_meas, B_field_I, R_mag)
  end
  if :sun_sensor in active && sun_valid
    R_sun = build_sun_noise(effector.config)
    ads_mekf_update_vector!(effector.mekf_state, s_meas, sun_dir_I, R_sun)
  end

  # Write output — adapt to actual NavigationEffect output struct fields
  valid = (mode != SAFE) && gyro_valid
  set_navigation_output!(model, sat_idx, effector.mekf_state.q_est, omega_meas - effector.mekf_state.gyro_bias, valid)
  effector.last_update_time = t
end

function register_ads_navigation!(sim_config, ads_config=ADSConfig())
  # Create effector, register with sim_config at ads_config.update_rate_hz
  # Follow pattern in navigation_hooks.jl and navigation_guidance_callbacks.jl
end

Look at src/gnc/navigation/navigation_hooks.jl to understand the exact calcNavigationEffect! signature and how to register. Adapt field names to match the actual SpaceAGORA model struct.