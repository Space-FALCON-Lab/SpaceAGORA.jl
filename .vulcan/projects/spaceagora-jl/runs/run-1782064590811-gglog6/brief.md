# Implement ADS sensor models (ads_sensor_models.jl)

Create src/gnc/navigation/ads_sensor_models.jl in SpaceAGORA.jl with four sensor measurement model structs.

Requirements: ADS-SNS-001 through ADS-SNS-004.

1. StarTrackerModel: quaternion measurement with 5 arcsec/axis Gaussian noise (std), 15 deg half-angle FOV, eclipse and sun-exclusion masks. measure() returns (q_meas::Quaternion, is_valid::Bool).

2. GyroModel: angular rate measurement with angle random walk (ARW=0.1 deg/sqrt(hr)) and bias instability (1 deg/hr). Maintains a bias state that random-walks. measure() returns (omega_meas::SVector{3}, is_valid::Bool). Bias is provided as internal state.

3. MagnetometerModel: calls get_magnetic_field_dipole(model, t, sat_idx) from perturbations.jl, rotates to body frame using q_true, adds 50 nT/axis Gaussian noise, normalizes. measure() returns (b_meas_normalized::SVector{3}, is_valid::Bool).

4. SunSensorModel: uses existing SpaceAGORA sun ephemeris to get inertial sun direction, rotates to body frame, adds 1 deg/axis noise, checks 70 deg half-angle FOV and eclipse occultation. measure() returns (s_meas::SVector{3}, is_valid::Bool).

All structs: use StaticArrays (SVector, etc.) if available; else plain Vector. Follow existing SpaceAGORA struct conventions (struct, keyword constructor). All measure() functions take (model, truth_state, env_state, t, sat_idx).

Reference files:
- src/dynamics/coupled/perturbations.jl: get_magnetic_field_dipole
- src/environment/ephemerides/: sun direction
- src/dynamics/rotational/rigid_body_dynamics.jl: truth quaternion