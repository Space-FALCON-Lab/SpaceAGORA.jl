# Implement ADS mode manager (ads_mode_manager.jl)

Create src/gnc/navigation/ads_mode_manager.jl implementing the ADS operational state machine.

Requirements: ADS-SYS-004, ADS-INT-001.

@enum ADSMode INIT COARSE_ACQ FINE_TRACK ECLIPSE SAFE

mutable struct ADSModeManager
  current_mode::ADSMode
  eclipse_flag::Bool
  sensor_health::Dict{Symbol, Bool}   # :star_tracker, :gyro, :magnetometer, :sun_sensor
  nees_history::Vector{Float64}
  nees_threshold::Float64
  estimated_error_deg::Float64
end

ads_mode_init(; nees_threshold=9.49) -> ADSModeManager

ads_mode_update!(mgr::ADSModeManager, sensor_health::Dict{Symbol,Bool}, eclipse_flag::Bool,
                  nees::Float64, estimated_error_deg::Float64, reset_command::Bool=false) -> ADSMode
  Implement all transitions from the ADS Mode State Machine:
  - INIT -> COARSE_ACQ: all sensors healthy
  - COARSE_ACQ -> FINE_TRACK: star_tracker healthy AND estimated_error < 2 deg
  - FINE_TRACK -> ECLIPSE: eclipse_flag becomes true
  - ECLIPSE -> FINE_TRACK: eclipse_flag false AND star_tracker healthy
  - ECLIPSE -> COARSE_ACQ: eclipse_flag false AND star_tracker NOT healthy
  - Any -> SAFE: !gyro_health OR nees > threshold
  - SAFE -> INIT: reset_command=true

ads_mode_get_active_sensors(mgr::ADSModeManager) -> Vector{Symbol}
  Returns which sensors to include in MEKF updates given current mode.
  FINE_TRACK: [:star_tracker, :gyro, :magnetometer, :sun_sensor]
  ECLIPSE: [:gyro, :magnetometer]
  COARSE_ACQ: [:magnetometer, :sun_sensor]
  SAFE/INIT: []