# Attitude Determination Subsystem — Operational Mode Manager
# DRM 1 | Requirements: ADS-SYS-004
# Implements mode state machine from ADS — Mode State Machine sheet

@enum ADSMode::Int8 INIT=1 COARSE_ACQ=2 FINE_TRACK=3 ECLIPSE=4 SAFE=5

mutable struct ADSModeManager
    current_mode::ADSMode
    nees_threshold::Float64
    fine_track_error_threshold_deg::Float64
end

function ads_mode_init(;
    nees_threshold::Float64 = 9.49,
    fine_track_threshold_deg::Float64 = 2.0)
    ADSModeManager(INIT, nees_threshold, fine_track_threshold_deg)
end

"""
    ads_mode_update!(mgr, sensor_health, eclipse_flag, nees, error_deg, reset=false)

Advance the ADS mode state machine. Returns updated mode.
"""
function ads_mode_update!(mgr::ADSModeManager,
                           sensor_health::Dict{Symbol,Bool},
                           eclipse_flag::Bool,
                           nees::Float64,
                           error_deg::Float64,
                           reset_command::Bool = false)
    gyro_ok = get(sensor_health, :gyro, false)
    st_ok   = get(sensor_health, :star_tracker, false)
    mag_ok  = get(sensor_health, :magnetometer, false)
    sun_ok  = get(sensor_health, :sun_sensor, false)

    # Priority: fault → SAFE (from any state except SAFE)
    if !gyro_ok && mgr.current_mode != SAFE
        mgr.current_mode = SAFE
        return mgr.current_mode
    end
    if nees > mgr.nees_threshold && mgr.current_mode ∈ (FINE_TRACK, ECLIPSE)
        mgr.current_mode = SAFE
        return mgr.current_mode
    end

    if mgr.current_mode == INIT
        (gyro_ok && (mag_ok || sun_ok)) && (mgr.current_mode = COARSE_ACQ)

    elseif mgr.current_mode == COARSE_ACQ
        if eclipse_flag
            mgr.current_mode = ECLIPSE
        elseif st_ok && error_deg < mgr.fine_track_error_threshold_deg
            mgr.current_mode = FINE_TRACK
        end

    elseif mgr.current_mode == FINE_TRACK
        eclipse_flag && (mgr.current_mode = ECLIPSE)

    elseif mgr.current_mode == ECLIPSE
        if !eclipse_flag
            mgr.current_mode = st_ok ? FINE_TRACK : COARSE_ACQ
        end

    elseif mgr.current_mode == SAFE
        reset_command && (mgr.current_mode = INIT)
    end

    return mgr.current_mode
end

"""Returns sensor symbols to use in MEKF update for the current mode."""
function ads_mode_get_active_sensors(mgr::ADSModeManager)
    mgr.current_mode == FINE_TRACK  && return [:star_tracker, :gyro]
    mgr.current_mode == ECLIPSE     && return [:gyro, :magnetometer]
    mgr.current_mode == COARSE_ACQ  && return [:magnetometer, :sun_sensor]
    return Symbol[]
end

