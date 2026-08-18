@kwdef struct AerobrakingCampaignPropulsiveManeuverGuidanceModel <: AbstractGuidanceModel
    maneuver_orbit_number::Vector{Int64} # The orbit number at which to perform the propulsive maneuver (e.g., 1 for the first apoapsis, 2 for the second, etc.)
    maneuver_Δv::Vector{Float64} # The desired Δv for the propulsive maneuver
    # Diagnostic replay mode: when non-empty (one entry per maneuver), each
    # burn's Δv is scaled by flight_apoapsis_radius / current osculating
    # apoapsis radius, so the burn delivers the flight's periapsis change
    # (Δr_p ∝ r_a·Δv) rather than the flight's Δv. Injects flight truth;
    # diagnostics only — empty (the default) replays Δv verbatim.
    maneuver_flight_apoapsis_radius_m::Vector{Float64} = Float64[]
end

@kwdef mutable struct ApoapsisTargetPeriapsisRaiseGuidanceModel <: AbstractGuidanceModel
    target_apoapsis_radius_m::Float64
    target_periapsis_altitude_m::Float64 = 200.0e3
    apoapsis_tolerance_m::Float64 = 0.0
    apoapsis_window_rad::Float64 = deg2rad(30.0)
    command_state::Vector{Int64} = Int64[]
end
