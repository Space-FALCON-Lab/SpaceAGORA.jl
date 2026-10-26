# The two propulsive maneuver guidance models, for the burn scheduler
# (BaseThrusterModel). Either acts on every spacecraft at configuration level, its
# `spacecraft_index` 0, or on one: declared on a spacecraft, it is bound to that
# spacecraft's run position (`bind_spacecraft`, guidance/hooks.jl), so each
# spacecraft can carry its own maneuver table or targets.
@kwdef struct AerobrakingCampaignPropulsiveManeuverGuidanceModel <: AbstractGuidanceModel
    maneuver_orbit_number::Vector{Int64} # The orbit number at which to perform the propulsive maneuver (e.g., 1 for the first apoapsis, 2 for the second, etc.)
    maneuver_Δv::Vector{Float64} # The desired Δv for the propulsive maneuver
    # Diagnostic replay mode: when non-empty (one entry per maneuver), each
    # burn's Δv is scaled by flight_apoapsis_radius / current osculating
    # apoapsis radius, so the burn delivers the flight's periapsis change
    # (Δr_p ∝ r_a·Δv) rather than the flight's Δv. Injects flight truth;
    # diagnostics only — empty (the default) replays Δv verbatim.
    maneuver_flight_apoapsis_radius_m::Vector{Float64} = Float64[]
    spacecraft_index::Int = 0 # The run position it acts on; 0, every spacecraft
end
# The positional form without the spacecraft, as before SpaceAGORA 0.7.0.
AerobrakingCampaignPropulsiveManeuverGuidanceModel(orbits, Δv, flight_apoapsis_radius_m) =
    AerobrakingCampaignPropulsiveManeuverGuidanceModel(orbits, Δv, flight_apoapsis_radius_m, 0)

@kwdef mutable struct ApoapsisTargetPeriapsisRaiseGuidanceModel <: AbstractGuidanceModel
    target_apoapsis_radius_m::Float64
    target_periapsis_altitude_m::Float64 = 200.0e3
    apoapsis_tolerance_m::Float64 = 0.0
    apoapsis_window_rad::Float64 = deg2rad(30.0)
    command_state::Vector{Int64} = Int64[] # Per spacecraft, by run position
    spacecraft_index::Int = 0 # The run position it acts on; 0, every spacecraft
end
# The positional form without the spacecraft, as before SpaceAGORA 0.7.0.
ApoapsisTargetPeriapsisRaiseGuidanceModel(r_a, h_p, tolerance, window, command_state) =
    ApoapsisTargetPeriapsisRaiseGuidanceModel(r_a, h_p, tolerance, window, command_state, 0)
