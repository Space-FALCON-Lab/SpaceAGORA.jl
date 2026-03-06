@kwdef struct AerobrakingCampaignPropulsiveManeuverGuidanceModel <: AbstractGuidanceModel
    maneuver_orbit_number::Vector{Int64} # The orbit number at which to perform the propulsive maneuver (e.g., 1 for the first apoapsis, 2 for the second, etc.)
    maneuver_Δv::Vector{Float64} # The desired Δv for the propulsive maneuver
end