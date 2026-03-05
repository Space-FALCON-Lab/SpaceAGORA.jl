using ComponentArrays
function calcGuidanceEffect!(guidanceAlg::AerobrakingCampaignPropulsiveManeuverGuidanceModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    # Placeholder function for calculating the guidance effect based on the guidance model and current state, and store it in the shared buffers for use in the dynamics calculations
    # This function can be implemented later to include specific guidance logic for aerobraking campaigns, such as adjusting the target periapsis altitude for each pass based on the achieved Δv from the previous pass, etc.
    thruster = p.args.control_model.control_effectors[1] # Assuming the first control effector is the thruster model for the propulsive maneuver
    orbit_counter = p.orbit_counter[i]
    if orbit_counter in guidanceAlg.maneuver_orbit_number
        thruster.Δv[i] = guidanceAlg.maneuver_Δv[findfirst(==(orbit_counter), guidanceAlg.maneuver_orbit_number)]
        if thruster.Δv[i] < 0.0
            thruster.direction[i] = π # Retrograde burn
            # println("Lower maneuver: Orbit $orbit_counter, Δv = $(thruster.Δv[i]) m/s (retrograde)")
        else            
            thruster.direction[i] = 0.0 # Prograde burn
            # println("Raise maneuver: Orbit $orbit_counter, Δv = $(thruster.Δv[i]) m/s (prograde)")
        end
        thruster.Δv[i] = abs(thruster.Δv[i])
    else
        thruster.Δv[i] = 0.0
    end
end