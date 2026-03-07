using ComponentArrays
function calcGuidanceEffect!(guidanceAlg::AerobrakingCampaignPropulsiveManeuverGuidanceModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    # Placeholder function for calculating the guidance effect based on the guidance model and current state, and store it in the shared buffers for use in the dynamics calculations
    # This function can be implemented later to include specific guidance logic for aerobraking campaigns, such as adjusting the target periapsis altitude for each pass based on the achieved Δv from the previous pass, etc.
    if i < 1 || i > length(p.shared_buffers.maneuver_commands)
        return nothing
    end

    orbit_counter = p.orbit_counter[i]
    maneuver_idx = findfirst(==(orbit_counter), guidanceAlg.maneuver_orbit_number)

    command = if maneuver_idx === nothing
        PropulsiveManeuverCommand(valid=true, source_orbit=orbit_counter)
    else
        delta_v_cmd = Float64(guidanceAlg.maneuver_Δv[maneuver_idx])
        PropulsiveManeuverCommand(
            valid=true,
            delta_v_mps=abs(delta_v_cmd),
            direction_rad=(delta_v_cmd > 0.0 ? π : 0.0),
            source_orbit=orbit_counter
        )
    end

    p.shared_buffers.maneuver_commands[i] = command
    return nothing
end
