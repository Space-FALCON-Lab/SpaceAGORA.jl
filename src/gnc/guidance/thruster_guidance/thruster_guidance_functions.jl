using ComponentArrays
function calcGuidanceEffect!(guidanceAlg::AerobrakingCampaignPropulsiveManeuverGuidanceModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    if i < 1 || i > length(p.shared_buffers.maneuver_commands)
        return nothing
    end

    orbit_counter = p.orbit_counter[i]
    maneuver_idx = findfirst(==(orbit_counter), guidanceAlg.maneuver_orbit_number)

    # Campaign guidance currently emits at most one burn command per orbit by
    # looking up the precomputed passage-to-delta-v table.
    command = if maneuver_idx === nothing
        PropulsiveManeuverCommand(valid=true, source_orbit=orbit_counter)
    else
        delta_v_cmd = Float64(guidanceAlg.maneuver_Δv[maneuver_idx])
        PropulsiveManeuverCommand(
            valid=true,
            delta_v_mps=abs(delta_v_cmd),
            direction_rad=(delta_v_cmd > 0.0 ? 0.0 : π),
            source_orbit=orbit_counter
        )
    end

    p.shared_buffers.maneuver_commands[i] = command
    return nothing
end
