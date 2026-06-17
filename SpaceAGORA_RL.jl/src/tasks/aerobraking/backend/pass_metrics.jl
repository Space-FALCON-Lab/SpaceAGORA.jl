Base.@kwdef struct AerobrakingPassMetrics
    pass_index::Int = 0
    apoapsis_radius_m::Float64 = NaN
    periapsis_altitude_m::Float64 = NaN
    argument_of_periapsis_rad::Float64 = NaN
    raan_rad::Float64 = NaN
    inclination_rad::Float64 = NaN
    drag_passage_time_s::Float64 = NaN
    date_ordinal::Float64 = NaN
    max_density_kg_m3::Float64 = NaN
    max_heat_rate_w_cm2::Float64 = NaN
    total_delta_v_mps::Float64 = 0.0
    maneuver_count::Int = 0
    mission_elapsed_s::Float64 = 0.0
    solver_converged::Bool = true
end

function pass_metrics_from_state(state; pass_index=state.pass_index,
                                 solver_converged::Bool=true)
    return AerobrakingPassMetrics(
        pass_index = pass_index,
        apoapsis_radius_m = state.apoapsis_radius_m,
        periapsis_altitude_m = state.periapsis_altitude_m,
        argument_of_periapsis_rad = state.argument_of_periapsis_rad,
        raan_rad = state.raan_rad,
        inclination_rad = state.inclination_rad,
        drag_passage_time_s = state.last_drag_passage_time_s,
        date_ordinal = python_ordinal(state.epoch),
        max_density_kg_m3 = state.previous_density_kg_m3,
        max_heat_rate_w_cm2 = state.previous_heat_rate_w_cm2,
        total_delta_v_mps = state.total_delta_v_mps,
        maneuver_count = state.maneuver_count,
        mission_elapsed_s = state.mission_elapsed_s,
        solver_converged = solver_converged,
    )
end
