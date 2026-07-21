function build_simulation_configuration(config, state, action::AerobrakingAction)
    return (
        backend_mode = config.backend_mode,
        phase = config.phase,
        pass_index = state.pass_index + 1,
        apoapsis_radius_m = state.apoapsis_radius_m,
        periapsis_altitude_m = state.periapsis_altitude_m,
        action_delta_v_mps = action.delta_v_mps,
        action_phi_deg = action.phi_deg,
        drag_coefficient_scale = state.drag_coefficient_scale,
        lift_coefficient_scale = state.lift_coefficient_scale,
        marsgram_seed = state.marsgram_seed,
        marsgram_prediction_seed = state.marsgram_prediction_seed,
    )
end
