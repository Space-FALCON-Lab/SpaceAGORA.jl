function compute_t_edg_guidance_window!(input::AerobrakingGuidanceInput)
    args = input.args
    cnf_state = _bridge_get_cnf(args; cnf=input.cnf)

    # T-EDG targeting is solved in targeting_solver/eom_predictor paths.
    # For runtime control-hook integration, we retain current switch state unless a
    # dedicated targeting schedule is provided by mission policy.
    return AerobrakingGuidanceOutput(
        time_switch_1=cnf_state.time_switch_1,
        time_switch_2=cnf_state.time_switch_2,
        security_mode=cnf_state.security_mode,
    )
end

function compute_aerobraking_guidance(::TEdgStrategy, input::AerobrakingGuidanceInput)
    return compute_t_edg_guidance_window!(input)
end
