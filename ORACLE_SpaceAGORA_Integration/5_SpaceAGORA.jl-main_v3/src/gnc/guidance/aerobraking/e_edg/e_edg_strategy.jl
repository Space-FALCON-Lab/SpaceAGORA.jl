function compute_e_edg_guidance_window!(input::AerobrakingGuidanceInput)
    ip = input.ip
    m = input.mission
    args = input.args
    index_ratio = input.index_ratio
    t = input.t
    position = input.position
    current_position = input.current_position
    gram_atmosphere = input.gram_atmosphere
    heat_rate_control = input.heat_rate_control

    cnf_state = _bridge_get_cnf(args; cnf=input.cnf)
    if args[:flash2_through_integration] == 1
        args[:security_mode] = false
    end

    start_reevaluation = cnf_state.ascending_phase ||
        (!cnf_state.ascending_phase && abs(cnf_state.time_switch_2 - t) < 0.2 * cnf_state.time_switch_2)

    if !cnf_state.evaluate_switch_heat_load
        if args[:flash2_through_integration] == 1
            if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 1
                cnf_state.time_switch_1, cnf_state.time_switch_2 =
                    switch_calculation_with_integration(ip, m, position, args, t, heat_rate_control, 1, gram_atmosphere, position; cnf=cnf_state)
            end
            if args[:heat_load_sol] == 2 || args[:heat_load_sol] == 3
                cnf_state.time_switch_1, cnf_state.time_switch_2 =
                    second_time_switch_recalc_with_integration(ip, m, position, args, t, heat_rate_control, 1, gram_atmosphere, position; cnf=cnf_state)
            end
        else
            cnf_state.time_switch_1, cnf_state.time_switch_2 =
                switch_calculation(ip, m, position, args, t, heat_rate_control, 1, position; cnf=cnf_state)
        end
        cnf_state.evaluate_switch_heat_load = true
    elseif cnf_state.evaluate_switch_heat_load && start_reevaluation &&
           ((t - cnf_state.time_switch_1 > 1 && t - cnf_state.timer_revaluation > 10 && t - cnf_state.time_switch_2 < 0) ||
            (3 < cnf_state.time_switch_2 - t < 50 && t - cnf_state.timer_revaluation > 3) ||
            (0 < cnf_state.time_switch_2 - t < 3 && t - cnf_state.timer_revaluation > 0.8)) &&
           !cnf_state.security_mode
        reevaluation_mode = (t - cnf_state.timer_revaluation) > 3 ? 1 : 2
        cnf_state.timer_revaluation = t

        if args[:second_switch_reevaluation] == 1
            if args[:flash2_through_integration] == 1
                cnf_state.time_switch_1, cnf_state.time_switch_2 =
                    second_time_switch_recalc_with_integration(ip, m, position, args, t, heat_rate_control, reevaluation_mode, gram_atmosphere, current_position; cnf=cnf_state)
            else
                cnf_state.time_switch_1, cnf_state.time_switch_2 =
                    second_time_switch_recalc(ip, m, position, args, t, heat_rate_control, current_position, reevaluation_mode; cnf=cnf_state)
            end
        end
    elseif !isempty(index_ratio) && length(index_ratio) >= 2 &&
           maximum(cnf_state.heat_load_past) > 0.98 * m.aerodynamics.heat_load_limit &&
           any(i -> i > 2, cnf_state.heat_load_past - cnf_state.heat_load_ppast) &&
           args[:security_mode] == 1 && !cnf_state.security_mode && index_ratio[2] == 1
        cnf_state.time_switch_1, cnf_state.time_switch_2 = security_mode(ip, m, position, args, t, false; cnf=cnf_state)
    end

    return AerobrakingGuidanceOutput(
        time_switch_1=cnf_state.time_switch_1,
        time_switch_2=cnf_state.time_switch_2,
        security_mode=cnf_state.security_mode,
    )
end

function compute_aerobraking_guidance(::EEdgStrategy, input::AerobrakingGuidanceInput)
    return compute_e_edg_guidance_window!(input)
end
