include("../../utils/Closed_form_solution.jl")
include("../../physical_models/Density_models.jl")
include("../Eoms.jl")

using Roots

 # import .config

function second_time_switch_recalc_with_integration(ip, m, position, args, t, heat_rate_control, reevaluation_mode, gram_atmosphere=nothing, current_position=0)
    time_switch = config.cnf.time_switch_2

    function func(t_s)
        # y = asim(ip, m, t, current_position, args, 0, heat_rate_control, false, t_s, reevaluation_mode)
        y, time_switch = asim_ctrl(ip, m, t, current_position, args, 0, heat_rate_control, false, gram_atmosphere, t_s, reevaluation_mode)

        Q = y[end,end]

        return Q - m.aerodynamics.heat_load_limit
    end

    # Evaluates max Q
    delta_Q_switch = func(time_switch)
    delta_Q_current_time = func(t)

    delta_Q_with_current_time = abs(delta_Q_current_time - delta_Q_switch) # if the difference between the heat load at the current time and the switch is small but the time is too large, recheck

    if abs(delta_Q_switch) <= 0.01
        if (args[:heat_load_sol] == 0 || args[:heat_load_sol] == 3) && !(delta_Q_with_current_time < 0.5 && abs(t - time_switch) > 20) && (delta_Q_switch < delta_Q_current_time)
            return config.cnf.time_switch_1, config.cnf.time_switch_2
        elseif (args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2)
            return config.cnf.time_switch_1, config.cnf.time_switch_2
        end
    elseif delta_Q_current_time < 0
        if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 3
            return config.cnf.time_switch_1, t
        end
    end

    if reevaluation_mode == 1
        x_tol = 0.01
    else
        x_tol = 0.01
    end

    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 1
        b = t + 200
    else
        b = t + 1500
    end

    # try
        time_switch = fzero(k -> func(k), [t, b], Roots.Brent())
    # catch
    #     nothing
    # end

    return config.cnf.time_switch_1, time_switch
end

function second_time_switch_recalc(ip, m, position, args, t, heat_rate_control, current_position=0, reevaluation_mode=0)
    # Evaluates past heat load
    aoa_past = config.cnf.α_list
    time_switch_1 = config.cnf.time_switch_1
    time_switch_2 = config.cnf.time_switch_2

    Q_past = config.cnf.heat_load_past

    function func(time_switch)
        # predict the rest part of the passage heat load
        T = m.planet.T #fixed temperature
        t_cf =closed_form(args, m, position, T, true, m.aerodynamics.α)[1] # closed-form solution only for length of t

        mask = (t_cf .< time_switch_1) .|| (t_cf .> time_switch)
        aoa_list = m.aerodynamics.α .* mask
        
        t_cf, h_cf, γ_cf, v_cf = closed_form(args, m, position, T, true, m.aerodynamics.α, aoa_list)

        RT = T*m.planet.R
        S = v_cf / sqrt(2*RT)

        index_tilltsw = (t_cf .> t) .&& (t_cf .<= time_switch)
        index_remaining = t_cf .> time_switch

        t_tilltsw = t_cf .* index_tilltsw
        t_tilltsw = filter(!iszero, t_tilltsw)

        t_remaining = t_cf .* index_remaining
        t_remaining = filter(!iszero, t_remaining)

        if length(t_tilltsw) + length(t_remaining) == 0
            return Q_past
        end

        ρ = density_exp(h_cf, m.planet)[1]

        ρ_tilltsw = ρ[index_tilltsw]
        ρ_remaining = ρ[index_remaining]

        aoa_cf = zeros(length(t_tilltsw))
        S_tilltsw = S[index_tilltsw]

        heat_rate_tilltsw = heat_rate_calc(args[:multiplicative_factor_heatload] * m.aerodynamics.thermal_accomodation_factor, ρ_tilltsw, T, T, m.planet.R, m.planet.γ, S_tilltsw, aoa_cf)

        if length(t_tilltsw) > 1
            Q_rate_tilltsw = sum(heat_rate_tilltsw) * (t_tilltsw[end] - t_tilltsw[1]) / length(t_tilltsw)
        else
            Q_rate_tilltsw = 0
        end

        S_remaining = S[index_remaining]
        aoa_cf = ones(length(t_remaining)) * m.aerodynamics.α
        heat_rate_remaining = heat_rate_calc(args[:multiplicative_factor_heatload] * m.aerodynamics.thermal_accomodation_factor, ρ_remaining, T, T, m.planet.R, m.planet.γ, S_remaining, aoa_cf)

        if args[:control_mode] == 3 # Account for max heat rate possible
            index_hr = heat_rate_remaining .> args[:max_heat_rate]
            index_hr_not = heat_rate_remaining .<= args[:max_heat_rate]
            temp = args[:max_heat_rate] .* index_hr
            heat_rate_remaining = heat_rate_remaining .* index_hr_not + temp
        end

        if length(t_remaining) > 1
            Q_rate_remaining = sum(heat_rate_remaining) * (t_remaining[end] - t_remaining[1]) / length(t_remaining)
        else
            Q_rate_remaining = 0
        end

        Q = Q_past + Q_rate_tilltsw + Q_rate_remaining

        return Q - m.aerodynamics.heat_load_limit
    end

    # Evaluates max Q
    delta_Q_switch = func(time_switch_2)
    delta_Q_current_time = func(t)

    if reevaluation_mode == 1
        x_tol = 0.1
    else
        x_tol = 0.01
    end

    delta_Q_current_time = abs(delta_Q_current_time - delta_Q_switch) # if the difference between the heat load at the current time and the switch is small but the time is too large, recheck
    if abs(delta_Q_switch) <= x_tol && !(delta_Q_current_time < 0.5 && abs(t - time_switch_2) > 20)
        return config.cnf.time_switch_1, config.cnf.time_switch_2
    elseif delta_Q_current_time < 0
        return config.cnf.time_switch_1, t
    elseif delta_Q_switch > 0
        mult = 1
    elseif delta_Q_switch < 0
        mult = -1
    end

    if reevaluation_mode == 1
        x_tol = 0.1
    else
        x_tol = 0.05
    end

    b = t + 200

    time_switch_2 = fzero(k -> func(k), [t, b], Roots.Brent())

    time_switch_2 -= time_switch_2*0.1

    return config.cnf.time_switch_1, time_switch_2
end