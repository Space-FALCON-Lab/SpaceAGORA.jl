include("../utils/Closed_form_solution.jl")
include("../physical_models/Density_models.jl")
include("../physical_models/MonteCarlo_pertrubations.jl")
include("heatload_control/Time_switch_calcs.jl")
# include("heatload_control/Second_tsw_calcs.jl")
# include("heatload_control/Security_mode.jl")

import .config

using SpecialFunctions
using Roots

function no_control(ip, m, args=0, index_ratio=0, state=0, t=0, position=0, current_position=0, heat_rate_control=true)
    α = m.aerodynamics.α

    return α
end

function control_solarpanels_heatrate(ip, m, args, index_ratio, state, t=0, position=0, current_position=0)
    # println(state)
    # α = nothing
    if index_ratio[1] == 1
        T_p = state[1]
        ρ = state[2]
        S = state[3]

        if args != 0
            if args[:montecarlo] == true
                ρ, T_p, S = monte_carlo_guidance_environment(ρ, T_p, S, args)
            end
        end

        T_w = T_p

        taf = m.aerodynamics.thermal_accomodation_factor
        R = m.planet.R
        γ = m.planet.γ

        max_α = m.aerodynamics.α
        min_α = 0.0001

        thermal_limit = m.aerodynamics.heat_rate_limit - 0.00001

        heat_rate_max = heat_rate_calc(taf, ρ, T_w, T_p, R, γ, S, max_α)
        heat_rate_min = heat_rate_calc(taf, ρ, T_w, T_p, R, γ, S, min_α)

        L = (taf * ρ * R * T_p) * (sqrt(R * T_p / (2*pi))) * 1e-4

        f(x) = L * ((S^2 + (γ) / (γ - 1) - (γ + 1) / (2 * (γ - 1)) * (T_w / T_p)) * (exp(-(S * sin(x))^2) + (pi^0.5) * (S * sin(x)) * (1 + erf(S * sin(x)))) - 0.5 * exp(-(S * sin(x))^2)) - thermal_limit

        if (heat_rate_max < thermal_limit)
            α = max_α
        elseif (heat_rate_min > thermal_limit)
            α = min_α
        elseif (heat_rate_max >= thermal_limit) && (heat_rate_min <= thermal_limit)
            x_0 = config.cnf.α_past
            try
                df(x) = L * S * cos(x) * ((pi^0.5) * (S^2 + γ / (γ - 1) + (γ + 1) / (2 * (γ - 1)) * T_w / T_p) * (1 + erf(S * sin(x))) + S * sin(x) * exp(-(S * sin(x))^2))
                α = find_zero((f, df), x_0, Roots.Newton())

                if α < 0 || α > pi/2
                    α = find_zero((f, df), 1e-1, Roots.Newton())
                end

                if α < 0 || α > pi/2
                    error("False Exception")
                end 
            catch
                try
                    if abs(heat_rate_max - thermal_limit) < abs(heat_rate_min - thermal_limit)   # Newton method is unable to find a solution since there are multiple ones. We need to provide a good initial guess
                        x_0 = 2 * max_α / 3
                    elseif abs(heat_rate_max - thermal_limit) > abs(heat_rate_min - thermal_limit)
                        x_0 = 2 * max_α / 6
                        df(x) = L * S * cos(x) * ((pi^0.5) * (S^2 + γ / (γ - 1) + (γ + 1) / (2 * (γ - 1)) * T_w / T_p) * (1 + erf(S * sin(x))) + S * sin(x) * exp(-(S * sin(x))^2))
                        α = find_zero((f, df), x_0, Roots.Newton())
                    end
                catch
                    println("Check - heat rate controller does not converge")
                end
            end
        else
            println("Check Controller - Second Check")

            if (α > max_α) || (α < 0)
                α = 0
            end
        end

        # if α === nothing
        #     println("alpha is cooked")
        # end

        return α
    else
        return config.cnf.α
    end

end

function heat_rate_calc(taf, ρ, T_w, T_p, R, γ, S, angle)
    first_term = ρ .* (1e-4 * taf * R * T_p  * sqrt(R * T_p / (2*pi)))

    term_a = exp(-(S .* sin.(angle)).^2) + (sqrt(pi) * S .* sin.(angle) .* (1 + erf.(S .* sin.(angle))))

    term_b = (S^2 + γ/(γ - 1) - (γ + 1)/(2 * (γ - 1)) * (T_w/T_p)) .* term_a

    return (term_b - 0.5 * exp(-(S .* sin.(angle)).^2)) .* first_term
end

function control_solarpanels_heatload(ip, m, args, index_ratio, state=0, t=0, position=0, current_position=0, heat_rate_control=false)
    # IF increase too much the 50 sec, the solution becomes a bit unstable because of Mars Gram
    start_reevaluation = false

    if config.cnf.ascending_phase || (config.cnf.ascending_phase == false && abs(config.cnf.time_switch_2 - t) < 0.2 * config.cnf.time_switch_2)  # in periapsis passed or the time switch is too close to the current time but periapsis not passed, start reevaluation
        start_reevaluation = true
    end

    if (config.cnf.evaluate_switch_heat_load == false)
        if args[:flash2_through_integration] == 1
            if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 1
                config.cnf.time_switch_1, config.cnf.time_switch_2 = switch_calculation_with_integration(ip, m, position, args, t, heat_rate_control, 1, position)
            end
            if args[:heat_load_sol] == 2 || args[:heat_load_sol] == 3
                config.cnf.time_switch_1, config.cnf.time_switch_2 = second_time_switch_recalc_with_integration(ip, m, position, args, t, heat_rate_control, 1, position)
            end

            args[:security_mode] = 0
        else
            config.cnf.time_switch_1, config.cnf.time_switch_2 = switch_calculation(ip, m, position, args, t, heat_rate_control, 1, position)
        end
        
        config.cnf.evaluate_switch_heat_load = true
    elseif config.cnf.evaluate_switch_heat_load == true && start_reevaluation && ((t - config.cnf.time_switch_1 > 1 && t - config.cnf.timer_revaluation > 10 && t - config.cnf.time_switch_2 < 0) || (3 < config.cnf.time_switch_2 - t < 50 && t - config.cnf.timer_revaluation > 3) || (0 < config.cnf.time_switch_2 - t < 3 && t - config.cnf.timer_revaluation > 0.8)) && config.cnf.security_mode == false
        if (t - config.cnf.timer_revaluation) > 3
            reevaluation_mode = 1
        else
            reevaluation_mode = 2
        end

        config.cnf.timer_revaluation = t

        if args[:second_switch_reevaluation] == 1
            if args[:flash2_through_integration] == 1
                config.cnf.time_switch_1, config.cnf.time_switch_2 = second_time_switch_recalc_with_integration(ip, m, position, args, t, heat_rate_control, reevaluation_mode, current_position)
    
                args[:security_mode] = 0
            else
                config.cnf.time_switch_1, config.cnf.time_switch_2 = second_time_switch_recalc(ip, m, position, args, t, heat_rate_control, reevaluation_mode, current_position)
            end

        end
    elseif config.cnf.heat_load_past > 0.98 && m.aerodynamics.heat_load_limit && (config.cnf.heat_load_past - config.cnf.heat_load_ppast) < 2 && args[:security_mode] == 1 && config.snf.security_mode == false && index_ratio[2] == 1
        config.cnf.time_switch_1, config.cnf.time_switch_2 = security_mode(ip, m, position, args, t, false)
    end

    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 3
        if (t > config.cnf.time_switch_1) && (t < config.cnf.time_switch_2)
            α = 0
        else
            α = m.aerodynamics.α
        end
    elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2
        if (t > config.cnf.time_switch_1) && (t < config.cnf.time_switch_2)
            α = m.aerodynamics.α
        else
            α = 0
        end
    end

    config.cnf.heat_load_ppast = config.cnf.heat_load_past

    return α
end

function control_solarpanels_openloop(ip, m, args, index_ratio, state, t=0, position=0, current_position=0, heat_rate_control=true)
    control_solarpanels_heatload(ip, m, args, index_ratio, t, position, current_position, heat_rate_control, state)

    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 3
        if t >= config.cnf.time_switch_1 && t <= config.cnf.time_switch_2
            α = 0
        else
            α = control_solarpanels_heatrate(ip, m, index_ratio, state, args, t)
        end
    elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2
        if t >= config.cnf.time_switch_1 && t <= config.cnf.time_switch_2
            α = control_solarpanels_heatrate(ip, m, index_ratio, state, args, t)
        else
            α = 0
        end
    end

    config.cnf.α_past = α

    return α
end