include("../utils/Closed_form_solution.jl")
include("../physical_models/Density_models.jl")
include("../physical_models/MonteCarlo_pertrubations.jl")

import .config
using SpecialFunctions

function no_control(ip, args=0, ρ=0, T_p=0, T_w=0, S=0, position=0, current_position=0, t=0, index_ratio=0, state=0)
    α = m.aerodynamics.α

    return α
end

function control_solarpanels_heatrate(ip, m, index_ratio, state, args=0, t=0, position=0, current_position=0)
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
            x0 = config.cnf.α_past
            try
                α = 

                if α < 0 || α > pi/2

                end

                if α < 0 || α > pi/2
                    
                end
                    
            catch
            end
        else
            println("Check Controller - Second Check")
        end

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

function control_solarpanels_heatload(ip, m, args, index_ratio, t=0, position=0, crrent_position=0, heat_rate_control=false, state=0)
    if args[:flash2_through_integration] == 1
        if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 1

        end
    end
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

    config.cnf.α = α

    return α
end