include("../../physical_models/Aerodynamic_models.jl")
include("Utils_timeswitch.jl")
include("../../utils/Closed_form_solution.jl")
include("../utils/Eoms.jl")

using Roots

function switch_calculation_with_integration(ip, m, position, args, t, heat_rate_control, reevaluation_mode, gram_atmosphere, current_position=0)
    # println("time_switch 1: ", time_switch)

    function func(k_cf)
        global time_switch, k

        k = k_cf/100
        y, time_switch = asim_ctrl(ip, m, t, position, args, k, heat_rate_control, true, gram_atmosphere)

        # println("time_switch 2: ", time_switch)

        Q = y[end,end]

        return Q - m.aerodynamics.heat_load_limit
    end

    delta_Q_max = func(5)
    delta_Q_min = func(0)

    # println("time_switch 3: ", time_switch)

    if m.planet.name == "mars"
        low_lim = 3.3
        high_lim = 4.5
    elseif m.planet.name == "venus"
        low_lim = 1.0
        high_lim = 10.0
    elseif m.planet.name == "earth"
        low_lim = 1.0
        high_lim = 10.0
    elseif m.planet.name == "titan"
        low_lim = 0.00001
        high_lim = 10.0 
    end

    if delta_Q_max * delta_Q_min < 0
        find_zero(k -> func(k), [low_lim, high_lim], Bisection(), rtol=1e-5)
        # println("time_switch 4: ", time_switch)
    elseif delta_Q_max < 0
        if args[:heat_load_sol] == 0 && args[:heat_load_sol] == 3
            return [0, 0]
        elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2
            return [0, 1000]
        end
    elseif delta_Q_min > 0
        if args[:heat_load_sol] == 0 && args[:heat_load_sol] == 3
            return [0, 1000]
        elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2
            return [0, 0]
        end
    end

    # println("time_switch 5: ", time_switch)

    return time_switch
end

function switch_calculation(ip, m, position, args, t, heat_rate_control, reevaluation_mode, current_position=0)

    # global t_cf, h_cf, γ_cf, v_cf

    T = m.planet.T

    t_cf, h_cf, γ_cf, v_cf = closed_form(args, m, position, T, true, m.aerodynamics.α)

    RT = T * m.planet.R
    S = (v_cf / sqrt(2 * RT))

    CL_90, CD_90 = aerodynamic_coefficient_fM(pi/2, m.body, T, S[1], m.aerodynamics, 0)
    CL_0, CD_0 = aerodynamic_coefficient_fM(0, m.body, T, S[1], m.aerodynamics, 0)

    CD_slope = (CD_90 - CD_0) / (pi/2)

    coeff = (CD_slope, CL_0, CD_0)

    # Evaluates max Q
    aoa_cf = aoa(m, 0.1, t_cf, h_cf, γ_cf, v_cf, coeff)[1]

    approx_sol = (t_cf, h_cf, γ_cf, v_cf)

    delta_Q_max = func(0.1, m, args, coeff, position, heat_rate_control, approx_sol, aoa_cf)
    delta_Q_min = func(0.0, m, args, coeff, position, heat_rate_control, approx_sol, zeros(length(aoa_cf)))

    if delta_Q_max * delta_Q_min < 0
        k_cf = fzero(k -> func(k, m, args, coeff, position, heat_rate_control, approx_sol, aoa_cf), [0.0, 0.1], Roots.Brent())
    elseif delta_Q_max < 0.0
        return [0.0, 0.0]
    elseif delta_Q_min > 0.0
        return [0.0, t_cf[end]/2]
    end

    t_cf, v_cf, γ_cf, h_cf = func(k_cf, m, args, coeff, position, heat_rate_control, approx_sol, aoa_cf, false, true)

    lambda_switch, lambdav = lambdas(m, aoa_cf, k_cf, t_cf, h_cf, γ_cf, v_cf, coeff)[1:2]

    index_array = lambdav .< lambda_switch

    temp = t_cf .* index_array

    temp = filter(!iszero, temp)

    t_switch = [temp[1], temp[end]]

    if abs(t_switch[end] - t_cf[end]) < 5
        t_switch[1] -= t_cf[end]*0.04
    elseif abs(t_switch[end] - t_cf[end]) < 60
        t_switch[1] -= 5
    end

    t_switch[2] -= t_switch[2]*0.1

    return t_switch
end