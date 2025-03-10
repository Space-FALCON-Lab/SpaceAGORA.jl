include("../../physical_models/Density_models.jl")
include("../../utils/Closed_form_solution.jl")

function lambdas(m, aoa, k, t, h, γ, v, coeff)
    CD_slope, CL_0, CD_0 = coeff
    S_ref = m.body.area_tot
    mass = m.body.mass
    Rp = m.planet.Rp_e
    μ = m.planet.μ
    g0 = m.planet.g_ref
    H = m.planet.H

    # evaluate lambda switch
    lambda_switch = (k * 2 * mass * v) ./ (S_ref * CD_slope * pi)

    lambdav = zeros(length(t))
    lambdag = zeros(length(t))
    lambdah = zeros(length(t))

    lambdav[end] = v[end]
    lambdag[end] = 0
    lambdah[end] = μ / (Rp + h[end])^2

    ρ = density_exp(h, m.planet)[1]

    g = g0 * Rp^2 ./ (Rp .+ h).^2

    for ii in collect(range(start=length(t), step=-1, stop=2))
        lambdav_dot = -3 * k * ρ[ii - 1] * v[ii - 1]^2 * aoa[ii - 1] / pi + lambdav[ii] * (ρ[ii - 1] * S_ref * (CD_0 + aoa[ii - 1] * CD_slope) * v[ii - 1]) / mass - lambdag[ii] * ((ρ[ii - 1] * S_ref * CL_0) / (2 * mass) + g[ii - 1] / v[ii - 1]^2 + 1 / (Rp + h[ii - 1])) - lambdah[ii] * γ[ii - 1]
        lambdag_dot = lambdav[ii] * g[ii - 1] - lambdah[ii] * v[ii - 1]
        lambdah_dot = k * ρ[ii - 1] * v[ii - 1]^3 * aoa[ii - 1] / (pi * H) - lambdav[ii] * ((ρ[ii - 1] * S_ref * (CD_0 + aoa[ii - 1] * CD_slope) * v[ii - 1]^2) / (2 * mass * H) + 2 * g[ii - 1] * γ[ii - 1] / (Rp + h[ii - 1])) + lambdag[ii] * (ρ[ii - 1] * S_ref * CL_0 * v[ii - 1] / (2 * mass * H) - 2 * g[ii - 1] / ((Rp + h[ii - 1]) * v[ii - 1]) + v[ii - 1] / (Rp + h[ii - 1])^2)

        lambdav[ii - 1] = lambdav[ii] - lambdav_dot * (t[ii] - t[ii - 1])
        lambdag[ii - 1] = lambdag[ii] - lambdag_dot * (t[ii] - t[ii - 1])
        lambdah[ii - 1] = lambdah[ii] - lambdah_dot * (t[ii] - t[ii - 1])
    end

    in_cond_lambda = [lambdav[2], lambdag[2], lambdah[2]]

    return lambda_switch, lambdav, in_cond_lambda
end

function aoa(m, k_cf, t_cf, h_cf, γ_cf, v_cf, coeff, aoa_cf=[])
    if length(aoa_cf) == 0
        aoa_cf = ones(length(t_cf)) * m.aerodynamics.α
    end

    lambda_switch, lambda_v, in_cond_lambda = lambdas(m, aoa_cf, k_cf, t_cf, h_cf, γ_cf, v_cf, coeff)

    aoa_cf = ones(length(lambda_v)) * m.aerodynamics.α
    index_array = lambda_v .>= lambda_switch
    aoa_cf = aoa_cf .* index_array

    return aoa_cf, in_cond_lambda
end

function func(k_cf, m, args, coeff, position, heat_rate_control, approx_sol, aoa_cf, initial_guess=false, approx_calc=false)
    t_cf, h_cf, γ_cf, v_cf = approx_sol

    Q = 0
    temp_Q = 1000
    Q_prec = 0

    count = 0
    
    while temp_Q > 1e-3 # stop if two following trajectory evaluation provides the same Q
        # define angle of attack lagrangian multipliers
        T = m.planet.T  # fixed temperature
        aoa_cf, in_cond_lambda = aoa(m, k_cf, t_cf, h_cf, γ_cf, v_cf, coeff, aoa_cf)  # update angle of attack profile with new k
        t_cf, h_cf, γ_cf, v_cf = closed_form(args, m, position, T, true, m.aerodynamics.α, aoa_cf)  # re-evaluate the closed form solution using previous angle of attack profile

        a = sqrt(m.planet.γ * m.planet.R * T)
        M = v_cf / a
        S = sqrt(m.planet.γ) * M
        ρ = density_exp(h_cf, m.planet)[1]  # density calculated through exponential density

        heat_rate = heat_rate_calc(args[:multiplicative_factor_heatload] * m.aerodynamics.thermal_accomodation_factor, ρ, T, T, m.planet.R, m.planet.γ, S, aoa_cf)

        if heat_rate_control == true # Account for max heat rate possible
            index_hr = heat_rate .> args[:max_heat_rate]
            heat_rate[index_hr] = args[:max_heat_rate]
        end

        Q = sum(heat_rate) * t_cf[end] / length(t_cf)

        # update error
        temp_Q = Q - Q_prec

        Q_prec = Q
        count += 1
    end

    if initial_guess == true
        return in_cond_lambda
    elseif approx_calc == true
        return t_cf, v_cf, γ_cf, h_cf
    end

    delta_Q = Q - m.aerodynamics.heat_load_limit

    return delta_Q
end