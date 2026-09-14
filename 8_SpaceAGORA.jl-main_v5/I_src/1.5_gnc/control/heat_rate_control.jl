using SpecialFunctions
using Roots

function _energy_depletion_heat_rate_calc(
    taf::Float64,
    rho::Float64,
    T_w::Float64,
    T_p::Float64,
    R::Float64,
    gamma::Float64,
    S::Float64,
    alpha::Float64,
)::Float64
    first_term = rho * (1e-4 * taf * R * T_p * sqrt(R * T_p / (2pi)))
    s_sin = S * sin(alpha)
    exp_term = exp(-(s_sin)^2)
    term_a = exp_term + sqrt(pi) * s_sin * (1.0 + erf(s_sin))
    term_b = (S^2 + gamma / (gamma - 1.0) - (gamma + 1.0) / (2.0 * (gamma - 1.0)) * (T_w / T_p)) * term_a
    qdot = (term_b - 0.5 * exp_term) * first_term
    return isfinite(qdot) ? qdot : Inf
end

function _energy_depletion_heatrate_root_alpha(;
    taf::Float64,
    rho::Float64,
    T_p::Float64,
    R::Float64,
    gamma::Float64,
    S::Float64,
    max_alpha::Float64,
    min_alpha::Float64,
    heat_rate_limit::Float64,
    alpha_past::Float64,
)::Float64
    if !(isfinite(heat_rate_limit) && heat_rate_limit > 0.0)
        return max_alpha
    end
    if !(isfinite(rho) && isfinite(T_p) && isfinite(S) && rho > 0.0 && T_p > 0.0 && S > 0.0)
        return max_alpha
    end

    T_w = T_p
    thermal_margin = max(1e-5, 5e-4 * heat_rate_limit)
    thermal_limit = heat_rate_limit - thermal_margin
    thermal_limit > 0.0 || return min_alpha

    heat_rate_max = _energy_depletion_heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, max_alpha)
    heat_rate_min = _energy_depletion_heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, min_alpha)

    if heat_rate_max < thermal_limit
        return max_alpha
    elseif heat_rate_min > thermal_limit
        return min_alpha
    end

    L = (taf * rho * R * T_p) * sqrt(R * T_p / (2pi)) * 1e-4
    f(alpha) = begin
        s_sin = S * sin(alpha)
        exp_term = exp(-(s_sin)^2)
        term_a = exp_term + sqrt(pi) * s_sin * (1.0 + erf(s_sin))
        term_b = (S^2 + gamma / (gamma - 1.0) - (gamma + 1.0) / (2.0 * (gamma - 1.0)) * (T_w / T_p)) * term_a
        (term_b - 0.5 * exp_term) * L - thermal_limit
    end
    df(alpha) = begin
        s_sin = S * sin(alpha)
        L * S * cos(alpha) *
        (
            sqrt(pi) * (S^2 + gamma / (gamma - 1.0) + (gamma + 1.0) / (2.0 * (gamma - 1.0)) * (T_w / T_p)) *
            (1.0 + erf(s_sin)) +
            s_sin * exp(-(s_sin)^2)
        )
    end

    x0 = isfinite(alpha_past) && min_alpha <= alpha_past <= max_alpha ? alpha_past : 0.5 * (min_alpha + max_alpha)
    alpha = try
        Roots.find_zero((f, df), x0, Roots.Newton())
    catch
        NaN
    end

    if !(isfinite(alpha) && min_alpha <= alpha <= max_alpha)
        alpha = try
            Roots.find_zero((f, df), 1e-1, Roots.Newton())
        catch
            NaN
        end
    end

    if !(isfinite(alpha) && min_alpha <= alpha <= max_alpha)
        alpha = try
            if abs(heat_rate_max - thermal_limit) < abs(heat_rate_min - thermal_limit)
                Roots.find_zero((f, df), 2.0 * max_alpha / 3.0, Roots.Newton())
            else
                Roots.find_zero((f, df), 2.0 * max_alpha / 6.0, Roots.Newton())
            end
        catch
            NaN
        end
    end

    if !(isfinite(alpha) && min_alpha <= alpha <= max_alpha)
        alpha = Roots.find_zero(f, (min_alpha, max_alpha), Roots.Bisection())
    end

    return clamp(alpha, min_alpha, max_alpha)
end

function _edg_maxwellian_heat_rate(p::ODEParams, env, alpha::Float64)::Float64
    thermal_model = p.args.environment_model.thermal_model
    taf = hasproperty(thermal_model, :thermal_accomodation_factor) ? Float64(thermal_model.thermal_accomodation_factor) : 1.0
    planet = p.args.environment_model.planet
    S = env.molecular_speed_ratio
    T_p = env.temperature
    rho = env.rho
    if !(isfinite(S) && isfinite(T_p) && isfinite(rho) && S > 0.0 && T_p > 0.0 && rho > 0.0)
        return 0.0
    end
    gamma = planet.γ
    T_w = T_p
    first_term = rho * (1e-4 * taf * planet.R * T_p * sqrt(planet.R * T_p / (2pi)))
    s_sin = S * sin(alpha)
    exp_term = exp(-(s_sin)^2)
    erf_term = 1.0 + erf(s_sin)
    term_a = exp_term + sqrt(pi) * s_sin * erf_term
    term_b = (S^2 + gamma / (gamma - 1.0) - (gamma + 1.0) / (2.0 * (gamma - 1.0)) * (T_w / T_p)) * term_a
    qdot = (term_b - 0.5 * exp_term) * first_term
    return (isfinite(qdot) && qdot > 0.0) ? qdot : 0.0
end

function _edg_heat_rate_alpha(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    env,
    base_alpha::Float64;
    limit_override::Float64=config.heat_rate_limit_w_cm2,
    alpha_past::Float64=base_alpha,
)
    limit = limit_override
    if !(isfinite(limit) && limit > 0.0)
        return base_alpha
    end
    thermal_model = p.args.environment_model.thermal_model
    taf = hasproperty(thermal_model, :thermal_accomodation_factor) ? Float64(thermal_model.thermal_accomodation_factor) : 1.0
    planet = p.args.environment_model.planet
    return _energy_depletion_heatrate_root_alpha(
        taf=taf,
        rho=env.rho,
        T_p=env.temperature,
        R=planet.R,
        gamma=planet.γ,
        S=env.molecular_speed_ratio,
        max_alpha=base_alpha,
        min_alpha=config.min_alpha_rad,
        heat_rate_limit=limit,
        alpha_past=alpha_past,
    )
end
