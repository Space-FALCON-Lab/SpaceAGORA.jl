using SpecialFunctions
using Roots

const _EDG_LEGACY_CONSTRAINT_MIN_ALPHA_RAD = 1e-4

@inline _edg_constraint_min_alpha(config) = max(config.min_alpha_rad, _EDG_LEGACY_CONSTRAINT_MIN_ALPHA_RAD)

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
    thermal_limit = heat_rate_limit - 1e-5
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
    alpha = try
        Roots.find_zero(f, (min_alpha, max_alpha), Roots.Brent())
    catch
        min_alpha
    end
    return (isfinite(alpha) && 0.0 <= alpha <= max_alpha) ? alpha : 0.0
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
    min_alpha = _edg_constraint_min_alpha(config)
    base_alpha < min_alpha && return base_alpha
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
        min_alpha=min_alpha,
        heat_rate_limit=limit,
        alpha_past=alpha_past,
    )
end
