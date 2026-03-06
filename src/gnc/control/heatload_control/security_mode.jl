if !isdefined(@__MODULE__, :__legacy_montecarlo_perturbations_included__)
    include(joinpath(@__DIR__, "..", "..", "..", "mission", "campaigns", "montecarlo_perturbations.jl"))
    const __legacy_montecarlo_perturbations_included__ = true
end
if !isdefined(@__MODULE__, :__legacy_density_models_included__)
    include(joinpath(@__DIR__, "..", "..", "..", "environment", "atmosphere", "density_models.jl"))
    const __legacy_density_models_included__ = true
end
if !isdefined(@__MODULE__, :__legacy_closed_form_solution_included__)
    include(joinpath(@__DIR__, "..", "..", "..", "analysis", "reports", "closed_form_solution.jl"))
    const __legacy_closed_form_solution_included__ = true
end

if !isdefined(@__MODULE__, :_legacy_get_cnf)
    @inline function _legacy_get_cnf(args=nothing; cnf=nothing)
        if cnf !== nothing
            return cnf
        end
        if args isa AbstractDict && haskey(args, :cnf)
            return args[:cnf]
        end
        if (@isdefined config) && isdefined(config, :cnf)
            return getproperty(config, :cnf)
        end
        throw(ArgumentError("Legacy control state `cnf` not found. Pass `cnf=` or args[:cnf]."))
    end
end

function security_mode(ip, m, position, args, t, heat_rate_control=false; cnf=nothing)
    cnf_state = _legacy_get_cnf(args; cnf=cnf)
    T = m.planet.T

    t_cf, h_cf, γ_cf, v_cf = closed_form(args, m, position, T, true, m.aerodynamics.α)

    RT = T * m.planet.R

    S = v_cf / sqrt(2*RT)

    ρ = density_polyfit(h_cf/1e3, m.planet)[1]

    # Security mode
    aoa_cf_min = zeros(length(t_cf))
    heat_rate_min = heat_rate_calc(args[:multiplicative_factor_heatload] * m.aerodynamics.thermal_accomodation_factor, ρ, T, T, m.planet.R, m.planet.γ, S, aoa_cf_min)

    index_future = t_cf .> T
    heat_rate_min = heat_rate_min .* index_future
    traj_rate = t_cf[2] - t_cf[1]

    if (sum(heat_rate_min) * traj_rate) + cnf_state.heat_load_past > m.aerodynamics.heat_load_limit
        cnf_state.security_mode = true

        return [0, t_cf[end] + 10000] # switch angle of attack to 0 for the rest # security made on the check
    end

    return [cnf_state.time_switch_1, cnf_state.time_switch_2]
end
