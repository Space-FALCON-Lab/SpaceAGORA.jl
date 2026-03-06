include(joinpath(@__DIR__, "..", "..", "internal", "bridge_helpers.jl"))
_bridge_include_heatload_security_deps!()

include(joinpath(@__DIR__, "..", "..", "internal", "bridge_helpers.jl"))

function security_mode(ip, m, position, args, t, heat_rate_control=false; cnf=nothing)
    cnf_state = _bridge_get_cnf(args; cnf=cnf)
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
