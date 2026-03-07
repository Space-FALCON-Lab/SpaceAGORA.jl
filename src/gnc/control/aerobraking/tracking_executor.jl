using SpecialFunctions
using Roots
using Logging

@inline function _compat_control_log_enabled(args)::Bool
    if get(ENV, "SPACEAGORA_DEBUG_LEGACY_CONTROL", "0") == "1"
        return true
    end
    if hasproperty(args, :simulation_settings) && hasproperty(args.simulation_settings, :verbose)
        return Bool(getproperty(args.simulation_settings, :verbose))
    end
    if hasproperty(args, :verbose)
        return Bool(getproperty(args, :verbose))
    end
    return false
end

@inline function _compat_control_strict_exceptions(args)::Bool
    if get(ENV, "SPACEAGORA_STRICT_LEGACY_CONTROL_EXCEPTIONS", "0") == "1"
        return true
    end
    if hasproperty(args, :strict_compat_control_exceptions)
        return Bool(getproperty(args, :strict_compat_control_exceptions))
    end
    return false
end

@inline function _compat_control_exception_fallback(args, location::AbstractString, err, bt, fallback)
    if _compat_control_log_enabled(args)
        @warn "Legacy control fallback in $(location)." exception=(err, bt)
    end
    if _compat_control_strict_exceptions(args)
        throw(err)
    end
    return fallback
end

function no_control(ip, m, args=0, index_ratio=0, state=0, t=0, position=0, current_position=0, heat_rate_control=true)
    α = m.aerodynamics.α

    return α
end

function control_struct_load(ip, m, args, S, T_p, q, MonteCarlo=false)

    max_α = m.aerodynamics.α
    min_α = 0.0001

    area_tot = config.get_spacecraft_reference_area(m.body)

    CL90, CD90 = aerodynamic_coefficient_fM(pi/2, m.body, T_p, S, m.aerodynamics, MonteCarlo)

    drag_max = q * aerodynamic_coefficient_fM(max_α, m.body, T_p, S, m.aerodynamics, MonteCarlo)[2] * area_tot
    drag_min = q * aerodynamic_coefficient_fM(min_α, m.body, T_p, S, m.aerodynamics, MonteCarlo)[2] * area_tot
    drag_limit = args[:max_dyn_press] * CD90 * area_tot

    f(x) = q * aerodynamic_coefficient_fM(x, m.body, T_p, S, m.aerodynamics, MonteCarlo)[2] * area_tot - drag_limit

    # α = cnf_state.α # pi/2

    if (drag_max < drag_limit)
        α = max_α
    elseif (drag_min > drag_limit)
        α = min_α
    elseif (drag_max >= drag_limit) && (drag_min <= drag_limit)
        try
            α = find_zero(f, (0, pi/2), Roots.Bisection())
        catch err
            α = _compat_control_exception_fallback(args, "control_struct_load.find_zero", err, catch_backtrace(), min_α)
        end
    else
        if _compat_control_log_enabled(args)
            println("Check Controller - Second Check")
        end
    end

    if (α > max_α) || (α < 0)
        α = 0
    end

    # Update solar panel angle
    # Assumes that the spacecraft is the standard 2 panels one bus
    # root = m.body.roots[1]
    # bodies, root_index = config.traverse_bodies(m.body, root)
    # for body in bodies
    #     if !body.root
    #         axis = SVector{3, Float64}(abs.(body.r))
    #         # Rotate the solar panel to the angle α
    #         config.rotate_link(body, axis, - α + m.body.roots[root_index].α)
    #     end
    # end

    return α

end

function control_solarpanels_heatrate(ip, m, args, index_ratio, state, t=0, position=0, current_position=0; cnf=nothing)
    cnf_state = _bridge_get_cnf(args; cnf=cnf)
    # println(state)
    # α = nothing
    if index_ratio[1] == 1
        T_p = state[1]
        ρ = state[2]
        S = state[3]

        if length(args) != 0
            if Bool(get(args, :montecarlo, false))
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

        # println("Heat Rate Max: ", heat_rate_max)
        # println("Heat Rate Min: ", heat_rate_min)

        L = (taf * ρ * R * T_p) .* (sqrt(R * T_p / (2*pi))) * 1e-4

        f(x) = L .* ((S.^2 .+ (γ) / (γ - 1) - (γ + 1) / (2 * (γ - 1)) * (T_w ./ T_p)) * (exp.(-(S .* sin.(x)).^2) + (pi^0.5) * (S .* sin.(x)) * (1 + erf.(S .* sin.(x)))) - 0.5 * exp.(-(S .* sin.(x)).^2)) .- thermal_limit

        α = cnf_state.α # pi/2

        if (heat_rate_max < thermal_limit)
            α = max_α
        elseif (heat_rate_min > thermal_limit)
            α = min_α
        elseif (heat_rate_max >= thermal_limit) && (heat_rate_min <= thermal_limit)
            x_0 = cnf_state.α_past
            df(x) = L .* S * cos.(x) * ((pi^0.5) * (S.^2 .+ γ / (γ - 1) + (γ + 1) / (2 * (γ - 1)) .* T_w ./ T_p) .* (1 + erf.(S .* sin.(x))) .+ S .* sin.(x) * exp.(-(S .* sin.(x)).^2))

            # println("try")
            try
                α = find_zero((f, df), x_0, Roots.Newton())
                # println(find_zero((f, df), x_0, Roots.Newton()))

                if α < 0 || α > pi/2
                    α = find_zero((f, df), 1e-1, Roots.Newton())
                end

            catch err
                if _compat_control_log_enabled(args)
                    @warn "Legacy control Newton solve failed; trying alternate initial guess." exception=(err, catch_backtrace())
                end

            # if α < 0 || α > pi/2
                try
                    if abs(heat_rate_max - thermal_limit) < abs(heat_rate_min - thermal_limit)   # Newton method is unable to find a solution since there are multiple ones. We need to provide a good initial guess
                        x_0 = 2 * max_α / 3
                    elseif abs(heat_rate_max - thermal_limit) > abs(heat_rate_min - thermal_limit)
                        x_0 = 2 * max_α / 6
                        α = find_zero((f, df), x_0, Roots.Newton())
                    end
                catch err_retry
                    α = _compat_control_exception_fallback(args, "control_solarpanels_heatrate.find_zero", err_retry, catch_backtrace(), min_α)
                # end
                end
            # end
            end

        else
            if _compat_control_log_enabled(args)
                println("Check Controller - Second Check")
            end
        end

        if (α > max_α) || (α < 0)
            α = 0
        end

        # println("control: ", α)
        # Update solar panel angle
        # Assumes that the spacecraft is the standard 2 panels one bus
        # root = m.body.roots[1]
        # bodies, root_index = config.traverse_bodies(m.body, root)
        # for body in bodies
        #     if !body.root
        #         axis = SVector{3, Float64}(abs.(body.r))
        #         # Rotate the solar panel to the angle α
        #         config.rotate_link(body, axis, - α + m.body.roots[root_index].α)
        #     end
        # end
        return α
    else
        return cnf_state.α
    end

end

function heat_rate_calc(taf, ρ, T_w, T_p, R, γ, S, angle)
    first_term = ρ .* (1e-4 * taf * R * T_p  * sqrt.(R * T_p / (2*pi)))

    term_a = exp.(-(S .* sin.(angle)).^2) .+ (sqrt(pi) * S .* sin.(angle) .* (1 .+ erf.(S .* sin.(angle))))

    term_b = (S.^2 .+ γ/(γ - 1) .- (γ + 1)/(2 * (γ - 1)) * (T_w./T_p)) .* term_a

    return (term_b .- 0.5 * exp.(-(S .* sin.(angle)).^2)) .* first_term
end

function control_solarpanels_heatload(ip, m, args, index_ratio, state=0, t=0, position=0, current_position=0, gram_atmosphere=nothing, heat_rate_control=false; cnf=nothing)
    lock(CONTROL_BRIDGE_STATE_LOCK)
    try
        return _control_solarpanels_heatload_impl(ip, m, args, index_ratio, state, t, position, current_position, gram_atmosphere, heat_rate_control; cnf=cnf)
    finally
        unlock(CONTROL_BRIDGE_STATE_LOCK)
    end
end

function _control_solarpanels_heatload_impl(ip, m, args, index_ratio, state=0, t=0, position=0, current_position=0, gram_atmosphere=nothing, heat_rate_control=false; cnf=nothing)
    cnf_state = _bridge_get_cnf(args; cnf=cnf)
    policy_selector = DefaultAerobrakingPolicySelector()
    policy_config = AerobrakingPolicyConfig()
    guidance_input = AerobrakingGuidanceInput(
        ip=ip,
        mission=m,
        args=args,
        index_ratio=index_ratio isa AbstractVector{<:Integer} ? collect(Int, index_ratio) : Int[],
        state=state,
        t=Float64(t),
        position=position,
        current_position=current_position,
        gram_atmosphere=gram_atmosphere,
        heat_rate_control=Bool(heat_rate_control),
        cnf=cnf_state,
    )
    guidance_output = dispatch_aerobraking_guidance(policy_selector, policy_config, guidance_input)
    cnf_state.time_switch_1 = guidance_output.time_switch_1
    cnf_state.time_switch_2 = guidance_output.time_switch_2
    cnf_state.security_mode = guidance_output.security_mode

    # println("Time Switches: ", cnf_state.time_switch_1, " , ", cnf_state.time_switch_2)

    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 3
        if (t > cnf_state.time_switch_1) && (t < cnf_state.time_switch_2)
            α = 0
        else
            α = m.aerodynamics.α
        end
    elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2
        if (t > cnf_state.time_switch_1) && (t < cnf_state.time_switch_2)
            α = m.aerodynamics.α
        else
            α = 0
        end
    end

    if isempty(cnf_state.heat_load_ppast)
        cnf_state.heat_load_ppast = zeros(length(m.body.links))
    end
    cnf_state.heat_load_ppast .= cnf_state.heat_load_past

    # Update solar panel angle
    # Assumes that the spacecraft is the standard 2 panels one bus
    root = m.body.roots[1]
    bodies, root_index = config.traverse_bodies(m.body, root)
    for body in bodies
        if !body.root
            axis = SVector{3, Float64}(abs.(body.r))
            # Rotate the solar panel to the angle α
            config.rotate_link(body, axis, - α + m.body.roots[root_index].α)
        end
    end

    return α
end

function control_solarpanels_openloop(ip, m, args, index_ratio, state, t=0, position=0, current_position=0, heat_rate_control=true, gram_atmosphere=nothing; cnf=nothing)
    lock(CONTROL_BRIDGE_STATE_LOCK)
    try
        return _control_solarpanels_openloop_impl(ip, m, args, index_ratio, state, t, position, current_position, heat_rate_control, gram_atmosphere; cnf=cnf)
    finally
        unlock(CONTROL_BRIDGE_STATE_LOCK)
    end
end

function _control_solarpanels_openloop_impl(ip, m, args, index_ratio, state, t=0, position=0, current_position=0, heat_rate_control=true, gram_atmosphere=nothing; cnf=nothing)
    cnf_state = _bridge_get_cnf(args; cnf=cnf)
    control_solarpanels_heatload(ip, m, args, index_ratio, 0, t, position, current_position, gram_atmosphere, heat_rate_control; cnf=cnf_state)

    if args[:heat_load_sol] == 0 || args[:heat_load_sol] == 3
        if t >= cnf_state.time_switch_1 && t <= cnf_state.time_switch_2
            α = 0
        else
            α = control_solarpanels_heatrate(ip, m, args, index_ratio, state, t; cnf=cnf_state)
        end
    elseif args[:heat_load_sol] == 1 || args[:heat_load_sol] == 2
        if t >= cnf_state.time_switch_1 && t <= cnf_state.time_switch_2
            α = control_solarpanels_heatrate(ip, m, args, index_ratio, state, t; cnf=cnf_state)
        else
            α = 0
        end
    end

    # cnf_state.α_past = α
    # # Update solar panel angle
    # # Assumes that the spacecraft is the standard 2 panels one bus
    # root = m.body.roots[1]
    # bodies, root_index = config.traverse_bodies(m.body, root)
    # for body in bodies
    #     if !body.root
    #         axis = SVector{3, Float64}(abs.(body.r))
    #         # Rotate the solar panel to the angle α
    #         config.rotate_link(body, axis, - α + m.body.roots[root_index].α)
    #     end
    # end
    return α
end

# function targeting_solarpanels(ip, m, args, index_ratio, state, t=0, position=0, current_position=0, heat_rate_control=true, gram_atmosphere=nothing)

#     t_switch = 

#     if t > t_switch 

#     else
#         α = control_solarpanels_heatrate(ip, m, args, index_ratio, state, t; cnf=cnf_state)
#     end
    
# end
