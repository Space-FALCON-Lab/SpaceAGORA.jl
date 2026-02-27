using .SimulationModel

# include("../utils/Define_mission.jl")
# include("../utils/MonteCarlo_set.jl")
# include("../utils/Initial_cond_calc.jl")
include("Set_and_run.jl")

function execute_orbital_elements_campaign(args::SimulationConfiguration; isolate_state::Bool=true)
    return execute_campaign(args; isolate_state=isolate_state)
end

function execute_vgamma_campaign(args::SimulationConfiguration; isolate_state::Bool=true)
    return execute_campaign(args; isolate_state=isolate_state)
end

@inline function _legacy_run_planet_id(selector)::Int
    if selector isa Integer
        return Int(selector)
    end
    key = lowercase(strip(string(selector)))
    key == "earth" && return 0
    key == "mars" && return 1
    key == "venus" && return 2
    key == "sun" && return 3
    key == "moon" && return 4
    key == "jupiter" && return 5
    key == "saturn" && return 6
    key == "titan" && return 7
    throw(ArgumentError("Unsupported legacy planet selector $(repr(selector))."))
end

function _legacy_typed_planet(selector, args)
    spice_path = get(args, :directory_spice, get(args, :spice_path, "GRAM Suite 2.0/SPICE"))
    pid = _legacy_run_planet_id(selector)
    pid == 0 && return Earth("", spice_path)
    pid == 1 && return Mars("", spice_path)
    pid == 2 && return Venus("", spice_path)
    pid == 7 && return Titan("", spice_path)
    # Bodies without typed constructors in current API: return minimal fields needed
    # by `ic_calculation_rptoae` / `ic_calculation_ae`.
    pid == 3 && return (Rp_e=6.9634e8, μ=1.3271244002331e20)
    pid == 4 && return (Rp_e=1.7381e6, μ=4.9028005821478e12)
    pid == 5 && return (Rp_e=7.1492e7, μ=1.26686534e17)
    pid == 6 && return (Rp_e=6.0268e7, μ=3.7931187e16)
    throw(ArgumentError("Unsupported legacy planet selector $(repr(selector))."))
end

function execute_orbital_elements_campaign(args)
    apoapsis, periapsis_alt, inclination, Ω, ω = collect(range(start=round(args[:ra_initial_a]), stop=round(args[:ra_initial_b]), step=round(args[:ra_step]))), 
                                                 collect(range(start=round(args[:hp_initial_a]), stop=round(args[:hp_initial_b]), step=round(args[:hp_step]))), 
                                                 args[:inclination], args[:Ω], args[:ω]
    if args[:initial_condition_type] == 2
        apoapsis, periapsis_alt = ic_calculation_ae(apoapsis, periapsis_alt, args)
    end
    final_apoapsis = args[:final_apoapsis]

    for periapsis_item in periapsis_alt
        for apoapsis_item in apoapsis

            if Bool(args[:print_res])
                println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
            end

            state = Dict()

            MC, count, args = MonteCarlo_setting(args)

            for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size], step=1)
                state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

                args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_ra=" * string(Float64(apoapsis_item/1e3)) * "_rp=" * string(Float64(periapsis_item/1e3)) * "_hl=" * string(args[:max_heat_rate]) * "_" * string(args[:α])

                if args[:montecarlo] == true
                    args = MonteCarlo_setting_passage(mc_index, args)
                end

                execute_campaign(args, state)
                count = MonteCarlo_append(MC, args, count)
            end

            if args[:montecarlo] == true
                MonteCarlo_save(args, state, MC)
            end
        end
    end
end

function execute_vgamma_campaign(args)
    γ_0, v_0, inclination, Ω, ω = collect(range(start=round(args[:γ_initial_a]*100), stop=round(args[:γ_initial_a]*100), step=round(args[:γ_step]*100))), 
                                  collect(range(start=round(args[:v_initial_a]), stop=round(args[:v_initial_b]), step=round(args[:v_step]))),
                                  args[:inclination], args[:Ω], args[:ω]
    final_apoapsis = args[:final_apoapsis]

    for γ in γ_0
        γ = -γ/100

        for v in v_0
            state = Dict()
            planet = _legacy_typed_planet(args[:planet], args)
            apoapsis, periapsis_alt = ic_calculation_rptoae(planet, γ, v, args)

            if Bool(args[:print_res])
                println("Velocity: " * string(v) * " m/s, Flight-Path Angle: " * string(γ) * " deg")
            end

            MC, count, args = MonteCarlo_setting(args)

            for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size])
                state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis, periapsis_alt * 1e-3, inclination, Ω, ω, final_apoapsis

                args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_v=" * string(v) * "_gamma=" * string(abs(γ)) * "_" * string(args[:α]) * "deg"

                if args[:montecarlo] == true
                    args = MonteCarlo_setting_passage(mc_index, args)
                end

                execute_campaign(args, state)
                count = MonteCarlo_append(MC, args, count)
            end

            if args[:montecarlo] == true
                MonteCarlo_save(args, state, MC)
            end
        end
    end
end

function execute_ae_campaign(args::SimulationConfiguration)
    # a, e, inclination, Ω, ω = collect(range(start=round(args[:a_initial_a]), stop=round(args[:a_initial_b]), step=round(args[:a_step]))), 
    #                                              collect(range(start=args[:e_initial_a], stop=args[:e_initial_b], step=args[:e_step])), 
    #                                              args[:inclination], args[:Ω], args[:ω]

    # final_apoapsis = args[:final_apoapsis]
    execute_campaign(args)
    # planet = planet_data(args.planet)
    # for e_item in e
    #     for a_item in a
    #         apoapsis_item, periapsis_item = ic_calculation_ae(planet, a_item, e_item, args)
    #         if Bool(args[:print_res])
    #             println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
    #         end

    #         state = Dict()

    #         MC, count, args = MonteCarlo_setting(args)

    #         for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size], step=1)
    #             state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

    #             args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_ra=" * string(Int64(round(apoapsis_item/1e3))) * "_rp=" * string(Float64(periapsis_item/1e3)) * "_hl=" * string(args[:max_heat_rate])

    #             if args[:montecarlo] == true
    #                 args = MonteCarlo_setting_passage(mc_index, args)
    #             end
    #             println("State: ", state)
    #             execute_campaign(args, state)
    #             # MonteCarlo_append(MC, args, count)
    #         end

    #         if args[:montecarlo] == true
    #             MonteCarlo_save(args, state, MC)
    #         end
    #     end
    # end
end

function execute_analysis(args::SimulationConfiguration; isolate_state::Bool=true)
    # args = def_miss(args)
    execute_campaign(args; isolate_state=isolate_state)
    # if args[:initial_condition_type] == 1 && (Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone")
    #     execute_vgamma_campaign(args)
    # elseif args[:initial_condition_type] == 0
    #     execute_orbital_elements_campaign(args)
    # elseif args[:initial_condition_type] == 2
    #     execute_ae_campaign(args)
    # end

    # if Bool(args[:passresults])
    #     return config.solution
    # else
    #     return true
    # end
end
