include("../utils/Define_mission.jl")
include("../utils/MonteCarlo_set.jl")
include("../utils/Initial_cond_calc.jl")
include("Set_and_run.jl")

function run_orbitalelements(args)
    apoapsis, periapsis_alt, inclination, Ω, ω = collect(range(start=round(args[:ra_initial_a]), stop=round(args[:ra_initial_b]), step=round(args[:ra_step]))), 
                                                 collect(range(start=round(args[:hp_initial_a]), stop=round(args[:hp_initial_b]), step=round(args[:hp_step]))), 
                                                 args[:inclination], args[:Ω], args[:ω]
    
    final_apoapsis = args[:final_apoapsis]

    for periapsis_item in periapsis_alt
        for apoapsis_item in apoapsis

            if Bool(args[:print_res])
                println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
            end

            state = Dict()

            MC, count, args = MonteCarlo_setting(args)

            for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size]-1, step=1)
                state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

                args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_ra=" * string(Int64(round(apoapsis_item/1e3))) * "_rp=" * string(Float64(periapsis_item/1e3)) * "_hl=" * string(args[:max_heat_rate]) * "_" * string(args[:α])

                if args[:montecarlo] == true
                    args = MonteCarlo_setting_passage(mc_index, args)
                end

                aerobraking_campaign(args, state)
                MonteCarlo_append(MC, args, count)
            end

            if args[:montecarlo] == true
                MonteCarlo_save(args, state, MC)
            end
        end
    end
end

function run_vgamma(args)
    γ_0, v_0, inclination, Ω, ω = collect(range(start=round(args[:γ_initial_a]*100), stop=round(args[:γ_initial_a]*100), step=round(args[:γ_step]*100))), 
                                  collect(range(start=round(args[:v_initial_a]), stop=round(args[:v_initial_b]), step=round(args[:v_step]))),
                                  args[:inclination], args[:Ω], args[:ω]
    final_apoapsis = args[:final_apoapsis]

    for γ in γ_0
        γ = -γ

        for v in v_0
            state = Dict()
            planet = planet_data(args[:planet])
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

                aerobraking_campaign(args, state)
                MonteCarlo_append(MC, args, count)
            end

            if args[:montecarlo] == true
                MonteCarlo_save(args, state, MC)
            end
        end
    end
end

function run_analysis(args)

    args = def_miss(args)

    if args[:initial_condition_type] == 1 && (Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone")
        run_vgamma(args)
    else
        run_orbitalelements(args)
    end

    if Bool(args[:passresults])
        return config.solution
    else
        return true
    end
end