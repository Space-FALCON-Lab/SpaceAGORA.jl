include("../utils/Define_mission.jl")
include("../utils/MonteCarlo_set.jl")


function run_orbitalelements(args)
    apoapsis, periapsis_alt, inclination, Ω, ω = [range(Int64(args[:ra_initial_a]), Int64(args[:ra_initial_b]), step=Int64(args[:ra_step]))], 
                                                 [range(Int64(args[:hp_initial_a]), Int64(args[:hp_initial_b]), step=Int64(args[:hp_step]))], 
                                                 args[:inclination], args[:Ω], args[:ω]
    
    final_apoapsis = args[:final_apoapsis]

    for periapsis_item in periapsis_alt
        for apoapsis_item in apoapsis
            if args[:print_res]
                println("Apoapsis Radius: " * string(apoapsis_item/10^3) * " km, Periapsis Altitude: " * string(periapsis_item/10^3) * " km")  
            end

            state = Dict()

            MC, count, args = MonteCarlo_setting(args)

            for mc_index in range(args[:initial_montecarlo_number], args[:montecarlo_size], step=1)
                state[:Apoapsis], state[:Periapsis], state[:Inclination], state[:Ω], state[:ω], state[:Final_sma] = apoapsis_item, Float64(periapsis_item*1e-3), inclination, Ω, ω, final_apoapsis

                args[:simulation_filename] = "Results_ctrl=" * string(args[:control_mode]) * "_ra=" * string(Int64(apoapsis_item/1e3)) * "_rp=" * string(Float64(periapsis_item/1e3)) * "_hl=" * string(args[:max_heat_rate]) * "_" * string(args[:angle_of_attack])

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
    
end

function run_analysis(args)

    args = def_miss(args)

    if args[:initial_condition_type] == 1 && (args[:drag_passage] || args[:body_shape] == "Blunted Cone")
        run_vgamma(args)
    else
        run_orbitalelements(args)
    end

    if args[:passresults]
        return config.solution
    else
        return Bool(1)
    end
end