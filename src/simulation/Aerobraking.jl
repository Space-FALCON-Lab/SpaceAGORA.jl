include("Complete_passage.jl")
include("../utils/Ref_system_conf.jl")
include("../utils/Closed_form_solution.jl")
include("../utils/Odyssey_maneuver_plan.jl")
include("../utils/Save_results.jl")

import .config

function aerobraking(ip, m, args)
    initial_state = m.initial_condition
    FinalState = true
    continue_campaign = true
    numberofpassage = 0
    config.cnf.count_numberofpassage = 0

    clean_results()

    config.cnf.time_OP = 0
    config.cnf.time_IP = 0

    # Aerobraking Campaign
    while continue_campaign && FinalState
        config.index_Mars_Gram_call = 0
        firing_orbit = 0
        numberofpassage = 1 + numberofpassage

        if args[:print_res]
            println("--> Start Passage #" * string(numberofpassage))
        end

        t = tok()

        if args[:Odyssey_sim]
            args = Odyssey_firing_plan(numberofpassage, args)
        end

        if args[:ip.tc] == 1
            if args[:delta_v] != 0.0
                if ip.tc == 1
                    initial_state = propulsion_ic_calcs(m, args, initial_state)
                end
            end
        elseif ip.tc == 2
            if round(rad2deg(args[:ϕ])) == 180
                println("DECELERATE DRAG FIRING!!")
            elseif round(rad2deg(args[:ϕ])) == 0
                println("ACCELERATE DRAG FIRING!!")
            end
        end

        if numberofpassage != 1
            # Orbtial Elements Results
            initial_state.a = config.solution.orientation.oe[1][end]
            initial_state.e = config.solution.orientation.oe[2][end]
            initial_state.i = config.solution.orientation.oe[3][end]
            initial_state.Ω = config.solution.orientation.oe[4][end]
            initial_state.ω = config.solution.orientation.oe[5][end]
            initial_state.m = config.solution.performance.mass[end]
            initial_state.vi = config.solution.orientation.oe[6][end]

            m.initial_condition.year = Int64(config.solution.orientation.year[end])
            m.initial_condition.month = Int64(config.solution.orientation.month[end])
            m.initial_condition.day = Int64(config.solution.orientation.day[end])
            m.initial_condition.hour = Int64(config.solution.orientation.hour[end])
            m.initial_condition.minute = Int64(config.solution.orientation.minute[end])
            m.initial_condition.second = config.solution.orientation.second[end]

            if (args[:drag_passage] || args[:body_shape] == "Blunted Cone") && continue_campaign
                r = m.planet.Rp_e + args[:EI]*1e3
                initial_state.vi = acos(1 / initial_state.e * (initial_state.a * (1 - initial_state.e^2) / r - 1))
            end
        end

        continue_campaign = asim(ip, m, initial_state, numberofpassage, args)

        r_a = config.solution.orientation.oe[1][end] * (1 + config.solution.orientation.oe[2][end])
        r_p = config.solution.orientation.oe[1][end] * (1 - config.solution.orientation.oe[2][end])
        elapsed = tok() - t

        if args[:print_res]
            println("Computational time: " * string(elapsed) * " seconds")
            println("--> PASSAGE #" * string(numberofpassage) * " COMPLETE")
        end

        if args[:number_of_orbits] == numberofpassage
            continue_campaign = false
        end

        if r_a <= args[:final_apoapsis]
            FinalState = false
            println("Reached FinalState! R_a = " * string(r_a*1e-3) * " km")
            println("Thermal Limit overcomed totally " * string(config.cnf.count_overcome_hr) * " times")
        end

        if r_p - m.planet.Rp_e >= 180*1e3
            FinalState = false
            println("Periapsis too high, final state unreachable! R_a = " * string(r_p*1e-3) * " km")
        end

        println(" ")
    end

    closed_form(args, m)
end