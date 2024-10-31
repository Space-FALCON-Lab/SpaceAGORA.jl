include("Complete_passage.jl")
include("../utils/Ref_system_conf.jl")
include("../utils/Closed_form_solution.jl")
include("../utils/Odyssey_maneuver_plan.jl")
include("../utils/Save_results.jl")
include("../config.jl")

function aerobraking(ip, m, args)
    inital_state = m.initialcondition
    FinalState = true
    continue_campaign = true
    numberofpassage = 0
    config.numberofpassage = 0

    clean_results()

    time_OP = 0
    time_IP = 0

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
        elseif ip.tc == 1
            if round(rad2deg(args[:ϕ])) == 180
                println("DECELERATE DRAG FIRING!!")
            elseif round(rad2deg(args[:ϕ])) == 0
                println("ACCELERATE DRAG FIRING!!")
            end
        end

        if numberofpassage != 1
            # Orbtial Elements Results
            initial_state.a = solution.orientation.oe[1][end]
            initial_state.e = solution.orientation.oe[2][end]
            initial_state.i = solution.orientation.oe[3][end]
            initial_state.Ω = solution.orientation.oe[4][end]
            initial_state.ω = solution.orientation.oe[5][end]
            initial_state.m = solution.performance.mass[end]
            initial_state.vi = solution.orientation.oe[6][end]

            m.initialcondition.year = Int64(solution.orientation.year[end])
            m.initialcondition.month = Int64(solution.orientation.month[end])
            m.initialcondition.day = Int64(solution.orientation.day[end])
            m.initialcondition.hour = Int64(solution.orientation.hour[end])
            m.initialcondition.minute = Int64(solution.orientation.minute[end])
            m.initialcondition.second = solution.orientation.second[end]

            if (args[:drag_passage] || args[:body_shape] == "Blunted Cone") && continue_campaign
                r = m.planet.Rp_e + args[:EI]*1e3
                initial_state.vi = acos(1 / initial_state.e * (initial_state.a * (1 - initial_state.e^2) / r - 1))
            end
        end

        continue_campaign = asim(ip, m, initial_state, numberofpassage, args)

        r_a = solution.orientation.oe[1][end] * (1 + solution.orientation.oe[2][end])
        r_p = solution.orientation.oe[1][end] * (1 - solution.orientation.oe[2][end])
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
            println("Thermal Limit overcomed totally " * string(config.count_oversome_hr) * " times")
        end

        if r_p - m.planet.Rp_e >= 180*1e3
            FinalState = false
            println("Periapsis too high, final state unreachable! R_a = " * string(r_p*1e-3) * " km")
        end

        println(" ")
    end

    closed_form(args, m)
end