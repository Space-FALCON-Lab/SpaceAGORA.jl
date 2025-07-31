include("Complete_passage.jl")
include("../utils/Ref_system_conf.jl")
include("../utils/Closed_form_solution.jl")
include("../utils/Save_results.jl")
include("../physical_models/Propulsive_maneuvers.jl")

using PythonCall

sys = pyimport("sys")

os = pyimport("os")

function aerobraking(ip, m, args, gram, gram_atmosphere, filename, temp_name)

    initial_state = m.initial_condition
    FinalState = true
    continue_campaign = true
    numberofpassage = 0
    config.cnf.count_numberofpassage = 0

    clean_results()

    config.cnf.time_OP = 1
    config.cnf.time_IP = 1
    
    # Aerobraking Campaign
    while continue_campaign && FinalState
        config.cnf.index_Mars_Gram_call = 0
        firing_orbit = 0
        numberofpassage += 1

        if args[:print_res] == true
            println("--> Start Passage #" * string(numberofpassage))
        end 

        t_el_ab = @elapsed begin
            # Define maneuver
            if uppercase(args[:thrust_control]) == "AEROBRAKING MANEUVER" && numberofpassage != 1
                r_a = config.solution.orientation.oe[1][end] * (1 + config.solution.orientation.oe[2][end])
                r_p = config.solution.orientation.oe[1][end] * (1 - config.solution.orientation.oe[2][end])
                args = args[:maneuver_plan](m.planet, r_a, r_p, numberofpassage, args)
            end
            
            if ip.tc == 1
                if numberofpassage == 1
                    args[:delta_v] = 0.0
                end

                if args[:delta_v] != 0.0
                    if round(rad2deg(args[:phi])) == 180
                        append!(config.cnf.raise_man_orbit, numberofpassage)
                    else
                        append!(config.cnf.lower_man_orbit, numberofpassage)
                    end

                    if ip.tc == 1
                        initial_state = propulsion_ic_calcs(m, args, initial_state)
                    end
                end
            elseif ip.tc == 2
                if round(rad2deg(args[:phi])) == 180
                    println("DECELERATE DRAG FIRING!!")
                elseif round(rad2deg(args[:phi])) == 0
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

                m.initial_condition.year = round(config.solution.orientation.year[end])
                m.initial_condition.month = round(config.solution.orientation.month[end])
                m.initial_condition.day = round(config.solution.orientation.day[end])
                m.initial_condition.hour = round(config.solution.orientation.hour[end])
                m.initial_condition.minute = round(config.solution.orientation.minute[end])
                m.initial_condition.second = round(config.solution.orientation.second[end])
                m.initial_condition.el_time = config.solution.orientation.time[end]
                if (Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone") && continue_campaign
                    r = m.planet.Rp_e + args[:EI]*1e3
                    initial_state.vi = -acos(1 / initial_state.e * (initial_state.a * (1 - initial_state.e^2) / r - 1))
                end
            end

            clean_results()
            if uppercase(args[:density_model]) == "GRAM"
                continue_campaign = asim(ip, m, initial_state, numberofpassage, args, gram_atmosphere, gram)
            else
                continue_campaign = asim(ip, m, initial_state, numberofpassage, args)
            end

            r_a = config.solution.orientation.oe[1][end] * (1 + config.solution.orientation.oe[2][end])
            r_p = config.solution.orientation.oe[1][end] * (1 - config.solution.orientation.oe[2][end])
        end

        if Bool(args[:print_res])
            println("Computational time: " * string(t_el_ab) * " seconds")
            println("--> PASSAGE #" * string(numberofpassage) * " COMPLETE")
        end

        if lowercase(args[:type_of_mission]) != "time" && args[:number_of_orbits] == numberofpassage
            continue_campaign = false
        end

        if r_a <= args[:final_apoapsis] && args[:keplerian] == false
            FinalState = false
            println("Reached FinalState! R_a = " * string(r_a*1e-3) * " km")
            println("Thermal Limit overcomed totally " * string(config.cnf.count_overcome_hr) * " times")
        end

        if r_p - m.planet.Rp_e >= args[:EI]*1e3 && args[:keplerian] == false
            FinalState = false
            println("Periapsis too high, final state unreachable! R_a = " * string(r_p*1e-3) * " km")
        end

        if Bool(args[:print_res])
            println(" ")
        end

        if args[:results] == 1
            # Save the current passage results
            save_csv(filename, args, temp_name)
            # Clear the config data buffer
            # clean_results()
        end
    end

    if args[:closed_form] == 1 && (m.planet.name == "mars" || m.planet.name == "venus" || m.planet.name == "earth" || m.planet.name == "titan")
        closed_form(args, m)
    else
        len_sol = length(config.solution.orientation.time)
        results(zeros(len_sol), (args[:EI] - 10)*1e3*ones(len_sol), zeros(len_sol), zeros(len_sol))
    end

end