include("Complete_passage.jl")
include("../utils/Ref_system_conf.jl")
include("../utils/Closed_form_solution.jl")
include("../utils/Save_results.jl")
include("../physical_models/Propulsive_maneuvers.jl")

using .SimulationModel
using PythonCall

sys = pyimport("sys")

os = pyimport("os")

function aerobraking(ip, args, gram, gram_atmosphere, filename, temp_name, params)
    cnf, m, solution = params
    initial_state = m.initial_condition
    FinalState = true
    continue_campaign = true
    numberofpassage = 0
    cnf.count_numberofpassage = 0

    # clean_results()
    solution = Solution()
    params = (cnf, m, solution)

    cnf.time_OP = 1
    cnf.time_IP = 1
    
    # Aerobraking Campaign
    while continue_campaign && FinalState
        cnf.index_Mars_Gram_call = 0
        firing_orbit = 0
        numberofpassage += 1

        if args[:print_res] == true
            println("--> Start Passage #" * string(numberofpassage))
        end 

        t_el_ab = @elapsed begin
            # Define maneuver
            if uppercase(args[:thrust_control]) == "AEROBRAKING MANEUVER" && numberofpassage != 1
                r_a = solution.orientation.oe[1][end] * (1 + solution.orientation.oe[2][end])
                r_p = solution.orientation.oe[1][end] * (1 - solution.orientation.oe[2][end])
                args = args[:maneuver_plan](m.planet, r_a, r_p, numberofpassage, args)
            end
            
            if ip.tc == 1
                if numberofpassage == 1
                    args[:delta_v] = 0.0
                end

                if args[:delta_v] != 0.0
                    if round(rad2deg(args[:phi])) == 180
                        append!(cnf.raise_man_orbit, numberofpassage)
                    else
                        append!(cnf.lower_man_orbit, numberofpassage)
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
                initial_state.a = solution.orientation.oe[1][end]
                initial_state.e = solution.orientation.oe[2][end]
                initial_state.i = solution.orientation.oe[3][end]
                initial_state.Ω = solution.orientation.oe[4][end]
                initial_state.ω = solution.orientation.oe[5][end]
                initial_state.m = solution.performance.mass[end]
                initial_state.vi = solution.orientation.oe[6][end]

                m.initial_condition.year = round(solution.orientation.year[end])
                m.initial_condition.month = round(solution.orientation.month[end])
                m.initial_condition.day = round(solution.orientation.day[end])
                m.initial_condition.hour = round(solution.orientation.hour[end])
                m.initial_condition.minute = round(solution.orientation.minute[end])
                m.initial_condition.second = solution.orientation.second[end]
                m.initial_condition.el_time = solution.orientation.time[end]
                println("Initial Date and Time of the Passage: " * string(m.initial_condition.year) * "-" * string(m.initial_condition.month) * "-" * string(m.initial_condition.day) * " " * string(m.initial_condition.hour) * ":" * string(m.initial_condition.minute) * ":" * string(round(m.initial_condition.second, digits=2)))
                if (Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone") && continue_campaign
                    r = m.planet.Rp_e + args[:EI]*1e3
                    initial_state.vi = -acos(1 / initial_state.e * (initial_state.a * (1 - initial_state.e^2) / r - 1))
                end
            end

            params = (cnf, m, Solution()) # Reset solution struct for new passage
            if uppercase(args[:density_model]) == "GRAM"
                continue_campaign = asim(ip, initial_state, numberofpassage, args, params, gram_atmosphere, gram)
            else
                continue_campaign = asim(ip, initial_state, numberofpassage, args, params)
            end
            solution = params[3]
            r_a = solution.orientation.oe[1][end] * (1 + solution.orientation.oe[2][end])
            r_p = solution.orientation.oe[1][end] * (1 - solution.orientation.oe[2][end])
        end

        if Bool(args[:print_res])
            println("Computational time: " * string(t_el_ab) * " seconds")
            println("--> PASSAGE #" * string(numberofpassage) * " COMPLETE")
        end

        if lowercase(args[:type_of_mission]) != "time" && args[:number_of_orbits] == numberofpassage
            continue_campaign = false
        end

        if (r_a <= args[:ra_fin_orbit] || config.cnf.targeting == 1) && args[:keplerian] == false
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

        if args[:closed_form] == 1 && (m.planet.name == "mars" || m.planet.name == "venus" || m.planet.name == "earth" || m.planet.name == "titan")
            closed_form(args, m)
        else
            len_sol = length(solution.orientation.time)
            results(solution, zeros(len_sol), (args[:EI] - 10)*1e3*ones(len_sol), zeros(len_sol), zeros(len_sol))
        end

        if args[:results] == 1
            # Save the current passage results
            save_csv(filename, args, temp_name, params)
            # Clear the config data buffer
            # clean_results()
        end
    end

    # if args[:closed_form] == 1 && (m.planet.name == "mars" || m.planet.name == "venus" || m.planet.name == "earth" || m.planet.name == "titan")
    #     println("Computing Closed-Form Solution...")
    #     closed_form(args, m)
    # else
    #     len_sol = length(config.solution.orientation.time)
    #     results(zeros(len_sol), (args[:EI] - 10)*1e3*ones(len_sol), zeros(len_sol), zeros(len_sol))
    # end

end