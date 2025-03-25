include("Complete_passage.jl")
include("../utils/Ref_system_conf.jl")
# include("../utils/Closed_form_solution.jl")
include("../utils/Odyssey_maneuver_plan.jl")
include("../utils/VEx_maneuver_plan.jl")
include("../utils/Magellan_maneuver_plan.jl")
include("../utils/Save_results.jl")
include("../physical_models/Propulsive_maneuvers.jl")

using PythonCall

 # import .config

sys = pyimport("sys")

os = pyimport("os")

# sys.path.append(os.path.join(os.path.dirname(os.path.abspath(@__FILE__)), "GRAMpy"))

function aerobraking(ip, m, args)

    sys.path.append(args[:directory_Gram])

    gram = pyimport("gram")

    initial_state = m.initial_condition
    FinalState = true
    continue_campaign = true
    numberofpassage = 0
    config.cnf.count_numberofpassage = 0

    clean_results()

    config.cnf.time_OP = 1
    config.cnf.time_IP = 1

    if args[:density_model] == "Gram" || args[:density_model] == "GRAM"
        inputParameters = Dict("earth" => gram.EarthInputParameters(),
                               "mars" => gram.MarsInputParameters(),
                               "venus" => gram.VenusInputParameters())
        
        namelistReaders = Dict("earth" => gram.EarthNamelistReader(),
                               "mars" => gram.MarsNamelistReader(),
                               "venus" => gram.VenusNamelistReader())
            
        atmospheres = Dict("earth" => gram.EarthAtmosphere(),
                           "mars" => gram.MarsAtmosphere(),
                           "venus" => gram.VenusAtmosphere())

        planet_name = m.planet.name
        input_parameters = inputParameters[planet_name]

        # Mars has some weird specific parameters, so this line is just to check to make sure the it doesn't do it for the other planets
        if planet_name == "mars"
            # input_parameters.dataPath = os.path.join(os.path.dirname(os.path.abspath(@__FILE__)),"..", "GRAM_Data", "Mars", "data", "")
            input_parameters.dataPath = args[:directory_Gram_data] * "/Mars/data/"
            if !Bool(os.path.exists(input_parameters.dataPath))
                throw(ArgumentError("GRAM data path not found: " * input_parameters.dataPath))
            end
        end

        reader = namelistReaders[planet_name]
        reader.tryGetSpicePath(input_parameters)

        gram_atmosphere = atmospheres[planet_name]
        gram_atmosphere.setInputParameters(input_parameters)

        gram_atmosphere.setPerturbationScales(1.5)
        gram_atmosphere.setMinRelativeStepSize(0.5)
        gram_atmosphere.setSeed(1001)
        if planet_name == "mars"
            gram_atmosphere.setMOLAHeights(false)
        end

        ttime = gram.GramTime()
        ttime.setStartTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs], gram.UTC, gram.PET)
        gram_atmosphere.setStartTime(ttime)
    end
    
    # Aerobraking Campaign
    while continue_campaign && FinalState
        config.cnf.index_Mars_Gram_call = 0
        firing_orbit = 0
        numberofpassage = 1 + numberofpassage

        if args[:print_res] == true
            println("--> Start Passage #" * string(numberofpassage))
        end

        t_el_ab = @elapsed begin

            if args[:Odyssey_sim] == true
                args = Odyssey_firing_plan(numberofpassage, args)
            end

            if args[:vex_sim] == true
                args = Venus_Express_firing_plan(numberofpassage, args)
            end

            if args[:magellan_sim] == true
                args = Magellan_firing_plan(numberofpassage, args)
            end

            if ip.tc == 1
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

                if (Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone") && continue_campaign
                    r = m.planet.Rp_e + args[:EI]*1e3
                    initial_state.vi = -acos(1 / initial_state.e * (initial_state.a * (1 - initial_state.e^2) / r - 1))
                end
            end

            if args[:density_model] == "Gram"
                continue_campaign = asim(ip, m, initial_state, numberofpassage, args, gram_atmosphere)
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

        if args[:number_of_orbits] == numberofpassage
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

        println(" ")
    end

    # closed_form(args, m)
    save_fitting_data(args, m)
end