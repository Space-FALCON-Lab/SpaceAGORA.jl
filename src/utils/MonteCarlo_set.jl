import .config

function MonteCarlo_setting(args)
    MC = Dict()
    MC[:N_passages], MC[:Duration], MC[:Median_heat], MC[:Periapsis_min], MC[:Periapsis_max], count = [], [], [], [], [], 0

    if args[:montecarlo] == false
        args[:montecarlo_size] = 2
        args[:intial_montecarlo_number] = 1
    end

    return MC, count, args
end

function MonteCarlo_setting_passage(mc_index, args)
    config.cnf.counter_random = 0

    heat_passage = []

    if args[:print_res]
        println("--> MC number: " * string(mc_index))
    end

    args[:simulation_filename] = args[:simulation_filename] * "_nMC=" * string(mc_index)

    # Initialization
    config.cnf.altitude_periapsis, config.cnf.max_heatrate = [], []
    config.cnf.counter_random = mc_index
    config.cnf.index_MonteCarlo = mc_index

    return args
end

function MonteCarlo_append(MC, args, count)
    # append results to MC
    config.cnf.index_MonteCarlo += 1

    # save results
    append!(MC[:N_passages], config.solution.orientation.number_of_passage[end])
    append!(MC[:Duration], config.solution.orientation.time[end])
    append!(MC[:Median_heat], median(config.cnf.max_heatrate))
    append!(MC[:Periapsis_min], minimum(config.cnf.altitude_periapsis))
    append!(MC[:Periapsis_max], maximum(config.cnf.altitude_periapsis))

    heat_rate_max = maximum(config.solution.performance.heat_rate) 
    if heat_rate_max > args[:max_heat_rate]
        count += 1
    end

    if args[:print_res]
        println("--> Count = " * count)
    end
end

function MonteCarlo_save(args, state, MC)
    if args[:save_results]
        folder_name = args[:simulation_filename][1:indexin("_nMC", args[:simulation_filename])]

        name = args[:directory_results] * folder_name * "/MC_results_control=" * string(args[:control_mode]) * "_ra=" * string(Int64(state[:Apoapsis]/1e3)) * "_rp=" * string(state[:Periapsis]) * "_hl=" * string(args[:max_heat_rate])
        filename = name * ".csv"

        # NEED TO FINISH

    end
end