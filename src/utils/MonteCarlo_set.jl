
if !isdefined(@__MODULE__, :_legacy_get_mc_runtime_state)
    @inline function _legacy_get_mc_runtime_state(args=nothing; cnf=nothing, solution=nothing)
        module_runtime = isdefined(@__MODULE__, :config) ? getfield(@__MODULE__, :config) : nothing

        cnf_state = if cnf !== nothing
            cnf
        elseif args isa AbstractDict && haskey(args, :cnf)
            args[:cnf]
        elseif module_runtime !== nothing && hasproperty(module_runtime, :cnf)
            getproperty(module_runtime, :cnf)
        else
            nothing
        end

        solution_state = if solution !== nothing
            solution
        elseif args isa AbstractDict && haskey(args, :solution)
            args[:solution]
        elseif module_runtime !== nothing && hasproperty(module_runtime, :solution)
            getproperty(module_runtime, :solution)
        else
            nothing
        end

        if cnf_state === nothing || solution_state === nothing
            throw(ArgumentError("MonteCarlo state missing; provide `cnf` and `solution` via args."))
        end

        return (cnf=cnf_state, solution=solution_state)
    end
end

function MonteCarlo_setting(args)
    MC = Dict()
    MC[:N_passages], MC[:Duration], MC[:Median_heat], MC[:Periapsis_min], MC[:Periapsis_max], count = [], [], [], [], [], 0

    if args[:montecarlo] == false
        args[:montecarlo_size] = 1
        args[:intial_montecarlo_number] = 1
    end

    return MC, count, args
end

function MonteCarlo_setting_passage(mc_index, args)
    runtime = _legacy_get_mc_runtime_state(args)
    cnf_state = runtime.cnf
    cnf_state.counter_random = 0

    heat_passage = []

    if Bool(args[:print_res])
        println("--> MC number: " * string(mc_index))
    end

    args[:simulation_filename] = args[:simulation_filename] * "_nMC=" * string(mc_index)

    # Initialization
    cnf_state.altitude_periapsis, cnf_state.max_heatrate = [], []
    cnf_state.counter_random = mc_index
    cnf_state.index_MonteCarlo = mc_index

    return args
end

function MonteCarlo_append(MC, args, count)
    runtime = _legacy_get_mc_runtime_state(args)
    cnf_state = runtime.cnf
    solution_state = runtime.solution

    # append results to MC
    cnf_state.index_MonteCarlo += 1

    # save results
    append!(MC[:N_passages], solution_state.orientation.number_of_passage[end])
    append!(MC[:Duration], solution_state.orientation.time[end])
    append!(MC[:Median_heat], median(cnf_state.max_heatrate))
    if args[:type_of_mission] != "Entry" && !isempty(cnf_state.altitude_periapsis)
        append!(MC[:Periapsis_min], minimum(cnf_state.altitude_periapsis))
        append!(MC[:Periapsis_max], maximum(cnf_state.altitude_periapsis))
    end
    heat_rate_max = maximum(maximum.(solution_state.performance.heat_rate))
    if heat_rate_max > args[:max_heat_rate]
        count += 1
    end

    if Bool(args[:print_res])
        println("--> Count = " * string(count))
    end

    return count
end

function MonteCarlo_save(args, state, MC)
    if Bool(args[:save_results])
        folder_name = args[:simulation_filename][1:indexin("_nMC", args[:simulation_filename])]

        name = args[:directory_results] * folder_name * "/MC_results_control=" * string(args[:control_mode]) * "_ra=" * string(Int64(state[:Apoapsis]/1e3)) * "_rp=" * string(state[:Periapsis]) * "_hl=" * string(args[:max_heat_rate])
        filename = name * ".csv"

        if !isdir(args[:directory_results])
            mkpath(args[:directory_results])
        end

        touch(filename)

        writer = open(filename, "w")

        data_push = DataFrame(MC_size = range(args[:montecarlo_size]),
                              N_passages = MC[:N_passages],
                              Duration = MC[:Duration],
                              Median_heat = MC[:Median_heat],
                              Periapsis_min = MC[:Periapsis_min],
                              Periapsis_max = MC[:Periapsis_max])
        
        CSV.write(filename, data_push)
    end
end
