function run_entry_duration_sweep(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    selected = _entry_duration_selected_cases(spec, cases)
    interface_counts = _entry_duration_interface_counts(spec)
    println(
        "[entry-duration-sweep] scenarios=$(length(selected)), interface_counts=$(join(interface_counts, ",")), " *
        "profile=$(spec.name)"
    )
    isempty(selected) && return DataFrame()

    rows = NamedTuple[]
    backend = perf_parallel_backend()

    # Always collect serial references first so event-time error can be reported for any backend path.
    println("[entry-duration-sweep] building serial event-time references")
    for (idx, base_case) in enumerate(selected)
        println("  reference scenario $(idx)/$(length(selected)) = $(base_case.name)")
        ref_rows, ref_logs = measure_entry_duration_scenario(
            base_case,
            spec,
            interface_counts;
            outer_route=:none,
            run_role="reference",
            warmup_override=min(spec.warmup, 1),
            repeats_override=1
        )
        append!(rows, ref_rows)
        for log_line in ref_logs
            println(log_line)
        end
    end

    if backend == :auto
        scenario_rows = Vector{Vector{NamedTuple}}(undef, length(selected))
        scenario_logs = Vector{Vector{String}}(undef, length(selected))
        chosen_routes = fill(:none, length(selected))
        process_tasks = Tuple{Int, BenchmarkCase}[]
        thread_indices = Int[]
        serial_indices = Int[]

        for (idx, base_case) in enumerate(selected)
            route = auto_backend_for_case(base_case; spec=spec)
            chosen_routes[idx] = route
            if route == :process
                push!(process_tasks, (idx, base_case))
            elseif route == :threads
                if _case_outer_threads_safe(base_case)
                    push!(thread_indices, idx)
                else
                    push!(serial_indices, idx)
                    chosen_routes[idx] = :none
                end
            else
                push!(serial_indices, idx)
            end
        end

        if !isempty(process_tasks)
            ensure_perf_workers!()
            process_results = pmap(process_tasks) do task
                _, base_case = task
                local_rows, local_logs = perf_worker_measure_entry_duration_scenario(
                    base_case,
                    spec,
                    interface_counts,
                    :process,
                    "measured"
                )
                return (rows=local_rows, logs=local_logs)
            end
            for (k, task) in enumerate(process_tasks)
                idx = task[1]
                scenario_rows[idx] = process_results[k].rows
                scenario_logs[idx] = process_results[k].logs
            end
        end

        if !isempty(thread_indices)
            for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
                withenv(env_pairs...) do
                    Threads.@threads for j in eachindex(payload)
                        idx, _ = payload[j]
                        local_rows, local_logs = measure_entry_duration_scenario(
                            selected[idx],
                            spec,
                            interface_counts;
                            outer_route=:threads,
                            apply_env=false,
                            run_role="measured"
                        )
                        scenario_rows[idx] = local_rows
                        scenario_logs[idx] = local_logs
                    end
                end
            end
        end

        for idx in serial_indices
            local_rows, local_logs = measure_entry_duration_scenario(
                selected[idx],
                spec,
                interface_counts;
                outer_route=:none,
                run_role="measured"
            )
            scenario_rows[idx] = local_rows
            scenario_logs[idx] = local_logs
        end

        for (idx, base_case) in enumerate(selected)
            println("  measured scenario $(idx)/$(length(selected)) = $(base_case.name)")
            _record_outer_route_feedback!(base_case, scenario_rows[idx]; route=chosen_routes[idx])
            append!(rows, scenario_rows[idx])
            for log_line in scenario_logs[idx]
                println(log_line)
            end
        end
    elseif backend == :process
        ensure_perf_workers!()
        tasks = collect(enumerate(selected))
        scenario_results = pmap(tasks) do task
            _, base_case = task
            local_rows, local_logs = perf_worker_measure_entry_duration_scenario(
                base_case,
                spec,
                interface_counts,
                :process,
                "measured"
            )
            return (rows=local_rows, logs=local_logs)
        end
        for (idx, base_case) in enumerate(selected)
            println("  measured scenario $(idx)/$(length(selected)) = $(base_case.name)")
            _record_outer_route_feedback!(base_case, scenario_results[idx].rows; route=:process)
            append!(rows, scenario_results[idx].rows)
            for log_line in scenario_results[idx].logs
                println(log_line)
            end
        end
    elseif backend == :threads
        scenario_rows = Vector{Vector{NamedTuple}}(undef, length(selected))
        scenario_logs = Vector{Vector{String}}(undef, length(selected))
        thread_indices, serial_indices = _split_threadsafe_indices(selected, collect(eachindex(selected)))

        for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
            withenv(env_pairs...) do
                Threads.@threads for j in eachindex(payload)
                    idx, _ = payload[j]
                    local_rows, local_logs = measure_entry_duration_scenario(
                        selected[idx],
                        spec,
                        interface_counts;
                        outer_route=:threads,
                        apply_env=false,
                        run_role="measured"
                    )
                    scenario_rows[idx] = local_rows
                    scenario_logs[idx] = local_logs
                end
            end
        end
        for idx in serial_indices
            local_rows, local_logs = measure_entry_duration_scenario(
                selected[idx],
                spec,
                interface_counts;
                outer_route=:none,
                run_role="measured"
            )
            scenario_rows[idx] = local_rows
            scenario_logs[idx] = local_logs
        end
        for (idx, base_case) in enumerate(selected)
            println("  measured scenario $(idx)/$(length(selected)) = $(base_case.name)")
            route = idx in serial_indices ? :none : :threads
            _record_outer_route_feedback!(base_case, scenario_rows[idx]; route=route)
            append!(rows, scenario_rows[idx])
            for log_line in scenario_logs[idx]
                println(log_line)
            end
        end
    else
        for (idx, base_case) in enumerate(selected)
            println("  measured scenario $(idx)/$(length(selected)) = $(base_case.name)")
            local_rows, local_logs = measure_entry_duration_scenario(
                base_case,
                spec,
                interface_counts;
                outer_route=:none,
                run_role="measured"
            )
            _record_outer_route_feedback!(base_case, local_rows; route=:none)
            append!(rows, local_rows)
            for log_line in local_logs
                println(log_line)
            end
        end
    end

    entry_raw_df = DataFrame(rows)
    if nrow(entry_raw_df) == 0
        return entry_raw_df
    end

    reference_rows = entry_raw_df[
        (entry_raw_df.entry_run_role .== "reference") .&
        (entry_raw_df.solve_success .== true) .&
        .!ismissing.(entry_raw_df.terminal_time_s), :
    ]
    reference_map = Dict{Tuple{String, Int}, Float64}()
    if nrow(reference_rows) > 0
        for grp in groupby(reference_rows, [:scenario, :entry_atmospheric_interface_count])
            key = (String(grp.scenario[1]), Int(grp.entry_atmospheric_interface_count[1]))
            reference_map[key] = mean(Float64.(grp.terminal_time_s))
        end
    end

    entry_raw_df[!, :entry_reference_terminal_time_s] = [
        get(reference_map, (String(sc), Int(ic)), missing)
        for (sc, ic) in zip(entry_raw_df.scenario, entry_raw_df.entry_atmospheric_interface_count)
    ]
    entry_raw_df[!, :entry_event_time_abs_error_s] = [
        (t isa Missing || ref isa Missing) ? missing : abs(Float64(t) - Float64(ref))
        for (t, ref) in zip(entry_raw_df.terminal_time_s, entry_raw_df.entry_reference_terminal_time_s)
    ]

    return entry_raw_df
end

