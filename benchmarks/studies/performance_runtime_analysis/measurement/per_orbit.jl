function measure_entry_duration_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    interface_counts::Vector{Int};
    outer_route::Symbol=:none,
    apply_env::Bool=true,
    run_role::String="measured",
    warmup_override::Union{Nothing, Int}=nothing,
    repeats_override::Union{Nothing, Int}=nothing
)
    rows = NamedTuple[]
    logs = String[]
    stream_logs = _perf_stream_orbit_logs()
    warmup_count = isnothing(warmup_override) ? spec.warmup : Int(warmup_override)
    repeat_count = isnothing(repeats_override) ? spec.repeats : Int(repeats_override)
    warmup_count >= 0 || throw(ArgumentError("warmup_override must be >= 0, got $warmup_count"))
    repeat_count > 0 || throw(ArgumentError("repeats_override must be > 0, got $repeat_count"))

    for interface_count in interface_counts
        case = _entry_sweep_case(base_case, interface_count)
        plan = parallel_priority_plan(case, outer_route)

        run_case = () -> begin
            if stream_logs
                println(
                    "    entry_interface_count=$(interface_count): warmup x$(warmup_count), " *
                    "repeats x$(repeat_count), role=$(run_role)"
                )
                flush(stdout)
            end
            run_warmup(case, warmup_count, spec.name)
            for rep in 1:repeat_count
                last_row = nothing
                for attempt in 1:spec.max_attempts
                    row = measure_case(case, spec.name, rep; attempt=attempt, plan=plan)
                    row_entry = merge(
                        row,
                        (
                            entry_atmospheric_interface_count=interface_count,
                            entry_run_role=run_role,
                            entry_passage_duration_s=row.terminal_time_s,
                            entry_wall_time_per_passage_s=_safe_div(row.total_time_s, Float64(interface_count))
                        )
                    )
                    last_row = row_entry
                    if row_entry.solve_success
                        push!(rows, row_entry)
                        line = "    entry_interface_count=$(interface_count) repeat $(rep)/$(repeat_count): total=$(round(row_entry.total_time_s; digits=3)) s"
                        push!(logs, line)
                        if stream_logs
                            println(line)
                            flush(stdout)
                        end
                        break
                    end
                end
                if !(last_row === nothing) && !last_row.solve_success
                    push!(rows, last_row)
                    line = "    entry_interface_count=$(interface_count) repeat $(rep)/$(repeat_count): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)"
                    push!(logs, line)
                    if stream_logs
                        println(line)
                        flush(stdout)
                    end
                end
            end
        end

        if apply_env
            env_pairs = parallel_priority_env_pairs(plan)
            withenv(env_pairs...) do
                run_case()
            end
        else
            run_case()
        end
    end
    return rows, logs
end

function measure_per_orbit_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    period_s::Float64,
    orbit_counts::Vector{Int};
    outer_route::Symbol=:none,
    apply_env::Bool=true
)
    rows = NamedTuple[]
    logs = String[]
    stream_logs = _perf_stream_orbit_logs()
    for orbit_count in orbit_counts
        mission_time = orbit_count * period_s
        args_template = deepcopy(base_case.args_template)
        args_template = SimulationConfiguration(
            file_paths=args_template.file_paths,
            simulation_settings=args_template.simulation_settings,
            mission_configuration=MissionConfiguration(
                mission_type=args_template.mission_configuration.mission_type,
                keplerian=args_template.mission_configuration.keplerian,
                number_of_orbits=args_template.mission_configuration.number_of_orbits,
                mission_time=mission_time,
                orientation_sim=args_template.mission_configuration.orientation_sim,
                num_steps_to_save=args_template.mission_configuration.num_steps_to_save
            ),
            environment_model=args_template.environment_model,
            dynamics_model=args_template.dynamics_model,
            guidance_model=args_template.guidance_model,
            navigation_model=args_template.navigation_model,
            control_model=args_template.control_model,
            initial_time=args_template.initial_time,
            integration_tolerances=args_template.integration_tolerances
        )

        case = BenchmarkCase(
            name=base_case.name,
            category=base_case.category,
            description=base_case.description,
            args_template=args_template,
            run_in_quick=base_case.run_in_quick,
            solver_mode_override=base_case.solver_mode_override,
            split_imex_solver_override=base_case.split_imex_solver_override,
            entry_target_count_override=base_case.entry_target_count_override
        )

        plan = parallel_priority_plan(case, outer_route)
        run_case = () -> begin
            if stream_logs
                println("    mission_time_multiplier=x$(orbit_count): warmup x$(spec.warmup), repeats x$(spec.repeats)")
                flush(stdout)
            end
            run_warmup(case, spec.warmup, spec.name)
            for rep in 1:spec.repeats
                last_row = nothing
                for attempt in 1:spec.max_attempts
                    row = measure_case(case, spec.name, rep; attempt=attempt, plan=plan)
                    row_orbit = merge(
                        row,
                        (
                            orbit_count=orbit_count,
                            mission_time_multiplier=orbit_count,
                            orbital_period_s=period_s
                        )
                    )
                    last_row = row_orbit
                    if row_orbit.solve_success
                        push!(rows, row_orbit)
                        line = "    mission_time_multiplier=x$(orbit_count) repeat $(rep)/$(spec.repeats): total=$(round(row_orbit.total_time_s; digits=3)) s"
                        push!(logs, line)
                        if stream_logs
                            println(line)
                            flush(stdout)
                        end
                        break
                    end
                end
                if !(last_row === nothing) && !last_row.solve_success
                    push!(rows, last_row)
                    line = "    mission_time_multiplier=x$(orbit_count) repeat $(rep)/$(spec.repeats): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)"
                    push!(logs, line)
                    if stream_logs
                        println(line)
                        flush(stdout)
                    end
                end
            end
        end
        if apply_env
            env_pairs = parallel_priority_env_pairs(plan)
            withenv(env_pairs...) do
                run_case()
            end
        else
            run_case()
        end
    end
    return rows, logs
end

function run_montecarlo_per_orbit!(
    rows::Vector{NamedTuple},
    spec::ProfileSpec,
    planet::Earth,
    period_s::Float64,
    orbit_counts::Vector{Int}
)
    seeds = collect(1:spec.montecarlo_samples)
    scenarios = _active_montecarlo_scenarios()
    if _perf_smoke_mode()
        scenarios = scenarios[1:min(end, 1)]
    end
    scenario_names = join([s.name for s in scenarios], ", ")
    println("  montecarlo scenarios (mission-time sweep, seeds=$(length(seeds))): $(scenario_names)")
    backend = perf_parallel_backend()
    mars = perf_worker_mars()
    workers_ready = false

    for scenario in scenarios
        variant = scenario.variant
        println("  scenario $(scenario.name)")
        for orbit_count in orbit_counts
            mission_time = orbit_count * period_s
            orbit_case = make_montecarlo_case(first(seeds), mission_time, variant, planet; mars=mars)
            mc_backend = backend == :auto ? auto_backend_for_case(orbit_case; spec=spec) : backend
            if mc_backend == :process && !workers_ready
                ensure_perf_workers!()
                workers_ready = true
            end
            println("    mission_time_multiplier=x$(orbit_count)")
            orbit_rows = Vector{NamedTuple}(undef, length(seeds))
            orbit_msgs = Vector{String}(undef, length(seeds))

            if mc_backend == :process
                seed_results = pmap(seed -> perf_worker_measure_montecarlo_seed(spec, mission_time, seed, variant, :process), seeds)
                for i in eachindex(seeds)
                    seed = seeds[i]
                    row, err = seed_results[i]
                    row_orbit = merge(row, (orbit_count=orbit_count, mission_time_multiplier=orbit_count, orbital_period_s=period_s))
                    orbit_rows[i] = row_orbit
                    if row_orbit.solve_success
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                    else
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                    end
                end
            elseif mc_backend == :threads
                threaded_plan = parallel_priority_plan(orbit_case, :threads)
                threaded_env = parallel_priority_env_pairs(threaded_plan)
                withenv(threaded_env...) do
                    Threads.@threads for i in eachindex(seeds)
                        seed = seeds[i]
                        row, err = measure_montecarlo_seed(
                            spec,
                            planet,
                            mission_time,
                            seed;
                            variant=variant,
                            mars=mars,
                            outer_route=:threads,
                            plan=threaded_plan,
                            apply_env=false
                        )
                        row_orbit = merge(row, (orbit_count=orbit_count, mission_time_multiplier=orbit_count, orbital_period_s=period_s))
                        orbit_rows[i] = row_orbit
                        if row_orbit.solve_success
                            orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                        else
                            orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                        end
                    end
                end
            else
                for i in eachindex(seeds)
                    seed = seeds[i]
                    row, err = measure_montecarlo_seed(
                        spec,
                        planet,
                        mission_time,
                        seed;
                        variant=variant,
                        mars=mars,
                        outer_route=:none
                    )
                    row_orbit = merge(row, (orbit_count=orbit_count, mission_time_multiplier=orbit_count, orbital_period_s=period_s))
                    orbit_rows[i] = row_orbit
                    if row_orbit.solve_success
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                    else
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                    end
                end
            end

            for i in eachindex(seeds)
                push!(rows, orbit_rows[i])
                println(orbit_msgs[i])
            end
            _record_outer_route_feedback!(orbit_case, orbit_rows; route=mc_backend)
        end
    end
    return nothing
end
function run_per_orbit_for_scenarios(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    baseline_sc = make_spacecraft(planet; id=1, with_panel=false)
    period_s = orbital_period_seconds(baseline_sc, planet)
    orbit_counts = _perf_smoke_mode() ? [1] : (spec.name == "full" ? collect(1:5) : collect(1:3))
    include_control_stress = _include_control_stress_per_orbit()
    selected_base = include_control_stress ? selected_cases(spec, cases) : [c for c in selected_cases(spec, cases) if c.category != "control_stress"]
    # Entry scenarios use entry-interface counting semantics, not baseline-period multipliers.
    selected = [c for c in selected_base if c.category != "entry"]
    if _perf_smoke_mode()
        selected = selected[1:min(end, 1)]
    end
    excluded_entry = length(selected_base) - length(selected)

    println(
        "[mission-time-sweep] scenarios=$(length(selected)), baseline period=$(round(period_s; digits=3)) s, " *
        "multipliers=x$(first(orbit_counts)):x$(last(orbit_counts)), include_control_stress=$(include_control_stress), " *
        "entry_excluded=$(excluded_entry)"
    )
    rows = NamedTuple[]
    backend = perf_parallel_backend()

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
                idx, base_case = task
                local_rows, local_logs = perf_worker_measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts, :process)
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
                        local_rows, local_logs = measure_per_orbit_scenario(
                            selected[idx],
                            spec,
                            period_s,
                            orbit_counts;
                            outer_route=:threads,
                            apply_env=false
                        )
                        scenario_rows[idx] = local_rows
                        scenario_logs[idx] = local_logs
                    end
                end
            end
        end

        for idx in serial_indices
            local_rows, local_logs = measure_per_orbit_scenario(selected[idx], spec, period_s, orbit_counts; outer_route=:none)
            scenario_rows[idx] = local_rows
            scenario_logs[idx] = local_logs
        end

        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
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
            idx, base_case = task
            local_rows, local_logs = perf_worker_measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts, :process)
            return (rows=local_rows, logs=local_logs)
        end
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
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
                    local_rows, local_logs = measure_per_orbit_scenario(
                        selected[idx],
                        spec,
                        period_s,
                        orbit_counts;
                        outer_route=:threads,
                        apply_env=false
                    )
                    scenario_rows[idx] = local_rows
                    scenario_logs[idx] = local_logs
                end
            end
        end
        for idx in serial_indices
            local_rows, local_logs = measure_per_orbit_scenario(
                selected[idx],
                spec,
                period_s,
                orbit_counts;
                outer_route=:none
            )
            scenario_rows[idx] = local_rows
            scenario_logs[idx] = local_logs
        end
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            route = idx in serial_indices ? :none : :threads
            _record_outer_route_feedback!(base_case, scenario_rows[idx]; route=route)
            append!(rows, scenario_rows[idx])
            for log_line in scenario_logs[idx]
                println(log_line)
            end
        end
    else
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            local_rows, local_logs = measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts; outer_route=:none)
            _record_outer_route_feedback!(base_case, local_rows; route=:none)
            append!(rows, local_rows)
            for log_line in local_logs
                println(log_line)
            end
        end
    end

    run_montecarlo_per_orbit!(rows, spec, planet, period_s, orbit_counts)
    return DataFrame(rows)
end
