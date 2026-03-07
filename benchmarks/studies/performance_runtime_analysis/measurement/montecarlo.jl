function measure_montecarlo_seed(
    spec::ProfileSpec,
    planet::Earth,
    mission_time_s::Float64,
    seed::Int;
    variant::Symbol=:high_accuracy,
    mars::Union{Nothing, Mars}=nothing,
    outer_route::Symbol=:none,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing,
    apply_env::Bool=true
)
    case = make_montecarlo_case(seed, mission_time_s, variant, planet; mars=mars)
    resolved_plan = isnothing(plan) ? parallel_priority_plan(case, outer_route) : plan
    run_seed = () -> begin
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, 1; seed=seed, attempt=attempt, plan=resolved_plan)
            last_row = row
            if row.solve_success
                return row, nothing
            end
        end
        return last_row, last_row === nothing ? "failed without attempt data" : "failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)"
    end
    if apply_env
        env_pairs = parallel_priority_env_pairs(resolved_plan)
        return withenv(env_pairs...) do
            run_seed()
        end
    end
    return run_seed()
end

function perf_worker_montecarlo_warmup(
    spec::ProfileSpec,
    mission_time_s::Float64,
    seed::Int,
    variant::Symbol,
    outer_route::Symbol=:process
)
    planet = perf_worker_planet()
    mars = perf_worker_mars()
    warmup_case = make_montecarlo_case(seed, mission_time_s, variant, planet; mars=mars)
    plan = parallel_priority_plan(warmup_case, outer_route)
    env_pairs = parallel_priority_env_pairs(plan)
    withenv(env_pairs...) do
        run_warmup(warmup_case, spec.warmup, spec.name)
    end
    return nothing
end

function perf_worker_measure_montecarlo_seed(
    spec::ProfileSpec,
    mission_time_s::Float64,
    seed::Int,
    variant::Symbol,
    outer_route::Symbol=:process
)
    planet = perf_worker_planet()
    mars = perf_worker_mars()
    return measure_montecarlo_seed(
        spec,
        planet,
        mission_time_s,
        seed;
        variant=variant,
        mars=mars,
        outer_route=outer_route
    )
end

function perf_worker_run_case_batch(
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int,
    outer_route::Symbol=:process
)
    local_rows = NamedTuple[]
    run_case_batch!(local_rows, case, spec, idx, total; outer_route=outer_route)
    return local_rows
end

function perf_worker_measure_per_orbit_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    period_s::Float64,
    orbit_counts::Vector{Int},
    outer_route::Symbol=:process
)
    return measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts; outer_route=outer_route)
end

function perf_worker_measure_entry_duration_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    interface_counts::Vector{Int},
    outer_route::Symbol=:process,
    run_role::String="measured";
    warmup_override::Union{Nothing, Int}=nothing,
    repeats_override::Union{Nothing, Int}=nothing
)
    return measure_entry_duration_scenario(
        base_case,
        spec,
        interface_counts;
        outer_route=outer_route,
        run_role=run_role,
        warmup_override=warmup_override,
        repeats_override=repeats_override
    )
end

function run_montecarlo_batch!(rows::Vector{NamedTuple}, spec::ProfileSpec, planet::Earth)
    seeds = collect(1001:(1000 + spec.montecarlo_samples))
    scenarios = _active_montecarlo_scenarios()
    if _perf_smoke_mode()
        scenarios = scenarios[1:min(end, 1)]
    end
    scenario_names = join([s.name for s in scenarios], ", ")
    println("[montecarlo] warmup x$(spec.warmup), seeds=$(length(seeds)), scenarios=$(scenario_names)")

    backend = perf_parallel_backend()
    mars = perf_worker_mars()

    for scenario in scenarios
        variant = scenario.variant
        mission_time_s = _montecarlo_batch_mission_time_s(spec, variant)
        warmup_case = make_montecarlo_case(first(seeds), mission_time_s, variant, planet; mars=mars)
        println("  scenario $(scenario.name) (mission_time=$(round(mission_time_s; digits=1)) s)")

        mc_backend = backend == :auto ? auto_backend_for_case(warmup_case; spec=spec) : backend
        if mc_backend == :threads && !_case_outer_threads_safe(warmup_case)
            mc_backend = :none
        end
        if mc_backend == :process
            ensure_perf_workers!()
            warmup_seed = first(seeds)
            for w in workers()
                remotecall_wait(perf_worker_montecarlo_warmup, w, spec, mission_time_s, warmup_seed, variant, :process)
            end
        else
            plan = parallel_priority_plan(warmup_case, mc_backend)
            env_pairs = parallel_priority_env_pairs(plan)
            withenv(env_pairs...) do
                run_warmup(warmup_case, spec.warmup, spec.name)
            end
        end

        seed_rows = Vector{NamedTuple}(undef, length(seeds))
        seed_msgs = Vector{String}(undef, length(seeds))

        if mc_backend == :process
            seed_results = pmap(seed -> perf_worker_measure_montecarlo_seed(spec, mission_time_s, seed, variant, :process), seeds)
            for i in eachindex(seeds)
                seed = seeds[i]
                row, err = seed_results[i]
                seed_rows[i] = row
                if row.solve_success
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
                else
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): $(err)"
                end
            end
        elseif mc_backend == :threads
            threaded_plan = parallel_priority_plan(warmup_case, :threads)
            threaded_env = parallel_priority_env_pairs(threaded_plan)
            withenv(threaded_env...) do
                Threads.@threads for i in eachindex(seeds)
                    seed = seeds[i]
                    row, err = measure_montecarlo_seed(
                        spec,
                        planet,
                        mission_time_s,
                        seed;
                        variant=variant,
                        mars=mars,
                        outer_route=:threads,
                        plan=threaded_plan,
                        apply_env=false
                    )
                    seed_rows[i] = row
                    if row.solve_success
                        seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
                    else
                        seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): $(err)"
                    end
                end
            end
        else
            for i in eachindex(seeds)
                seed = seeds[i]
                row, err = measure_montecarlo_seed(
                    spec,
                    planet,
                    mission_time_s,
                    seed;
                    variant=variant,
                    mars=mars,
                    outer_route=:none
                )
                seed_rows[i] = row
                if row.solve_success
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
                else
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): $(err)"
                end
            end
        end

        for i in eachindex(seeds)
            push!(rows, seed_rows[i])
            println(seed_msgs[i])
        end
        _record_outer_route_feedback!(warmup_case, seed_rows; route=mc_backend)
    end
    return nothing
end

function run_benchmarks(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    selected = selected_cases(spec, cases)
    selected = _split_rollout_benchmark_cases(selected)
    selected = _multirate_rollout_benchmark_cases(selected)
    rows = NamedTuple[]
    total = length(selected)
    backend = perf_parallel_backend()
    if backend == :auto
        case_rows = Vector{Vector{NamedTuple}}(undef, total)
        chosen_routes = fill(:none, total)
        process_tasks = Tuple{Int, BenchmarkCase}[]
        thread_indices = Int[]
        serial_indices = Int[]

        for (idx, case) in enumerate(selected)
            route = auto_backend_for_case(case; spec=spec)
            chosen_routes[idx] = route
            if route == :process
                push!(process_tasks, (idx, case))
            elseif route == :threads
                if _case_outer_threads_safe(case)
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
            process_rows = pmap(process_tasks) do task
                idx, case = task
                perf_worker_run_case_batch(case, spec, idx, total, :process)
            end
            for (k, task) in enumerate(process_tasks)
                idx = task[1]
                case_rows[idx] = process_rows[k]
            end
        end

        if !isempty(thread_indices)
            for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
                withenv(env_pairs...) do
                    Threads.@threads for j in eachindex(payload)
                        idx, plan = payload[j]
                        local_rows = NamedTuple[]
                        run_case_batch!(
                            local_rows,
                            selected[idx],
                            spec,
                            idx,
                            total;
                            outer_route=:threads,
                            plan=plan,
                            apply_env=false
                        )
                        case_rows[idx] = local_rows
                    end
                end
            end
        end

        for idx in serial_indices
            local_rows = NamedTuple[]
            run_case_batch!(local_rows, selected[idx], spec, idx, total; outer_route=:none)
            case_rows[idx] = local_rows
        end

        for idx in eachindex(case_rows)
            _record_outer_route_feedback!(selected[idx], case_rows[idx]; route=chosen_routes[idx])
            append!(rows, case_rows[idx])
        end
    elseif backend == :process
        ensure_perf_workers!()
        tasks = collect(enumerate(selected))
        case_rows = pmap(tasks) do task
            idx, case = task
            perf_worker_run_case_batch(case, spec, idx, total, :process)
        end
        for idx in eachindex(case_rows)
            _record_outer_route_feedback!(selected[idx], case_rows[idx]; route=:process)
            append!(rows, case_rows[idx])
        end
    elseif backend == :threads
        case_rows = Vector{Vector{NamedTuple}}(undef, total)
        thread_indices, serial_indices = _split_threadsafe_indices(selected, collect(eachindex(selected)))
        for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
            withenv(env_pairs...) do
                Threads.@threads for j in eachindex(payload)
                    idx, plan = payload[j]
                    local_rows = NamedTuple[]
                    run_case_batch!(
                        local_rows,
                        selected[idx],
                        spec,
                        idx,
                        total;
                        outer_route=:threads,
                        plan=plan,
                        apply_env=false
                    )
                    case_rows[idx] = local_rows
                end
            end
        end
        for idx in serial_indices
            local_rows = NamedTuple[]
            run_case_batch!(local_rows, selected[idx], spec, idx, total; outer_route=:none)
            case_rows[idx] = local_rows
        end
        for idx in eachindex(case_rows)
            route = idx in serial_indices ? :none : :threads
            _record_outer_route_feedback!(selected[idx], case_rows[idx]; route=route)
            append!(rows, case_rows[idx])
        end
    else
        for (idx, case) in enumerate(selected)
            local_rows = NamedTuple[]
            run_case_batch!(local_rows, case, spec, idx, total; outer_route=:none)
            _record_outer_route_feedback!(case, local_rows; route=:none)
            append!(rows, local_rows)
        end
    end

    run_montecarlo_batch!(rows, spec, planet)
    return DataFrame(rows)
end
@inline function selected_cases(spec::ProfileSpec, cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    selected = spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
    if _perf_smoke_mode()
        return selected[1:min(end, 3)]
    end
    return selected
end

@inline function _entry_duration_interface_counts(spec::ProfileSpec)::Vector{Int}
    if _perf_smoke_mode()
        return [1]
    end
    return spec.name == "full" ? collect(1:3) : collect(1:2)
end

@inline function _entry_duration_selected_cases(spec::ProfileSpec, cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    selected = [c for c in (spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]) if c.category == "entry"]
    if _perf_smoke_mode()
        return selected[1:min(end, 1)]
    end
    return selected
end

@inline function _entry_sweep_case(base_case::BenchmarkCase, interface_count::Int)::BenchmarkCase
    interface_count > 0 || throw(ArgumentError("entry interface count must be > 0, got $interface_count"))
    return BenchmarkCase(
        name=base_case.name,
        category=base_case.category,
        description=base_case.description,
        args_template=base_case.args_template,
        run_in_quick=base_case.run_in_quick,
        solver_mode_override=base_case.solver_mode_override,
        split_imex_solver_override=base_case.split_imex_solver_override,
        entry_target_count_override=interface_count
    )
end
