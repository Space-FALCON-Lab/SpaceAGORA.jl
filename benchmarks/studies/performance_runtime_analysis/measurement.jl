function measure_case(
    case::BenchmarkCase,
    profile_name::String,
    repeat_idx::Int;
    seed::Union{Missing, Int}=missing,
    attempt::Int=1,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing
)
    timestamp_utc = string(now(UTC))
    args_meta = case.args_template
    n_sats = length(args_meta.dynamics_model.spacecraft)
    mission_time_s = args_meta.mission_configuration.mission_time
    density_model = args_meta.environment_model.density_model
    density_family = _density_model_family(density_model)
    gram_surrogate_enabled = _gram_surrogate_flag(density_model)
    gram_static_grid_enabled = (
        density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate
    ) && (
        (_env_bool_token(get(ENV, "SPACEAGORA_GRAM_STATIC_GRID", "off")) === true) ||
        _gram_static_grid_flag_for_case(case, density_model)
    )
    gram_track_cache_mode = _gram_track_cache_mode_for_case(case)
    density_backend_bucket = _density_backend_bucket(
        density_family,
        gram_surrogate_enabled,
        gram_static_grid_enabled,
        gram_track_cache_mode
    )
    resolved_plan = isnothing(plan) ? ParallelPriorityPlan() : plan
    hardware = _runtime_hardware_snapshot()

    GC.gc()
    copy_timed = @timed deepcopy(case.args_template)
    args_run = copy_timed.value

    GC.gc()
    solve_timed = @timed begin
        try
            result = _run_perf_simulation(
                args_run;
                return_solution=true,
                return_solver_metadata=true,
                profile_name=profile_name,
                solver_mode_override=case.solver_mode_override,
                split_imex_solver_override=case.split_imex_solver_override,
                entry_target_count_override=case.entry_target_count_override
            )
            (ok=true, result=result, err=nothing, bt=nothing)
        catch err
            if err isa InterruptException
                rethrow()
            end
            if _perf_strict_errors()
                rethrow(err)
            end
            (ok=false, result=nothing, err=err, bt=catch_backtrace())
        end
    end

    copy_time_s = copy_timed.time
    solve_time_s = solve_timed.time
    total_time_s = copy_time_s + solve_time_s

    copy_bytes_mb = copy_timed.bytes / 1024^2
    solve_bytes_mb = solve_timed.bytes / 1024^2
    total_bytes_mb = copy_bytes_mb + solve_bytes_mb

    copy_alloc_calls = _alloc_calls(copy_timed.gcstats)
    solve_alloc_calls = _alloc_calls(solve_timed.gcstats)

    solver_mode = missing
    solver_sequence = missing
    solver_fallback_used = missing
    solver_fallback_count = missing
    solver_fallback_trigger = missing
    solve_retcode = missing
    solve_success = false
    saved_points = missing
    accepted_steps = missing
    rejected_steps = missing
    sim_seconds_per_wall_second = missing
    satellite_sim_seconds_per_wall_second = missing
    policy_decisions_total = missing
    policy_threads_enabled_total = missing
    policy_density_threads_enabled = missing
    policy_control_threads_enabled = missing
    policy_multibody_threads_enabled = missing
    policy_other_threads_enabled = missing
    nbody_spkpos_runtime_calls = missing
    nbody_spkpos_cache_build_calls = missing
    nbody_spkpos_total_calls = missing
    srp_spkpos_runtime_calls = missing
    srp_spkpos_cache_build_calls = missing
    srp_spkpos_total_calls = missing
    planet_pxform_runtime_calls = missing
    planet_pxform_cache_build_calls = missing
    planet_pxform_total_calls = missing
    final_primary_pos_norm_m = missing
    final_primary_vel_norm_mps = missing
    final_primary_mass_kg = missing
    terminal_time_s = missing
    entry_target_count = isnothing(case.entry_target_count_override) ? missing : Int(case.entry_target_count_override)

    solve_payload = solve_timed.value
    if solve_payload.ok
        solve_result = solve_payload.result
        sol = solve_result.solution
        solver_mode = solve_result.solver_mode
        solver_trace = solve_result.solver_trace
        solver_sequence = isempty(solver_trace) ? missing : join([meta.solver for meta in solver_trace], "->")
        solver_fallback_count = count(meta -> meta.fallback_used, solver_trace)
        solver_fallback_used = solver_fallback_count > 0
        fallback_triggers = [meta.trigger_retcode for meta in solver_trace if !(meta.trigger_retcode isa Missing)]
        solver_fallback_trigger = isempty(fallback_triggers) ? missing : _safe_unique_join(fallback_triggers; delimiter="|")
        solve_retcode = string(sol.retcode)
        solve_success = _solve_success_for_case(sol, case)
        saved_points = length(sol.t)
        accepted_steps = _destat_int(sol, :naccept)
        rejected_steps = _destat_int(sol, :nreject)
        sim_seconds_per_wall_second = _safe_div(mission_time_s, total_time_s)
        satellite_sim_seconds_per_wall_second = _safe_div(mission_time_s * n_sats, total_time_s)
        terminal = _primary_terminal_state_metrics(sol)
        final_primary_pos_norm_m = terminal.pos_norm_m
        final_primary_vel_norm_mps = terminal.vel_norm_mps
        final_primary_mass_kg = terminal.mass_kg
        if !isempty(sol.t)
            terminal_time_s = Float64(sol.t[end])
        end
        if hasproperty(solve_result, :parallel_policy)
            snapshot = solve_result.parallel_policy
            if !(snapshot isa Nothing)
                policy_decisions_total = getproperty(snapshot, :decisions_total)
                policy_threads_enabled_total = getproperty(snapshot, :threads_enabled_total)
                policy_density_threads_enabled = getproperty(snapshot, :density_threads_enabled)
                policy_control_threads_enabled = getproperty(snapshot, :control_threads_enabled)
                policy_multibody_threads_enabled = getproperty(snapshot, :multibody_threads_enabled)
                policy_other_threads_enabled = getproperty(snapshot, :other_threads_enabled)
            end
        end
        if hasproperty(solve_result, :spice_counters)
            counters = solve_result.spice_counters
            if !(counters isa Nothing)
                nbody_spkpos_runtime_calls = getproperty(counters, :nbody_spkpos_runtime_calls)
                nbody_spkpos_cache_build_calls = getproperty(counters, :nbody_spkpos_cache_build_calls)
                nbody_spkpos_total_calls = getproperty(counters, :nbody_spkpos_total_calls)
                srp_spkpos_runtime_calls = getproperty(counters, :srp_spkpos_runtime_calls)
                srp_spkpos_cache_build_calls = getproperty(counters, :srp_spkpos_cache_build_calls)
                srp_spkpos_total_calls = getproperty(counters, :srp_spkpos_total_calls)
                planet_pxform_runtime_calls = getproperty(counters, :planet_pxform_runtime_calls)
                planet_pxform_cache_build_calls = getproperty(counters, :planet_pxform_cache_build_calls)
                planet_pxform_total_calls = getproperty(counters, :planet_pxform_total_calls)
            end
        end
    else
        solve_err = solve_payload.err
        solve_bt = solve_payload.bt
        solve_retcode = _solve_retcode_from_error(solve_err)
        if ismissing(solve_retcode)
            solve_retcode = "Exception"
            @warn "[perf] case=$(case.name) repeat=$(repeat_idx) attempt=$(attempt) threw $(typeof(solve_err)): $(_perf_error_text(solve_err)) @ $(_perf_stack_head(solve_bt))"
        end
    end

    return (
        profile=profile_name,
        hardware_class=hardware.hardware_class,
        machine_label=hardware.machine_label,
        host_name=hardware.host_name,
        cpu_name=hardware.cpu_name,
        cpu_threads=hardware.cpu_threads,
        julia_threads=hardware.julia_threads,
        os=hardware.os,
        arch=hardware.arch,
        category=case.category,
        scenario=case.name,
        description=case.description,
        density_family=density_family,
        gram_surrogate_enabled=gram_surrogate_enabled,
        gram_static_grid_enabled=gram_static_grid_enabled,
        gram_track_cache_mode=gram_track_cache_mode,
        density_backend_bucket=density_backend_bucket,
        seed=seed,
        repeat=repeat_idx,
        attempt=attempt,
        satellites=n_sats,
        orientation=args_meta.mission_configuration.orientation_sim,
        mission_time_s=mission_time_s,
        outer_route=string(resolved_plan.outer_route),
        outer_threads_safe=_case_outer_threads_safe(case),
        density_parallel_mode=resolved_plan.density_mode,
        control_parallel_mode=resolved_plan.control_mode,
        multibody_parallel_mode=resolved_plan.multibody_mode,
        dt_max_orbit_s=args_meta.integration_tolerances.dt_max_orbit,
        dynamic_effectors=_effector_signature(args_meta.dynamics_model.dynamic_effectors),
        control_effectors=_effector_signature(args_meta.control_model.control_effectors),
        copy_time_s=copy_time_s,
        solve_time_s=solve_time_s,
        total_time_s=total_time_s,
        copy_compile_time_s=copy_timed.compile_time,
        solve_compile_time_s=solve_timed.compile_time,
        copy_gctime_s=copy_timed.gctime,
        solve_gctime_s=solve_timed.gctime,
        solver_mode=solver_mode,
        split_imex_solver_override=isnothing(case.split_imex_solver_override) ? missing : String(case.split_imex_solver_override),
        solver_sequence=solver_sequence,
        solver_fallback_used=solver_fallback_used,
        solver_fallback_count=solver_fallback_count,
        solver_fallback_trigger=solver_fallback_trigger,
        solve_retcode=solve_retcode,
        solve_success=solve_success,
        copy_bytes_mb=copy_bytes_mb,
        solve_bytes_mb=solve_bytes_mb,
        total_bytes_mb=total_bytes_mb,
        copy_alloc_calls=copy_alloc_calls,
        solve_alloc_calls=solve_alloc_calls,
        saved_points=saved_points,
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
        policy_decisions_total=policy_decisions_total,
        policy_threads_enabled_total=policy_threads_enabled_total,
        policy_density_threads_enabled=policy_density_threads_enabled,
        policy_control_threads_enabled=policy_control_threads_enabled,
        policy_multibody_threads_enabled=policy_multibody_threads_enabled,
        policy_other_threads_enabled=policy_other_threads_enabled,
        nbody_spkpos_runtime_calls=nbody_spkpos_runtime_calls,
        nbody_spkpos_cache_build_calls=nbody_spkpos_cache_build_calls,
        nbody_spkpos_total_calls=nbody_spkpos_total_calls,
        srp_spkpos_runtime_calls=srp_spkpos_runtime_calls,
        srp_spkpos_cache_build_calls=srp_spkpos_cache_build_calls,
        srp_spkpos_total_calls=srp_spkpos_total_calls,
        planet_pxform_runtime_calls=planet_pxform_runtime_calls,
        planet_pxform_cache_build_calls=planet_pxform_cache_build_calls,
        planet_pxform_total_calls=planet_pxform_total_calls,
        final_primary_pos_norm_m=final_primary_pos_norm_m,
        final_primary_vel_norm_mps=final_primary_vel_norm_mps,
        final_primary_mass_kg=final_primary_mass_kg,
        terminal_time_s=terminal_time_s,
        entry_target_count=entry_target_count,
        sim_seconds_per_wall_second=sim_seconds_per_wall_second,
        satellite_sim_seconds_per_wall_second=satellite_sim_seconds_per_wall_second,
        timestamp_utc=timestamp_utc
    )
end

function run_warmup(case::BenchmarkCase, warmup::Int, profile_name::String)
    log_warmup = _perf_warmup_logs()
    for i in 1:warmup
        args_run = deepcopy(case.args_template)
        if log_warmup
            println("  warmup $(i)/$(warmup): start")
            flush(stdout)
        end
        warmup_started_ns = time_ns()
        try
            _run_perf_simulation(
                args_run;
                return_solution=false,
                profile_name=profile_name,
                solver_mode_override=case.solver_mode_override,
                split_imex_solver_override=case.split_imex_solver_override,
                entry_target_count_override=case.entry_target_count_override
            )
            warmup_elapsed_s = (time_ns() - warmup_started_ns) / 1e9
            if log_warmup
                println("  warmup $(i)/$(warmup): done total=$(round(warmup_elapsed_s; digits=3)) s")
                flush(stdout)
            end
        catch err
            if err isa InterruptException
                rethrow()
            end
            if _perf_strict_errors()
                rethrow(err)
            end
            err_bt = catch_backtrace()
            solve_retcode = _solve_retcode_from_error(err)
            warmup_elapsed_s = (time_ns() - warmup_started_ns) / 1e9
            if ismissing(solve_retcode)
                @warn "[warmup] $(case.name) $(i)/$(warmup) threw $(typeof(err)): $(_perf_error_text(err)) @ $(_perf_stack_head(err_bt)); continuing (elapsed=$(round(warmup_elapsed_s; digits=3)) s)"
            elseif log_warmup
                println("  warmup $(i)/$(warmup): failed retcode=$(solve_retcode), continuing (elapsed=$(round(warmup_elapsed_s; digits=3)) s)")
                flush(stdout)
            end
        end
    end
    return nothing
end

@inline function _case_sample_schedule(case::BenchmarkCase, spec::ProfileSpec)::Tuple{Int, Int}
    warmup = spec.warmup
    repeats = spec.repeats
    if case.category == "control_stress"
        if spec.name == "full"
            warmup = _control_stress_warmup_full()
            repeats = _control_stress_repeats_full()
        else
            warmup = min(warmup, 1)
            repeats = min(repeats, 2)
        end
    end
    return warmup, repeats
end

function _run_case_batch_core!(
    rows::Vector{NamedTuple},
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int,
    plan::ParallelPriorityPlan
)
    warmup_count, repeat_count = _case_sample_schedule(case, spec)
    println(
        "[$idx/$total] $(case.name) :: warmup x$(warmup_count), repeats x$(repeat_count), " *
        "outer=$(plan.outer_route), density=$(plan.density_mode), control=$(plan.control_mode), multibody=$(plan.multibody_mode), " *
        "solver_override=$(isnothing(case.solver_mode_override) ? "none" : case.solver_mode_override), " *
        "split_override=$(isnothing(case.split_imex_solver_override) ? "none" : case.split_imex_solver_override)"
    )
    run_warmup(case, warmup_count, spec.name)
    for rep in 1:repeat_count
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, rep; attempt=attempt, plan=plan)
            last_row = row
            if row.solve_success
                push!(rows, row)
                println("  repeat $(rep)/$(repeat_count) attempt $(attempt)/$(spec.max_attempts): total=$(round(row.total_time_s; digits=3)) s, solve=$(round(row.solve_time_s; digits=3)) s")
                flush(stdout)
                break
            end
            println("  repeat $(rep)/$(repeat_count) attempt $(attempt)/$(spec.max_attempts): failed retcode=$(row.solve_retcode), retrying")
            flush(stdout)
        end
        if !(last_row === nothing) && !last_row.solve_success
            push!(rows, last_row)
            println("  repeat $(rep)/$(repeat_count): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)")
            flush(stdout)
        end
    end
    return nothing
end

function run_case_batch!(
    rows::Vector{NamedTuple},
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int;
    outer_route::Symbol=:none,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing,
    apply_env::Bool=true
)
    resolved_plan = isnothing(plan) ? parallel_priority_plan(case, outer_route) : plan
    if apply_env
        env_pairs = case_env_pairs(case, resolved_plan)
        return withenv(env_pairs...) do
            _run_case_batch_core!(rows, case, spec, idx, total, resolved_plan)
        end
    end
    return _run_case_batch_core!(rows, case, spec, idx, total, resolved_plan)
end

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
    selected = spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
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
    return spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
end

@inline function _entry_duration_interface_counts(spec::ProfileSpec)::Vector{Int}
    return spec.name == "full" ? collect(1:3) : collect(1:2)
end

@inline function _entry_duration_selected_cases(spec::ProfileSpec, cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    return [c for c in selected_cases(spec, cases) if c.category == "entry"]
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
    orbit_counts = spec.name == "full" ? collect(1:5) : collect(1:3)
    include_control_stress = _include_control_stress_per_orbit()
    selected_base = include_control_stress ? selected_cases(spec, cases) : [c for c in selected_cases(spec, cases) if c.category != "control_stress"]
    # Entry scenarios use entry-interface counting semantics, not baseline-period multipliers.
    selected = [c for c in selected_base if c.category != "entry"]
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

