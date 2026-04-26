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

