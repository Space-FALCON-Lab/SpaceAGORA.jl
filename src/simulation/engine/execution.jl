function _build_block_diagonal_jac_prototype(u0::ComponentVector)::SparseMatrixCSC{Float64, Int}
    n_sats = length(u0.sc)
    n_total = length(u0)
    rows = Int[]
    cols = Int[]
    # Probe each satellite block: set its elements to 1.0 in a zero copy of u0,
    # read back the flat indices, then declare all pairs within that block as
    # structurally non-zero. This handles heterogeneous state sizes (different
    # n_bodies / heat_loads per satellite) without assuming a uniform block size.
    for i in 1:n_sats
        probe = zero(u0)
        probe.sc[i] .= 1.0
        flat = Vector{Float64}(probe)
        sat_idx = findall(!iszero, flat)
        for k in sat_idx, j in sat_idx
            push!(rows, k)
            push!(cols, j)
        end
    end
    return sparse(rows, cols, ones(Float64, length(rows)), n_total, n_total)
end

@inline function _build_typed_solver_problem(u0, tspan, p, callbacks, solver_mode::Symbol,
    jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int}}=nothing)
    mode = solver_mode
    if mode == :gravity_backbone_split
        q0, dq0 = _gravity_backbone_initial_states(u0, p.args)
        return SecondOrderODEProblem(
            spacecraft_dynamics_gravity_backbone!,
            dq0,
            q0,
            tspan,
            p,
            callback=callbacks
        )
    elseif mode == :split_imex
        return SplitODEProblem(
            spacecraft_dynamics_implicit_atmosphere!,
            spacecraft_dynamics_explicit_remainder!,
            u0,
            tspan,
            p,
            callback=callbacks
        )
    elseif mode == :multirate
        return SplitODEProblem(
            spacecraft_dynamics_slow!,
            spacecraft_dynamics_fast_control!,
            u0,
            tspan,
            p,
            callback=callbacks
        )
    end
    if jac_prototype !== nothing
        f = ODEFunction(spacecraft_dynamics!; jac_prototype=jac_prototype)
        return ODEProblem(f, u0, tspan, p; callback=callbacks)
    end
    return ODEProblem(spacecraft_dynamics!, u0, tspan, p, callback=callbacks)
end

@inline function _append_backbone_saved_segment!(
    times_acc::Vector{Float64},
    data_acc::Vector{SimulationModel.SaveData},
    sol,
    save_fields,
    p,
)
    length(sol.t) <= 1 && return nothing
    integrator_view = (p=p,)
    @inbounds for idx in 2:length(sol.t)
        t_sample = Float64(sol.t[idx])
        push!(times_acc, t_sample)
        push!(data_acc, SimulationModel.SimulationCallbacks._save_snapshot(save_fields, sol.u[idx], t_sample, integrator_view))
    end
    return nothing
end

@inline function _append_checkpoint_saved_segment!(
    times_acc::Vector{Float64},
    data_acc::Vector,
    saved_values,
    sol,
    save_fields,
    p,
)
    _append_saved_segment!(times_acc, data_acc, saved_values)
    isempty(sol.t) && return nothing

    t_final = Float64(sol.t[end])
    if isempty(times_acc) || !isapprox(times_acc[end], t_final; atol=0.0, rtol=0.0)
        integrator_view = (p=p,)
        push!(times_acc, t_final)
        push!(data_acc, SimulationModel.SimulationCallbacks._save_snapshot(save_fields, sol.u[end], t_final, integrator_view))
    end
    return nothing
end

function _save_simulation_results_if_enabled!(
    args,
    solver_mode::Symbol,
    checkpoint_active::Bool,
    save_fields_resolved,
    saved_values,
    checkpoint_saved_times,
    checkpoint_saved_data,
    backbone_saved_times,
    backbone_saved_data,
)
    args.simulation_settings.results || return nothing
    results_times = if solver_mode == :gravity_backbone_split
        checkpoint_active ? checkpoint_saved_times : backbone_saved_times
    else
        checkpoint_active ? checkpoint_saved_times : saved_values.t
    end
    results_data = if solver_mode == :gravity_backbone_split
        checkpoint_active ? checkpoint_saved_data : backbone_saved_data
    else
        checkpoint_active ? checkpoint_saved_data : saved_values.saveval
    end
    results_df = _build_results_dataframe(results_times, results_data, save_fields_resolved, args)
    csv_path = _write_results_csv!(results_df, args)
    if _typed_save_bundle_enabled()
        _write_results_bundle!(results_df, results_times, args; csv_path=csv_path)
    end
    return csv_path
end

function _try_save_simulation_results_if_enabled!(args...)
    try
        return _save_simulation_results_if_enabled!(args...)
    catch err
        @warn "Failed to save partial simulation results after solve failure." reason=sprint(showerror, err)
        return nothing
    end
end

function run_simulation(
    args::SimulationConfiguration;
    isolate_state::Bool=true,
    return_solution::Bool=false,
    return_solver_metadata::Bool=false,
    save_fields=nothing,
    extra_callbacks=(),
    solver_cache::Union{Nothing, SolverIntegratorCache}=nothing
)
    return SimulationModel.ParallelPolicy.with_policy_context() do
    # Isolate mutable campaign/model state by default so repeated/concurrent runs
    # do not alias shared in-memory objects.
    args = isolate_state ? deepcopy(args) : args
    # Resolve effective solver config: explicit field on args wins; otherwise read from env
    # (which respects any active SimulationEngineConfig overrides via _engine_env_get).
    solver_cfg = isnothing(args.solver_config) ? _solver_config_from_env() : args.solver_config
    solver_mode = solver_cfg.solver_mode

    # Typed pipeline is SI-native (meters, seconds, kilograms). The
    # `simulation_settings.normalize` field is legacy-only and rejected by default.
    _enforce_typed_normalize_policy!(args)
    _validate_orientation_inertia!(args)
    _validate_thermal_model_support!(args)
    _validate_ephemerides_support!(args)
    try
        SimulationModel.ParallelPolicy.reset_policy_telemetry!()
        if SimulationModel.ParallelPolicy.persistent_hints_state_reset_requested()
            SimulationModel.ParallelPolicy.reset_persistent_hint_state!()
        end
    catch
    end

    # Set up the model and initial conditions
    initial_conditions = build_initial_conditions(args)
    if args.simulation_settings.verbose
        println("Initial conditions:")
        println(initial_conditions)
    end

    # Define the ODE parameters and callbacks
    p = SimulationModel.ODEParams{length(args.dynamics_model.spacecraft)}(args=args) # Define the parameters for the ODE problem, including the shared buffers for the callbacks
    _initialize_heat_rate_buffers!(p)
    _initialize_density_model_instances!(p)
    _initialize_density_cache_buffers!(p)
    _initialize_gram_isolated_pool_buffers!(p)
    _initialize_harmonics_workspace_buffers!(p)
    _initialize_nbody_workspace_buffers!(p)
    _initialize_aero_workspace_buffers!(p)
    _initialize_nbody_ephemeris_cache_buffer!(p)
    _initialize_srp_sun_cache_buffer!(p)
    _initialize_planet_frame_cache_buffer!(p)
    _initialize_spice_rhs_memo_mode!(p)
    _reset_spice_runtime_counters!(p)
    _reset_spice_rhs_memo!(p)
    p.shared_buffers.debug_control[] = _engine_env_get("SPACEAGORA_DEBUG_CONTROL", "0") == "1"
    p.shared_buffers.debug_initial_derivative[] = _engine_env_get("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE", "0") == "1"
    save_fields_resolved = isnothing(save_fields) ? SimulationModel.default_save_fields(args) : collect(save_fields)
    save_field_names = Symbol[field.name for field in save_fields_resolved]
    length(unique(save_field_names)) == length(save_field_names) || throw(ArgumentError("save_fields names must be unique. Got $(save_field_names)."))
    saved_values = SavedValues(Float64, SimulationModel.SaveData)
    callbacks = SimulationModel.get_callbacks(
        length(args.dynamics_model.spacecraft),
        args.dynamics_model.dynamic_effectors,
        args;
        saved_values=saved_values,
        save_fields=save_fields_resolved,
        extra_callbacks=extra_callbacks
    ) # Get the callbacks based on the number of satellites and the dynamic effectors being used in the simulation
    ephemerides_model = args.environment_model.ephemerides_model
    et_start = SimulationModel.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    p.shared_buffers.et_start[] = et_start
    args.environment_model.planet.L_PI .= SimulationModel.planet_frame_lpi(args.environment_model.planet, et_start, ephemerides_model)
    if SimulationModel.ephemerides_requires_spice(ephemerides_model)
        Base.Threads.atomic_add!(p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls, 1)
    end
    mission_end = args.mission_configuration.mission_time
    _initialize_nbody_ephemeris_cache!(p, et_start, mission_end)
    _initialize_srp_sun_ephemeris_cache!(p, et_start, mission_end)
    _initialize_planet_frame_ephemeris_cache!(p, et_start, mission_end)
    checkpoint_active = _typed_checkpoint_enabled(args)
    if checkpoint_active && args.simulation_settings.checkpoint_interval_s <= 0.0
        throw(ArgumentError("SimulationSettings.checkpoint_interval_s must be > 0 when checkpointing is enabled."))
    end

    t_start = 0.0
    u_start = initial_conditions
    if args.simulation_settings.resume_from_checkpoint
        ckpt = _load_checkpoint(args)
        if ckpt === nothing
            @warn "resume_from_checkpoint=true but no checkpoint file was found; starting from initial conditions."
        else
            if solver_mode == :gravity_backbone_split
                ckpt.solver_mode == "gravity_backbone_split" || throw(ArgumentError(
                    "gravity_backbone_split can only resume from checkpoints written by SPACEAGORA_SOLVER_MODE=gravity_backbone_split."
                ))
            elseif ckpt.solver_mode == "gravity_backbone_split"
                throw(ArgumentError(
                    "Checkpoint was written by gravity_backbone_split and cannot be resumed with solver mode $(solver_mode)."
                ))
            end
            t_start = ckpt.t
            u_start = ckpt.u
            if args.simulation_settings.verbose
                println("Resuming simulation from checkpoint at t=$(round(t_start, digits=6)) s")
            end
        end
    end

    # println("Initial conditions:")
    # println(initial_conditions)
    # println("ODE parameters:")
    # println(p)
    # println("args.mission_configuration.mission_time: $(args.mission_configuration.mission_time)")
    p.shared_buffers.solve_segment_end_time[] = mission_end
    prob_debug_state = solver_mode == :gravity_backbone_split ? initial_conditions : u_start
    prob_debug = ODEProblem(spacecraft_dynamics!, prob_debug_state, (t_start, mission_end), p, callback=callbacks)
    if p.shared_buffers.debug_initial_derivative[] && solver_mode != :gravity_backbone_split
        # 1. Manually evaluate the derivative at the start
        du_test = copy(prob_debug.u0)
        try
            prob_debug.f(du_test, prob_debug.u0, prob_debug.p, prob_debug.tspan[1])
        catch e
            @error "The derivative function itself crashed!" exception=e
            throw(ErrorException("Initial derivative evaluation failed; aborting solve in debug mode."))
        end

        # 2. Check for NaNs and print exactly where they are
        if any(isnan, du_test)
            println("--- INITIAL NaN DETECTED ---")

            # Check global parameters in p
            _debug_print_nan_parameter_paths!(prob_debug.p)

            # Check the state vector (u)
            # Assuming your u has a .sc field for satellites
            for (i, sat) in enumerate(du_test.sc)
                if any(isnan, sat.pos) || any(isnan, sat.vel)
                    println("NaN found in Satellite $i derivative!")
                    println("  Pos: $(sat.pos)")
                    println("  Vel: $(sat.vel)")
                end
            end

            throw(ErrorException("Initial derivative contains NaN; aborting solve in debug mode."))
        end
    elseif p.shared_buffers.debug_initial_derivative[] && solver_mode == :gravity_backbone_split
        @warn "Initial derivative debug is not supported for gravity_backbone_split; skipping first-order RHS debug probe."
    end

    reltol_tol, abstol_tol = if solver_mode == :gravity_backbone_split
        args.integration_tolerances.reltol_orbit, args.integration_tolerances.abstol_orbit
    else
        _build_solver_tolerances(u_start, args)
    end
    jac_prototype = (
        solver_mode ∉ (:gravity_backbone_split, :split_imex, :multirate) &&
        u_start isa ComponentVector &&
        length(u_start.sc) > 1
    ) ? _build_block_diagonal_jac_prototype(u_start) : nothing

    # Auto-calibration: time candidate execution plans before the solve and pin the
    # fastest one for the duration.  No-ops when budget <= 1, SPACEAGORA_RHS_CALIBRATE=off,
    # or a cached result already exists for this machine + scenario signature.
    _calibrate_rhs_plan_if_needed!(p, u_start, args)

    last_sol = nothing
    solver_trace = NamedTuple[]
    checkpoint_saved_times = Float64[]
    checkpoint_saved_data = SimulationModel.SaveData[]
    backbone_saved_times = Float64[]
    backbone_saved_data = SimulationModel.SaveData[]

    if t_start < mission_end && checkpoint_active
        interval = args.simulation_settings.checkpoint_interval_s
        t_cursor = t_start
        u_cursor = deepcopy(u_start)

        while t_cursor < mission_end
            t_next = min(t_cursor + interval, mission_end)
            empty!(saved_values.t)
            empty!(saved_values.saveval)
            p.shared_buffers.solve_segment_end_time[] = t_next
            prob = _build_typed_solver_problem(u_cursor, (t_cursor, t_next), p, callbacks, solver_mode, jac_prototype)
            seg_sol, solve_meta = try
                _solve_with_solver_policy(prob, solver_cfg, args, reltol_tol, abstol_tol)
            catch err
                _try_save_simulation_results_if_enabled!(
                    args,
                    solver_mode,
                    checkpoint_active,
                    save_fields_resolved,
                    saved_values,
                    checkpoint_saved_times,
                    checkpoint_saved_data,
                    backbone_saved_times,
                    backbone_saved_data,
                )
                rethrow()
            end
            push!(solver_trace, solve_meta)
            if solver_mode == :gravity_backbone_split
                _append_backbone_saved_segment!(checkpoint_saved_times, checkpoint_saved_data, seg_sol, save_fields_resolved, p)
            else
                _append_checkpoint_saved_segment!(
                    checkpoint_saved_times,
                    checkpoint_saved_data,
                    saved_values,
                    seg_sol,
                    save_fields_resolved,
                    p,
                )
            end
            last_sol = seg_sol
            if !SciMLBase.successful_retcode(seg_sol.retcode)
                _try_save_simulation_results_if_enabled!(
                    args,
                    solver_mode,
                    checkpoint_active,
                    save_fields_resolved,
                    saved_values,
                    checkpoint_saved_times,
                    checkpoint_saved_data,
                    backbone_saved_times,
                    backbone_saved_data,
                )
                throw(ErrorException("Checkpointed solve failed with retcode=$(seg_sol.retcode)."))
            end
            t_cursor = Float64(seg_sol.t[end])
            u_cursor = deepcopy(seg_sol.u[end])
            _write_checkpoint!(args, t_cursor, u_cursor, string(solver_mode))
            if string(seg_sol.retcode) != "Success" || !_gravity_backbone_time_reached(t_cursor, t_next)
                break
            end
        end

    elseif t_start < mission_end
        p.shared_buffers.solve_segment_end_time[] = mission_end
        prob = _build_typed_solver_problem(u_start, (t_start, mission_end), p, callbacks, solver_mode, jac_prototype)
        sol, solve_meta = try
            _solve_with_solver_policy(prob, solver_cfg, args, reltol_tol, abstol_tol; solver_cache=solver_cache)
        catch err
            _try_save_simulation_results_if_enabled!(
                args,
                solver_mode,
                checkpoint_active,
                save_fields_resolved,
                saved_values,
                checkpoint_saved_times,
                checkpoint_saved_data,
                backbone_saved_times,
                backbone_saved_data,
            )
            rethrow()
        end
        push!(solver_trace, solve_meta)
        last_sol = sol
        if solver_mode == :gravity_backbone_split
            _append_backbone_saved_segment!(backbone_saved_times, backbone_saved_data, sol, save_fields_resolved, p)
        end
        if !SciMLBase.successful_retcode(sol.retcode)
            _try_save_simulation_results_if_enabled!(
                args,
                solver_mode,
                checkpoint_active,
                save_fields_resolved,
                saved_values,
                checkpoint_saved_times,
                checkpoint_saved_data,
                backbone_saved_times,
                backbone_saved_data,
            )
            throw(ErrorException("Solve failed with retcode=$(sol.retcode)."))
        end
    end

    # Process and save results
    _save_simulation_results_if_enabled!(
        args,
        solver_mode,
        checkpoint_active,
        save_fields_resolved,
        saved_values,
        checkpoint_saved_times,
        checkpoint_saved_data,
        backbone_saved_times,
        backbone_saved_data,
    )

    if return_solution
        if checkpoint_active && args.simulation_settings.checkpoint_interval_s < mission_end
            @warn "return_solution=true with checkpointed integration returns the final segment ODESolution, not a stitched full-history ODESolution."
        end
        if return_solver_metadata
            parallel_policy = try
                SimulationModel.ParallelPolicy.policy_telemetry_snapshot()
            catch
                nothing
            end
            return (
                solution=last_sol,
                solver_mode=string(solver_mode),
                solver_trace=solver_trace,
                parallel_policy=parallel_policy,
                spice_counters=_spice_runtime_counters_snapshot(p)
            )
        end
        return last_sol
    end
    return nothing
    end
end
