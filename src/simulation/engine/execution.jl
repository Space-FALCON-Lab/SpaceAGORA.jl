@inline function _build_typed_solver_problem(u0, tspan, p, callbacks)
    mode = _solver_policy_mode()
    if mode == :split_imex || mode == :multirate
        return SplitODEProblem(
            spacecraft_dynamics_slow!,
            spacecraft_dynamics_fast_control!,
            u0,
            tspan,
            p,
            callback=callbacks
        )
    end
    return ODEProblem(spacecraft_dynamics!, u0, tspan, p, callback=callbacks)
end

function run_simulation(
    args;
    isolate_state::Bool=true,
    return_solution::Bool=false,
    return_solver_metadata::Bool=false,
    save_fields=nothing
)
    return SimulationModel.ParallelPolicy.with_policy_context() do
    # Isolate mutable campaign/model state by default so repeated/concurrent runs
    # do not alias shared in-memory objects.
    args = isolate_state ? deepcopy(args) : args

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
        save_fields=save_fields_resolved
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
    prob_debug = ODEProblem(spacecraft_dynamics!, u_start, (t_start, mission_end), p, callback=callbacks)
    if p.shared_buffers.debug_initial_derivative[]
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
    end

    reltol_tol, abstol_tol = _build_solver_tolerances(u_start, args)
    last_sol = nothing
    solver_trace = NamedTuple[]
    checkpoint_saved_times = Float64[]
    checkpoint_saved_data = SimulationModel.SaveData[]

    if t_start < mission_end && checkpoint_active
        interval = args.simulation_settings.checkpoint_interval_s
        t_cursor = t_start
        u_cursor = deepcopy(u_start)

        while t_cursor < mission_end
            t_next = min(t_cursor + interval, mission_end)
            empty!(saved_values.t)
            empty!(saved_values.saveval)
            prob = _build_typed_solver_problem(u_cursor, (t_cursor, t_next), p, callbacks)
            seg_sol, solve_meta = _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
            push!(solver_trace, solve_meta)
            if !SciMLBase.successful_retcode(seg_sol.retcode)
                throw(ErrorException("Checkpointed solve failed with retcode=$(seg_sol.retcode)."))
            end
            _append_saved_segment!(checkpoint_saved_times, checkpoint_saved_data, saved_values)
            last_sol = seg_sol
            t_cursor = Float64(seg_sol.t[end])
            u_cursor = deepcopy(seg_sol.u[end])
            _write_checkpoint!(args, t_cursor, u_cursor)
        end

    elseif t_start < mission_end
        prob = _build_typed_solver_problem(u_start, (t_start, mission_end), p, callbacks)
        sol, solve_meta = _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
        push!(solver_trace, solve_meta)
        if !SciMLBase.successful_retcode(sol.retcode)
            throw(ErrorException("Solve failed with retcode=$(sol.retcode)."))
        end
        last_sol = sol
    end

    # Process and save results
    if args.simulation_settings.results
        results_times = checkpoint_active ? checkpoint_saved_times : saved_values.t
        results_data = checkpoint_active ? checkpoint_saved_data : saved_values.saveval
        results_df = _build_results_dataframe(results_times, results_data, save_fields_resolved, args)
        # Keep backwards-compatible CSV contract used by existing scripts/tests.
        csv_path = _write_results_csv!(results_df, args)
        if _typed_save_bundle_enabled()
            _write_results_bundle!(results_df, results_times, args; csv_path=csv_path)
        end
    end

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
                solver_mode=string(_solver_policy_mode()),
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
