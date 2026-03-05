function get_control_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    # Perform the control effects' calculations at specific rates given by the control_rates field in the ControlModel
    control_models = args.control_model.control_effectors
    control_rates = args.control_model.control_rates
    use_invokelatest = callback_use_invokelatest()
    callbacks = Vector{DiscreteCallback}(undef, length(control_models))
    for i in eachindex(control_models)
        control_model = control_models[i]
        control_rate = control_rates[i]
        if control_model isa BaseThrusterModel
            n_slots = length(control_model.thrust)
            if n_slots != num_sats
                throw(ArgumentError(
                    "BaseThrusterModel vector length ($n_slots) must match number of spacecraft ($num_sats). " *
                    "Use one shared model with per-spacecraft vectors."
                ))
            end
        end
        # Each control effector callback runs at its own rate and updates
        # all spacecraft states. The spacecraft index is passed explicitly
        # to avoid conflating effector-index with spacecraft-index.
        apply_control! = if use_invokelatest
            # Dev mode: keep Revise/hot-reload workflows free of world-age errors.
            (integrator, sat_idx) -> Base.invokelatest(calcControlEffect!, control_model, integrator.u, integrator.p, integrator.t, sat_idx)
        else
            # Production mode: direct dispatch avoids invokelatest overhead.
            (integrator, sat_idx) -> calcControlEffect!(control_model, integrator.u, integrator.p, integrator.t, sat_idx)
        end
        control_func = (integrator) -> begin
            decision = _control_callback_thread_decision(control_model, num_sats, use_invokelatest)
            use_threads = decision.use_threads
            started_ns = time_ns()
            if use_threads
                ParallelPolicy.threaded_foreach_persistent(:control_callback, num_sats, decision.allotment) do sat_idx
                    @inbounds apply_control!(integrator, sat_idx)
                end
            else
                @inbounds for sat_idx in 1:num_sats
                    apply_control!(integrator, sat_idx)
                end
            end
            if decision.policy_applied
                ParallelPolicy.record_policy_observation!(
                    :control_callback;
                    mode=decision.mode,
                    num_items=num_sats,
                    use_threads=use_threads,
                    elapsed_ns=(time_ns() - started_ns)
                )
            end
        end
        callbacks[i] = PeriodicCallback(control_func, control_rate)
    end
    return callbacks
end

