@inline function _maybe_add_control_tstop!(integrator, tstop::Float64)
    if !(isfinite(tstop) && tstop > integrator.t)
        return nothing
    end
    current_max = try
        DiffEqBase.get_tstops_max(integrator)
    catch
        Inf
    end
    if tstop <= current_max + 1e-9 && applicable(DiffEqBase.add_tstop!, integrator, tstop)
        DiffEqBase.add_tstop!(integrator, tstop)
    end
    return nothing
end

@inline function _register_control_tstops!(integrator, control_model::BaseThrusterModel, sat_idx::Int)
    if sat_idx < 1 || sat_idx > length(control_model.start_burn_time)
        return nothing
    end
    _maybe_add_control_tstop!(integrator, control_model.start_burn_time[sat_idx])
    _maybe_add_control_tstop!(integrator, control_model.stop_burn_time[sat_idx])
    return nothing
end

@inline control_requires_periodic_callback(::Any) = true
@inline control_requires_periodic_callback(::BaseThrusterModel) = false

@inline function _run_guidance_for_thruster_schedule!(integrator, sat_idx::Int)
    @inbounds for guidance_model in integrator.p.args.guidance_model.guidance_effectors
        calcGuidanceEffect!(guidance_model, integrator.u, integrator.p, integrator.t, sat_idx)
    end
    return nothing
end

@inline function _schedule_thruster_control!(integrator, control_model::BaseThrusterModel, sat_idx::Int)
    _run_guidance_for_thruster_schedule!(integrator, sat_idx)
    calcControlEffect!(control_model, integrator.u, integrator.p, integrator.t, sat_idx)
    _register_control_tstops!(integrator, control_model, sat_idx)
    return nothing
end

function schedule_event_driven_thruster_controls!(integrator, sat_idx::Int)
    @inbounds for control_model in integrator.p.args.control_model.control_effectors
        if control_model isa BaseThrusterModel
            _schedule_thruster_control!(integrator, control_model, sat_idx)
        end
    end
    return nothing
end

function _thruster_schedule_callbacks(control_model::BaseThrusterModel, num_sats::Int, args::SimulationConfiguration)
    n_slots = length(control_model.thrust)
    if n_slots != num_sats
        throw(ArgumentError(
            "BaseThrusterModel vector length ($n_slots) must match number of spacecraft ($num_sats). " *
            "Use one shared model with per-spacecraft vectors."
        ))
    end

    function schedule_all!(integrator)
        @inbounds for sat_idx in 1:num_sats
            _schedule_thruster_control!(integrator, control_model, sat_idx)
        end
        return nothing
    end

    init_callback = DiscreteCallback(
        (u, t, integrator) -> false,
        schedule_all!;
        initialize=(cb, u, t, integrator) -> schedule_all!(integrator)
    )
    return (init_callback,)
end

function get_control_callbacks(num_sats::Int, args::SimulationConfiguration)
    # Perform the control effects' calculations at specific rates given by the control_rates field in the ControlModel
    control_models = args.control_model.control_effectors
    control_rates = args.control_model.control_rates
    callbacks = Any[]
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
        if !control_requires_periodic_callback(control_model)
            append!(callbacks, _thruster_schedule_callbacks(control_model, num_sats, args))
            continue
        end
        # Each control effector callback runs at its own rate and updates
        # all spacecraft states. The spacecraft index is passed explicitly
        # to avoid conflating effector-index with spacecraft-index.
        apply_control! = (integrator, sat_idx) -> begin
            calcControlEffect!(control_model, integrator.u, integrator.p, integrator.t, sat_idx)
            if control_model isa BaseThrusterModel
                _register_control_tstops!(integrator, control_model, sat_idx)
            end
        end
        control_func = (integrator) -> begin
            decision = _control_callback_thread_decision(integrator.p, control_model, num_sats)
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
        push!(callbacks, PeriodicCallback(control_func, control_rate))
    end
    return callbacks
end
