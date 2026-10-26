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
    # The effective window: the scheduler's own, or else the plan installed in
    # its burn-plan buffer (by the burn-plan family).
    start_s, stop_s = _effective_burn_window(control_model, integrator.p, sat_idx)
    _maybe_add_control_tstop!(integrator, start_s)
    _maybe_add_control_tstop!(integrator, stop_s)
    return nothing
end

# After the control tick of an effector in `:mounted` actuation, the engine moves
# the spacecraft's mounted hinges and thruster gimbals toward the levels the
# effector holds.
@inline function _apply_mounted_hinges!(control_model, p, sat_idx::Int, t)
    actuation_mode(control_model) === :mounted || return nothing
    spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
    (mounted_hinge_count(spacecraft) == 0 && mounted_gimbal_count(spacecraft) == 0) && return nothing
    levels = actuator_levels(control_model, sat_idx)
    levels === nothing && return nothing
    applied = faulted_levels(spacecraft, levels)
    apply_hinge_levels!(spacecraft, applied, t)
    apply_gimbal_levels!(spacecraft, applied, t)
    return nothing
end

# After the control tick of an effector in `:mounted` actuation, each mounted
# thruster with firing dynamics schedules the average thrust it delivers until
# the next tick: the command as an on-time within the control period, through
# its minimum on-time, Schmitt trigger and thrust ramp. The schedule follows
# the command; the thruster's faults act on the thrust it delivers
# (`ActuatorMounting._delivered_thrust`), so they switch at their own times.
function _apply_mounted_thruster_firing!(control_model, p, sat_idx::Int, t, period::Float64)
    actuation_mode(control_model) === :mounted || return nothing
    spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
    any(link -> any(th -> th.firing_dynamics, link.thrusters), spacecraft.links) || return nothing
    levels = actuator_levels(control_model, sat_idx)
    (levels === nothing || length(levels) != mounted_actuator_count(spacecraft)) && return nothing
    k = 0
    for link in spacecraft.links, thruster in link.thrusters
        k += 1
        thruster.firing_dynamics || continue
        thrust_calculation_schmitt_trigger!(link, thruster, actuator_saturate(thruster, levels[k]), Float64(t); period=period)
    end
    return nothing
end

@inline function _apply_mounted_tick!(control_model, p, sat_idx::Int, t, period::Float64)
    _apply_mounted_hinges!(control_model, p, sat_idx, t)
    _apply_mounted_thruster_firing!(control_model, p, sat_idx, t, period)
    # What the tick wrote, which the next run's start compares with (`get_run_state_callbacks`).
    actuation_mode(control_model) === :mounted && _note_run_writes!(p.args.dynamics_model.spacecraft[sat_idx])
    return nothing
end

"""
    get_run_state_callbacks(num_sats, args) -> Tuple

For runs whose spacecraft carry hinges, thruster gimbals or thrusters with
firing dynamics: one callback that, at the run's start, starts their run state
(the angles, and the ramp factor and thrust) as a battery's flags and an
actuator's faults are started. A value the previous run's control ticks wrote
returns to its value at that run's start, and a value changed after the last
tick (a user's edit between two runs with `isolate_state=false`) is kept and
becomes the start of later runs (`ActuatorMounting._start_run_state!`). A run
resumed from a checkpoint is not a run start: its run state comes from the
checkpoint, and is recorded as the latest write. Empty without such parts.
"""
function get_run_state_callbacks(num_sats::Int, args::SimulationConfiguration)
    sats = [i for (i, sc) in enumerate(args.dynamics_model.spacecraft) if _has_run_state_parts(sc)]
    isempty(sats) && return ()
    # A fresh callback is built for every run. Its initialization also runs at
    # checkpoint-segment boundaries, which keep the run state.
    initialized = Ref(false)
    function initialize!(cb, u, t, integrator)
        if !initialized[]
            resumed = integrator.p.shared_buffers.callback_resume[]
            for i in sats
                sc = integrator.p.args.dynamics_model.spacecraft[i]
                resumed ? _note_run_writes!(sc) : _start_run_state!(sc)
            end
            initialized[] = true
        end
        _mark_discontinuity!(integrator, false)
        return nothing
    end
    return (DiscreteCallback((u, t, integrator) -> false, integrator -> nothing; initialize=initialize!,
        save_positions=(false, false)),)
end

# The burn's bookkeeping for spacecraft `sat_idx` (its window, burn flag and trace)
# goes on while the spacecraft is inactive, but a released command stays released
# (`_released`): no scheduling pass or burn edge turns its mounted engine back on.
@inline function _refresh_burn_keeping_release!(control_model::BaseThrusterModel, p, sat_idx::Int, t)
    refresh_burn!(control_model, p, sat_idx, t)
    _released(control_model, p, sat_idx) && release_command!(control_model, t, sat_idx)
    return nothing
end

@inline function _schedule_thruster_control!(integrator, control_model::BaseThrusterModel, sat_idx::Int)
    _gnc_tick!(update_control!, control_model, integrator.u, integrator.p, integrator.t, sat_idx)
    _register_control_tstops!(integrator, control_model, sat_idx)
    _refresh_burn_keeping_release!(control_model, integrator.p, sat_idx, integrator.t)
    return nothing
end

# Planning is a separate late GNC phase: all guidance due at this time has
# published its request. Both ideal and mounted schedulers use the run clock;
# mounted outputs are applied after planning on that same tick.
function get_maneuver_planning_callbacks(num_sats::Int, args::SimulationConfiguration)
    callbacks = Any[]
    for (j, model) in pairs(args.control_model.control_effectors)
        model isa BaseThrusterModel || continue
        period = Float64(args.control_model.control_rates[j])
        tick! = integrator -> begin
            slot = _run_scope_slot(integrator.p, :control, j)
            scope = spacecraft_scope(model)
            _check_run_scope(slot, scope, num_sats)
            for sat_idx in _run_tick_positions(slot, scope, num_sats)
                _schedule_thruster_control!(integrator, model, sat_idx)
                if actuation_mode(model) === :mounted
                    _gnc_tick_in_scope!(mounted_burn_tick!, model, integrator.u, integrator.p, integrator.t, sat_idx)
                    _apply_mounted_tick!(model, integrator.p, sat_idx, integrator.t, period)
                end
            end
            return nothing
        end
        # Planning updates commands, not the integrated state. Output rows have
        # their own saving callback; avoid two extra states at every control tick.
        push!(callbacks, _run_clock_periodic_callback(tick!, period; save_positions=(false, false)))
        # Explicit external plans can arrive off the control grid. Register
        # their boundaries after guidance without processing a new request.
        push!(callbacks, _burn_edge_callback(model, num_sats))
    end
    return callbacks
end

# A burn switches only between integration steps: the solver stops at the
# window's ends (registered as stops when the burn is set) and this callback
# switches the burn there. A discrete callback counts as a derivative
# discontinuity unless it clears the flag, so the solver evaluates the
# derivative again, with the new thrust, before the next step; the callback
# must not call `derivative_discontinuity!(integrator, false)`, which would keep
# the stale one. The callback also takes up a plan installed in the burn-plan
# buffer since the scheduler's last pass: it adds the plan's window ends as stops
# and records it, at the end of the step it was installed in.
function _burn_edge_callback(control_model::BaseThrusterModel, num_sats::Int)
    condition = (u, t, integrator) -> burn_edge(control_model, t) || _burn_window_pending(control_model, integrator.p)
    affect! = integrator -> begin
        @inbounds for sat_idx in 1:num_sats
            _burn_window_changed(control_model, integrator.p, sat_idx) &&
                _register_control_tstops!(integrator, control_model, sat_idx)
            _refresh_burn_keeping_release!(control_model, integrator.p, sat_idx, integrator.t)
        end
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions=(false, false))
end

# An atmosphere exit may restore boundary bookkeeping, but new requests wait
# for the next defined control tick; the event never calls guidance or planning.
_reschedule_event_driven!(integrator, control_model, sat_idx::Int) = nothing
function _reschedule_event_driven!(integrator, control_model::BaseThrusterModel, sat_idx::Int)
    _register_control_tstops!(integrator, control_model, sat_idx)
    _refresh_burn_keeping_release!(control_model, integrator.p, sat_idx, integrator.t)
    return nothing
end

function schedule_event_driven_thruster_controls!(integrator, sat_idx::Int)
    control_models = integrator.p.args.control_model.control_effectors
    @inbounds for j in _scope_members(integrator.p, :control, control_models, sat_idx)
        control_model = control_models[j]
        if tick_schedule(control_model) === :event_driven
            _reschedule_event_driven!(integrator, control_model, sat_idx)
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
            _register_control_tstops!(integrator, control_model, sat_idx)
            _refresh_burn_keeping_release!(control_model, integrator.p, sat_idx, integrator.t)
        end
        return nothing
    end

    init_callback = DiscreteCallback(
        (u, t, integrator) -> false,
        schedule_all!;
        initialize=(cb, u, t, integrator) -> begin
            schedule_all!(integrator)
            # A burn already on at the start: the first derivative must thrust.
            if any(control_model.burning) && applicable(DiffEqBase.derivative_discontinuity!, integrator, true)
                DiffEqBase.derivative_discontinuity!(integrator, true)
            end
        end,
        save_positions=(false, false)
    )
    return (init_callback,)
end

# The run-time declarations of the built-in control models (GNCRuntime). The
# burn scheduler keeps exact event-driven window boundaries. Its request
# processing and mounted actuation have a late periodic phase after guidance.
GNCRuntime.tick_schedule(::BaseThrusterModel) = :event_driven
function GNCRuntime.gnc_callbacks(control_model::BaseThrusterModel, num_sats, args, period_s)
    callbacks = Any[_thruster_schedule_callbacks(control_model, num_sats, args)...]
    return Tuple(callbacks)
end
GNCRuntime.register_tstops!(integrator, control_model::BaseThrusterModel, spacecraft_index) =
    _register_control_tstops!(integrator, control_model, spacecraft_index)
# The EDG controller keeps the clock of its own ticks in the EDG state, where its
# predictions read it; recording the actual start matches the periodic callback's
# phase on resumes.
function GNCRuntime.gnc_clock_started!(control_model::_simulation_model_module.AerobrakingEnergyDepletionControlModel,
                                       period_s, origin_s)
    control_model.state.control_clock_period_s = Float64(period_s)
    control_model.state.control_clock_origin_s = Float64(origin_s)
    return nothing
end
# The two model-predictive controllers tick at their control step.
GNCRuntime.gnc_period_s(control_model::_simulation_model_module.RPOMPCControlModel) = control_model.control_dt_s
GNCRuntime.gnc_period_s(control_model::_simulation_model_module.RobotArmControlEffector) = control_model.control_dt_s

# The step a control model's controller was discretized with, when it has one.
_discretization_step(control_model) = nothing
_discretization_step(control_model::_simulation_model_module.RPOMPCControlModel) =
    control_model.controller isa _simulation_model_module.RpoLQMPCController ? control_model.controller.dt_s : nothing
_discretization_step(control_model::_simulation_model_module.RobotArmControlEffector) =
    control_model.controller isa _simulation_model_module.RobotArmJointMPCController ? control_model.controller.dt_s : nothing

"""
    check_control_periods(args) -> nothing

Before a run: refuse it with an `ArgumentError` when a periodic control model's
declared period (`gnc_period_s`) differs from the rate of its callback, or when
its controller was discretized with a different step. Equality is exact, since
builders set both from one value. (SpaceAGORA 0.3 warned instead.)
"""
function check_control_periods(args)
    models, rates = args.control_model.control_effectors, args.control_model.control_rates
    for (i, control_model) in enumerate(models)
        tick_schedule(control_model) === :periodic || continue
        T = typeof(control_model)
        period = gnc_period_s(control_model)
        period !== nothing && Float64(period) != Float64(rates[i]) && throw(ArgumentError(
            "$(nameof(T)) declares a $(period) s period (gnc_period_s) but its callback runs every " *
            "$(rates[i]) s (control_rates[$(i)]). Set both from one value."))
        step = _discretization_step(control_model)
        period !== nothing && step !== nothing && step != Float64(period) && throw(ArgumentError(
            "$(nameof(T)) ticks every $(period) s but its controller was discretized with a $(step) s step. " *
            "Build the controller with the model's control step."))
    end
    return nothing
end

# Links of the built-in models (GNCRuntime). The RPO plan buffer is written by the
# planner (its method lives with the planner) and read by the RPO MPC controller.
# The energy-depletion state is shared by the EDG guidance and control; the link
# initializes it once per run.
GNCRuntime.gnc_links(c::_simulation_model_module.RPOMPCControlModel) = (Consumes(c.plan_buffer),)
GNCRuntime.gnc_links(g::_simulation_model_module.AerobrakingEnergyDepletionGuidanceModel) = (Produces(g.state),)
GNCRuntime.gnc_links(c::_simulation_model_module.AerobrakingEnergyDepletionControlModel) = (Consumes(c.state),)
GNCRuntime.link_initialize!(s::_simulation_model_module.AerobrakingEnergyDepletionState, producers, consumers, u, p, t) =
    _simulation_model_module.GuidanceHooks._edg_initialize_state!(s)

"""
    gnc_link_groups(args) -> Vector

The run's links, in configuration order: for each shared object that a guidance,
navigation or control model declares in `gnc_links`, the models that produce it
and those that consume it, each in configuration order. Links are grouped by
identity.
"""
function gnc_link_groups(args)
    order = Any[]
    members = IdDict()                       # link => (producers, consumers), by identity
    models = (args.guidance_model.guidance_effectors..., args.navigation_model.navigation_effectors...,
              args.control_model.control_effectors...)
    for model in models, declared in gnc_links(model)
        declared isa Union{Produces, Consumes} || throw(ArgumentError(
            "gnc_links($(nameof(typeof(model)))) must return Produces and Consumes wrappers, not $(typeof(declared))."))
        link = declared.link
        haskey(members, link) || (members[link] = (Any[], Any[]); push!(order, link))
        push!(declared isa Produces ? members[link][1] : members[link][2], model)
    end
    return [(link=link, producers=members[link][1], consumers=members[link][2]) for link in order]
end

"""
    check_gnc_links(args) -> nothing

Before a run: refuse a link whose numbers of producers and consumers fall outside
its `link_cardinality`, then run each link's `check_link`.
"""
function check_gnc_links(args)
    for group in gnc_link_groups(args)
        limits = link_cardinality(group.link)
        T = nameof(typeof(group.link))
        for (role, models, allowed) in (("produce", group.producers, limits.producers),
                                        ("consume", group.consumers, limits.consumers))
            length(models) in allowed || throw(ArgumentError(
                "$(length(models)) models $(role) one $(T); $(allowed) may: " *
                join(unique(string.(nameof.(typeof.(models)))), ", ") * "."))
        end
        check_link(group.link, group.producers, group.consumers, args)
    end
    return nothing
end

# At the run's start, before the models' own initializers: each periodic control
# effector learns the clock its callback reports at the solve's start (its configured
# rate, from the grid's origin `origin`), so an initializer that predicts before the
# solve has it (the EDG guidance family's first request).
function _report_gnc_clocks!(args, origin)
    model = args.control_model
    for i in eachindex(model.control_effectors)
        control = model.control_effectors[i]
        tick_schedule(control) === :periodic || continue
        gnc_clock_started!(control, Float64(model.control_rates[i]), Float64(origin))
    end
    return nothing
end

"""
    initialize_gnc_links!(args, u, p, t) -> nothing

At a run's start, before the models' own initializers: each link's
`link_initialize!`, once.
"""
function initialize_gnc_links!(args, u, p, t)
    for group in gnc_link_groups(args)
        link_initialize!(group.link, group.producers, group.consumers, u, p, t)
    end
    return nothing
end

# The models that produce a link `control_model` consumes: the run's links in the
# order models first declare them (`gnc_link_groups`), and within each link its
# producers in configuration order, a model already listed for an earlier link not
# repeated. Not configuration order overall: with A producing X and Y, B producing
# Y and C producing X, a consumer of X and Y gets A, C, B.
function _linked_producers(control_model, args)
    consumed = [declared.link for declared in gnc_links(control_model) if declared isa Consumes]
    isempty(consumed) && return ()
    producers = Any[]
    for group in gnc_link_groups(args), producer in group.producers
        any(link -> link === group.link, consumed) || continue
        any(model -> model === producer, producers) || push!(producers, producer)
    end
    return Tuple(producers)
end

# One periodic control tick of `control_model`, for every spacecraft in its scope,
# in three steps (see "Callback order" in docs/contributing/architecture.md):
#
# 1. the models that produce a link it consumes prepare it;
# 2. each spacecraft in the model's scope gets the model's tick and then its
#    mounted actuators' update, on several threads when the parallel policy
#    decides so;
# 3. the model's solver stops are registered, one spacecraft at a time.
#
# The scope is read after step 1, and its positions (`_tick_positions`) serve steps
# 2 and 3; the thread decision counts them. A model that consumes a link also has
# its scope read before step 1, to decide whether its producers prepare it
# (`_ticks_any`). That is two reads per tick for a linked consumer and one for any
# other model. Each read is compared with the scope the run took when it started
# (`slot`, the model's entry in the run's GNC scopes; `_check_run_scope`), so a
# producer that changes its consumer's scope stops the run at the second read
# (SpaceAGORA 0.6 followed the change). The parallel policy's observation of the
# tick covers steps 2 and 3.
function _periodic_control_tick!(integrator, control_model, producers, num_sats::Int, period::Float64, slot=nothing)
    _prepare_linked_producers!(integrator, control_model, producers, num_sats, slot)
    scope = spacecraft_scope(control_model)
    _check_run_scope(slot, scope, num_sats)
    positions = _run_tick_positions(slot, scope, num_sats)
    decision = _control_callback_thread_decision(integrator.p, control_model, length(positions))
    started_ns = time_ns()
    _apply_control_ticks!(integrator, control_model, positions, period, decision)
    _register_control_stops!(integrator, control_model, positions)
    decision.policy_applied && _record_control_tick_policy!(integrator, decision, length(positions), started_ns)
    return nothing
end

# Step 1. The models that write what this one reads (its links) prepare it first,
# each once per tick, in the order `_linked_producers` lists them, on this thread:
# the RPO planner checks its reference. Not when the tick releases or skips every
# spacecraft (`_ticks_any`).
function _prepare_linked_producers!(integrator, control_model, producers, num_sats::Int, slot=nothing)
    (isempty(producers) || !_ticks_any(control_model, integrator.p, num_sats, slot)) && return nothing
    for producer in producers
        SimulationLifecycle.before_reference_control!(producer, control_model, integrator.u, integrator.p,
            integrator.t)
    end
    return nothing
end

# Step 2, for spacecraft `sat_idx` in the model's scope: the model's tick, then,
# in `:mounted` actuation, the update of the mounted actuators it drives. An
# inactive spacecraft gets the model's `release_command!` in place of the tick
# (`_gnc_tick_in_scope!`).
@inline function _apply_control_and_mounted!(integrator, control_model, sat_idx::Int, period::Float64)
    _gnc_tick_in_scope!(update_control!, control_model, integrator.u, integrator.p, integrator.t, sat_idx)
    _apply_mounted_tick!(control_model, integrator.p, sat_idx, integrator.t, period)
    return nothing
end

# Step 2, for every spacecraft at `positions` (the model's scope). The spacecraft
# index is passed explicitly, apart from the effector's index.
function _apply_control_ticks!(integrator, control_model, positions, period::Float64, decision)
    if decision.use_threads
        ParallelPolicy.threaded_foreach_persistent(:control_callback, length(positions), decision.allotment) do k
            @inbounds _apply_control_and_mounted!(integrator, control_model, positions[k], period)
        end
    else
        @inbounds for sat_idx in positions
            _apply_control_and_mounted!(integrator, control_model, sat_idx, period)
        end
    end
    return nothing
end

# Step 3. The solver's stop list belongs to the integrator that all spacecraft
# share, so stops are registered one spacecraft at a time, after the ticks, even
# when the ticks ran on several threads: for each spacecraft at `positions`.
function _register_control_stops!(integrator, control_model, positions)
    @inbounds for sat_idx in positions
        register_tstops!(integrator, control_model, sat_idx)
    end
    return nothing
end

function _record_control_tick_policy!(integrator, decision, num_items::Int, started_ns)
    ParallelPolicy.record_policy_observation!(
        :control_callback;
        mode=decision.mode,
        num_items=num_items,
        use_threads=decision.use_threads,
        elapsed_ns=(time_ns() - started_ns),
        env=_policy_env_config(integrator.p),
        ctx=ParallelPolicy.policy_context_hint(integrator.p)
    )
    return nothing
end

# A periodic control model's callback: its tick every `control_rate` seconds on
# the run's tick grid (`_run_clock_periodic_callback`). At each solve's start the
# model learns its clock's period and origin: the solve's start, or the run's in
# a checkpointed or resumed run.
function _periodic_control_callback(control_model, control_rate, num_sats::Int, args, j::Int)
    period = Float64(control_rate)
    producers = _linked_producers(control_model, args)
    tick! = integrator -> _periodic_control_tick!(integrator, control_model, producers, num_sats, period,
        _run_scope_slot(integrator.p, :control, j))
    start_clock = (origin, integrator) -> gnc_clock_started!(control_model, period, origin)
    return _run_clock_periodic_callback(tick!, period; on_clock=start_clock)
end

# The control models' callbacks, in configuration order: an event-driven model's
# own callbacks (`gnc_callbacks`), added as they are, or a periodic model's tick at
# its rate in `control_rates`. With event-driven models, one more callback, last,
# checks their scopes once per accepted step (`_event_driven_scope_check`).
function get_control_callbacks(num_sats::Int, args::SimulationConfiguration)
    control_models = args.control_model.control_effectors
    control_rates = args.control_model.control_rates
    callbacks = Any[]
    event_driven = Int[]
    for i in eachindex(control_models)
        control_model = control_models[i]
        control_rate = control_rates[i]
        schedule_kind = tick_schedule(control_model)
        schedule_kind in (:periodic, :event_driven) || throw(ArgumentError(
            "tick_schedule($(nameof(typeof(control_model)))) must return :periodic or :event_driven, not $(repr(schedule_kind))."))
        if schedule_kind === :event_driven
            append!(callbacks, gnc_callbacks(control_model, num_sats, args, control_rate))
            push!(event_driven, i)
        else
            push!(callbacks, _periodic_control_callback(control_model, control_rate, num_sats, args, i))
        end
    end
    isempty(event_driven) || push!(callbacks, _event_driven_scope_check(event_driven, num_sats))
    return callbacks
end

# An event-driven control model has no periodic tick to compare its scope at, and
# its own callbacks may read or change the scope at any step. So once per accepted
# step, after its callbacks (this callback comes after them), each event-driven
# model's scope is compared with the one the run took when it started
# (`_check_run_scope`). The condition does the comparison and never fires, so the
# callback changes neither the state nor the saved points, and its initialization
# leaves the solver's flags alone.
function _event_driven_scope_check(models::Vector{Int}, num_sats::Int)
    condition = (u, t, integrator) -> begin
        _check_event_driven_scopes(integrator.p, models, num_sats)
        false
    end
    return DiscreteCallback(condition, integrator -> nothing; save_positions=(false, false),
        initialize=(c, u, t, integrator) -> nothing)
end

function _check_event_driven_scopes(p, models::Vector{Int}, num_sats::Int)
    _run_scopes(p) === nothing && return nothing
    control_models = p.args.control_model.control_effectors
    for j in models
        _check_run_scope(_run_scope_slot(p, :control, j), spacecraft_scope(control_models[j]), num_sats)
    end
    return nothing
end
