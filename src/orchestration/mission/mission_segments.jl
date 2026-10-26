# Operational segments run inside one solver invocation. The callbacks own only
# mission sequencing; physical state and model lifecycle remain engine-owned.

"""
    MissionDeadline(time_s)

End a segment at an absolute elapsed mission time in seconds. The deadline
must be finite, positive and later than the time this segment is entered.
`MissionTime` remains the existing mission-duration enum value.
"""
struct MissionDeadline
    time_s::Float64
    function MissionDeadline(time_s::Real)
        value = Float64(time_s)
        !(time_s isa Bool) && isfinite(value) && value > 0 || throw(ArgumentError(
            "MissionDeadline needs a finite positive elapsed time."))
        new(value)
    end
end

"""
    MissionDuration(duration_s)

End a segment this many seconds after its actual entry time. The duration must
be finite and positive. The solver is stopped exactly at the resulting time.
"""
struct MissionDuration
    duration_s::Float64
    function MissionDuration(duration_s::Real)
        value = Float64(duration_s)
        !(duration_s isa Bool) && isfinite(value) && value > 0 || throw(ArgumentError(
            "MissionDuration needs a finite positive duration."))
        new(value)
    end
end

"""
    MissionCondition(f; direction=:rising)

End a segment at a zero crossing of the finite real scalar `f(context)`.
`:rising` crosses from negative to positive; `:falling` crosses from positive to
negative. At segment entry the value must be strictly on the starting side:
already-satisfied conditions, including zero, are rejected. Conditions must be
pure and must not mutate configuration, controller state or external objects.

The context contains `time_s`, the isolated `configuration`, immutable numerical
`spacecraft` samples (`pos_ii`, `vel_ii`, `mass_kg`, `q_ib`, `ω_body`), and a copy
of the full solver `state`. Neither numerical view aliases live physical state.
"""
struct MissionCondition{F}
    f::F
    direction::Symbol
    function MissionCondition(f; direction::Symbol=:rising)
        direction in (:rising, :falling) || throw(ArgumentError(
            "MissionCondition direction must be :rising or :falling."))
        new{typeof(f)}(f, direction)
    end
end

_no_segment_action(context) = nothing
_default_segment_actions(configuration) = (; on_enter=_no_segment_action, on_exit=_no_segment_action)

"""
    MissionSegment(name; until=nothing, build_actions=configuration -> (;))

An operational part of one continuous mission. `until` is a `MissionDeadline`,
`MissionDuration` or `MissionCondition`; only the final segment may omit it.
`build_actions` is called once with the isolated run configuration and returns a
named tuple with optional `on_enter` and `on_exit` functions. Each function takes
the same context as `MissionCondition`. Capture the run's controller, reference
or actuator handles in the factory and change their commands in these actions.
The factory itself must avoid mutation and I/O.

Actions do not receive an integrator or writable physical state. They must not
change vehicle topology, physical models, solver settings or ownership. They do
not automatically disable controllers: an action that starts a held command
must also stop or replace it explicitly. Initializers are not rerun at boundaries.
Supported solver modes are `:tsit5`, `:auto_stiff`, `:rodas5p` and `:dp8`.
External propagators and checkpoint workflows are not supported yet.
On horizon or earlier simulation termination, the active segment's exit action
runs once at the actual final time. Existing `MissionPhase` windows remain
independent observation metadata.
"""
struct MissionSegment{U, F}
    name::Symbol
    until::U
    build_actions::F
    function MissionSegment(name::Union{Symbol, AbstractString}; until=nothing,
        build_actions=_default_segment_actions)
        until isa Union{Nothing, MissionDeadline, MissionDuration, MissionCondition} ||
            throw(ArgumentError("MissionSegment.until must be a MissionDeadline, MissionDuration, MissionCondition or nothing."))
        new{typeof(until), typeof(build_actions)}(Symbol(_name(name, "MissionSegment")), until, build_actions)
    end
end

_empty_segment_history() = (; segments=NamedTuple[], transitions=NamedTuple[])

mutable struct _SegmentRuntime{S, A, C}
    segments::S
    actions::A
    configuration::C
    active::Int
    deadline::Float64
    history::typeof(_empty_segment_history())
end

function _segment_context(u, t, configuration)
    samples = map(eachindex(configuration.dynamics_model.spacecraft)) do i
        sample = SimulationEngine.build_state_sample(u.sc[i], nothing,
            configuration.mission_configuration.orientation_sim, i)
        sample
    end
    return (; time_s=Float64(t), configuration, spacecraft=Tuple(samples), state=deepcopy(u))
end

function _condition_value(condition::MissionCondition, context)
    applicable(condition.f, context) || throw(ArgumentError("MissionCondition function must accept one mission context."))
    value = condition.f(context)
    value isa Real && !(value isa Bool) || throw(ArgumentError("MissionCondition must return a finite real scalar, not $(typeof(value))."))
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("MissionCondition must return a finite real scalar; got $result."))
    return result
end

function _segment_actions(segment, configuration)
    applicable(segment.build_actions, configuration) || throw(ArgumentError(
        "MissionSegment $(segment.name) build_actions must accept the isolated configuration."))
    actions = segment.build_actions(configuration)
    actions isa NamedTuple || throw(ArgumentError("MissionSegment $(segment.name) build_actions must return a named tuple."))
    all(key -> key in (:on_enter, :on_exit), keys(actions)) || throw(ArgumentError(
        "MissionSegment $(segment.name) actions may contain only on_enter and on_exit."))
    return merge((; on_enter=_no_segment_action, on_exit=_no_segment_action), actions)
end

function _call_segment_action(action, context, name, kind)
    applicable(action, context) || throw(ArgumentError("MissionSegment $name $kind must accept one mission context."))
    action(context)
    return nothing
end

function _enter_segment!(runtime, index, integrator)
    segment = runtime.segments[index]
    t = Float64(integrator.t)
    horizon = runtime.configuration.mission_configuration.mission_time
    context = _segment_context(integrator.u, t, runtime.configuration)
    exit = segment.until
    deadline = exit isa MissionDeadline ? exit.time_s : exit isa MissionDuration ? t + exit.duration_s : Inf
    # Relative durations can round one ulp beyond the horizon (0.1 + 0.2).
    # Snap only that arithmetic boundary, never an explicit absolute deadline.
    if exit isa MissionDuration && isfinite(deadline) && horizon < deadline <= nextfloat(Float64(horizon))
        deadline = Float64(horizon)
    end
    if exit isa Union{MissionDeadline, MissionDuration}
        isfinite(deadline) && t < deadline <= horizon || throw(ArgumentError(
            "MissionSegment $(segment.name) deadline must be later than entry at $t s and no later than mission_time=$horizon s; got $deadline s."))
    end
    runtime.active = index
    runtime.deadline = deadline
    push!(runtime.history.segments, (; name=String(segment.name), start_s=t, end_s=t))
    _call_segment_action(runtime.actions[index].on_enter, context, segment.name, :on_enter)
    if exit isa MissionCondition
        value = _condition_value(exit, _segment_context(integrator.u, t, runtime.configuration))
        starting_side = exit.direction == :rising ? value < 0 : value > 0
        starting_side || throw(ArgumentError(
            "MissionSegment $(segment.name) condition is already satisfied at entry; its $(exit.direction) crossing must lie ahead."))
    end
    if isfinite(deadline)
        SciMLBase.add_tstop!(integrator, deadline)
    end
    return nothing
end

function _close_segment!(runtime, t)
    index = length(runtime.history.segments)
    previous = runtime.history.segments[index]
    runtime.history.segments[index] = (; name=previous.name, start_s=previous.start_s, end_s=Float64(t))
    return nothing
end

function _advance_segment!(runtime, integrator, reason)
    index = runtime.active
    index == 0 && return nothing
    segment = runtime.segments[index]
    t = Float64(integrator.t)
    _call_segment_action(runtime.actions[index].on_exit,
        _segment_context(integrator.u, t, runtime.configuration), segment.name, :on_exit)
    _close_segment!(runtime, t)
    next_index = index + 1
    next_name = next_index <= length(runtime.segments) ? String(runtime.segments[next_index].name) : nothing
    push!(runtime.history.transitions, (; time_s=t, from=String(segment.name), to=next_name, reason))
    if next_name === nothing
        runtime.active = 0
        runtime.deadline = Inf
        SciMLBase.terminate!(integrator)
    else
        _enter_segment!(runtime, next_index, integrator)
    end
    # Actions may change the RHS through held commands without editing u. Drop
    # the cached derivative so the next step uses the newly commanded model.
    SciMLBase.derivative_discontinuity!(integrator, true)
    return nothing
end

function _finalize_segments!(runtime, u, t)
    index = runtime.active
    index == 0 && return nothing
    segment = runtime.segments[index]
    _call_segment_action(runtime.actions[index].on_exit,
        _segment_context(u, t, runtime.configuration), segment.name, :on_exit)
    _close_segment!(runtime, t)
    horizon = runtime.configuration.mission_configuration.mission_time
    reason = Float64(t) >= horizon ? "horizon" : "simulation_stop"
    push!(runtime.history.transitions, (; time_s=Float64(t), from=String(segment.name), to=nothing, reason))
    runtime.active = 0
    runtime.deadline = Inf
    return nothing
end

function _segment_callbacks(segments, configuration, runtime_ref)
    actions = map(segment -> _segment_actions(segment, configuration), segments)
    runtime = _SegmentRuntime(segments, actions, configuration, 0, Inf, _empty_segment_history())
    runtime_ref[] = runtime
    function initialize(cb, u, t, integrator)
        _enter_segment!(runtime, 1, integrator)
        SciMLBase.derivative_discontinuity!(integrator, true)
        return nothing
    end
    timed_condition(u, t, integrator) = runtime.active != 0 && isfinite(runtime.deadline) && t >= runtime.deadline
    timed_affect!(integrator) = _advance_segment!(runtime, integrator,
        runtime.segments[runtime.active].until isa MissionDeadline ? "deadline" : "duration")
    finalize(cb, u, t, integrator) = _finalize_segments!(runtime, u, t)
    timed = SciMLBase.DiscreteCallback(timed_condition, timed_affect!;
        initialize=initialize, finalize=finalize, save_positions=(false, false))
    continuous = map(eachindex(segments)) do index
        exit = segments[index].until
        exit isa MissionCondition || return nothing
        inactive_value = exit.direction == :rising ? -1.0 : 1.0
        condition(u, t, integrator) = runtime.active == index ?
            _condition_value(exit, _segment_context(u, t, configuration)) : inactive_value
        affect!(integrator) = runtime.active == index ? _advance_segment!(runtime, integrator, "condition") : nothing
        positive = exit.direction == :rising ? affect! : nothing
        negative = exit.direction == :falling ? affect! : nothing
        SciMLBase.ContinuousCallback(condition, positive, negative; save_positions=(false, false))
    end
    return SciMLBase.CallbackSet(filter(!isnothing, continuous)..., timed)
end

function _validate_segments(segments, configuration)
    isempty(segments) && return nothing
    isempty(configuration.external_propagators) || throw(ArgumentError(
        "Operational mission segments do not yet support external propagators; keep ownership fixed in an unsegmented run."))
    solver = configuration.solver_config === nothing ? SimulationEngine._solver_config_from_env() : configuration.solver_config
    solver.solver_mode in (:tsit5, :auto_stiff, :rodas5p, :dp8) || throw(ArgumentError(
        "Operational mission segments require :tsit5, :auto_stiff, :rodas5p or :dp8; solver mode $(solver.solver_mode) is not supported yet."))
    horizon = configuration.mission_configuration.mission_time
    for segment in segments
        exit = segment.until
        exit isa MissionDeadline && exit.time_s > horizon && throw(ArgumentError(
            "MissionSegment $(segment.name) deadline exceeds mission_time=$horizon s."))
        exit isa MissionDuration && exit.duration_s > horizon && throw(ArgumentError(
            "MissionSegment $(segment.name) duration exceeds mission_time=$horizon s."))
    end
    return nothing
end

function _write_mission_history(prefix, name, history)
    path = prefix * "_mission.json"
    SimulationModel.IOSerialization._atomic_write_file(path, tmp -> open(tmp, "w") do io
        JSON.print(io, (; schema=1, name, segments=history.segments, transitions=history.transitions), 2)
        println(io)
    end)
    return path
end
