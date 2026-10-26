# Periodic navigation, guidance and control ticks.
#
# DiffEqCallbacks' `PeriodicCallback(f, Δt)` ticks at t0 + k*Δt from the start
# t0 of each solve and never at the solve's final time. A checkpointed run
# solves consecutive segments, so under it the tick at a segment boundary on the
# grid was lost, and a boundary off the grid shifted every later tick. The
# callback below keeps one grid for the whole run instead. The engine sets the
# run clock in the shared buffers for a checkpointed or resumed run: the grid's
# origin is the run's start, and a tick fires at any grid time before the run's
# end, including a segment's final time. There it fires after the step to that
# time and before the samples saved there, as in one solve, and before the
# checkpoint written there; the next segment, or a run resumed from that
# checkpoint, starts with the following tick. A run does not tick at its own end,
# so a checkpoint written there marks the ticks due then as pending, and a run
# resumed from it takes them when its first solve starts, before the samples
# saved there (`initial_affect` of PeriodicCallback).
#
# Without the run clock (a run in one solve, or a solve outside the engine) it
# is `PeriodicCallback(f, Δt)`: the origin is the solve's start and the tick
# times, condition and stops are the same. Under the run clock the run's end
# replaces the solve's final time in the condition: a cached integrator that
# `reinit!` moves to the next segment keeps the first segment's
# `integrator.sol.prob.tspan`.
#
# Under the run clock a tick and a solve's start or final time that differ only
# by rounding are one time (`_same_run_time`). A checkpoint boundary and the grid
# time meant to fall on it can differ by an ulp: `3 * 0.1` is 0.30000000000000004,
# while a 0.3 s checkpoint interval ends the first segment at 0.3. A stop that
# close to a solve's start leaves a step of an ulp or two, which the solver can
# reject and then abort as unstable. So a tick that close to a solve's final time
# fires there, before the checkpoint written there, and adds no stop of its own;
# at a solve's start it is the tick taken at the previous solve's end, or the
# pending one. At the run's end it is pending, as a tick exactly there. Elsewhere
# in a solve, and without the run clock, tick times are compared exactly, as
# PeriodicCallback compares them.
#
# What each clock did is what the checkpoint records: the run lists its clocks in
# the shared buffers, and a checkpoint names those whose next tick, not yet taken,
# is due at its time (`SimulationEngine._periodic_tick_checkpoint_state`). That
# covers a tick at the run's end, one just before a horizon where the run stopped,
# and one a later callback's `terminate!` cut off. A run resumed from it takes those
# clocks' ticks at its start, once, and no other clock's. A checkpoint written
# before that list is read as SpaceAGORA 0.4.0 left it
# (`SimulationEngine._pending_periodic_ticks`).

const _PeriodicTickClock = NamedTuple{(:origin, :index, :horizon),
    Tuple{Base.RefValue{Float64}, Base.RefValue{Int}, Base.RefValue{Float64}}}

# Each Task running the callback keeps its own clock, as PeriodicCallback keeps
# its cache, so concurrent solves of one callback do not share a phase.
@inline _periodic_tick_clock(key::Symbol) = task_local_storage(key)::_PeriodicTickClock

@inline _periodic_tick_time(clock::_PeriodicTickClock, Δt::Float64, k::Int) = clock.origin[] + k * Δt

# Two run times within 100 units in the last place of the larger are one time:
# the tolerance OrdinaryDiffEqCore's `modify_dt_for_tstops!` gives a step that
# reaches a stop. It covers the rounding of `k * Δt` and of a checkpoint boundary,
# an ulp or two each, and stays far below a tick period: 100 ulps of a year in
# seconds is under a microsecond.
@inline _same_run_time(a::Float64, b::Float64) = abs(a - b) <= 100 * eps(max(abs(a), abs(b)))

# Tick `k` falls at time `t`: exactly, or under the run clock within `_same_run_time`.
@inline function _periodic_tick_at(clock::_PeriodicTickClock, Δt::Float64, k::Int, t)
    tick = _periodic_tick_time(clock, Δt, k)
    return tick == t || (!isnan(clock.horizon[]) && _same_run_time(tick, Float64(t)))
end

# The solve ends at this step: at its final time, the largest stop, or where an
# event or an earlier discrete callback terminated it (`terminate!` empties the stops).
@inline function _at_solve_end(integrator, t)
    isempty(DiffEqBase.get_tstops_array(integrator)) && return true
    return integrator.tdir * t == DiffEqBase.get_tstops_max(integrator)
end

@inline function _periodic_run_clock(p)
    hasproperty(p, :shared_buffers) || return (NaN, NaN, NaN)
    buffers = p.shared_buffers
    return (buffers.periodic_clock_origin_s[], buffers.periodic_clock_horizon_s[],
            buffers.periodic_clock_pending_s[])
end

# The clock's position in the run's list of clocks, adding it at its first
# initialization (0 outside the engine's parameters).
function _register_periodic_clock!(p, key::Symbol, Δt::Float64, clock::_PeriodicTickClock)
    hasproperty(p, :shared_buffers) || return 0
    clocks = p.shared_buffers.periodic_clocks
    entry = (key=key, period_s=Δt, origin=clock.origin, index=clock.index)
    for i in eachindex(clocks)
        clocks[i].key === key || continue
        clocks[i] = entry
        return i
    end
    push!(clocks, entry)
    return length(clocks)
end

# The tick at a resumed run's start that the run which wrote its checkpoint did
# not take: for the clock at `position` with period `Δt`, when the checkpoint lists
# its clocks; otherwise any clock's tick at the start from the earliest pending
# tick time on.
function _start_tick_pending(p, position::Int, Δt::Float64, tick::Float64)
    hasproperty(p, :shared_buffers) || return false
    buffers = p.shared_buffers
    clocks = buffers.periodic_clock_pending_clocks[]
    clocks === nothing && return tick >= buffers.periodic_clock_pending_from_s[]
    return (position, Δt) in clocks
end

# Stop at the next tick unless the solve ends first, as PeriodicCallback does. A
# solve terminated in this step (`terminate!` empties its stops; an event's or an
# earlier discrete callback's) has no next stop, and the tick it is taking still runs.
# Under the run clock a tick at the solve's final time (`_same_run_time`) fires at
# that stop rather than at one of its own.
@inline function _schedule_next_periodic_tick!(integrator, clock::_PeriodicTickClock, Δt::Float64)
    clock.index[] += 1
    isempty(DiffEqBase.get_tstops_array(integrator)) && return nothing
    tnew = _periodic_tick_time(clock, Δt, clock.index[])
    tend = DiffEqBase.get_tstops_max(integrator)
    isnan(clock.horizon[]) || !_same_run_time(tnew, integrator.tdir * tend) || return nothing
    if integrator.tdir * tnew < tend
        DiffEqBase.add_tstop!(integrator, tnew)
    end
    return nothing
end

# The callback's affect, shaped as DiffEqCallbacks' `PeriodicCallbackAffect`:
# `affect!` is the tick function itself, which tests call without scheduling the
# next stop.
struct _RunClockTickAffect{F}
    affect!::F
    Δt::Float64
    key::Symbol
end
function (tick::_RunClockTickAffect)(integrator)
    _schedule_next_periodic_tick!(integrator, _periodic_tick_clock(tick.key), tick.Δt)
    return tick.affect!(integrator)
end

"""
    _run_clock_periodic_callback(f, Δt; on_clock=(origin, integrator) -> nothing)

A `DiscreteCallback` applying `f(integrator)` every `Δt` seconds on the run's
tick grid (see the comment at the top of this file). `on_clock` receives the
grid's origin when the callback initializes, before the first tick is scheduled.
"""
function _run_clock_periodic_callback(f, Δt::Float64; on_clock=(origin, integrator) -> nothing,
        save_positions=(true, true))
    key = gensym(:run_clock_periodic_callback)

    condition = function (u, t, integrator)
        clock = _periodic_tick_clock(key)
        tick = _periodic_tick_time(clock, Δt, clock.index[])
        horizon = clock.horizon[]
        t == tick || (!isnan(horizon) && _same_run_time(t, tick) && _at_solve_end(integrator, t)) || return false
        hasfield(typeof(integrator), :iter) && integrator.iter == integrator.opts.maxiters && return false
        return isnan(horizon) ? t != last(integrator.sol.prob.tspan) : tick < horizon && !_same_run_time(tick, horizon)
    end

    affect! = _RunClockTickAffect(f, Δt, key)

    initialize = function (cb, u, t, integrator)
        @assert integrator.tdir == sign(Δt)
        DiffEqBase.derivative_discontinuity!(integrator, false)
        storage = task_local_storage()
        clock = get(storage, key, nothing)
        if clock === nothing
            clock = (origin=Ref(NaN), index=Ref(0), horizon=Ref(NaN))
            storage[key] = clock
        end
        clock = clock::_PeriodicTickClock
        run_origin, run_horizon, pending = _periodic_run_clock(integrator.p)
        clock.origin[] = isnan(run_origin) ? Float64(t) : run_origin
        clock.horizon[] = run_horizon
        position = _register_periodic_clock!(integrator.p, key, Δt, clock)
        # The first tick after the solve's start. A tick at the start itself
        # fired before the checkpoint the solve starts from, unless it is pending.
        k = max(1, floor(Int, (t - clock.origin[]) / Δt))
        while _periodic_tick_time(clock, Δt, k) <= t || _periodic_tick_at(clock, Δt, k, t)
            k += 1
        end
        while k > 1 && _periodic_tick_time(clock, Δt, k - 1) > t && !_periodic_tick_at(clock, Δt, k - 1, t)
            k -= 1
        end
        clock.index[] = k - 1
        on_clock(clock.origin[], integrator)
        if k > 1 && t == pending && _periodic_tick_at(clock, Δt, k - 1, t) &&
                _start_tick_pending(integrator.p, position, Δt, _periodic_tick_time(clock, Δt, k - 1))
            # The tick at the start, which the run that wrote the checkpoint
            # did not take. It changes the derivative.
            DiffEqBase.derivative_discontinuity!(integrator, true)
            affect!(integrator)
            return nothing
        end
        _schedule_next_periodic_tick!(integrator, clock, Δt)
        return nothing
    end

    return DiscreteCallback(condition, affect!; initialize=initialize, save_positions=save_positions)
end
