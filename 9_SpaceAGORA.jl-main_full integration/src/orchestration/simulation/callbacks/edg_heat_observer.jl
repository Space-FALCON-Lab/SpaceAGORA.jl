"""
    EDGGeodeticBoundaryCondition(condition, spacecraft_index, altitude_m)

Register a scalar continuous condition as the declared spacecraft's geodetic
EDG entry-interface surface. `condition` is called unchanged and must describe
that surface (a positive scaling is allowed). This declaration lets the observer
capture pre-effect heat at the original solver's event time, including a left
representable root. It does not add a root or alter callback direction/options.
Index and EI altitude are checked at initialization; event geometry must agree
within frame-transform roundoff. Generic callbacks need no registration to keep
their original behavior, but receive only raw sign-bracketed EDG crossings.
"""
struct EDGGeodeticBoundaryCondition{F}
    condition::F
    spacecraft_index::Int
    altitude_m::Float64
end
(c::EDGGeodeticBoundaryCondition)(u,t,i) = c.condition(u,t,i)

function _edg_observer_distance(u, p, t, idx)
    value = _edg_heat_boundary_distance(u, p, Float64(t), idx)
    isfinite(value) || throw(ArgumentError("EDG heat boundary distance must be finite."))
    return value
end
function _edg_observed_boundary_side(u, p, t, idx)
    d = _edg_observer_distance(u, p, t, idx)
    v = d == 0.0 ? _edg_heat_boundary_vertical_speed(u, p, Float64(t), idx) : 0.0
    return (inside=d == 0.0 ? v <= 0.0 : d < 0.0,
            at_boundary=d == 0.0, vertical_speed=v)
end
function _edg_validate_boundary_condition(c::EDGGeodeticBoundaryCondition, p, num_sats)
    1 <= c.spacecraft_index <= num_sats || throw(ArgumentError("EDG boundary spacecraft index is out of range."))
    c.altitude_m == 1e3 * Float64(p.args.environment_model.EI) ||
        throw(ArgumentError("Registered EDG geodetic boundary must equal the configured EI altitude in metres."))
    return nothing
end
function _edg_registered_boundary_side(c::EDGGeodeticBoundaryCondition, u, p, t)
    idx = c.spacecraft_index
    distance = _edg_observer_distance(u, p, t, idx)
    position = CallbackServices._state_position_ii(u, idx)
    roundoff = 64eps(Float64) * max(1.0, norm(position), abs(p.args.environment_model.planet.Rp_e))
    abs(distance) <= roundoff || throw(ArgumentError("Registered callback event is not on the declared EDG geodetic boundary."))
    v = _edg_heat_boundary_vertical_speed(u, p, Float64(t), idx)
    return v == 0.0 ? nothing : v < 0.0
end
function _edg_reconcile_observed_heat!(u, integrator, root_t, observed_t, idx, inside)
    _edg_reconcile_heat_boundary!(integrator.p.args, u, idx, inside;
        event=true, crossing_time_s=Float64(root_t), observed_at_s=Float64(observed_t))
    return nothing
end

# Observe passage heat without adding a continuous root that can consume an
# existing drag, termination, or user callback. The original callback set still
# owns root location and simultaneous-effect priority.
function _edg_observe_heat_callbacks(callbacks::CallbackSet, num_sats::Int)
    processed_t = Ref(-Inf)
    initializing = Ref(true)
    remaining_initializers = Ref(0)
    # A declared event can be a left-representable root. Preserve its identity
    # across later effects at the exact same time and unchanged position only.
    declared_events = Any[nothing for _ in 1:num_sats]
    function effect_side!(u, p, t, idx)
        event = declared_events[idx]
        if event !== nothing
            position = CallbackServices._state_position_ii(u, idx)
            if Float64(t) == event.time && position == event.position
                v = _edg_heat_boundary_vertical_speed(u, p, Float64(t), idx)
                inside = v == 0.0 ? event.inside : v < 0.0
                declared_events[idx] = merge(event, (; inside))
                return inside
            end
            declared_events[idx] = nothing
        end
        return _edg_observed_boundary_side(u, p, t, idx).inside
    end
    function record_initial_side!(u, t, integrator)
        for idx in 1:num_sats
            initial = _edg_observed_boundary_side(u, integrator.p, t, idx)
            # Preserve the reviewed startup convention for a supplied EI state:
            # transform-level altitude residue uses normal velocity to choose
            # its initial side. This never captures heat or resets a plan, and
            # the band is not used for subsequent generic event crossings.
            position = CallbackServices._state_position_ii(u, idx)
            roundoff = 64eps(Float64) * max(1.0, norm(position), abs(integrator.p.args.environment_model.planet.Rp_e))
            distance = _edg_observer_distance(u, integrator.p, t, idx)
            initial_inside = abs(distance) <= roundoff ?
                _edg_heat_boundary_vertical_speed(u, integrator.p, Float64(t), idx) <= 0.0 : initial.inside
            for state in _edg_heat_states(integrator.p.args)
                # This records the initial side, without inventing an entry or
                # altering a guidance plan already constructed at solve start.
                state.heat_boundary_inside[idx] = initial_inside
            end
        end
        return nothing
    end
    function initialize_observer!(u,t,integrator)
        processed_t[] = Float64(t)
        fill!(declared_events, nothing)
        initializing[] = true
        remaining_initializers[] = length(callbacks.continuous_callbacks) + length(callbacks.discrete_callbacks)
        record_initial_side!(u,t,integrator)
        remaining_initializers[] == 0 && (initializing[] = false)
        return nothing
    end
    function observe!(integrator)
        now = Float64(integrator.t)
        now <= processed_t[] && return nothing
        for idx in 1:num_sats
            event = declared_events[idx]
            event !== nothing && event.time != now && (declared_events[idx] = nothing)
        end
        lo = max(Float64(integrator.tprev), processed_t[])
        if !isfinite(lo) || now <= lo
            initialize_observer!(integrator.u, now, integrator)
            return nothing
        end
        # Ten subintervals are a bounded dense-output scan, not a guarantee for
        # arbitrarily short excursions between nodes. Resolve every detected
        # crossing, including an entry and exit within the same accepted step.
        times = range(lo, now; length=11)
        for idx in 1:num_sats
            left = first(times)
            left_u = integrator(left)
            # Already processed effects are reconciled by their wrappers.
            # Re-reading the raw left side would undo a registered left root
            # or replay the pre-mutation interpolation as a second transition.
            fl = _edg_observer_distance(left_u, integrator.p, left, idx)
            for right in Iterators.drop(times, 1)
                right_u = right == now ? integrator.u : integrator(right)
                fr = _edg_observer_distance(right_u, integrator.p, right, idx)
                crossed = (fl < 0.0 && fr >= 0.0) || (fl > 0.0 && fr <= 0.0) ||
                    (fl == 0.0 && fr != 0.0)
                if crossed
                    l, r, fleft = left, right, fl
                    root = if fl == 0.0
                        left
                    elseif fr == 0.0
                        right
                    else
                        for _ in 1:64
                            mid = l + (r-l)/2
                            (mid == l || mid == r) && break
                            fm = _edg_observer_distance(integrator(mid), integrator.p, mid, idx)
                            if fm == 0.0
                                l = r = mid
                                break
                            elseif signbit(fm) == signbit(fleft)
                                l, fleft = mid, fm
                            else
                                r = mid
                            end
                        end
                        l + (r-l)/2
                    end
                    root_u = root == now ? integrator.u : integrator(root)
                    root_side = _edg_observed_boundary_side(root_u, integrator.p, root, idx)
                    # An exact tangent is not a crossing. A later nonzero
                    # sampled sign resolves whether a boundary was traversed.
                    if fr != 0.0 || root_side.vertical_speed != 0.0
                        inside = fr == 0.0 ? root_side.inside : fr < 0.0
                        _edg_reconcile_observed_heat!(root_u, integrator, root, now, idx, inside)
                    end
                end
                left, fl = right, fr
            end
        end
        processed_t[] = now
        return nothing
    end
    function observe_mutated_side!(integrator, before_sides)
        if initializing[]
            record_initial_side!(integrator.u, integrator.t, integrator)
            return nothing
        end
        # The original effect may terminate or mutate u. Do not rescan its old
        # interpolation. Only reconcile a zero-time change of side using the
        # post-effect state; the original effect owns derivative invalidation.
        for idx in 1:num_sats
            inside = effect_side!(integrator.u, integrator.p, integrator.t, idx)
            # A no-op or a velocity-only mutation away from a declared/zero
            # boundary must not invent a side change.
            inside == before_sides[idx] && continue
            _edg_reconcile_observed_heat!(integrator.u, integrator, integrator.t, integrator.t, idx, inside)
        end
        return nothing
    end
    wrap_effect(::Nothing, marker=nothing) = nothing
    function wrap_effect(effect!, marker=nothing)
        return function (integrator, args...)
            observe!(integrator)
            tagged_side = marker === nothing || initializing[] ? nothing :
                _edg_registered_boundary_side(marker, integrator.u, integrator.p, integrator.t)
            if tagged_side !== nothing
                _edg_reconcile_observed_heat!(integrator.u, integrator, integrator.t,
                    integrator.t, marker.spacecraft_index, tagged_side)
                declared_events[marker.spacecraft_index] = (time=Float64(integrator.t),
                    position=copy(CallbackServices._state_position_ii(integrator.u, marker.spacecraft_index)),
                    inside=tagged_side)
            end
            before_sides = [effect_side!(integrator.u, integrator.p, integrator.t, idx)
                for idx in 1:num_sats]
            result = effect!(integrator, args...)
            observe_mutated_side!(integrator, before_sides)
            return result
        end
    end
    function wrap_initialize(callback, first_callback)
        return function (cb,u,t,i)
            first_callback && initialize_observer!(u,t,i)
            callback.condition isa EDGGeodeticBoundaryCondition &&
                _edg_validate_boundary_condition(callback.condition, i.p, num_sats)
            # Lifecycle hooks belong to the original callback. SavingCallback and
            # user-defined hooks inspect or mutate fields of their own affect
            # object; the observer's effect closure is not that object.
            result = callback.initialize(callback,u,t,i)
            # Initializers can mutate u or invoke their callback's effect.
            # All such changes still define the initial side, not a passage.
            record_initial_side!(i.u,t,i)
            remaining_initializers[] -= 1
            remaining_initializers[] == 0 && (initializing[] = false)
            return result
        end
    end
    wrap_finalize(callback) = (cb,u,t,i) -> callback.finalize(callback,u,t,i)
    function wrap_callback(callback::ContinuousCallback, first_callback)
        options = (; (name => getfield(callback, name) for name in fieldnames(typeof(callback))
            if name ∉ (:condition, :affect!, :affect_neg!, :initialize, :finalize))...)
        initialize = wrap_initialize(callback, first_callback)
        marker = callback.condition isa EDGGeodeticBoundaryCondition ? callback.condition : nothing
        return ContinuousCallback(callback.condition, wrap_effect(callback.affect!, marker),
            wrap_effect(callback.affect_neg!, marker); initialize=initialize, finalize=wrap_finalize(callback), options...)
    end
    function wrap_callback(callback::VectorContinuousCallback, first_callback)
        options = (; (name => getfield(callback, name) for name in fieldnames(typeof(callback))
            if name ∉ (:condition, :affect!, :len, :initialize, :finalize))...)
        initialize = wrap_initialize(callback, first_callback)
        return VectorContinuousCallback(callback.condition, wrap_effect(callback.affect!),
            callback.len; initialize=initialize, finalize=wrap_finalize(callback), options...)
    end
    function wrap_callback(callback::DiscreteCallback, first_callback)
        options = (; (name => getfield(callback, name) for name in fieldnames(typeof(callback))
            if name ∉ (:condition, :affect!, :initialize, :finalize))...)
        return DiscreteCallback(callback.condition, wrap_effect(callback.affect!);
            initialize=wrap_initialize(callback, first_callback), finalize=wrap_finalize(callback), options...)
    end
    continuous = map(enumerate(callbacks.continuous_callbacks)) do (idx, callback)
        wrap_callback(callback, idx == 1)
    end
    initialize = isempty(continuous) ? ((cb,u,t,i) -> initialize_observer!(u,t,i)) :
        ((cb,u,t,i) -> nothing)
    function observe_without_derivative_change!(integrator)
        observe!(integrator)
        # This observer changes only EDG bookkeeping, not u or force inputs.
        # Pre/post-effect wrappers never clear the original effect's decision.
        DiffEqBase.derivative_discontinuity!(integrator, false)
        return nothing
    end
    discrete = map(cb -> wrap_callback(cb, false), callbacks.discrete_callbacks)
    observer = DiscreteCallback((u,t,i) -> true, observe_without_derivative_change!;
        initialize=initialize, save_positions=(false,false))
    return CallbackSet(continuous..., observer, discrete...)
end
