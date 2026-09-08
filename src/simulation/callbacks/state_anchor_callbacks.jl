# Scheduled state anchors.
#
# At given elapsed times a satellite's inertial position and velocity are
# replaced by supplied values: a reconstructed trajectory, a navigation
# solution, a checkpoint from another run. Between anchors the propagation is
# the engine's own. The telemetry regression uses this to re-anchor the Odyssey
# aerobraking replay to the NAV reconstruction after each of the mission's trim
# burns, so that a 340-orbit comparison is a sequence of short open-loop
# segments rather than one chaotic arc (see the manifest's state_anchors block).
#
# Each anchor is a tstop, so the solver lands on it exactly; the affect
# overwrites the state and reports how far the propagated state was from the
# anchor, which is the open-loop drift of the segment just ended.

"""
    StateAnchor(elapsed_s, sat_idx, state; orbit_count=nothing)

A scheduled state anchor: at `elapsed_s` seconds after the initial time,
satellite `sat_idx`'s inertial position and velocity are replaced by `state`
(six components: position in m, velocity in m/s, in the run's planet-centred
inertial axes). With `orbit_count` given, the satellite's orbit counter is set
to it as well: the value the counter holds between the apoapsis the anchored
trajectory has just passed and the next one (the counter starts at 1 and
increments at every apoapsis). That keeps orbit-keyed logic such as the
campaign burn replay aligned with the anchored trajectory when the propagated
one had drifted across an apoapsis crossing. Build the callback with
[`get_state_anchor_callback`](@ref).
"""
struct StateAnchor
    elapsed_s::Float64
    sat_idx::Int
    pos_ii_m::SVector{3, Float64}
    vel_ii_mps::SVector{3, Float64}
    orbit_count::Union{Nothing, Int}
end

function StateAnchor(elapsed_s::Real, sat_idx::Integer, state::AbstractVector{<:Real}; orbit_count::Union{Nothing, Integer}=nothing)
    length(state) == 6 || throw(ArgumentError("StateAnchor state must have six components (position m, velocity m/s); got $(length(state))."))
    isfinite(elapsed_s) || throw(ArgumentError("StateAnchor elapsed_s must be finite; got $elapsed_s."))
    all(isfinite, state) || throw(ArgumentError("StateAnchor state must be finite."))
    sat_idx >= 1 || throw(ArgumentError("StateAnchor sat_idx must be >= 1; got $sat_idx."))
    orbit_count === nothing || orbit_count >= 1 || throw(ArgumentError("StateAnchor orbit_count must be >= 1; got $orbit_count."))
    return StateAnchor(
        Float64(elapsed_s),
        Int(sat_idx),
        SVector{3, Float64}(state[1], state[2], state[3]),
        SVector{3, Float64}(state[4], state[5], state[6]),
        orbit_count === nothing ? nothing : Int(orbit_count),
    )
end

const _STATE_ANCHOR_TIME_TOL_S = 1.0e-6

# The staged density callback runs before the extra callbacks and records an
# atmosphere sample for the pre-anchor state at this same time; the buffered
# and freeze-per-step readers would hand that sample to the next RHS
# evaluation. Clearing the sample times makes every reader resample at the
# anchored position (see _buffered_atmosphere_valid).
@inline function _invalidate_environment_samples!(p, sat_idx::Int)::Nothing
    buffers = p.shared_buffers
    if hasproperty(buffers, :density_sample_t) && sat_idx <= length(buffers.density_sample_t)
        buffers.density_sample_t[sat_idx] = NaN
    end
    if hasproperty(buffers, :in_atmosphere_sample_t) && sat_idx <= length(buffers.in_atmosphere_sample_t)
        buffers.in_atmosphere_sample_t[sat_idx] = NaN
    end
    return nothing
end

"""
    get_state_anchor_callback(anchors; verbose=true) -> DiscreteCallback

Discrete callback that applies `anchors` (sorted by time) as the solve reaches
each anchor's elapsed time. Every anchor time is registered as a tstop at
initialisation, so the integrator steps onto it exactly; anchors at or before
the initial time are skipped. With `verbose` on, one line per anchor reports
the propagated-minus-anchor position and velocity differences, i.e. the
open-loop drift accumulated since the previous anchor. Supported on the
first-order solver paths; under the gravity-backbone split the callback
throws at initialisation rather than truncating the solve.
"""
function get_state_anchor_callback(anchors::AbstractVector{StateAnchor}; verbose::Bool=true)
    isempty(anchors) && throw(ArgumentError("get_state_anchor_callback needs at least one anchor."))
    sorted = sort(collect(anchors); by=a -> a.elapsed_s)
    times = Float64[a.elapsed_s for a in sorted]
    next = Ref(1)

    function condition(u, t, integrator)
        k = next[]
        k <= length(times) || return false
        return abs(t - times[k]) <= _STATE_ANCHOR_TIME_TOL_S
    end

    function affect!(integrator)
        engine = _simulation_engine_module()
        u = integrator.u
        p = integrator.p
        t = integrator.t
        # Every anchor at this time is applied in this one invocation: the
        # callback fires once per solver time, so anchors for several
        # satellites at the same instant would otherwise be skipped along
        # with everything after them.
        k = next[]
        @inbounds while k <= length(times) && abs(t - times[k]) <= _STATE_ANCHOR_TIME_TOL_S
            anchor = sorted[k]
            r_before = engine._state_position_ii(u, anchor.sat_idx)
            v_before = engine._state_velocity_ii(u, anchor.sat_idx)
            engine._set_state_position_velocity_ii!(u, anchor.sat_idx, anchor.pos_ii_m, anchor.vel_ii_mps)
            _invalidate_environment_samples!(p, anchor.sat_idx)
            count_before = anchor.sat_idx <= length(p.orbit_counter) ? p.orbit_counter[anchor.sat_idx] : 0
            if anchor.orbit_count !== nothing && anchor.sat_idx <= length(p.orbit_counter)
                p.orbit_counter[anchor.sat_idx] = anchor.orbit_count
            end
            if verbose
                dr = norm(r_before - anchor.pos_ii_m)
                dv = norm(v_before - anchor.vel_ii_mps)
                count_note = anchor.orbit_count === nothing ? "" : " orbit_count=$(anchor.orbit_count) was=$(count_before)"
                println("state_anchor sat=$(anchor.sat_idx) index=$k t_s=$t drift_pos_m=$(round(dr, digits=3)) drift_vel_mps=$(round(dv, digits=6))$(count_note)")
            end
            k += 1
        end
        next[] = k
        if applicable(DiffEqBase.u_modified!, integrator, true)
            DiffEqBase.u_modified!(integrator, true)
        end
        return nothing
    end

    function initialize(cb, u, t, integrator)
        _simulation_engine_module()._is_gravity_backbone_state(u) && throw(ArgumentError(
            "Scheduled state anchors are supported on the first-order solver paths only; the gravity-backbone split (SPACEAGORA_SOLVER_MODE=gravity_backbone_split) cannot be re-anchored."
        ))
        first_future = findfirst(tk -> tk > t + _STATE_ANCHOR_TIME_TOL_S, times)
        next[] = first_future === nothing ? length(times) + 1 : first_future
        @inbounds for k in next[]:length(times)
            _maybe_add_control_tstop!(integrator, times[k])
        end
        return nothing
    end

    return DiscreteCallback(condition, affect!; initialize=initialize, save_positions=(false, false))
end

export StateAnchor, get_state_anchor_callback
