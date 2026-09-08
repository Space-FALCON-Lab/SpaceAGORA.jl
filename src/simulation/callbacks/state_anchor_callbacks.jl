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

struct StateAnchor
    elapsed_s::Float64
    sat_idx::Int
    pos_ii_m::SVector{3, Float64}
    vel_ii_mps::SVector{3, Float64}
end

function StateAnchor(elapsed_s::Real, sat_idx::Integer, state::AbstractVector{<:Real})
    length(state) == 6 || throw(ArgumentError("StateAnchor state must have six components (position m, velocity m/s); got $(length(state))."))
    isfinite(elapsed_s) || throw(ArgumentError("StateAnchor elapsed_s must be finite; got $elapsed_s."))
    all(isfinite, state) || throw(ArgumentError("StateAnchor state must be finite."))
    sat_idx >= 1 || throw(ArgumentError("StateAnchor sat_idx must be >= 1; got $sat_idx."))
    return StateAnchor(
        Float64(elapsed_s),
        Int(sat_idx),
        SVector{3, Float64}(state[1], state[2], state[3]),
        SVector{3, Float64}(state[4], state[5], state[6]),
    )
end

const _STATE_ANCHOR_TIME_TOL_S = 1.0e-6

"""
    get_state_anchor_callback(anchors; verbose=true) -> DiscreteCallback

Discrete callback that applies `anchors` (sorted by time) as the solve reaches
each anchor's elapsed time. Every anchor time is registered as a tstop at
initialisation, so the integrator steps onto it exactly; anchors at or before
the initial time are skipped. With `verbose` on, one line per anchor reports
the propagated-minus-anchor position and velocity differences, i.e. the
open-loop drift accumulated since the previous anchor.
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
        k = next[]
        anchor = sorted[k]
        engine = _simulation_engine_module()
        u = integrator.u
        r_before = engine._state_position_ii(u, anchor.sat_idx)
        v_before = engine._state_velocity_ii(u, anchor.sat_idx)
        engine._set_state_position_velocity_ii!(u, anchor.sat_idx, anchor.pos_ii_m, anchor.vel_ii_mps)
        if applicable(DiffEqBase.u_modified!, integrator, true)
            DiffEqBase.u_modified!(integrator, true)
        end
        next[] = k + 1
        if verbose
            dr = norm(r_before - anchor.pos_ii_m)
            dv = norm(v_before - anchor.vel_ii_mps)
            println("state_anchor sat=$(anchor.sat_idx) index=$k t_s=$(integrator.t) drift_pos_m=$(round(dr, digits=3)) drift_vel_mps=$(round(dv, digits=6))")
        end
        return nothing
    end

    function initialize(cb, u, t, integrator)
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
