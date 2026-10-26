# Control moment gyro gimbal stops (see `ControlMomentGyro`): an event at each
# finite stop puts the gimbal angle exactly on the stop when the integration
# reaches it, so the accepted state never passes a stop. The attitude equation
# then holds the gimbal there against a command that would drive it further
# (`cmg_angle_saturate`); a command back into the range moves it off again.

# The gimbal at position `j` of the `cmg_angles` state: a CMG's (outer) gimbal,
# then the double-gimbal CMGs' inner gimbals.
@inline function _cmg_gimbal(assembly, j::Int)
    n = length(assembly.cmgs)
    return j <= n ? assembly.cmgs[j] : assembly.cmg_inner[j - n]
end

"""
    get_cmg_callbacks(num_sats, args) -> Tuple

The gimbal-stop event of the run's control moment gyros with a finite stop: a
`VectorContinuousCallback` with two conditions per such CMG, its distance to
each stop (positive inside the range), which at a crossing onto a stop sets
the angle to the stop and marks the derivative discontinuous. Empty for runs
without such CMGs.
"""
function get_cmg_callbacks(num_sats::Int, args::SimulationConfiguration)
    stopped = Tuple{Int, Int}[]
    for (i, sc) in enumerate(args.dynamics_model.spacecraft)
        assembly = sc.root.rw_assembly
        for j in 1:(length(assembly.cmgs) + length(assembly.cmg_inner))
            g = _cmg_gimbal(assembly, j)
            (isfinite(g.min_gimbal_angle_rad) || isfinite(g.max_gimbal_angle_rad)) && push!(stopped, (i, j))
        end
    end
    isempty(stopped) && return ()
    function condition!(out, u, t, integrator)
        spacecraft = integrator.p.args.dynamics_model.spacecraft
        @inbounds for (k, (i, j)) in enumerate(stopped)
            cmg = _cmg_gimbal(spacecraft[i].root.rw_assembly, j)
            angle = u.sc[i].cmg_angles[j]
            out[2k - 1] = isfinite(cmg.max_gimbal_angle_rad) ? cmg.max_gimbal_angle_rad - angle : 1.0
            out[2k] = isfinite(cmg.min_gimbal_angle_rad) ? angle - cmg.min_gimbal_angle_rad : 1.0
        end
        return nothing
    end
    # `events[e]` is -1 when condition e crossed from positive to negative: the
    # gimbal reached that stop from inside its range.
    function stop!(integrator, events)
        spacecraft = integrator.p.args.dynamics_model.spacecraft
        changed = false
        for e in eachindex(events)
            events[e] < 0 || continue
            i, j = stopped[(e + 1) ÷ 2]
            cmg = _cmg_gimbal(spacecraft[i].root.rw_assembly, j)
            integrator.u.sc[i].cmg_angles[j] = isodd(e) ? cmg.max_gimbal_angle_rad : cmg.min_gimbal_angle_rad
            changed = true
        end
        _mark_discontinuity!(integrator, changed)
        return nothing
    end
    return (VectorContinuousCallback(condition!, stop!, 2 * length(stopped); save_positions=(false, false)),)
end
