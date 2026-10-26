# Actuator faults during a run (see ActuatorFaults): the faults switch exactly
# at their window ends, which are registered as solver stops, and a stuck
# actuator captures its level at onset; the fault columns.

_fault_spacecraft(args) = [i for i in eachindex(args.dynamics_model.spacecraft) if has_faults(args.dynamics_model.spacecraft[i])]

# The command the first `:mounted` effector acting on spacecraft `i` holds for it,
# or nothing.
function _fault_command_levels(p, i::Int)
    sc = p.args.dynamics_model.spacecraft[i]
    n = mounted_actuator_count(sc)
    n == 0 && return nothing
    effectors = p.args.control_model.control_effectors
    return _mounted_effector_levels(effectors, _scope_members(p, :control, effectors, i), i, n)
end

"""
    get_fault_callbacks(num_sats, args) -> Tuple

For runs with actuator faults: one callback at the faults' window ends, where
the solver stops, that switches them and declares a derivative discontinuity
when one switched;
its initialization resolves each fault's mounting position, registers the
window ends as stops and sets the faults active at the start. Empty without
faults.

The faults' state carries over from one checkpoint segment to the next, and a
run resumed from a checkpoint file takes it from the file: whether each fault
is active and the level a stuck actuator holds.
"""
function get_fault_callbacks(num_sats::Int, args::SimulationConfiguration)
    sats = _fault_spacecraft(args)
    isempty(sats) && return ()
    function switch!(integrator)
        changed = false
        for i in sats
            sc = integrator.p.args.dynamics_model.spacecraft[i]
            changed |= refresh_faults!(sc, _fault_command_levels(integrator.p, i), integrator.t)
        end
        _mark_discontinuity!(integrator, changed)   # declared only when a fault switched
        return nothing
    end
    condition(u, t, integrator) = any(i -> fault_edge(integrator.p.args.dynamics_model.spacecraft[i], t), sats)
    # A fresh callback is built for every run. Its initialization also runs at
    # checkpoint-segment boundaries, which keep the faults' state (a stuck
    # actuator keeps the level it had at onset) and register the window ends as
    # stops of the next segment's solve. A resumed run's faults come from its
    # checkpoint, active and held as they were.
    initialized = Ref(false)
    function initialize!(cb, u, t, integrator)
        resumed = integrator.p.shared_buffers.callback_resume[]
        for i in sats
            sc = integrator.p.args.dynamics_model.spacecraft[i]
            initialized[] || resolve_faults!(sc)
            for f in sc.root.faults.faults
                initialized[] || resumed || (f.active = false)
                _maybe_add_control_tstop!(integrator, f.start_s)
                _maybe_add_control_tstop!(integrator, f.stop_s)
            end
        end
        initialized[] = true
        switch!(integrator)
        _mark_discontinuity!(integrator, true)
        return nothing
    end
    return (DiscreteCallback(condition, switch!; initialize=initialize!, save_positions=(false, false)),)
end

"""
    fault_save_fields(args) -> Vector{SaveField}

The `fault_active` field when a spacecraft carries faults: per spacecraft one
column per fault, in the order they were added, 1 while it is active and 0
otherwise (`sc{i}_fault_active_{j}`).
"""
function fault_save_fields(args::SimulationConfiguration)
    any(has_faults, args.dynamics_model.spacecraft) || return SaveField[]
    num_sats = length(args.dynamics_model.spacecraft)
    getter = (u, t, integrator) -> [Float64[f.active ? 1.0 : 0.0 for f in integrator.p.args.dynamics_model.spacecraft[i].root.faults.faults]
                                    for i in 1:num_sats]
    return SaveField[SaveField(:fault_active, getter; per_satellite=true, column_prefix="fault_active")]
end
