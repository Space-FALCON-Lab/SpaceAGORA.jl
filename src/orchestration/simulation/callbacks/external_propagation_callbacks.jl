# Sync callback for externally propagated spacecraft (SimulationModel.ExternalPropagation).
#
# After every accepted solver step each owner is advanced by whole owner steps while its time stays within
# the integrator time (the owner keeps an integer step counter, so the solver is never forced to take a
# step of the owner's size), and the shadow entries it owns in `u.sc` are overwritten from its absolute
# state. The callback is installed before every other discrete callback so guidance, navigation, control
# and saving callbacks see the synced state.

"""Boolean mask over the run's spacecraft that an external propagator owns; `nothing` when none does."""
function _external_owned_mask(args::SimulationConfiguration, num_sats::Int)
    isempty(args.external_propagators) && return nothing
    mask = fill(false, num_sats)
    for ep in args.external_propagators
        for i in ExternalPropagation.external_spacecraft(ep)
            1 <= i <= num_sats && (mask[i] = true)
        end
    end
    return mask
end

function _synchronize_external_owners!(integrator, num_sats::Int)
    buffers = integrator.p.shared_buffers
    u = integrator.u
    t = Float64(integrator.t)
    modified = false
    for (r, rt) in enumerate(buffers.external_runtimes)
        advanced = ExternalPropagation.external_sync!(rt, t) > 0
        configured = ExternalPropagation.external_capabilities(integrator.p.args.external_propagators[r]).configured_forces
        (advanced || configured) || continue
        for i in 1:num_sats
            buffers.external_owner[i] == r || continue
            st = ExternalPropagation.external_state(rt, t, buffers.external_local[i])
            sc = u.sc[i]
            sc.pos .= st.pos
            sc.vel .= st.vel
            if hasproperty(sc, :q)
                sc.q .= st.q
                sc.ω .= st.ω
            end
            modified = true
        end
    end
    modified && DiffEqBase.derivative_discontinuity!(integrator, true)
    return nothing
end

function _external_before_control_tick!(integrator, num_sats::Int)
    any(ep -> ExternalPropagation.external_capabilities(ep).configured_forces, integrator.p.args.external_propagators) || return nothing
    _synchronize_external_owners!(integrator, num_sats)
    for (ep, rt) in zip(integrator.p.args.external_propagators, integrator.p.shared_buffers.external_runtimes)
        ExternalPropagation.external_capabilities(ep).configured_forces || continue
        t = Float64(integrator.t)
        abs(ExternalPropagation.external_time(rt)-t) <= 16eps(max(abs(t),1.0)) ||
            throw(ArgumentError("External owner must reach the physical control boundary before the command changes."))
    end
    return nothing
end

# After the control models of a tick have updated: a control update acts on an owned spacecraft
# only through forces, so its position and velocity must still be the owner's (the next sync
# would silently overwrite a write). Then the owner learns that the held command changed.
function _external_after_control_tick!(integrator)
    buffers = integrator.p.shared_buffers
    t = Float64(integrator.t)
    for (r, (ep, rt)) in enumerate(zip(integrator.p.args.external_propagators, buffers.external_runtimes))
        ExternalPropagation.external_capabilities(ep).configured_forces || continue
        for i in eachindex(buffers.external_owner)
            buffers.external_owner[i] == r || continue
            st = ExternalPropagation.external_state(rt, t, buffers.external_local[i])
            sc = integrator.u.sc[i]
            sc.pos == st.pos && sc.vel == st.vel || throw(ArgumentError(
                "A control update changed the position or velocity of spacecraft $i, which an external owner propagates; " *
                "control must act on it through forces."))
        end
        ExternalPropagation.external_control_changed!(rt, t)
    end
    return nothing
end

function get_external_propagation_callback(num_sats::Int)
    condition(u, t, integrator) = true
    affect!(integrator) = _synchronize_external_owners!(integrator, num_sats)
    return DiscreteCallback(condition, affect!; save_positions=(false, false))
end

# Saved output epochs of configured-force owners are solver stops.
#
# The data-saving SavingCallback records `integrator.u` only when a requested save time equals
# `integrator.t`; a save time between two stops is read from the outer solver's interpolant. For an
# external owner that interpolant belongs to the shadow integration, not to the owner, so an
# uninterrupted run saved interpolated shadow rows between stops while a checkpointed run saved the
# owner state at its segment endpoints. With every save time a stop, the owner sync callback (which
# precedes the SavingCallback) installs the owner state there and the saved row is that state. A
# configured-force owner is synchronized into the shadow at every accepted step; other owners only when
# they advance, so the stops are limited to runs with a configured-force owner.
#
# The stops are the saver's own pending requests. The run's saver is the run-clock wrapper
# (`_run_clock_saving_callback`), whose initialization builds this solve's inner SavingCallback on the
# run's output clock and initializes it (DiffEqCallbacks `saving_initialize` fills its request queue).
# This callback follows the saver, so its initialization runs after that and reads the live inner
# saver (`_current_saving_affect`), never one from an earlier segment. Each queued time
# strictly inside the current solve becomes a stop with the identical Float64 value, and a step that
# reaches a stop lands on that value, so each saved time equals `integrator.t`. The solve's final time
# is already its last stop. Runs without a configured-force owner are unchanged.
function _external_output_stops_required(args::SimulationConfiguration)
    return any(ep -> ExternalPropagation.external_capabilities(ep).configured_forces, args.external_propagators)
end

function _register_external_output_stops!(integrator, saver)
    integrator.tdir > 0 || throw(ArgumentError("External output stops support forward integration only."))
    affect = _current_saving_affect(saver.affect!)
    hasproperty(affect, :saveat) ||
        throw(ArgumentError("External output stops need the SavingCallback request queue (`saveat`)."))
    isempty(DiffEqBase.get_tstops_array(integrator)) && return nothing
    t_now = integrator.t
    t_last = DiffEqBase.get_tstops_max(integrator)
    pending = deepcopy(affect.saveat)
    while !isempty(pending)
        s = pop!(pending)
        t_now < s < t_last && DiffEqBase.add_tstop!(integrator, s)
    end
    return nothing
end

"""Initialization-only callback that makes each pending save time of `saver` (the run's data-saving
SavingCallback) a solver stop. Install it after `saver`. `nothing` without a configured-force owner."""
function get_external_output_stop_callback(args::SimulationConfiguration, saver)
    _external_output_stops_required(args) || return nothing
    condition(u, t, integrator) = false
    affect!(integrator) = nothing
    initialize(c, u, t, integrator) = _register_external_output_stops!(integrator, saver)
    return DiscreteCallback(condition, affect!; initialize=initialize, save_positions=(false, false))
end
