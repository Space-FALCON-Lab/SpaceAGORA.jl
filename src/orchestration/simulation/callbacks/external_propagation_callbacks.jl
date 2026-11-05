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

function _external_after_control_tick!(integrator)
    for (ep, rt) in zip(integrator.p.args.external_propagators, integrator.p.shared_buffers.external_runtimes)
        ExternalPropagation.external_capabilities(ep).configured_forces || continue
        ExternalPropagation.external_control_changed!(rt, Float64(integrator.t))
    end
    return nothing
end

function get_external_propagation_callback(num_sats::Int)
    condition(u, t, integrator) = true
    affect!(integrator) = _synchronize_external_owners!(integrator, num_sats)
    return DiscreteCallback(condition, affect!; save_positions=(false, false))
end
