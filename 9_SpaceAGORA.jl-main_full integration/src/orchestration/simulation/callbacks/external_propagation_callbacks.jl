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

function get_external_propagation_callback(num_sats::Int)
    condition(u, t, integrator) = true
    function affect!(integrator)
        buffers = integrator.p.shared_buffers
        u = integrator.u
        t = Float64(integrator.t)
        modified = false
        for (r, rt) in enumerate(buffers.external_runtimes)
            ExternalPropagation.external_sync!(rt, t) > 0 || continue
            for i in 1:num_sats
                buffers.external_owner[i] == r || continue
                st = ExternalPropagation.external_state(rt, buffers.external_local[i], t)
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
        modified && u_modified!(integrator, true)
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions=(false, false))
end
