# Onboard data during a run (see DataStorage): the events that mark a data store
# empty or full and the step-end update of its flags, and the storage columns.

# Spacecraft with a data store, whose stored data is a state.
_data_store_spacecraft(args) = [i for i in eachindex(args.dynamics_model.spacecraft)
                                if has_data_store(args.dynamics_model.spacecraft[i])]

# Step-end update of every store's flags.
function _update_data_states!(integrator, sats::Vector{Int})
    spacecraft = integrator.p.args.dynamics_model.spacecraft
    changed = false
    for i in sats
        changed |= _update_data_state!(spacecraft[i], integrator.u.sc[i].data_stored_bits)
    end
    _mark_discontinuity!(integrator, changed)
    return changed
end

"""
    get_data_callbacks(num_sats, args) -> Tuple

For runs with a data store: an event on each store emptying or reaching its
capacity, which marks it empty or full, and a step-end update of its flags. An
empty store sends no more than comes in, and a full one takes in no more than it
sends (see `InterlinkTerminals`). Both mark the derivative discontinuous when
they change anything, so no integration step sees a store's flags change. Empty
for runs without a store.

The flags carry over from one checkpoint segment to the next, and a run resumed
from a checkpoint file takes them from the file, as a battery's.
"""
function get_data_callbacks(num_sats::Int, args::SimulationConfiguration)
    sats = _data_store_spacecraft(args)
    isempty(sats) && return ()
    spacecraft = args.dynamics_model.spacecraft
    function condition!(out, u, t, integrator)
        @inbounds for (k, i) in enumerate(sats)
            bits = u.sc[i].data_stored_bits
            out[2k - 1] = bits
            out[2k] = data_store(spacecraft[i]).capacity_bits - bits
        end
        return nothing
    end
    # `events` holds one crossing direction per condition (0: no crossing).
    function boundary!(integrator, events)
        hit = false
        for e in eachindex(events)
            events[e] == 0 && continue
            i = sats[(e + 1) ÷ 2]
            hit |= _data_store_boundary!(integrator.p.args.dynamics_model.spacecraft[i], integrator.u.sc[i].data_stored_bits)
        end
        hit || (_mark_discontinuity!(integrator, false); return nothing)
        _update_data_states!(integrator, sats)
        _mark_discontinuity!(integrator, true)
        return nothing
    end
    boundary = VectorContinuousCallback(condition!, boundary!, 2 * length(sats); save_positions=(false, false))   # both crossing directions
    # A fresh callback is built for every run. Its initialization also runs at
    # checkpoint-segment boundaries, which keep the flags; a resumed run's come
    # from its checkpoint.
    initialized = Ref(false)
    function initialize!(cb, u, t, integrator)
        if !initialized[] && !integrator.p.shared_buffers.callback_resume[]
            for i in sats
                _start_data_state!(integrator.p.args.dynamics_model.spacecraft[i], u.sc[i].data_stored_bits)
            end
            initialized[] = true
        end
        _update_data_states!(integrator, sats)
        _mark_discontinuity!(integrator, true)
        return nothing
    end
    step_end = DiscreteCallback((u, t, integrator) -> true, integrator -> _update_data_states!(integrator, sats);
        initialize=initialize!, save_positions=(false, false))
    return (boundary, step_end)
end

function _save_data_storage(num_sats::Int, u, t, integrator, which::Symbol)
    spacecraft = integrator.p.args.dynamics_model.spacecraft
    out = fill(NaN, num_sats)
    @inbounds for i in 1:num_sats
        store = data_store(spacecraft[i])
        store === nothing && continue
        out[i] = which === :data_store_fill ? Float64(u.sc[i].data_stored_bits) / store.capacity_bits :
                 Float64(getproperty(u.sc[i], which))
    end
    return out
end

"""
    data_storage_save_fields(args) -> Vector{SaveField}

The storage columns, when a spacecraft carries a data store: per spacecraft
`data_stored_bits`, the data its store holds; `data_store_fill`, that as a
fraction of the capacity; `data_produced_bits`, the data its sources have
produced since the run started; and `data_lost_bits`, the production a full
store could not keep. NaN for a spacecraft without a store. Empty otherwise.
"""
function data_storage_save_fields(args::SimulationConfiguration)
    any(has_data_store, args.dynamics_model.spacecraft) || return SaveField[]
    num_sats = length(args.dynamics_model.spacecraft)
    return SaveField[SaveField(name, (u, t, integrator) -> _save_data_storage(num_sats, u, t, integrator, name); per_satellite=true)
                     for name in (:data_stored_bits, :data_store_fill, :data_produced_bits, :data_lost_bits)]
end
