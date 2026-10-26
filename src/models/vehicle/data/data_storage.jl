"""
    DataStorage

Onboard data: a spacecraft's data store and its data sources, mounted on links
with [`add_data!`](@ref). When a store is mounted, the engine integrates the data
it holds, the data its sources produce and the production a full store loses,
and the data links carry data only from a store that holds some and only into a
store with room (see `InterlinkTerminals`).
"""
module DataStorage

using ..Components: DataStore, DataSource, AbstractDataSource
using ..SpacecraftModels: Link, SpacecraftModel

export add_data!, has_data_store, data_store, data_production_bps

_check_own_link(model::SpacecraftModel, link::Link) = any(l -> l === link, model.links) || throw(ArgumentError(
    "Data parts mount on one of the spacecraft's own links; this link is not in its link list."))

"""
    add_data!(model::SpacecraftModel, link::Link, part) -> part

Mount an onboard-data part on `link`: the `DataStore` (the root link, at most
one) or a data source, an `AbstractDataSource` such as a `DataSource` (any link,
any number). Returns the part. A spacecraft with a data source needs a store.
"""
function add_data!(model::SpacecraftModel, link::Link, store::DataStore)
    _check_own_link(model, link)
    link === model.root || throw(ArgumentError("The data store mounts on the spacecraft's root link."))
    isempty(link.data_stores) || throw(ArgumentError("The spacecraft already has a data store; mount at most one."))
    push!(link.data_stores, store)
    return store
end

function add_data!(model::SpacecraftModel, link::Link, source::AbstractDataSource)
    _check_own_link(model, link)
    push!(link.data_sources, source)
    return source
end

"""Whether `model` carries a data store, whose stored data the engine integrates."""
has_data_store(model::SpacecraftModel)::Bool = !isempty(model.root.data_stores)
has_data_store(model)::Bool = false

"""The spacecraft's data store, or `nothing`."""
data_store(model::SpacecraftModel) = has_data_store(model) ? model.root.data_stores[1] : nothing
data_store(model) = nothing

# Whether `model` carries a data source.
_has_data_sources(model::SpacecraftModel) = any(link -> !isempty(link.data_sources), model.links)

"""
    data_production_bps(source, x::StateSample, ctx::LinkContext, t::Float64) -> Float64

The rate (bit/s) at which the data source `source`, on the spacecraft whose
state sample is `x`, produces data at `t` seconds into the run; `ctx` gives the
time and the planet. The engine asks for it at every right-hand-side evaluation,
for each source of each active spacecraft, and integrates it into the
spacecraft's `data_produced_bits`. It also evaluates the hook when saving data
rates and when `DataBacklogPolicy` scores a candidate spacecraft transmitter
with a data store at a scheduling step. These calls can revisit the same state and time; return the
rate without advancing model state or relying on a particular call count.
A source kind of your own implements it, for example imaging only over a
region, a duty cycle, or a rate that depends on the attitude. The rate must be finite and nonnegative. There is no default: a run
refuses at setup a source kind without a method.
"""
function data_production_bps end

data_production_bps(source::DataSource, x, ctx, t::Float64) = source.rate_bps

# The production (bit/s) of every source on `model`, at state sample `x`.
function _production_bps(model::SpacecraftModel, x, ctx, t::Float64)::Float64
    total = 0.0
    for link in model.links, source in link.data_sources
        rate = Float64(data_production_bps(source, x, ctx, t))
        isfinite(rate) && rate >= 0.0 || throw(ArgumentError(
            "A $(nameof(typeof(source))) data source produced $(rate) bit/s; a production rate is finite and nonnegative."))
        total += rate
    end
    return total
end

# A store counts as full when it holds within 1e-9 of its capacity, and as empty
# when it holds within 1e-9 of its capacity of nothing, as a battery does.
_within_store_full(store::DataStore, bits::Real) = bits >= store.capacity_bits * (1.0 - 1.0e-9)
_within_store_empty(store::DataStore, bits::Real) = bits <= 1.0e-9 * store.capacity_bits

# At a store boundary event: mark the store empty when the data it holds has
# fallen to nothing, full when it has reached its capacity. Returns whether a
# flag changed.
function _data_store_boundary!(model::SpacecraftModel, stored_bits::Real)::Bool
    store = data_store(model)
    store === nothing && return false
    if _within_store_empty(store, stored_bits) && !store.empty
        store.empty = true
        return true
    elseif _within_store_full(store, stored_bits) && !store.full
        store.full = true
        return true
    end
    return false
end

# At the end of an integration step: an empty store that holds data again, or a
# full one with room again, takes and gives data as usual. Returns whether a flag
# changed.
function _update_data_state!(model::SpacecraftModel, stored_bits::Real)::Bool
    store = data_store(model)
    store === nothing && return false
    changed = false
    if store.empty && !_within_store_empty(store, stored_bits)
        store.empty = false
        changed = true
    end
    if store.full && !_within_store_full(store, stored_bits)
        store.full = false
        changed = true
    end
    return changed
end

# The flags at the start of a run, from the data the store holds then.
function _start_data_state!(model::SpacecraftModel, stored_bits::Real)
    store = data_store(model)
    store === nothing && return nothing
    store.full = _within_store_full(store, stored_bits)
    store.empty = _within_store_empty(store, stored_bits)
    return nothing
end

end # module DataStorage
