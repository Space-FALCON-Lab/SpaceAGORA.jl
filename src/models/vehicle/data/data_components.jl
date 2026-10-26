# The onboard-data records: the data store and the data sources, mounted on
# spacecraft links with `add_data!` (`DataStorage`). Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# records, as it does the power records.

"""
    AbstractDataSource

Supertype of the onboard data sources: a payload, an instrument or a
housekeeping stream that produces data into its spacecraft's `DataStore`. A
source kind implements `data_production_bps(source, x, ctx, t)`, its production
rate in bits per second. The built-in kind is `DataSource`, a constant rate.
"""
abstract type AbstractDataSource end

"""
    DataStore(; capacity_bits, initial_bits=0.0)

The spacecraft's onboard data store, mounted on the root link (at most one).
The data it holds, in bits, is a state the engine integrates
(`data_stored_bits`): it grows with the spacecraft's production and with what
its incoming data links carry, and shrinks with what its outgoing links carry.
A link carries data only from a store that holds some and only into a store
with room, so links lose nothing: an empty store passes what comes in straight
on, and production that finds the store full is lost and counted
(`data_lost_bits`). `initial_bits` is the data held at the start of a run.
"""
mutable struct DataStore
    capacity_bits::Float64
    initial_bits::Float64
    # Run state, switched only at events and between integration steps (see
    # DataStorage): full (incoming links limited, production lost) and empty
    # (outgoing links limited to what comes in).
    full::Bool
    empty::Bool

    function DataStore(capacity_bits::Real, initial_bits::Real)
        (isfinite(capacity_bits) && capacity_bits > 0.0) || throw(ArgumentError(
            "DataStore capacity_bits must be finite and positive, got $(capacity_bits)."))
        (isfinite(initial_bits) && 0.0 <= initial_bits <= capacity_bits) || throw(ArgumentError(
            "DataStore initial_bits must be in [0, capacity_bits], got $(initial_bits)."))
        return new(Float64(capacity_bits), Float64(initial_bits), false, false)
    end
end

DataStore(; capacity_bits::Real, initial_bits::Real=0.0) = DataStore(capacity_bits, initial_bits)

"""
    DataSource(; rate_bps, name="")

A data source that produces `rate_bps` bits per second into its spacecraft's
`DataStore` while the spacecraft is active. A spacecraft may carry any number of
sources, on any of its links, and needs a `DataStore` for them.
"""
struct DataSource <: AbstractDataSource
    rate_bps::Float64
    name::String

    function DataSource(rate_bps::Real, name::AbstractString)
        (isfinite(rate_bps) && rate_bps >= 0.0) || throw(ArgumentError(
            "DataSource rate_bps must be finite and nonnegative, got $(rate_bps)."))
        return new(Float64(rate_bps), String(name))
    end
end

DataSource(; rate_bps::Real, name::AbstractString="") = DataSource(rate_bps, name)
