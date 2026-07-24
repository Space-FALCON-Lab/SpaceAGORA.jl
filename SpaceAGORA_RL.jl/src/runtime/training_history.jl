mutable struct DiskBackedHistory{T} <: AbstractVector{T}
    path::String
    offsets::Vector{Int64}
    writer::Union{Nothing,IOStream}
end

function DiskBackedHistory(::Type{T}, path::AbstractString) where {T}
    resolved = abspath(String(path))
    mkpath(dirname(resolved))
    history = DiskBackedHistory{T}(resolved, Int64[], open(resolved, "w+"))
    finalizer(history) do value
        close_history!(value)
    end
    return history
end

Base.IndexStyle(::Type{<:DiskBackedHistory}) = IndexLinear()
Base.size(history::DiskBackedHistory) = (length(history.offsets),)
Base.length(history::DiskBackedHistory) = length(history.offsets)

function Base.push!(history::DiskBackedHistory{T}, value::T) where {T}
    writer = history.writer
    writer === nothing && throw(ArgumentError("cannot append to a closed training history"))
    push!(history.offsets, Int64(position(writer)))
    serialize(writer, value)
    return history
end

function close_history!(history::DiskBackedHistory)
    writer = history.writer
    if writer !== nothing
        flush(writer)
        close(writer)
        history.writer = nothing
    end
    return history
end

function _flush_history_writer!(history::DiskBackedHistory)
    history.writer === nothing || flush(history.writer)
    return nothing
end

function Base.getindex(history::DiskBackedHistory{T}, index::Int) where {T}
    checkbounds(history, index)
    _flush_history_writer!(history)
    return open(history.path, "r") do io
        seek(io, history.offsets[index])
        deserialize(io)::T
    end
end

function Base.getindex(history::DiskBackedHistory{T}, indices::AbstractUnitRange{<:Integer}) where {T}
    checkbounds(history, indices)
    isempty(indices) && return T[]
    _flush_history_writer!(history)
    values = Vector{T}(undef, length(indices))
    open(history.path, "r") do io
        for (destination, index) in enumerate(indices)
            seek(io, history.offsets[index])
            values[destination] = deserialize(io)::T
        end
    end
    return values
end

function _iterate_disk_history(history::DiskBackedHistory{T}, io::IOStream, index::Int) where {T}
    if index > length(history)
        close(io)
        return nothing
    end
    try
        value = deserialize(io)::T
        return value, (io, index + 1)
    catch
        close(io)
        rethrow()
    end
end

function Base.iterate(history::DiskBackedHistory)
    isempty(history) && return nothing
    _flush_history_writer!(history)
    io = open(history.path, "r")
    seek(io, history.offsets[1])
    return _iterate_disk_history(history, io, 1)
end

Base.iterate(history::DiskBackedHistory, state::Tuple{IOStream,Int}) =
    _iterate_disk_history(history, state[1], state[2])

struct MappedHistory{S,F} <: AbstractVector{Any}
    source::S
    transform::F
end

Base.IndexStyle(::Type{<:MappedHistory}) = IndexLinear()
Base.size(history::MappedHistory) = size(history.source)
Base.length(history::MappedHistory) = length(history.source)
Base.getindex(history::MappedHistory, index::Int) =
    history.transform(history.source[index])

function Base.iterate(history::MappedHistory)
    next = iterate(history.source)
    next === nothing && return nothing
    value, state = next
    return history.transform(value), state
end

function Base.iterate(history::MappedHistory, state)
    next = iterate(history.source, state)
    next === nothing && return nothing
    value, next_state = next
    return history.transform(value), next_state
end
