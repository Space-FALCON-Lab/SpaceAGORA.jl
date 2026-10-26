module IOSerialization

using SHA
using TOML
using Serialization
using Dates

using ..IOConfig

function _atomic_write_file(path::String, writer::Function; force::Bool=true)
    dir = dirname(path)
    mkpath(dir)
    tmp = joinpath(
        dir,
        string(".", basename(path), ".tmp.", getpid(), ".", Threads.threadid(), ".", rand(UInt))
    )
    try
        writer(tmp)
        mv(tmp, path; force=force)
    finally
        if isfile(tmp)
            rm(tmp; force=true)
        end
    end
    return path
end

function _sha256_hex(path::String)::String
    open(path, "r") do io
        return bytes2hex(SHA.sha256(read(io)))
    end
end

function _write_checkpoint!(args, t::Float64, u_state, checkpoint_schema_version;
    solver_mode::Union{Nothing, String}=nothing, runtime_state=nothing)
    paths = IOConfig._checkpoint_paths(args)
    payload = (
        schema_version=checkpoint_schema_version,
        created_utc=string(now(UTC)),
        t=t,
        solver_mode=solver_mode,
        u=deepcopy(u_state)
    )
    # Ordinary checkpoints retain the established payload schema.
    runtime_state === nothing || (payload = merge(payload, (runtime_state=deepcopy(runtime_state),)))
    _atomic_write_file(paths.data, tmp -> open(tmp, "w") do io
        serialize(io, payload)
    end)

    manifest = Dict{String, Any}(
        "schema_version" => checkpoint_schema_version,
        "created_utc" => string(now(UTC)),
        "time_s" => t,
        "solver_mode" => solver_mode,
        "data_path" => paths.data,
        "data_size_bytes" => filesize(paths.data),
        "data_sha256" => _sha256_hex(paths.data)
    )
    _atomic_write_file(paths.manifest, tmp -> open(tmp, "w") do io
        TOML.print(io, manifest)
    end)
    return nothing
end

function _load_checkpoint(args)
    paths = IOConfig._checkpoint_paths(args)
    if !isfile(paths.data)
        return nothing
    end
    payload = open(paths.data, "r") do io
        deserialize(io)
    end
    if !haskey(payload, :t) || !haskey(payload, :u)
        throw(ArgumentError("Checkpoint payload missing required keys (:t, :u)."))
    end
    solver_mode = haskey(payload, :solver_mode) ? payload[:solver_mode] : nothing
    checkpoint = (
        t=Float64(payload[:t]),
        u=payload[:u],
        solver_mode=solver_mode === nothing ? nothing : String(solver_mode),
        data_path=paths.data,
        manifest_path=paths.manifest
    )
    runtime_state = get(payload, :runtime_state, nothing)
    return runtime_state === nothing ? checkpoint : merge(checkpoint, (runtime_state=runtime_state,))
end

function _clear_checkpoint!(args)
    paths = IOConfig._checkpoint_paths(args)
    isfile(paths.data) && rm(paths.data; force=true)
    isfile(paths.manifest) && rm(paths.manifest; force=true)
    isfile(paths.rows) && rm(paths.rows; force=true)
    rm(paths.interlinks; force=true)
    rm(paths.stream_spacecraft; force=true)
    rm(paths.stream_run; force=true)
    return nothing
end

# The results rows a checkpointed run has saved, kept beside its checkpoint so a
# resumed run can include them: an append-only stream of serialized records, a
# header first, then one record of rows (`times`, `data`) per checkpoint. A record
# may first drop rows from the end of those before it (`drop`), so no record is
# ever rewritten. Each checkpoint records the stream's length at its time, so a
# resume reads exactly the rows up to that checkpoint and cuts away any a crash
# wrote after it.

# Append one record of rows, which first drops `drop` rows; `header` starts a new
# stream in place of any old one. Returns the stream's length in bytes.
function _append_checkpoint_rows!(path::String, times::AbstractVector{Float64}, data::AbstractVector;
                                  drop::Integer=0, header=nothing)::Int
    mkpath(dirname(path))
    open(path, header === nothing ? "a" : "w") do io
        header === nothing || serialize(io, header)
        serialize(io, (drop=Int(drop), times=Vector{Float64}(times), data=collect(data)))
        flush(io)
    end
    return filesize(path)
end

# Hash only the requested byte range. Writers extend a live context with newly
# appended bytes; readers hash the checkpoint-owned prefix, excluding crash tails.
# Reading in bounded chunks avoids retaining a second copy of a large history.
function _hash_checkpoint_rows!(context, path::String, start_bytes::Integer, end_bytes::Integer)
    0 <= start_bytes <= end_bytes || throw(ArgumentError("Invalid checkpoint history byte range."))
    open(path, "r") do io
        filesize(io) >= end_bytes || throw(EOFError())
        seek(io, start_bytes)
        buffer = Vector{UInt8}(undef, min(64 * 1024, end_bytes - start_bytes))
        remaining = end_bytes - start_bytes
        while remaining > 0
            wanted = min(length(buffer), remaining)
            count = readbytes!(io, buffer, wanted)
            count == wanted || throw(EOFError())
            SHA.update!(context, buffer, count)
            remaining -= count
        end
    end
    return context
end

# The header and the rows of the first `length_bytes` bytes of the stream, or
# `nothing` when the file is missing, shorter, or does not end a record there.
function _read_checkpoint_rows(path::String, length_bytes::Integer)
    isfile(path) && filesize(path) >= length_bytes || return nothing
    return try
        open(path, "r") do io
            header = deserialize(io)
            times = Float64[]
            data = Any[]
            while position(io) < length_bytes
                record = deserialize(io)
                for _ in 1:record.drop
                    pop!(times)
                    pop!(data)
                end
                append!(times, record.times)
                append!(data, record.data)
            end
            position(io) == length_bytes ? (header=header, times=times, data=data) : nothing
        end
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
end

# Cut the stream back to `length_bytes`, dropping rows written after a checkpoint.
function _truncate_checkpoint_rows!(path::String, length_bytes::Integer)::Nothing
    filesize(path) > length_bytes || return nothing
    open(path, "r+") do io
        truncate(io, length_bytes)
    end
    return nothing
end

export _atomic_write_file
export _sha256_hex
export _write_checkpoint!
export _load_checkpoint
export _clear_checkpoint!
export _append_checkpoint_rows!
export _read_checkpoint_rows
export _truncate_checkpoint_rows!

# Resolve a constellation member's explicitly recorded results destination.
# Both disk merging and viewer discovery must make the same recovery choice.
function _member_results_directory(folder::String, destination)
    if ismissing(destination)
        # An old summary cannot identify which recovery bundle supplied its
        # rows. Never silently substitute the obsolete top-level bundle.
        if isdir(folder) && any(name -> startswith(name, "resumed_at_") && isdir(joinpath(folder, name)), readdir(folder))
            throw(ArgumentError("Member $(basename(folder)) has recovery results but no recorded results_subdirectory; " *
                "select its intended bundle explicitly before merging."))
        end
        return folder
    end
    relative = normpath(String(destination))
    (isabspath(relative) || relative == ".." || startswith(relative, "../") || startswith(relative, "..\\")) &&
        throw(ArgumentError("Member $(basename(folder))'s results_subdirectory must stay inside its member folder."))
    return relative == "." ? folder : joinpath(folder, rstrip(relative, ('/', '\\')))
end

end # module IOSerialization
