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

function _checkpoint_payload(t::Float64, u_state, checkpoint_schema_version; solver_mode::Union{Nothing, String}=nothing, runtime_state=nothing)
    base = (
        schema_version=checkpoint_schema_version,
        created_utc=string(now(UTC)),
        t=t,
        solver_mode=solver_mode,
        u=deepcopy(u_state)
    )
    runtime_state === nothing && return base
    return merge(base, (runtime_state=deepcopy(runtime_state),))
end

function _write_checkpoint_payload!(args, payload)
    paths = IOConfig._checkpoint_paths(args)
    _atomic_write_file(paths.data, tmp -> open(tmp, "w") do io
        serialize(io, payload)
    end)

    schema_version = haskey(payload, :schema_version) ? payload[:schema_version] : "unknown"
    t = haskey(payload, :t) ? Float64(payload[:t]) : NaN
    manifest = Dict{String, Any}(
        "schema_version" => schema_version,
        "created_utc" => haskey(payload, :created_utc) ? string(payload[:created_utc]) : string(now(UTC)),
        "time_s" => t,
        "solver_mode" => haskey(payload, :solver_mode) && payload[:solver_mode] !== nothing ? string(payload[:solver_mode]) : "",
        "data_path" => paths.data,
        "data_size_bytes" => filesize(paths.data),
        "data_sha256" => _sha256_hex(paths.data)
    )
    _atomic_write_file(paths.manifest, tmp -> open(tmp, "w") do io
        TOML.print(io, manifest)
    end)
    return nothing
end

function _write_checkpoint!(args, t::Float64, u_state, checkpoint_schema_version; solver_mode::Union{Nothing, String}=nothing, runtime_state=nothing)
    payload = _checkpoint_payload(
        t,
        u_state,
        checkpoint_schema_version;
        solver_mode=solver_mode,
        runtime_state=runtime_state
    )
    return _write_checkpoint_payload!(args, payload)
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
    return (
        t=Float64(payload[:t]),
        u=payload[:u],
        solver_mode=solver_mode === nothing ? nothing : String(solver_mode),
        runtime_state=haskey(payload, :runtime_state) ? payload[:runtime_state] : nothing,
        data_path=paths.data,
        manifest_path=paths.manifest
    )
end

function _clear_checkpoint!(args)
    paths = IOConfig._checkpoint_paths(args)
    isfile(paths.data) && rm(paths.data; force=true)
    isfile(paths.manifest) && rm(paths.manifest; force=true)
    return nothing
end

export _atomic_write_file
export _sha256_hex
export _write_checkpoint!
export _write_checkpoint_payload!
export _load_checkpoint
export _clear_checkpoint!

end # module IOSerialization
