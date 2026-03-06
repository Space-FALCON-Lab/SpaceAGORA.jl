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

function _write_checkpoint!(args, t::Float64, u_state, checkpoint_schema_version)
    paths = IOConfig._checkpoint_paths(args)
    payload = (
        schema_version=checkpoint_schema_version,
        created_utc=string(now(UTC)),
        t=t,
        u=deepcopy(u_state)
    )
    _atomic_write_file(paths.data, tmp -> open(tmp, "w") do io
        serialize(io, payload)
    end)

    manifest = Dict{String, Any}(
        "schema_version" => checkpoint_schema_version,
        "created_utc" => string(now(UTC)),
        "time_s" => t,
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
    return (t=Float64(payload[:t]), u=payload[:u], data_path=paths.data, manifest_path=paths.manifest)
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
export _load_checkpoint
export _clear_checkpoint!

end # module IOSerialization
