module ConstellationUtils

using YAML
using Dates
using SHA
using JLD2
using Tables
using Arrow
using DataFrames
using ProgressMeter
using Base.Threads

# ============================================================================
# Configuration Ingestion (from io/Config.jl)
# ============================================================================

"""
    ingest_yaml(path::AbstractString) -> Dict{String,Any}

Read a YAML file at `path` and return its contents as a Dict with String keys.
If the top-level YAML value is not a mapping, it will be returned as Dict("value" => parsed).
Throws ArgumentError if the file does not exist or if parsing fails.
"""
function ingest_yaml(path::AbstractString)::Dict{String,Any}
    if !isfile(path)
        throw(ArgumentError("YAML file not found at: $path"))
    end

    parsed = try
        YAML.load_file(path)
    catch err
        try
            open(path) do io
                YAML.load(io)
            end
        catch
            rethrow(err)
        end
    end

    converted = _to_string_keys(parsed)
    return converted isa AbstractDict ? Dict(converted) : Dict("value" => converted)
end

# Recursively convert Dict keys to String, preserve arrays and scalar values.
function _to_string_keys(x)
    if x === nothing
        return nothing
    elseif x isa AbstractDict
        out = Dict{String,Any}()
        for (k,v) in x
            out[string(k)] = _to_string_keys(v)
        end
        return out
    elseif x isa AbstractArray
        return [_to_string_keys(el) for el in x]
    else
        return x
    end
end

# ============================================================================
# Logging (from utils/Logging.jl)
# ============================================================================

const CONSTELLATION_LOG_IO = Ref{Union{Nothing, IO}}(nothing)
const CONSTELLATION_LOG_MIRROR = Ref{Bool}(true)

_constellation_ts() = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS.sss")

_bytes_to_gib(bytes::Integer) = round(Float64(bytes) / 1024^3; digits=3)

function constellation_memory_snapshot()
    try
        total_bytes = Sys.total_memory()
        free_bytes = Sys.free_memory()
        used_bytes = max(total_bytes - free_bytes, 0)
        return (
            sys_ram_used_gib = _bytes_to_gib(used_bytes),
            sys_ram_available_gib = _bytes_to_gib(free_bytes),
            sys_ram_total_gib = _bytes_to_gib(total_bytes),
        )
    catch
        return (
            sys_ram_used_gib = NaN,
            sys_ram_available_gib = NaN,
            sys_ram_total_gib = NaN,
        )
    end
end

function _constellation_memory_kv_pairs()
    snapshot = constellation_memory_snapshot()
    return (
        :sys_ram_used_gib => snapshot.sys_ram_used_gib,
        :sys_ram_available_gib => snapshot.sys_ram_available_gib,
        :sys_ram_total_gib => snapshot.sys_ram_total_gib,
    )
end

function _constellation_kv_string(kwargs)
    isempty(kwargs) && return ""
    parts = String[]
    for (k, v) in kwargs
        push!(parts, string(k, "=", repr(v)))
    end
    return " | " * join(parts, " ")
end

function _constellation_emit_line(line::String)
    io = CONSTELLATION_LOG_IO[]
    if io !== nothing
        try
            write(io, line * "\n")
            flush(io)
        catch
        end
    end
    if CONSTELLATION_LOG_MIRROR[]
        try
            println(line)
        catch
        end
    end
end

function _constellation_write(level::String, message::String; kwargs...)
    memory_kwargs = _constellation_memory_kv_pairs()
    line = "[$(_constellation_ts())] [$level] $message$(_constellation_kv_string((memory_kwargs..., kwargs...)))"
    Base.invokelatest(_constellation_emit_line, line)
end

function constellation_log(stage::AbstractString, message::AbstractString; kwargs...)
    _constellation_write("INFO", "[$stage] $message"; kwargs...)
end

function constellation_log_warn(stage::AbstractString, message::AbstractString; kwargs...)
    _constellation_write("WARN", "[$stage] $message"; kwargs...)
end

function constellation_log_error(stage::AbstractString, message::AbstractString; kwargs...)
    _constellation_write("ERROR", "[$stage] $message"; kwargs...)
end

function constellation_log_exception(stage::AbstractString, err; bt=nothing)
    constellation_log_error(stage, "exception"; error=sprint(showerror, err))
    if bt !== nothing
        _constellation_emit_line(sprint(show, bt))
    end
end

function make_constellation_log_path(config_dict; context::AbstractString="constellation_design")
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    custom_path = get(opt_params, "constellation_log_file", nothing)
    if custom_path !== nothing
        custom_str = String(custom_path)
        if !isempty(strip(custom_str))
            mkpath(dirname(custom_str))
            return custom_str
        end
    end
    log_dir = String(get(opt_params, "constellation_log_dir", "logs"))
    mkpath(log_dir)
    stamp = Dates.format(now(), dateformat"yyyymmdd_HHMMSS")
    return joinpath(log_dir, "$(context)_$(stamp).log")
end

function constellation_log_init!(config_dict; context::AbstractString="constellation_design")
    opt_params = get(config_dict, "optimizer_params", Dict{String,Any}())
    CONSTELLATION_LOG_MIRROR[] = Bool(get(opt_params, "constellation_log_console_mirror", true))
    path = make_constellation_log_path(config_dict; context=context)
    io = open(path, "a")
    CONSTELLATION_LOG_IO[] = io
    config_dict["constellation_log_path"] = path
    _constellation_write("INFO", "[logger] initialized";
        path=abspath(path),
        mirror_console=CONSTELLATION_LOG_MIRROR[],
        solver_trace=Bool(get(opt_params, "constellation_log_solver_trace", true)))
    return path
end

function constellation_log_close!()
    io = CONSTELLATION_LOG_IO[]
    if io !== nothing
        _constellation_write("INFO", "[logger] closing")
        close(io)
        CONSTELLATION_LOG_IO[] = nothing
    end
end

# ============================================================================
# Basic Utilities (from utils/utils.jl)
# ============================================================================

function summarize_index_vector(indices; max_items::Integer=50)
    idxs = Int[i for i in indices]
    total = length(idxs)
    total == 0 && return "[] (count=0)"
    return "(count=$(total))"
end

function hash_bin(pos, bin_size)
    return Tuple(floor.(Int, pos ./ bin_size))
end

function get_sat_positions(sat_positions, sc_id, timestep, axis=-1)
    if axis == -1
        x_out = sat_positions[Symbol("$(sc_id)_x")][timestep]
        y_out = sat_positions[Symbol("$(sc_id)_y")][timestep]
        z_out = sat_positions[Symbol("$(sc_id)_z")][timestep]
        return [x_out, y_out, z_out]
    else
        des_output = sat_positions[Symbol("$(sc_id)_$(axis)")][timestep]
        return des_output
    end
end

export ingest_yaml, constellation_log, constellation_log_warn, constellation_log_error
export constellation_log_exception, make_constellation_log_path
export constellation_log_init!, constellation_log_close!
export summarize_index_vector, hash_bin, get_sat_positions

end # module ConstellationUtils
