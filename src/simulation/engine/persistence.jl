@inline function _results_bundle_prefix(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results")
end

@inline function _legacy_results_csv_path(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
end

@inline function _collision_results_csv_path(args)::String
    stamp = Dates.format(now(UTC), dateformat"yyyymmddTHHMMSSsss")
    token = string(stamp, ".", getpid(), ".", rand(UInt))
    return joinpath(args.simulation_settings.results_directory, "simulation_results.$token.csv")
end

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

function _write_legacy_results_csv!(results_df::DataFrame, args)::String
    primary_path = _legacy_results_csv_path(args)
    started_s = time()
    existed_before = isfile(primary_path)
    try
        return _atomic_write_file(primary_path, tmp -> CSV.write(tmp, results_df); force=false)
    catch err
        if err isa ArgumentError && isfile(primary_path)
            mtime_s = try
                stat(primary_path).mtime
            catch
                0.0
            end
            # Keep per-run collision artifacts only when a concurrent writer is likely.
            concurrent_writer = (!existed_before) || (mtime_s >= started_s)
            if concurrent_writer
                collision_path = _collision_results_csv_path(args)
                _atomic_write_file(collision_path, tmp -> CSV.write(tmp, results_df); force=false)
            end
            return _atomic_write_file(primary_path, tmp -> CSV.write(tmp, results_df); force=true)
        end
        rethrow(err)
    end
end

function _sha256_hex(path::String)::String
    open(path, "r") do io
        return bytes2hex(SHA.sha256(read(io)))
    end
end


function _append_saved_segment!(times_acc::Vector{Float64}, data_acc::Vector, saved_values)
    seg_len = length(saved_values.t)
    seg_len == 0 && return nothing
    start_idx = 1
    if !isempty(times_acc) && isapprox(times_acc[end], saved_values.t[1]; atol=0.0, rtol=0.0)
        start_idx = 2
    end
    if start_idx <= seg_len
        append!(times_acc, @view saved_values.t[start_idx:seg_len])
        append!(data_acc, @view saved_values.saveval[start_idx:seg_len])
    end
    return nothing
end

@inline function _is_flat_scalar(value)::Bool
    return value === missing || value === nothing || value isa Number || value isa AbstractString || value isa Symbol || value isa Bool
end

function _find_sample_value(series)
    for value in series
        if value !== nothing
            return value
        end
    end
    return nothing
end

function _append_series_columns!(results_df::DataFrame, prefix::String, series)
    sample = _find_sample_value(series)
    if sample === nothing || _is_flat_scalar(sample)
        results_df[!, prefix] = collect(series)
        return nothing
    end

    if sample isa NamedTuple
        for key in keys(sample)
            child_series = [value === nothing ? nothing : getproperty(value, key) for value in series]
            _append_series_columns!(results_df, string(prefix, "_", key), child_series)
        end
        return nothing
    end

    if sample isa AbstractDict
        for key in sort!(collect(keys(sample)); by=string)
            child_series = [value === nothing ? nothing : value[key] for value in series]
            _append_series_columns!(results_df, string(prefix, "_", key), child_series)
        end
        return nothing
    end

    if sample isa Tuple || sample isa AbstractArray
        for idx in eachindex(sample)
            child_series = [value === nothing ? nothing : value[idx] for value in series]
            _append_series_columns!(results_df, string(prefix, "_", idx), child_series)
        end
        return nothing
    end

    results_df[!, prefix] = collect(series)
    return nothing
end

function _append_save_field_columns!(results_df::DataFrame, field, saved_data::Vector, num_sats::Int)
    field_series = [snapshot[field.name] for snapshot in saved_data]
    if field.per_satellite
        for sat_idx in 1:num_sats
            sat_series = [value[sat_idx] for value in field_series]
            _append_series_columns!(results_df, "sc$(sat_idx)_$(field.column_prefix)", sat_series)
        end
        return nothing
    end
    _append_series_columns!(results_df, field.column_prefix, field_series)
    return nothing
end

function _build_results_dataframe(times::Vector{Float64}, saved_data::Vector, save_fields, args)::DataFrame
    results_df = DataFrame(time=times)
    num_sats = length(args.dynamics_model.spacecraft)
    for field in save_fields
        _append_save_field_columns!(results_df, field, saved_data, num_sats)
    end
    return results_df
end

function _write_results_bundle!(
    results_df::DataFrame,
    times::Vector{Float64},
    args;
    csv_path::Union{Nothing, String}=nothing
)
    prefix = _results_bundle_prefix(args)
    feather_path = prefix * ".feather"
    manifest_path = prefix * ".manifest.toml"

    _atomic_write_file(feather_path, tmp -> Arrow.write(tmp, results_df))

    files = Dict{String, Any}()
    files["feather"] = Dict(
        "path" => feather_path,
        "size_bytes" => filesize(feather_path),
        "sha256" => _sha256_hex(feather_path)
    )

    if args.simulation_settings.save_csv
        csv_file_path = csv_path === nothing ? (prefix * ".csv") : csv_path
        if csv_path === nothing
            _atomic_write_file(csv_file_path, tmp -> CSV.write(tmp, results_df))
        end
        files["csv"] = Dict(
            "path" => csv_file_path,
            "size_bytes" => filesize(csv_file_path),
            "sha256" => _sha256_hex(csv_file_path)
        )
    end

    manifest = Dict{String, Any}(
        "schema_version" => RESULTS_BUNDLE_SCHEMA_VERSION,
        "created_utc" => string(now(UTC)),
        "mission_time_s" => args.mission_configuration.mission_time,
        "steps" => length(times),
        "spacecraft_count" => length(args.dynamics_model.spacecraft),
        "orientation_sim" => args.mission_configuration.orientation_sim,
        "files" => files
    )

    _atomic_write_file(manifest_path, tmp -> begin
        open(tmp, "w") do io
            TOML.print(io, manifest)
        end
    end)

    return nothing
end

