module IOOutputs

using Arrow
using CSV
using DataFrames
using Dates
using TOML

using ..IOConfig
using ..IOSerialization

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

# Direct assembly for the shape every built-in per-satellite save field has:
# each snapshot holds one `Vector` of `num_sats` entries, and each entry is
# either a real or a fixed-width vector of reals.
#
# The generic path reaches those columns through three layers of throwaway
# intermediates -- one `sat_series` per satellite, then one `child_series` per
# component, and `collect` copies that last one again. For a 3-vector field that
# is 72 bytes of allocation per (satellite, row) to produce 24 bytes of column.
# Here each final column is allocated once and filled in one pass.
#
# The reason this is split across a function barrier rather than written as one
# loop: `saved_data` is a `Vector{SaveData}` and `SaveData` is
# `Dict{Symbol, Any}`, so a snapshot lookup is `Any`. Doing the fills directly
# off that makes every element access a dynamic dispatch, which measured *three
# times worse* than the generic path it was meant to replace -- the generic
# comprehensions are fast precisely because `[snapshot[name] for ...]` narrows
# to a concrete element type. So the narrowing comprehension is kept, and the
# concrete `values` is handed to a method that dispatches on its element type.
#
# Column order and column element types must match the generic path exactly,
# because the feather schema is part of the output. On types, note that the
# generic component path builds `value === nothing ? nothing : value[idx]` but
# still yields a `Vector{T}`, not a `Vector{Union{Nothing, T}}` -- the
# comprehension narrows to the element type actually produced.
function _direct_assembly_columns!(results_df::DataFrame, field, saved_data::Vector, num_sats::Int)::Bool
    (isempty(saved_data) || num_sats <= 0) && return false
    values = [snapshot[field.name] for snapshot in saved_data]
    return _fill_per_satellite_columns!(results_df, field.column_prefix, values, num_sats)
end

# Anything the two concrete methods below do not claim stays with the generic path.
_fill_per_satellite_columns!(::DataFrame, ::String, values, ::Int)::Bool = false

# One scalar column per satellite.
function _fill_per_satellite_columns!(
    results_df::DataFrame, prefix::String, values::Vector{V}, num_sats::Int
)::Bool where {T <: Real, V <: AbstractVector{T}}
    n_rows = length(values)
    all(row -> length(row) == num_sats, values) || return false
    for sat_idx in 1:num_sats
        column = Vector{T}(undef, n_rows)
        @inbounds for row in 1:n_rows
            column[row] = values[row][sat_idx]
        end
        results_df[!, "sc$(sat_idx)_$(prefix)"] = column
    end
    return true
end

# `n_comp` columns per satellite, named `_1`.._n` to match the generic path's
# `eachindex` walk.
#
# The column set has to be a property of the field, so every entry must carry
# the same number of components and be indexed from 1. That is checked here
# rather than inferred from the element type. `isbitstype(S)` was used for it
# once and is not proof of anything of the kind: it says the element type has a
# fixed bit layout, not a fixed length, and `UnitRange{Int}` -- two `Int`s, so a
# bitstype -- carries its length as runtime data. Both ways of getting it wrong
# are silent. A ragged field whose first entry is short drops a real column
# (`[[1:2, 10:12]]` wrote four columns where the generic path writes five,
# losing `sc2_x_3`); one whose first entry is long fabricates a value, because
# the write loop is `@inbounds` and reading past the end of a `UnitRange`
# computes an element that was never there.
function _fill_per_satellite_columns!(
    results_df::DataFrame, prefix::String, values::Vector{V}, num_sats::Int
)::Bool where {T <: Real, S <: AbstractVector{T}, V <: AbstractVector{S}}
    n_rows = length(values)
    all(row -> length(row) == num_sats, values) || return false
    n_comp = length(first(first(values)))
    n_comp > 0 || return false
    for row in values, entry in row
        (length(entry) == n_comp && firstindex(entry) == 1) || return false
    end
    for sat_idx in 1:num_sats
        columns = [Vector{T}(undef, n_rows) for _ in 1:n_comp]
        @inbounds for row in 1:n_rows
            entry = values[row][sat_idx]
            for comp in 1:n_comp
                columns[comp][row] = entry[comp]
            end
        end
        for comp in 1:n_comp
            results_df[!, "sc$(sat_idx)_$(prefix)_$(comp)"] = columns[comp]
        end
    end
    return true
end

function _append_save_field_columns!(results_df::DataFrame, field, saved_data::Vector, num_sats::Int)
    if field.per_satellite && _direct_assembly_columns!(results_df, field, saved_data, num_sats)
        return nothing
    end
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

function _write_results_csv!(results_df::DataFrame, args)::String
    primary_path = IOConfig._results_csv_path(args)
    started_s = time()
    existed_before = isfile(primary_path)
    try
        return IOSerialization._atomic_write_file(primary_path, tmp -> CSV.write(tmp, results_df); force=false)
    catch err
        if err isa ArgumentError && isfile(primary_path)
            mtime_s = try
                stat(primary_path).mtime
            catch
                0.0
            end
            concurrent_writer = (!existed_before) || (mtime_s >= started_s)
            if concurrent_writer
                collision_path = IOConfig._collision_results_csv_path(args)
                IOSerialization._atomic_write_file(collision_path, tmp -> CSV.write(tmp, results_df); force=false)
            end
            return IOSerialization._atomic_write_file(primary_path, tmp -> CSV.write(tmp, results_df); force=true)
        end
        rethrow(err)
    end
end

function _write_results_bundle!(
    results_df::DataFrame,
    times::Vector{Float64},
    args,
    results_bundle_schema_version;
    csv_path::Union{Nothing, String}=nothing
)
    prefix = IOConfig._results_bundle_prefix(args)
    feather_path = prefix * ".feather"
    manifest_path = prefix * ".manifest.toml"

    IOSerialization._atomic_write_file(feather_path, tmp -> Arrow.write(tmp, results_df))

    files = Dict{String, Any}()
    files["feather"] = Dict(
        "path" => feather_path,
        "size_bytes" => filesize(feather_path),
        "sha256" => IOSerialization._sha256_hex(feather_path)
    )

    if args.simulation_settings.save_csv
        csv_file_path = csv_path === nothing ? (prefix * ".csv") : csv_path
        if csv_path === nothing
            IOSerialization._atomic_write_file(csv_file_path, tmp -> CSV.write(tmp, results_df))
        end
        files["csv"] = Dict(
            "path" => csv_file_path,
            "size_bytes" => filesize(csv_file_path),
            "sha256" => IOSerialization._sha256_hex(csv_file_path)
        )
    end

    manifest = Dict{String, Any}(
        "schema_version" => results_bundle_schema_version,
        "created_utc" => string(now(UTC)),
        "mission_time_s" => args.mission_configuration.mission_time,
        "steps" => length(times),
        "spacecraft_count" => length(args.dynamics_model.spacecraft),
        "orientation_sim" => args.mission_configuration.orientation_sim,
        "files" => files
    )

    IOSerialization._atomic_write_file(manifest_path, tmp -> begin
        open(tmp, "w") do io
            TOML.print(io, manifest)
        end
    end)

    return nothing
end

export _append_saved_segment!
export _append_series_columns!
export _direct_assembly_columns!
export _build_results_dataframe
export _write_results_csv!
export _write_results_bundle!

end # module IOOutputs
