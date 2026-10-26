# Saved simulation-result tables and position/velocity samples.
# Include in a caller that already binds SimulationConfiguration, as common.jl
# does for the mission examples. This file loads no plotting or SPICE helpers.
using CSV
using DataFrames

# Non-SPICE result consumers share a loading boundary; event tables keep their
# separate schemas and SPICE readers retain their kernel-furnishing order.
function _read_simulation_results(csv_path::String)::DataFrame
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))
    return CSV.read(csv_path, DataFrame)
end

function _require_float_column(df::DataFrame, name::Symbol)::Vector{Float64}
    return Float64.(df[!, name])
end

# Keep each request limited to its own three columns and convert metres to
# kilometres (or metres/second to kilometres/second) after Float64 conversion.
function _simulation_vector_samples(args::SimulationConfiguration, columns::NTuple{3,Symbol})
    csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
    df = _read_simulation_results(csv_path)
    time_s = _require_float_column(df, :time)
    x = _require_float_column(df, columns[1]) ./ 1e3
    y = _require_float_column(df, columns[2]) ./ 1e3
    z = _require_float_column(df, columns[3]) ./ 1e3
    return time_s, x, y, z
end

function _simulation_position_samples(args::SimulationConfiguration)
    return _simulation_vector_samples(args, (:sc1_pos_1, :sc1_pos_2, :sc1_pos_3))
end

function _simulation_velocity_samples(args::SimulationConfiguration)
    return _simulation_vector_samples(args, (:sc1_vel_1, :sc1_vel_2, :sc1_vel_3))
end

"""
    _quaternion_history_samples(df::DataFrame; spacecraft_index)

Prepare the saved raw quaternion columns for the one-based spacecraft position
in the output table. Preserve component values and row order, including signed,
non-unit and zero quaternions; convert only time from seconds to minutes.
"""
function _quaternion_history_samples(df::DataFrame; spacecraft_index)
    spacecraft_index isa Integer && !(spacecraft_index isa Bool) && spacecraft_index > 0 ||
        throw(ArgumentError("spacecraft_index must be a positive one-based integer (the spacecraft position, not model.id)."))
    columns = ntuple(i -> Symbol("sc$(spacecraft_index)_q_$(i)"), 4)
    for column in (:time, columns...)
        column in propertynames(df) || throw(ArgumentError(
            "Quaternion history requires column $(column). Save all four attitude components for the selected spacecraft (orientation_sim=true with default save fields)."))
    end
    nrow(df) > 0 || throw(ArgumentError("Quaternion history requires at least one saved row."))
    for column in (:time, columns...)
        for (row, value) in enumerate(df[!, column])
            value isa Real && !(value isa Bool) && isfinite(value) || throw(ArgumentError(
                "Quaternion history column $(column), row $(row) must contain a finite real number; got $(repr(value))."))
        end
    end
    # Copy without Float64 conversion: preparation must retain the saved raw
    # values, including signed zero and integer values above Float64 precision.
    components = ntuple(i -> copy(df[!, columns[i]]), 4)
    return (time_min=copy(df[!, :time]) ./ 60,
            components=components,
            labels=("q1", "q2", "q3", "q4"),
            colors=(:red, :green, :blue, :orange),
            spacecraft_index=spacecraft_index)
end


"""
    _body_rate_history_samples(df::DataFrame; spacecraft_index)

Prepare explicitly saved `scN_body_rate_1..3` body angular rates in rad/s.
Copy the raw rates and preserve row order; convert time to minutes and displayed
rates to deg/min with the legacy `rad2deg.(rates) .* 60` operation. These columns
are custom save fields, not part of the default orientation output.
"""
function _body_rate_history_samples(df::DataFrame; spacecraft_index)
    spacecraft_index isa Integer && !(spacecraft_index isa Bool) && spacecraft_index > 0 ||
        throw(ArgumentError("spacecraft_index must be a positive one-based integer (the spacecraft position, not model.id)."))
    columns = ntuple(i -> Symbol("sc$(spacecraft_index)_body_rate_$(i)"), 3)
    for column in (:time, columns...)
        column in propertynames(df) || throw(ArgumentError(
            "Body-rate history requires column $(column). Explicitly save all three body angular rates in rad/s with column_prefix=\"body_rate\"; default save fields do not include them."))
    end
    nrow(df) > 0 || throw(ArgumentError("Body-rate history requires at least one saved row."))
    for column in (:time, columns...)
        for (row, value) in enumerate(df[!, column])
            value isa Real && !(value isa Bool) && isfinite(value) || throw(ArgumentError(
                "Body-rate history column $(column), row $(row) must contain a finite real number; got $(repr(value))."))
        end
    end
    time_min = copy(df[!, :time]) ./ 60
    for (row, value) in enumerate(time_min)
        isfinite(value) || throw(ArgumentError(
            "Body-rate history column time, row $(row) is nonfinite after conversion to minutes."))
    end
    rates_rad_s = ntuple(i -> copy(df[!, columns[i]]), 3)
    # Keep the legacy operation order, including its ordinary numeric rounding.
    # These are stored body rates, never derivatives inferred from quaternion data.
    components = ntuple(i -> rad2deg.(rates_rad_s[i]) .* 60, 3)
    for i in 1:3
        for (row, value) in enumerate(components[i])
            isfinite(value) || throw(ArgumentError(
                "Body-rate history column $(columns[i]), row $(row) is nonfinite after conversion to deg/min."))
        end
    end
    return (time_min=time_min,
            rates_rad_s=rates_rad_s,
            components=components,
            labels=("ω1", "ω2", "ω3"),
            colors=(:red, :green, :blue),
            spacecraft_index=spacecraft_index)
end
