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
