@inline _results_bundle_prefix(args) = SimulationModel.IOConfig._results_bundle_prefix(args)
@inline _results_csv_path(args) = SimulationModel.IOConfig._results_csv_path(args)
@inline _collision_results_csv_path(args) = SimulationModel.IOConfig._collision_results_csv_path(args)

@inline _atomic_write_file(path::String, writer::Function; force::Bool=true) =
    SimulationModel.IOSerialization._atomic_write_file(path, writer; force=force)

@inline _sha256_hex(path::String) = SimulationModel.IOSerialization._sha256_hex(path)

@inline _append_saved_segment!(times_acc::Vector{Float64}, data_acc::Vector, saved_values) =
    SimulationModel.IOOutputs._append_saved_segment!(times_acc, data_acc, saved_values)

@inline _append_series_columns!(results_df::DataFrame, prefix::String, series) =
    SimulationModel.IOOutputs._append_series_columns!(results_df, prefix, series)

@inline _find_sample_value(series) = SimulationModel.IOOutputs._find_sample_value(series)

@inline _build_results_dataframe(times::Vector{Float64}, saved_data::Vector, save_fields, args)::DataFrame =
    SimulationModel.IOOutputs._build_results_dataframe(times, saved_data, save_fields, args)

@inline _write_results_csv!(results_df::DataFrame, args)::String =
    SimulationModel.IOOutputs._write_results_csv!(results_df, args)

# Viewer scene sidecar (docs/architecture/interactive_visualization_plan.md).
# Opt-in: nothing is written unless `save_visualization_scene` is set, and it
# never runs when results themselves are disabled.
@inline function _write_visualization_scene_if_enabled!(args; density_params=nothing)::Union{Nothing, String}
    args.simulation_settings.results || return nothing
    args.simulation_settings.save_visualization_scene || return nothing
    return SimulationModel.SceneVisualization.write_visualization_scene!(args; density_params=density_params)
end

@inline function _write_results_bundle!(
    results_df::DataFrame,
    times::Vector{Float64},
    args;
    csv_path::Union{Nothing, String}=nothing,
)
    return SimulationModel.IOOutputs._write_results_bundle!(
        results_df,
        times,
        args,
        RESULTS_BUNDLE_SCHEMA_VERSION;
        csv_path=csv_path,
    )
end
