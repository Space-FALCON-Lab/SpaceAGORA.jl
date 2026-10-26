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

@inline _write_link_schedule!(schedule::DataFrame, args; csv::Bool, feather::Bool) =
    SimulationModel.IOOutputs._write_link_schedule!(schedule, args; csv=csv, feather=feather)
@inline _build_long_results_tables(times::Vector{Float64}, saved_data::Vector, save_fields, args) =
    SimulationModel.IOOutputs._build_long_results_tables(times, saved_data, save_fields, args)

@inline _write_long_results!(tables, times::Vector{Float64}, args; write_bundle::Bool=true,
                             resumed_from_s::Union{Nothing, Float64}=nothing, link_schedule=nothing) =
    SimulationModel.IOOutputs._write_long_results!(tables, times, args, LONG_RESULTS_BUNDLE_SCHEMA_VERSION;
        write_bundle=write_bundle, resumed_from_s=resumed_from_s, link_schedule=link_schedule)

@inline _append_long_stream!(path::String, csv_path::Union{Nothing, String}, df::DataFrame, token::String) =
    SimulationModel.IOOutputs._append_long_stream!(path, csv_path, df, token)

@inline _long_stream_token(path::String) = SimulationModel.IOOutputs._long_stream_token(path)

@inline _finish_long_stream!(streams, args; steps::Integer, first_time_s, resumed_from_s, link_schedule=nothing) =
    SimulationModel.IOOutputs._finish_long_stream!(streams, args, LONG_RESULTS_BUNDLE_SCHEMA_VERSION;
        steps=steps, first_time_s=first_time_s, resumed_from_s=resumed_from_s, link_schedule=link_schedule)

@inline _write_long_manifest!(args, files, steps::Integer, first_time_s, resumed_from_s; extra=()) =
    SimulationModel.IOOutputs._write_long_manifest!(args, LONG_RESULTS_BUNDLE_SCHEMA_VERSION, files, steps,
        first_time_s, resumed_from_s; extra=extra)

# Viewer scene sidecar (docs/architecture/interactive_visualization_plan.md).
# Opt-in: nothing is written unless `save_visualization_scene` is set, and it
# never runs when results themselves are disabled.
@inline function _write_visualization_scene_if_enabled!(args; density_params=nothing)::Union{Nothing, String}
    args.simulation_settings.results || return nothing
    args.simulation_settings.save_visualization_scene || return nothing
    SimulationModel.IOOutputs._invalidate_mission_history!(args)
    return SimulationModel.SceneVisualization.write_visualization_scene!(args; density_params=density_params)
end

@inline function _write_results_bundle!(
    results_df::DataFrame,
    times::Vector{Float64},
    args;
    csv_path::Union{Nothing, String}=nothing,
    resumed_from_s::Union{Nothing, Float64}=nothing,
    link_schedule=nothing,
)
    return SimulationModel.IOOutputs._write_results_bundle!(
        results_df,
        times,
        args,
        RESULTS_BUNDLE_SCHEMA_VERSION;
        csv_path=csv_path,
        resumed_from_s=resumed_from_s,
        link_schedule=link_schedule,
    )
end
