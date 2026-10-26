# Plotly figures of a completed RPO planner comparison, part of RPOComparisonReports.
# Core keeps every public method. The figures come from the PlotlyJS package
# extension (ext/SpaceAGORAPlotlyJSExt.jl): when PlotlyJS loads, it registers its
# renderer here, and each figure function hands the renderer the data the
# preparation functions in rpo_comparison_reports.jl compute. Without PlotlyJS,
# every figure function, and each writer asked for figures, throws
# PlotlyJSUnavailableError.

"""
    PlotlyJSUnavailableError()

An RPO Plotly figure was requested, or a planner comparison was asked to write
its Plotly outputs, without the PlotlyJS package loaded. The figures of
`rpo_path_plot`, `rpo_tracking_plot` and the planner-comparison plots and
reports come from SpaceAGORA's PlotlyJS package extension: add PlotlyJS to the
environment (`Pkg.add("PlotlyJS")`) and load it with `using PlotlyJS`. A planner
comparison writes its CSV without figures when `write_plotly_outputs` and
`write_failed_path_outputs` are `false`.
"""
struct PlotlyJSUnavailableError <: Exception end

const _PLOTLYJS_UNAVAILABLE_MESSAGE = "The RPO Plotly figures (`rpo_path_plot`, `rpo_tracking_plot` " *
    "and the planner-comparison plots and reports) come from SpaceAGORA's PlotlyJS package extension: " *
    "add PlotlyJS to the environment that runs them (`Pkg.add(\"PlotlyJS\")`) and load it with " *
    "`using PlotlyJS`. To write a planner comparison's CSV without figures, set " *
    "`write_plotly_outputs=false` and `write_failed_path_outputs=false`."
Base.showerror(io::IO, ::PlotlyJSUnavailableError) = print(io, _PLOTLYJS_UNAVAILABLE_MESSAGE)

# The extension's renderer, registered by its `__init__`; `nothing` until PlotlyJS loads.
const _PLOTLY_RENDERER = Ref{Any}(nothing)
_register_plotly_renderer!(renderer) = (_PLOTLY_RENDERER[] = renderer; nothing)

function _plotly_renderer()
    renderer = _PLOTLY_RENDERER[]
    renderer === nothing && throw(PlotlyJSUnavailableError())
    return renderer
end

# Refuse before any work when one of the given plot flags is set and PlotlyJS is not loaded.
function _require_plotly_renderer(flags::Bool...)
    any(flags) && _plotly_renderer()
    return nothing
end

# Builders and the figure writer, with methods only for the extension's renderer.
function _metric_summary_figure end
function _station_mesh_trace end
function _station_trace end
function _path_family_figure end
function _single_path_figure end
function _failed_paths_figure end
function _cost_iteration_figure end
function _save_figure end

"""Create summary bar plots for planner comparison metrics."""
rpo_comparison_metric_summary_plot(batch) = _metric_summary_figure(_plotly_renderer(), batch)

"""Build a Plotly mesh trace for the station keepout geometry."""
rpo_comparison_station_mesh_trace(triangles) = _station_mesh_trace(_plotly_renderer(), triangles)

"""Build a downsampled Plotly station point-cloud trace."""
rpo_comparison_station_trace(batch; max_points::Integer=2_000) =
    _station_trace(_plotly_renderer(), batch; max_points=max_points)

"""Create a 3D Plotly view of planner paths for all comparison cases."""
rpo_comparison_path_family_plot(batch; planner=nothing) =
    _path_family_figure(_plotly_renderer(), batch; planner=planner)

"""Create a 3D Plotly view for one planner/case path."""
rpo_comparison_single_path_plot(batch, plan; planner) =
    _single_path_figure(_plotly_renderer(), batch, plan; planner=planner)

"""Create a 3D Plotly view containing all failed paths for one planner."""
rpo_comparison_failed_paths_plot(batch; planner, failed_indices) =
    _failed_paths_figure(_plotly_renderer(), batch; planner=planner, failed_indices=failed_indices)

"""Create cost-versus-iteration traces for planners that report histories."""
rpo_comparison_cost_iteration_plot(batch; planner) =
    _cost_iteration_figure(_plotly_renderer(), batch; planner=planner)

"""Write one path HTML per planner containing all failed planner/case results."""
function rpo_write_failed_path_outputs(batch, output_dir::AbstractString)
    renderer = _plotly_renderer()
    failed_dir = joinpath(output_dir, "rpo_failed_path_plots")
    failed_path_plots = Dict{Symbol, String}()
    for planner in batch.planner_types
        failed_indices = _rpo_comparison_failed_indices(batch.results_by_planner[planner])
        isempty(failed_indices) && continue
        mkpath(failed_dir)
        path = joinpath(failed_dir, "rpo_failed_paths_$(planner).html")
        _save_figure(renderer, rpo_comparison_failed_paths_plot(batch; planner=planner, failed_indices=failed_indices), path)
        failed_path_plots[planner] = path
    end
    return failed_path_plots
end

"""Write the comparison's metric, path-family, per-planner path and cost figures."""
function _rpo_write_comparison_plots(batch, output_dir::AbstractString)
    renderer = _plotly_renderer()
    metrics_plot = joinpath(output_dir, "rpo_planner_comparison_metrics.html")
    _save_figure(renderer, rpo_comparison_metric_summary_plot(batch), metrics_plot)

    path_plot = joinpath(output_dir, "rpo_planner_comparison_paths.html")
    _save_figure(renderer, rpo_comparison_path_family_plot(batch), path_plot)

    method_path_plots = Dict{Symbol, String}()
    for planner in batch.planner_types
        method_path = joinpath(output_dir, "rpo_planner_paths_$(planner).html")
        _save_figure(renderer, rpo_comparison_path_family_plot(batch; planner=planner), method_path)
        method_path_plots[planner] = method_path
    end

    cost_iteration_plots = Dict{Symbol, String}()
    for planner in batch.planner_types
        cost_path = joinpath(output_dir, "rpo_planner_cost_vs_iteration_$(planner).html")
        _save_figure(renderer, rpo_comparison_cost_iteration_plot(batch; planner=planner), cost_path)
        cost_iteration_plots[planner] = cost_path
    end
    return (; metrics_plot, path_plot, method_path_plots, cost_iteration_plots)
end
