"""Internal reporting implementation for completed RPO planner comparisons.

The imported function identities keep historical GuidanceHooks and
SimulationModel access compatible. Only presentation and output methods live
here; planner selection remains with guidance and tracking/acceptance are
supplied by the SimulationModel.RPOComparisonTracking implementation.

This file prepares the report data (result rows and their CSV, metric values,
station geometry, path and cost series) and writes the outputs. The Plotly
figures come from the PlotlyJS package extension through
`rpo_comparison_plots.jl`; the writer asks for them only when one of its two
plot flags is set, and then refuses before writing anything if PlotlyJS is not
loaded.
"""
module RPOComparisonReports
import ..GuidanceHooks: normalize_rpo_comparison_planner_type, rpo_comparison_planner_label
import ..GuidanceHooks: rpo_flatten_planner_results,
    rpo_write_namedtuple_csv,
    rpo_group_metric_mean,
    rpo_comparison_metric_specs,
    rpo_comparison_metric_value,
    rpo_comparison_axis_name,
    rpo_comparison_trace_axis_name,
    rpo_comparison_metric_summary_plot,
    rpo_comparison_station_mesh_trace,
    rpo_comparison_station_trace,
    rpo_comparison_path_family_plot,
    rpo_comparison_single_path_plot,
    rpo_comparison_failed_paths_plot,
    rpo_comparison_cost_iteration_xvalues,
    rpo_comparison_cost_iteration_plot,
    rpo_comparison_artifact_slug,
    rpo_write_failed_path_outputs,
    rpo_write_planner_comparison_outputs,
    rpo_comparison_metric_includes_result

"""Flatten nested planner-comparison output into one row per planner result."""
function rpo_flatten_planner_results(batch)
    rows = NamedTuple[]
    for planner in batch.planner_types
        append!(rows, batch.results_by_planner[planner])
    end
    return rows
end

"""Write named-tuple rows to a CSV file with stable column ordering."""
function rpo_write_namedtuple_csv(path::AbstractString, rows)
    isempty(rows) && return path
    mkpath(dirname(path))
    names = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for row in rows
            vals = String[]
            for name in names
                value = getproperty(row, name)
                push!(vals, replace(string(value), "," => ";"))
            end
            println(io, join(vals, ","))
        end
    end
    return path
end

"""Compute mean metric values grouped by planner label."""
function rpo_group_metric_mean(results, metric::Symbol)
    values = Float64[]
    for row in results
        value = getproperty(row, metric)
        value isa Bool && (value = value ? 1.0 : 0.0)
        f = Float64(value)
        isfinite(f) && push!(values, f)
    end
    isempty(values) && return NaN
    return sum(values) / length(values)
end

"""Return plot/export metadata for RPO planner comparison metrics."""
function rpo_comparison_metric_specs()
    return [
        (:success_pct, "Success", "success (%)"),
        (:planner_compute_time, "Planner Runtime", "s"),
        (:fuel_used_pct, "Fuel", "%"),
        (:control_effort_total, "Control Effort", "m/s"),
        (:thrust_saturation_fraction, "Thruster Sat.", "fraction"),
        (:min_clearance, "Min Clearance", "m"),
        (:final_pos_error, "Goal Error", "m"),
        (:actual_travel_duration, "Travel Duration", "s"),
        (:keepout_violations, "Keepout Violations", "count"),
    ]
end

"""Read one named comparison metric from a flattened result row."""
function rpo_comparison_metric_value(row, metric::Symbol)
    metric == :success_pct && return row.success ? 100.0 : 0.0
    value = getproperty(row, metric)
    value isa Bool && return value ? 1.0 : 0.0
    return Float64(value)
end

"""Decide whether a planner result contributes to a given metric aggregate."""
@inline rpo_comparison_metric_includes_result(row, metric::Symbol) = metric in (:success_pct, :keepout_violations) || row.success

"""Return the Plotly axis name for a metric subplot index."""
function rpo_comparison_axis_name(prefix::AbstractString, idx::Integer)
    idx == 1 && return Symbol(prefix)
    return Symbol(prefix, idx)
end

"""Return the Plotly trace axis reference for a metric subplot index."""
function rpo_comparison_trace_axis_name(prefix::AbstractString, idx::Integer)
    idx == 1 && return prefix
    return string(prefix, idx)
end

"""Choose x-axis values for a planner cost-history trace."""
function rpo_comparison_cost_iteration_xvalues(plan, costs)
    n = length(costs)
    n == 0 && return Int[]
    iterations = hasproperty(plan, :planner_iteration_count) ? Int(plan.planner_iteration_count) :
        (hasproperty(plan, :iterations) ? Int(plan.iterations) : n)
    if n == 1
        return [max(iterations, 0)]
    end
    return collect(1:n)
end

"""Return a filesystem-safe token for generated planner-comparison artifact names."""
function rpo_comparison_artifact_slug(value)
    token = replace(lowercase(String(value)), r"[^a-z0-9_=-]+" => "_")
    token = strip(token, '_')
    return isempty(token) ? "case" : token
end

"""Metric values per planner result, in `rpo_comparison_metric_specs` order: one
series per metric, with the planner label and case label of every included result."""
function _rpo_comparison_metric_series(batch)
    return map(rpo_comparison_metric_specs()) do (metric, label, units)
        planners = String[]
        values = Float64[]
        case_labels = String[]
        for planner in batch.planner_types
            planner_label = rpo_comparison_planner_label(planner)
            for result in batch.results_by_planner[planner]
                rpo_comparison_metric_includes_result(result, metric) || continue
                value = rpo_comparison_metric_value(result, metric)
                push!(planners, planner_label)
                push!(values, value)
                push!(case_labels, "case $(result.case_label)")
            end
        end
        (; metric, label, units, planners, values, case_labels)
    end
end

"""Unshared vertices and 0-based triangle indices of the station keepout mesh."""
function _rpo_comparison_station_mesh(triangles)
    tris = Matrix{Float64}(triangles)
    ntri = size(tris, 2) ÷ 3
    vertices = zeros(Float64, 3, 3 * ntri)
    i = zeros(Int, ntri)
    j = zeros(Int, ntri)
    k = zeros(Int, ntri)
    for tri_idx in 1:ntri
        src = 3 * (tri_idx - 1)
        vertices[:, src + 1] .= tris[:, src + 1]
        vertices[:, src + 2] .= tris[:, src + 2]
        vertices[:, src + 3] .= tris[:, src + 3]
        i[tri_idx] = src
        j[tri_idx] = src + 1
        k[tri_idx] = src + 2
    end
    return (; vertices, i, j, k)
end

"""Whether the comparison carries a station keepout mesh to draw instead of its point cloud."""
_rpo_comparison_has_station_mesh(batch) = hasproperty(batch, :station_triangles) && batch.station_triangles !== nothing

"""Station point-cloud columns kept when the cloud is downsampled to at most `max_points`."""
function _rpo_comparison_station_point_indices(points, max_points::Integer)
    stride = max(1, Int(ceil(size(points, 2) / max(1, max_points))))
    return 1:stride:size(points, 2)
end

"""One case's desired and tracked paths, its start and goal, and its raw planner path if it has one."""
function _rpo_comparison_case_path(plan)
    has_raw = hasproperty(plan, :raw_path)
    return (;
        case_label=plan.case.label,
        desired=plan.tracking.r_ref_rtn,
        tracked=plan.tracking.x_hist[1:3, :],
        start=plan.case.start_rtn,
        goal=plan.case.goal_rtn,
        has_raw,
        raw=has_raw ? plan.raw_path : nothing,
    )
end

"""Indices of a planner's failed results."""
_rpo_comparison_failed_indices(results) = [idx for idx in eachindex(results) if !results[idx].success]

"""Finite cost histories of a planner's cases, with their iteration values."""
function _rpo_comparison_cost_series(batch, planner_type)
    series = NamedTuple[]
    for plan in batch.plans_by_planner[planner_type]
        hasproperty(plan, :cost_history) || continue
        costs = [Float64(c) for c in plan.cost_history if isfinite(Float64(c))]
        isempty(costs) && continue
        iterations = rpo_comparison_cost_iteration_xvalues(plan, costs)
        push!(series, (; case_label=plan.case.label, iterations, costs))
    end
    return series
end

include("rpo_comparison_plots.jl")

"""Write planner comparison CSVs and plots to the requested output directory."""
function rpo_write_planner_comparison_outputs(
    batch;
    output_dir::AbstractString=batch.config.output_dir,
    write_plotly_outputs::Bool=batch.config.write_plotly_outputs,
    write_failed_path_outputs::Bool=batch.config.write_failed_path_outputs,
)
    _require_plotly_renderer(write_plotly_outputs, write_failed_path_outputs)
    mkpath(output_dir)
    rows = rpo_flatten_planner_results(batch)
    csv_path = rpo_write_namedtuple_csv(joinpath(output_dir, "rpo_planner_comparison_results.csv"), rows)

    plots = write_plotly_outputs ? _rpo_write_comparison_plots(batch, output_dir) :
        (metrics_plot=nothing, path_plot=nothing,
         method_path_plots=Dict{Symbol, String}(), cost_iteration_plots=Dict{Symbol, String}())
    failed_path_plots = write_failed_path_outputs ? rpo_write_failed_path_outputs(batch, output_dir) :
        Dict{Symbol, String}()

    return (
        csv=csv_path,
        metrics_plot=plots.metrics_plot,
        path_plot=plots.path_plot,
        method_path_plots=plots.method_path_plots,
        cost_iteration_plots=plots.cost_iteration_plots,
        failed_path_plots=failed_path_plots,
        plotly_outputs=write_plotly_outputs,
        failed_path_outputs=write_failed_path_outputs,
    )
end

end # module RPOComparisonReports
