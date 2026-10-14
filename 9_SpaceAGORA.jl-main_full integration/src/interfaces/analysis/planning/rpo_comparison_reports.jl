"""Internal reporting implementation for completed RPO planner comparisons.

The imported function identities keep historical GuidanceHooks and
SimulationModel access compatible. Only presentation and output methods live
here; planning, tracking and acceptance are supplied by the GNC owners.
"""
module RPOComparisonReports
import PlotlyJS
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

"""Create summary bar plots for planner comparison metrics."""
function rpo_comparison_metric_summary_plot(batch)
    specs = rpo_comparison_metric_specs()
    traces = PlotlyJS.GenericTrace[]
    rows = 3
    cols = 3
    layout_kwargs = Dict{Symbol, Any}(
        :title => "RPO Planner Comparison Metrics",
        :height => 920,
        :width => 1200,
        :showlegend => false,
        :margin => PlotlyJS.attr(l=70, r=30, t=80, b=50),
        :plot_bgcolor => "rgb(250,250,250)",
    )
    xgap = 0.055
    ygap = 0.095
    xspan = (1.0 - xgap * (cols - 1)) / cols
    yspan = (1.0 - ygap * (rows - 1)) / rows

    for (idx, (metric, label, units)) in enumerate(specs)
        row_idx = div(idx - 1, cols) + 1
        col_idx = mod(idx - 1, cols) + 1
        x0 = (col_idx - 1) * (xspan + xgap)
        x1 = x0 + xspan
        y1 = 1.0 - (row_idx - 1) * (yspan + ygap)
        y0 = y1 - yspan
        xaxis = rpo_comparison_trace_axis_name("x", idx)
        yaxis = rpo_comparison_trace_axis_name("y", idx)

        xvals = String[]
        yvals = Float64[]
        case_labels = String[]
        for planner in batch.planner_types
            planner_label = rpo_comparison_planner_label(planner)
            for result in batch.results_by_planner[planner]
                rpo_comparison_metric_includes_result(result, metric) || continue
                value = rpo_comparison_metric_value(result, metric)
                push!(xvals, planner_label)
                push!(yvals, value)
                push!(case_labels, "case $(result.case_label)")
            end
        end

        push!(traces, PlotlyJS.box(
            x=xvals,
            y=yvals,
            name=label,
            boxpoints=false,
            boxmean=true,
            fillcolor="rgba(64,120,165,0.28)",
            line=PlotlyJS.attr(color="rgb(42,86,125)", width=2),
            marker=PlotlyJS.attr(color="rgb(42,86,125)"),
            xaxis=xaxis,
            yaxis=yaxis,
            hovertemplate="$(label)<br>%{x}: %{y:.5g} $(units)<extra></extra>",
        ))
        push!(traces, PlotlyJS.scatter(
            x=xvals,
            y=yvals,
            mode="markers",
            name="$(label) cases",
            marker=PlotlyJS.attr(size=8, color="rgb(205,92,55)", line=PlotlyJS.attr(color="white", width=1)),
            text=case_labels,
            xaxis=xaxis,
            yaxis=yaxis,
            hovertemplate="%{text}<br>%{x}<br>$(label): %{y:.5g} $(units)<extra></extra>",
        ))

        layout_kwargs[rpo_comparison_axis_name("xaxis", idx)] = PlotlyJS.attr(
            domain=[x0, x1],
            anchor=yaxis,
            title=label,
        )
        layout_kwargs[rpo_comparison_axis_name("yaxis", idx)] = PlotlyJS.attr(
            domain=[y0, y1],
            anchor=xaxis,
            title=units,
            zeroline=true,
            gridcolor="rgb(225,225,225)",
        )
    end

    return PlotlyJS.Plot(traces, PlotlyJS.Layout(; layout_kwargs...))
end

"""Build a Plotly mesh trace for the station keepout geometry."""
function rpo_comparison_station_mesh_trace(triangles)
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
    return PlotlyJS.mesh3d(
        x=vertices[1, :],
        y=vertices[2, :],
        z=vertices[3, :],
        i=i,
        j=j,
        k=k,
        color="rgb(150,160,170)",
        opacity=0.55,
        flatshading=true,
        name="Gateway mesh",
        hoverinfo="skip",
    )
end

"""Build a downsampled Plotly station point-cloud trace."""
function rpo_comparison_station_trace(batch; max_points::Integer=2_000)
    if hasproperty(batch, :station_triangles) && batch.station_triangles !== nothing
        return rpo_comparison_station_mesh_trace(batch.station_triangles)
    end
    pts = batch.geometry.station.points_body
    stride = max(1, Int(ceil(size(pts, 2) / max(1, max_points))))
    idxs = 1:stride:size(pts, 2)
    return PlotlyJS.scatter3d(
        x=pts[1, idxs],
        y=pts[2, idxs],
        z=pts[3, idxs],
        mode="markers",
        marker=PlotlyJS.attr(size=2, color="rgb(150,160,170)", opacity=0.30),
        name="station point cloud",
        hoverinfo="skip",
    )
end

"""Create a 3D Plotly view of planner paths for all comparison cases."""
function rpo_comparison_path_family_plot(batch; planner=nothing)
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    planners = planner === nothing ? batch.planner_types : [normalize_rpo_comparison_planner_type(planner)]
    for planner_type in planners
        label = rpo_comparison_planner_label(planner_type)
        for plan in batch.plans_by_planner[planner_type]
            case_label = plan.case.label
            ref = plan.tracking.r_ref_rtn
            actual = plan.tracking.x_hist[1:3, :]
            push!(traces, PlotlyJS.scatter3d(
                x=ref[1, :],
                y=ref[2, :],
                z=ref[3, :],
                mode="lines",
                line=PlotlyJS.attr(width=4, dash="dot"),
                name="$(label) case $(case_label) desired",
            ))
            push!(traces, PlotlyJS.scatter3d(
                x=actual[1, :],
                y=actual[2, :],
                z=actual[3, :],
                mode="lines",
                line=PlotlyJS.attr(width=4),
                name="$(label) case $(case_label) tracked",
            ))
            push!(traces, PlotlyJS.scatter3d(
                x=[plan.case.start_rtn[1]],
                y=[plan.case.start_rtn[2]],
                z=[plan.case.start_rtn[3]],
                mode="markers",
                marker=PlotlyJS.attr(size=6, color="rgb(45,150,90)", symbol="circle"),
                name="$(label) case $(case_label) start",
            ))
            push!(traces, PlotlyJS.scatter3d(
                x=[plan.case.goal_rtn[1]],
                y=[plan.case.goal_rtn[2]],
                z=[plan.case.goal_rtn[3]],
                mode="markers",
                marker=PlotlyJS.attr(size=7, color="rgb(190,60,55)", symbol="diamond"),
                name="$(label) case $(case_label) goal",
            ))
        end
    end
    title = planner === nothing ? "RPO Planner Comparison Path Family" :
        "$(rpo_comparison_planner_label(planner)) RPO Path Family"
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title=title,
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis_title="radial (m)",
                yaxis_title="along-track (m)",
                zaxis_title="cross-track (m)",
            ),
        ),
    )
end

"""Create a 3D Plotly view for one planner/case path."""
function rpo_comparison_single_path_plot(batch, plan; planner)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    label = rpo_comparison_planner_label(planner_type)
    case_label = plan.case.label
    ref = plan.tracking.r_ref_rtn
    actual = plan.tracking.x_hist[1:3, :]
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    push!(traces, PlotlyJS.scatter3d(
        x=ref[1, :],
        y=ref[2, :],
        z=ref[3, :],
        mode="lines",
        line=PlotlyJS.attr(width=5, dash="dot", color="rgb(45,95,170)"),
        name="desired",
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=actual[1, :],
        y=actual[2, :],
        z=actual[3, :],
        mode="lines",
        line=PlotlyJS.attr(width=5, color="rgb(205,92,55)"),
        name="tracked",
    ))
    if hasproperty(plan, :raw_path)
        raw = Matrix{Float64}(plan.raw_path)
        push!(traces, PlotlyJS.scatter3d(
            x=raw[1, :],
            y=raw[2, :],
            z=raw[3, :],
            mode="lines+markers",
            line=PlotlyJS.attr(width=2, dash="dash", color="rgb(110,110,110)"),
            marker=PlotlyJS.attr(size=3, color="rgb(110,110,110)"),
            name="raw planner path",
        ))
    end
    push!(traces, PlotlyJS.scatter3d(
        x=[plan.case.start_rtn[1]],
        y=[plan.case.start_rtn[2]],
        z=[plan.case.start_rtn[3]],
        mode="markers",
        marker=PlotlyJS.attr(size=7, color="rgb(45,150,90)", symbol="circle"),
        name="start",
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=[plan.case.goal_rtn[1]],
        y=[plan.case.goal_rtn[2]],
        z=[plan.case.goal_rtn[3]],
        mode="markers",
        marker=PlotlyJS.attr(size=8, color="rgb(190,60,55)", symbol="diamond"),
        name="goal",
    ))
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="$(label) Failed Case $(case_label) Path",
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis_title="radial (m)",
                yaxis_title="along-track (m)",
                zaxis_title="cross-track (m)",
            ),
        ),
    )
end

"""Create a 3D Plotly view containing all failed paths for one planner."""
function rpo_comparison_failed_paths_plot(batch; planner, failed_indices)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    label = rpo_comparison_planner_label(planner_type)
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    plans = batch.plans_by_planner[planner_type]
    for idx in failed_indices
        plan = plans[idx]
        case_label = plan.case.label
        ref = plan.tracking.r_ref_rtn
        actual = plan.tracking.x_hist[1:3, :]
        push!(traces, PlotlyJS.scatter3d(
            x=ref[1, :],
            y=ref[2, :],
            z=ref[3, :],
            mode="lines",
            line=PlotlyJS.attr(width=4, dash="dot"),
            name="case $(case_label) desired",
        ))
        push!(traces, PlotlyJS.scatter3d(
            x=actual[1, :],
            y=actual[2, :],
            z=actual[3, :],
            mode="lines",
            line=PlotlyJS.attr(width=4),
            name="case $(case_label) tracked",
        ))
        if hasproperty(plan, :raw_path)
            raw = Matrix{Float64}(plan.raw_path)
            push!(traces, PlotlyJS.scatter3d(
                x=raw[1, :],
                y=raw[2, :],
                z=raw[3, :],
                mode="lines+markers",
                line=PlotlyJS.attr(width=2, dash="dash"),
                marker=PlotlyJS.attr(size=3),
                name="case $(case_label) raw planner path",
            ))
        end
        push!(traces, PlotlyJS.scatter3d(
            x=[plan.case.start_rtn[1]],
            y=[plan.case.start_rtn[2]],
            z=[plan.case.start_rtn[3]],
            mode="markers",
            marker=PlotlyJS.attr(size=6, color="rgb(45,150,90)", symbol="circle"),
            name="case $(case_label) start",
        ))
        push!(traces, PlotlyJS.scatter3d(
            x=[plan.case.goal_rtn[1]],
            y=[plan.case.goal_rtn[2]],
            z=[plan.case.goal_rtn[3]],
            mode="markers",
            marker=PlotlyJS.attr(size=7, color="rgb(190,60,55)", symbol="diamond"),
            name="case $(case_label) goal",
        ))
    end
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="$(label) Failed Case Paths ($(length(failed_indices)))",
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis_title="radial (m)",
                yaxis_title="along-track (m)",
                zaxis_title="cross-track (m)",
            ),
        ),
    )
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

"""Create cost-versus-iteration traces for planners that report histories."""
function rpo_comparison_cost_iteration_plot(batch; planner)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    haskey(batch.plans_by_planner, planner_type) ||
        throw(ArgumentError("No planner results found for $(planner_type)."))

    label = rpo_comparison_planner_label(planner_type)
    traces = PlotlyJS.GenericTrace[]
    for plan in batch.plans_by_planner[planner_type]
        hasproperty(plan, :cost_history) || continue
        costs = [Float64(c) for c in plan.cost_history if isfinite(Float64(c))]
        isempty(costs) && continue
        iters = rpo_comparison_cost_iteration_xvalues(plan, costs)
        push!(traces, PlotlyJS.scatter(
            x=iters,
            y=costs,
            mode=length(costs) == 1 ? "markers" : "lines+markers",
            name="case $(plan.case.label)",
            line=PlotlyJS.attr(width=2),
            marker=PlotlyJS.attr(size=7),
            hovertemplate="case $(plan.case.label)<br>iteration %{x}<br>cost %{y:.6g}<extra></extra>",
        ))
    end

    isempty(traces) && push!(traces, PlotlyJS.scatter(
        x=Int[],
        y=Float64[],
        mode="markers",
        name="no finite cost history",
    ))

    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="$(label) Cost vs Iteration",
            height=620,
            width=920,
            hovermode="closest",
            margin=PlotlyJS.attr(l=75, r=35, t=75, b=65),
            plot_bgcolor="rgb(250,250,250)",
            xaxis=PlotlyJS.attr(
                title="Iteration",
                zeroline=true,
                gridcolor="rgb(225,225,225)",
            ),
            yaxis=PlotlyJS.attr(
                title="Planner Cost",
                rangemode="tozero",
                gridcolor="rgb(225,225,225)",
            ),
        ),
    )
end

"""Return a filesystem-safe token for generated planner-comparison artifact names."""
function rpo_comparison_artifact_slug(value)
    token = replace(lowercase(String(value)), r"[^a-z0-9_=-]+" => "_")
    token = strip(token, '_')
    return isempty(token) ? "case" : token
end

"""Write one path HTML per planner containing all failed planner/case results."""
function rpo_write_failed_path_outputs(batch, output_dir::AbstractString)
    failed_dir = joinpath(output_dir, "rpo_failed_path_plots")
    failed_path_plots = Dict{Symbol, String}()
    for planner in batch.planner_types
        results = batch.results_by_planner[planner]
        failed_indices = [idx for idx in eachindex(results) if !results[idx].success]
        isempty(failed_indices) && continue
        mkpath(failed_dir)
        path = joinpath(failed_dir, "rpo_failed_paths_$(planner).html")
        PlotlyJS.savefig(rpo_comparison_failed_paths_plot(batch; planner=planner, failed_indices=failed_indices), path)
        failed_path_plots[planner] = path
    end
    return failed_path_plots
end

"""Write planner comparison CSVs and plots to the requested output directory."""
function rpo_write_planner_comparison_outputs(
    batch;
    output_dir::AbstractString=batch.config.output_dir,
    write_plotly_outputs::Bool=batch.config.write_plotly_outputs,
    write_failed_path_outputs::Bool=batch.config.write_failed_path_outputs,
)
    mkpath(output_dir)
    rows = rpo_flatten_planner_results(batch)
    csv_path = rpo_write_namedtuple_csv(joinpath(output_dir, "rpo_planner_comparison_results.csv"), rows)

    metric_plot_path = nothing
    path_plot_path = nothing
    method_path_plots = Dict{Symbol, String}()
    cost_iteration_plots = Dict{Symbol, String}()
    failed_path_plots = Dict{Symbol, String}()

    if write_plotly_outputs
        metric_plot_path = joinpath(output_dir, "rpo_planner_comparison_metrics.html")
        PlotlyJS.savefig(rpo_comparison_metric_summary_plot(batch), metric_plot_path)

        path_plot_path = joinpath(output_dir, "rpo_planner_comparison_paths.html")
        PlotlyJS.savefig(rpo_comparison_path_family_plot(batch), path_plot_path)

        for planner in batch.planner_types
            method_path = joinpath(output_dir, "rpo_planner_paths_$(planner).html")
            PlotlyJS.savefig(rpo_comparison_path_family_plot(batch; planner=planner), method_path)
            method_path_plots[planner] = method_path
        end

        for planner in batch.planner_types
            cost_path = joinpath(output_dir, "rpo_planner_cost_vs_iteration_$(planner).html")
            PlotlyJS.savefig(rpo_comparison_cost_iteration_plot(batch; planner=planner), cost_path)
            cost_iteration_plots[planner] = cost_path
        end
    end
    if write_failed_path_outputs
        failed_path_plots = rpo_write_failed_path_outputs(batch, output_dir)
    end

    return (
        csv=csv_path,
        metrics_plot=metric_plot_path,
        path_plot=path_plot_path,
        method_path_plots=method_path_plots,
        cost_iteration_plots=cost_iteration_plots,
        failed_path_plots=failed_path_plots,
        plotly_outputs=write_plotly_outputs,
        failed_path_outputs=write_failed_path_outputs,
    )
end

end # module RPOComparisonReports
