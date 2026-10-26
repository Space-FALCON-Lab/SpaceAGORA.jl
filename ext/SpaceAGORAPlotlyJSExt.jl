module SpaceAGORAPlotlyJSExt

# The RPO Plotly figures: the planner-comparison plots and reports
# (RPOComparisonReports) and `rpo_path_plot` and `rpo_tracking_plot`
# (RPOVisualization). Core keeps every public method and calls the builders
# below through the renderer this extension registers when PlotlyJS loads.

using SpaceAGORA
import PlotlyJS

const Reports = SpaceAGORA.SimulationModel.GuidanceHooks.RPOComparisonReports
const RPOVisualization = SpaceAGORA.RPOVisualization

using .Reports: normalize_rpo_comparison_planner_type, rpo_comparison_planner_label,
    rpo_comparison_axis_name, rpo_comparison_trace_axis_name,
    rpo_comparison_station_mesh_trace, rpo_comparison_station_trace,
    _rpo_comparison_metric_series, _rpo_comparison_station_mesh, _rpo_comparison_has_station_mesh,
    _rpo_comparison_station_point_indices, _rpo_comparison_case_path, _rpo_comparison_cost_series

struct PlotlyJSRenderer end

Reports._save_figure(::PlotlyJSRenderer, figure, path::AbstractString) = PlotlyJS.savefig(figure, path)

# Create summary bar plots for planner comparison metrics.
function Reports._metric_summary_figure(::PlotlyJSRenderer, batch)
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

    for (idx, series) in enumerate(_rpo_comparison_metric_series(batch))
        label = series.label
        units = series.units
        row_idx = div(idx - 1, cols) + 1
        col_idx = mod(idx - 1, cols) + 1
        x0 = (col_idx - 1) * (xspan + xgap)
        x1 = x0 + xspan
        y1 = 1.0 - (row_idx - 1) * (yspan + ygap)
        y0 = y1 - yspan
        xaxis = rpo_comparison_trace_axis_name("x", idx)
        yaxis = rpo_comparison_trace_axis_name("y", idx)

        push!(traces, PlotlyJS.box(
            x=series.planners,
            y=series.values,
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
            x=series.planners,
            y=series.values,
            mode="markers",
            name="$(label) cases",
            marker=PlotlyJS.attr(size=8, color="rgb(205,92,55)", line=PlotlyJS.attr(color="white", width=1)),
            text=series.case_labels,
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

# Build a Plotly mesh trace for the station keepout geometry.
function Reports._station_mesh_trace(::PlotlyJSRenderer, triangles)
    mesh = _rpo_comparison_station_mesh(triangles)
    return PlotlyJS.mesh3d(
        x=mesh.vertices[1, :],
        y=mesh.vertices[2, :],
        z=mesh.vertices[3, :],
        i=mesh.i,
        j=mesh.j,
        k=mesh.k,
        color="rgb(150,160,170)",
        opacity=0.55,
        flatshading=true,
        name="Gateway mesh",
        hoverinfo="skip",
    )
end

# Build a downsampled Plotly station point-cloud trace.
function Reports._station_trace(::PlotlyJSRenderer, batch; max_points::Integer)
    if _rpo_comparison_has_station_mesh(batch)
        return rpo_comparison_station_mesh_trace(batch.station_triangles)
    end
    pts = batch.geometry.station.points_body
    idxs = _rpo_comparison_station_point_indices(pts, max_points)
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

# Create a 3D Plotly view of planner paths for all comparison cases.
function Reports._path_family_figure(::PlotlyJSRenderer, batch; planner)
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    planners = planner === nothing ? batch.planner_types : [normalize_rpo_comparison_planner_type(planner)]
    for planner_type in planners
        label = rpo_comparison_planner_label(planner_type)
        for plan in batch.plans_by_planner[planner_type]
            path = _rpo_comparison_case_path(plan)
            case_label = path.case_label
            ref = path.desired
            actual = path.tracked
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
                x=[path.start[1]],
                y=[path.start[2]],
                z=[path.start[3]],
                mode="markers",
                marker=PlotlyJS.attr(size=6, color="rgb(45,150,90)", symbol="circle"),
                name="$(label) case $(case_label) start",
            ))
            push!(traces, PlotlyJS.scatter3d(
                x=[path.goal[1]],
                y=[path.goal[2]],
                z=[path.goal[3]],
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

# Create a 3D Plotly view for one planner/case path.
function Reports._single_path_figure(::PlotlyJSRenderer, batch, plan; planner)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    label = rpo_comparison_planner_label(planner_type)
    path = _rpo_comparison_case_path(plan)
    case_label = path.case_label
    ref = path.desired
    actual = path.tracked
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
    if path.has_raw
        raw = Matrix{Float64}(path.raw)
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
        x=[path.start[1]],
        y=[path.start[2]],
        z=[path.start[3]],
        mode="markers",
        marker=PlotlyJS.attr(size=7, color="rgb(45,150,90)", symbol="circle"),
        name="start",
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=[path.goal[1]],
        y=[path.goal[2]],
        z=[path.goal[3]],
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

# Create a 3D Plotly view containing all failed paths for one planner.
function Reports._failed_paths_figure(::PlotlyJSRenderer, batch; planner, failed_indices)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    label = rpo_comparison_planner_label(planner_type)
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    plans = batch.plans_by_planner[planner_type]
    for idx in failed_indices
        path = _rpo_comparison_case_path(plans[idx])
        case_label = path.case_label
        ref = path.desired
        actual = path.tracked
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
        if path.has_raw
            raw = Matrix{Float64}(path.raw)
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
            x=[path.start[1]],
            y=[path.start[2]],
            z=[path.start[3]],
            mode="markers",
            marker=PlotlyJS.attr(size=6, color="rgb(45,150,90)", symbol="circle"),
            name="case $(case_label) start",
        ))
        push!(traces, PlotlyJS.scatter3d(
            x=[path.goal[1]],
            y=[path.goal[2]],
            z=[path.goal[3]],
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

# Create cost-versus-iteration traces for planners that report histories.
function Reports._cost_iteration_figure(::PlotlyJSRenderer, batch; planner)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    haskey(batch.plans_by_planner, planner_type) ||
        throw(ArgumentError("No planner results found for $(planner_type)."))

    label = rpo_comparison_planner_label(planner_type)
    traces = PlotlyJS.GenericTrace[]
    for series in _rpo_comparison_cost_series(batch, planner_type)
        push!(traces, PlotlyJS.scatter(
            x=series.iterations,
            y=series.costs,
            mode=length(series.costs) == 1 ? "markers" : "lines+markers",
            name="case $(series.case_label)",
            line=PlotlyJS.attr(width=2),
            marker=PlotlyJS.attr(size=7),
            hovertemplate="case $(series.case_label)<br>iteration %{x}<br>cost %{y:.6g}<extra></extra>",
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

# RPOVisualization's figures.
function RPOVisualization._rpo_path_figure(::PlotlyJSRenderer, path_rtn, station_points; title::AbstractString)
    traces = PlotlyJS.GenericTrace[]
    path = Matrix{Float64}(path_rtn)
    push!(traces, PlotlyJS.scatter3d(
        x=path[1, :],
        y=path[2, :],
        z=path[3, :],
        mode="lines+markers",
        name="planned path",
    ))
    if station_points !== nothing
        pts = Matrix{Float64}(station_points)
        push!(traces, PlotlyJS.scatter3d(
            x=pts[1, :],
            y=pts[2, :],
            z=pts[3, :],
            mode="markers",
            marker=PlotlyJS.attr(size=3),
            name="station geometry",
        ))
    end
    return PlotlyJS.Plot(traces, PlotlyJS.Layout(title=title, scene=PlotlyJS.attr(aspectmode="data")))
end

function RPOVisualization._rpo_tracking_figure(::PlotlyJSRenderer, t, err_norm; title::AbstractString)
    return PlotlyJS.Plot(
        [PlotlyJS.scatter(x=t, y=err_norm, mode="lines", name="tracking error")],
        PlotlyJS.Layout(title=title, xaxis_title="time (s)", yaxis_title="error norm (m)"),
    )
end

__init__() = Reports._register_plotly_renderer!(PlotlyJSRenderer())

end # module SpaceAGORAPlotlyJSExt
