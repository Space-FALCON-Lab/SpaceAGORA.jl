module RPOVisualization

using PlotlyJS

export rpo_path_plot, rpo_tracking_plot

function rpo_path_plot(path_rtn, station_points=nothing; title::AbstractString="RPO Path")
    traces = PlotlyJS.GenericTrace[]
    path = Matrix{Float64}(path_rtn)
    push!(traces, scatter3d(
        x=path[1, :],
        y=path[2, :],
        z=path[3, :],
        mode="lines+markers",
        name="planned path",
    ))
    if station_points !== nothing
        pts = Matrix{Float64}(station_points)
        push!(traces, scatter3d(
            x=pts[1, :],
            y=pts[2, :],
            z=pts[3, :],
            mode="markers",
            marker=attr(size=3),
            name="station geometry",
        ))
    end
    return Plot(traces, Layout(title=title, scene=attr(aspectmode="data")))
end

function rpo_tracking_plot(t, err_norm; title::AbstractString="RPO Tracking Error")
    return Plot(
        [scatter(x=t, y=err_norm, mode="lines", name="tracking error")],
        Layout(title=title, xaxis_title="time (s)", yaxis_title="error norm (m)"),
    )
end

end # module RPOVisualization

