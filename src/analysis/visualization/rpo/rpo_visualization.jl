module RPOVisualization

export rpo_path_plot, rpo_tracking_plot

const _PLOTLYJS_MODULE = Ref{Union{Nothing,Module}}(nothing)
const _PLOTLYJS_LOAD_LOCK = ReentrantLock()

function _plotlyjs()
    return lock(_PLOTLYJS_LOAD_LOCK) do
        if _PLOTLYJS_MODULE[] === nothing
            _PLOTLYJS_MODULE[] = Base.require(
                Base.PkgId(
                    Base.UUID("f0f68f2c-4968-5e81-91da-67840de0976a"),
                    "PlotlyJS",
                ),
            )
        end
        return _PLOTLYJS_MODULE[]::Module
    end
end

function rpo_path_plot(path_rtn, station_points=nothing; title::AbstractString="RPO Path")
    PlotlyJS = _plotlyjs()
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
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title=title,
            scene=PlotlyJS.attr(aspectmode="data"),
        ),
    )
end

function rpo_tracking_plot(t, err_norm; title::AbstractString="RPO Tracking Error")
    PlotlyJS = _plotlyjs()
    return PlotlyJS.Plot(
        [
            PlotlyJS.scatter(
                x=t,
                y=err_norm,
                mode="lines",
                name="tracking error",
            ),
        ],
        PlotlyJS.Layout(
            title=title,
            xaxis_title="time (s)",
            yaxis_title="error norm (m)",
        ),
    )
end

end # module RPOVisualization
