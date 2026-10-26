module RPOVisualization

import ..SimulationModel.GuidanceHooks.RPOComparisonReports: _plotly_renderer

export rpo_path_plot, rpo_tracking_plot

# The PlotlyJS package extension builds both figures for its renderer. Without
# PlotlyJS loaded, both functions throw PlotlyJSUnavailableError.
function _rpo_path_figure end
function _rpo_tracking_figure end

rpo_path_plot(path_rtn, station_points=nothing; title::AbstractString="RPO Path") =
    _rpo_path_figure(_plotly_renderer(), path_rtn, station_points; title=title)

rpo_tracking_plot(t, err_norm; title::AbstractString="RPO Tracking Error") =
    _rpo_tracking_figure(_plotly_renderer(), t, err_norm; title=title)

end # module RPOVisualization
