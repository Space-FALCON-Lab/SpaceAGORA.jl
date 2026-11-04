#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

include("plot_all_results.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    options = AbsoluteRNPlots.OE.parse_options(vcat(["--output", AbsoluteRNPlots.OUTPUT], ARGS))
    rows, thresholds = AbsoluteRNPlots.OE.load_grid(options.input, options.windows; include_rtn=true)
    AbsoluteRNPlots.render_selected_altitudes(rows, thresholds, options.windows, options.output)
    println("Saved $(4 * length(options.windows)) selected-altitude figures in PNG and PDF, with matching CSVs.")
end
