include(normpath(joinpath(@__DIR__, "..", "..",
    "Output_3_rangeOnly_vs_positiveAlongTrack", "Code", "ComparisonResults.jl")))

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: plot_results.jl [OUTPUT_ROOT]")
    root = isempty(ARGS) ? ComparisonResults.CampaignPaths.NO_LASER_OUTPUT : abspath(ARGS[1])
    ComparisonResults.render_campaign("no_laser", root)
end
