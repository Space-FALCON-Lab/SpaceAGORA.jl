include(normpath(joinpath(@__DIR__, "..", "..",
    "Output_3_rangeOnly_vs_positiveAlongTrack", "Code", "ComparisonResults.jl")))

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (0, 3) ||
        error("Usage: plot_results.jl [LASER_OUTPUT_ROOT NO_LASER_OUTPUT_ROOT DIFFERENCE_OUTPUT_ROOT]")
    paths = ComparisonResults.CampaignPaths
    roots = isempty(ARGS) ? (paths.OUTPUT, paths.NO_LASER_OUTPUT, paths.DIFFERENCE_OUTPUT) :
        Tuple(abspath.(ARGS))
    ComparisonResults.render_differences(roots...)
end
