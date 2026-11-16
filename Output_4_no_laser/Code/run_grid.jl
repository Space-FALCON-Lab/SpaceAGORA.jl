include(normpath(joinpath(@__DIR__, "..", "..",
    "Output_3_rangeOnly_vs_positiveAlongTrack", "Code", "run_grid.jl")))

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: run_grid.jl [OUTPUT_ROOT]")
    run_campaign("no_laser", isempty(ARGS) ? CampaignPaths.NO_LASER_OUTPUT : abspath(ARGS[1]))
end
