include("ComparisonResults.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: plot_campaign.jl [OUTPUT_ROOT]")
    root = isempty(ARGS) ? ComparisonResults.CampaignPaths.OUTPUT : abspath(ARGS[1])
    ComparisonResults.render_campaign("range_along", root)
end
