include("CampaignPaths.jl")

module ProductionBatch
include(normpath(joinpath(@__DIR__, "..", "..",
    "Output_2_Simulation_absolute_RN", "Code", "run_grid.jl")))
end

function run_campaign(campaign, root)
    campaign in ("range_along", "no_laser") || error("Unknown comparison campaign: $campaign")
    code = campaign == "range_along" ? CampaignPaths.CODE_ROOT : CampaignPaths.NO_LASER_CODE
    ProductionBatch.run_grid(root; campaign,
        source_hashes_fn=() -> CampaignPaths.campaign_sources(ProductionBatch.StudyConfig.source_hashes),
        simulation_script=joinpath(code, "run_case.jl"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: run_grid.jl [OUTPUT_ROOT]")
    run_campaign("range_along", isempty(ARGS) ? CampaignPaths.OUTPUT : abspath(ARGS[1]))
end
