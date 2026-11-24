module CampaignPaths

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const LIBRARY_ROOT = joinpath(PROJECT_ROOT, "9_SpaceAGORA.jl-main_full integration")
const BASE_CODE = joinpath(PROJECT_ROOT, "Output_2_Simulation_absolute_RN", "Code")
const CODE_ROOT = @__DIR__
const OUTPUT = normpath(joinpath(CODE_ROOT, "..", "Output"))
const NO_LASER_CODE = joinpath(PROJECT_ROOT, "Output_4_no_laser", "Code")
const NO_LASER_OUTPUT = normpath(joinpath(NO_LASER_CODE, "..", "Output"))
const DIFFERENCE_CODE = joinpath(PROJECT_ROOT, "Output_5_Enable_Disable_Laser_compare", "Code")
const DIFFERENCE_OUTPUT = normpath(joinpath(DIFFERENCE_CODE, "..", "Output"))
const PAPER_CODE = joinpath(PROJECT_ROOT, "Results_For_Paper")
const JOB_DIRECTORY = joinpath(OUTPUT, "_comparison_job")

campaign_sources(hash_function) = hash_function(; additional_roots=(
    CODE_ROOT, NO_LASER_CODE, DIFFERENCE_CODE, PAPER_CODE))

end
