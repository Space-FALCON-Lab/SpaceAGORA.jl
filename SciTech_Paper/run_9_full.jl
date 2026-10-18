const INTEGRATION_ROOT = normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration"))
const EXAMPLES_DIR = "examples"
const PACKAGE_ENTRY = joinpath("src", "SpaceAGORA.jl")
const RESULT_CODE = "9_full"
include(joinpath(@__DIR__, "native_interlink_comparison.jl"))

abspath(PROGRAM_FILE) == (@__FILE__) && run_case()
