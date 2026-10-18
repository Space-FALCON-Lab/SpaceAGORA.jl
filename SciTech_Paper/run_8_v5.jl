const INTEGRATION_ROOT = normpath(joinpath(@__DIR__, "..", "8_SpaceAGORA.jl-main_v5"))
const EXAMPLES_DIR = "II_examples"
const PACKAGE_ENTRY = joinpath("I_src", "6_SpaceAGORA.jl")
const RESULT_CODE = "8_v5"
include(joinpath(@__DIR__, "native_interlink_comparison.jl"))

abspath(PROGRAM_FILE) == (@__FILE__) && run_case()
