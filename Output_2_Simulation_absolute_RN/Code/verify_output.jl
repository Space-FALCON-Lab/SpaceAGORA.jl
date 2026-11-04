include("VerifyOutput.jl")
using .StudyVerification

length(ARGS) == 1 || error("Usage: julia verify_output.jl CASE_DIRECTORY")
verify_output(abspath(only(ARGS)))
