const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "analysis", "verification", "TelemetryVerification.jl"))
using .TelemetryVerification

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_verification_cli(copy(ARGS))
end
