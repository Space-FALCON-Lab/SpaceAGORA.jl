const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using SpaceAGORA
using SpaceAGORA.TelemetryVerification

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_cygnss_estimator_validation_cli(copy(ARGS))
end
