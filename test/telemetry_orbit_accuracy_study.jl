const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "telemetry_orbit_accuracy_study.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_verification_cli(copy(ARGS))
end
