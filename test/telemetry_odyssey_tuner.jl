const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "telemetry_odyssey_tuner.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_tuner()
end
