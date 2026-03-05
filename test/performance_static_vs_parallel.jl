const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "performance_static_vs_parallel.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_static_vs_parallel()
end
