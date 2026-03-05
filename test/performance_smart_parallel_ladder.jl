const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "performance_smart_parallel_ladder.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_smart_parallel_ladder()
end
