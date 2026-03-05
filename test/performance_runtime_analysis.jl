const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "performance_runtime_analysis.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
