const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "performance_smart_parallel_ladder_cross_machine.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_cross_machine_replay()
end
