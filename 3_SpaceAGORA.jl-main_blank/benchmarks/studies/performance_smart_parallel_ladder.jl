include(joinpath(@__DIR__, "performance_smart_parallel_ladder", "cli.jl"))
include(joinpath(@__DIR__, "performance_smart_parallel_ladder", "rungs.jl"))
include(joinpath(@__DIR__, "performance_smart_parallel_ladder", "execution.jl"))
include(joinpath(@__DIR__, "performance_smart_parallel_ladder", "statistics.jl"))
include(joinpath(@__DIR__, "performance_smart_parallel_ladder", "accuracy_analysis.jl"))
include(joinpath(@__DIR__, "performance_smart_parallel_ladder", "reporting.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_smart_parallel_ladder()
end
