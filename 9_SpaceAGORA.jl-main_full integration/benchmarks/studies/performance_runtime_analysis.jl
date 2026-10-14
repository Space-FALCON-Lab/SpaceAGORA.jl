# Canonical launcher for runtime-analysis study; implementation is split by responsibility.
include(joinpath(@__DIR__, "performance_runtime_analysis", "main.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
