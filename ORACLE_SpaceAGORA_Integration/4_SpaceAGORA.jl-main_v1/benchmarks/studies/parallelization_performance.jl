include(joinpath(@__DIR__, "parallelization_performance", "cli.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "modes.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "cases.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "trajectory_parity.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "reporting.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "execution.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_parallelization_performance()
end
