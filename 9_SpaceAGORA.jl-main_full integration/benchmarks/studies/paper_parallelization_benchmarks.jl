include(joinpath(@__DIR__, "parallelization_performance", "cli.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "modes.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "cases.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "trajectory_parity.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "reporting.jl"))
include(joinpath(@__DIR__, "parallelization_performance", "execution.jl"))
include(joinpath(@__DIR__, "paper_parallelization_benchmarks", "cli.jl"))
include(joinpath(@__DIR__, "paper_parallelization_benchmarks", "reporting.jl"))
include(joinpath(@__DIR__, "paper_parallelization_benchmarks", "main.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    try
        main_paper_benchmarks()
    catch err
        err isa ArgumentError || rethrow()
        println(stderr, "Error: ", err.msg)
        exit(1)
    end
end
