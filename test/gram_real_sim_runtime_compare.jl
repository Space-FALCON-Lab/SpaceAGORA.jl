const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "gram_real_sim_runtime_compare.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_benchmark()
end
