const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "gram_real_sim_surrogate_matrix.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_matrix()
end
