const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "gram_real_sim_surrogate_decision_table.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end
