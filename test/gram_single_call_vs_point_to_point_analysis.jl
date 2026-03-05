const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "gram_single_call_vs_point_to_point_analysis.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_single_call_analysis()
end
