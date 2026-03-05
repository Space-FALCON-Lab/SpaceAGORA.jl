const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "benchmarks", "studies", "performance_effector_reduction_microbench.jl"))
