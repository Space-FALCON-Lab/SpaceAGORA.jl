# Simulation-dependent adapters live outside the importable report module.
include(joinpath(@__DIR__, "measurement", "rollout_solution.jl"))
include(joinpath(@__DIR__, "reporting", "plots_and_reporting.jl"))

include(joinpath(@__DIR__, "reporting", "runtime_adapters.jl"))
