include(joinpath(@__DIR__, "measurement", "routing.jl"))
using .BenchmarkRouting: partition_sweep_routes

include(joinpath(@__DIR__, "measurement", "warmup_and_batch.jl"))
include(joinpath(@__DIR__, "measurement", "montecarlo.jl"))
include(joinpath(@__DIR__, "measurement", "per_orbit.jl"))
include(joinpath(@__DIR__, "measurement", "entry_duration.jl"))
