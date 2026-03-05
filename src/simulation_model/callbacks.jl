# Compatibility wrapper: legacy callbacks path forwards to canonical split callbacks module.
include(joinpath(@__DIR__, "..", "simulation", "callbacks", "callbacks.jl"))
