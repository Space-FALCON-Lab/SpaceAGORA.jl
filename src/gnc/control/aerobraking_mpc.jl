#=
"""
    File loading order for the aerobraking MPC.

    Keep this order because the later files depend on the type definitions,
    constraint flags, KS utilities, and condensed QP helpers.
"""
=#
include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "constraints.jl"))
include(joinpath(@__DIR__, "ks_transform.jl"))
include(joinpath(@__DIR__, "condensed_core.jl"))
include(joinpath(@__DIR__, "objectives.jl"))
include(joinpath(@__DIR__, "scenario_config.jl"))
include(joinpath(@__DIR__, "ks_dynamics.jl"))
include(joinpath(@__DIR__, "trajectory_reference.jl"))
include(joinpath(@__DIR__, "ks_linearization.jl"))
include(joinpath(@__DIR__, "qp_solver.jl"))
include(joinpath(@__DIR__, "area_mapping.jl"))
include(joinpath(@__DIR__, "spaceagora_control_model.jl"))
