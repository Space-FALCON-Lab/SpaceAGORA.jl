# Data contracts and controller-independent algebra.
include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "constraints.jl"))
include(joinpath(@__DIR__, "condensed_core.jl"))
include(joinpath(@__DIR__, "objectives.jl"))

# SpaceAGORA environment adapters and scenario construction.
include(joinpath(@__DIR__, "density_adapter.jl"))
include(joinpath(@__DIR__, "scenario_config.jl"))

# Reference propagation, local prediction model, and QP solution.
include(joinpath(@__DIR__, "trajectory_reference.jl"))
include(joinpath(@__DIR__, "prediction_model.jl"))
include(joinpath(@__DIR__, "qp_solver.jl"))

# Command conversion, validation, and simulation callbacks.
include(joinpath(@__DIR__, "area_mapping.jl"))
include(joinpath(@__DIR__, "plan_validation.jl"))
include(joinpath(@__DIR__, "spaceagora_control_model.jl"))
include(joinpath(@__DIR__, "campaign_control_model.jl"))
