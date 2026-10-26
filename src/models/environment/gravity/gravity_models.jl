# Compatibility composition entry point for analytic gravity.
# DynamicEffectors.GravityEffectors owns every definition in these includes.
# Field kernels live here; spacecraft/runtime adapters live in dynamics.
# The ordered includes preserve the original definition and initialization order.
include(joinpath(@__DIR__, "..", "..", "dynamics", "effectors", "gravity", "models.jl"))
include(joinpath(@__DIR__, "encke_and_j2.jl"))
include(joinpath(@__DIR__, "..", "..", "dynamics", "effectors", "gravity", "runtime_force_helpers.jl"))
include(joinpath(@__DIR__, "central_acceleration.jl"))
include(joinpath(@__DIR__, "..", "..", "dynamics", "effectors", "gravity", "sampled_gravity_gradient.jl"))
include(joinpath(@__DIR__, "secular_rates.jl"))
include(joinpath(@__DIR__, "..", "..", "dynamics", "effectors", "gravity", "wrenches.jl"))
include(joinpath(@__DIR__, "gravity_gradient.jl"))
include(joinpath(@__DIR__, "..", "..", "dynamics", "effectors", "gravity", "gravity_gradient_effector.jl"))
