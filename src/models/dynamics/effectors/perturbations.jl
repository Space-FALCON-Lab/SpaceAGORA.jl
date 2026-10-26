# Ordered composition of perturbing force and torque models.
# Every include is evaluated in DynamicEffectors.PerturbationEffectors.
# Preserve this order: constants, constructors and shared caches retain their
# original initialization sequence. The LVLH controller remains a continuous
# torque effector, not a sampled guidance callback.
include(joinpath(@__DIR__, "perturbations", "dependencies.jl"))
include(joinpath(@__DIR__, "perturbations", "magnetics", "dipole_constants.jl"))
include(joinpath(@__DIR__, "perturbations", "gravity", "models_and_workspaces.jl"))
include(joinpath(@__DIR__, "perturbations", "radiation", "model.jl"))
include(joinpath(@__DIR__, "perturbations", "gravity", "third_body.jl"))
include(joinpath(@__DIR__, "perturbations", "radiation", "wrenches.jl"))
include(joinpath(@__DIR__, "perturbations", "ephemeris_interpolation.jl"))
include(joinpath(@__DIR__, "perturbations", "gravity", "harmonics.jl"))
include(joinpath(@__DIR__, "perturbations", "magnetics", "torques.jl"))
include(joinpath(@__DIR__, "perturbations", "control", "lvlh_cascade.jl"))
include(joinpath(@__DIR__, "perturbations", "radiation", "eclipse.jl"))
