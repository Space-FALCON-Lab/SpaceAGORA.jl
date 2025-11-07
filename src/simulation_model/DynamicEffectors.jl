"""
    Wrapper module  for all dynamic effector models (all forces/torques)
"""
module DynamicEffectors
    using ..ConfigTypes: Planet, ODEParams # Get the Planet struct
    using ..AbstractTypes       # Assuming this is also in types.jl now
    using ..LinearAlgebra       # Get deps from parent
    using ..StaticArrays        # Get deps from parent
    
    # Public members to export
    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel # Gravity models
    export NBodyGravityModel, GravitationalHarmonicsModel # N-body gravity model
    export calcForceTorque

    include("../physical_models/Gravity_models.jl")
    include("../physical_models/Perturbations.jl")
end