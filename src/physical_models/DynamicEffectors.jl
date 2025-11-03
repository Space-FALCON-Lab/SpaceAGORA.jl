"""
    Wrapper module  for all dynamic effector models (all forces/torques)
"""
module DynamicEffectors
    # Public members
    # export AbstractForceTorqueModel # Abstract type for force/torque models
    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel # Gravity models
    export NBodyGravityModel # N-body gravity model
    export calcForceTorque

    include("../simulation_model/abstract_types.jl")
    using .AbstractTypes
    include("Gravity_models.jl")
    # include("Perturbations.jl")
end