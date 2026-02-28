"""
    Wrapper module  for all dynamic effector models (all forces/torques)
"""
module DynamicEffectors
    using ..Analysis

    using ..ConfigTypes: ODEParams # Get the Planet struct
    using ..AbstractTypes: AbstractPlanet, AbstractForceTorqueModel, AbstractThrusterModel, AbstractGuidanceModel
    using ..ParallelPolicy
    using ..LinearAlgebra       # Get deps from parent
    using ..StaticArrays        # Get deps from parent
    using ..Kinematics

    # Public members to export
    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel # Gravity models
    export NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel # N-body gravity + SRP models
    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight # Aerodynamic models
    export calcForceTorque
    
    include("../utils/Reference_system.jl")
    include("../physical_models/Gravity_models.jl")
    include("../physical_models/Aerodynamic_models.jl")
    include("../physical_models/Perturbations.jl")
    
    # Control forces/torques
    export BaseThrusterModel
    include("../physical_models/Thruster_models.jl")

    # Guidance effectors
    export AerobrakingCampaignPropulsiveManeuverGuidanceModel
    include("../guidance/Thruster_guidance_models.jl")
end
