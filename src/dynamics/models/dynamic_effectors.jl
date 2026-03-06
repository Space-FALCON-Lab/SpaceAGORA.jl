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
    export srp, srp_cannonball_accel
    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight # Aerodynamic models
    export calcForceTorque
    
    include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
    include(joinpath(@__DIR__, "..", "..", "environment", "gravity", "gravity_models.jl"))
    include(joinpath(@__DIR__, "..", "..", "dynamics", "translational", "aerodynamic_models.jl"))
    include(joinpath(@__DIR__, "..", "..", "dynamics", "coupled", "perturbations.jl"))
    
    # Control forces/torques
    export BaseThrusterModel
    include(joinpath(@__DIR__, "..", "..", "dynamics", "models", "thruster_models.jl"))

    # Guidance effectors
    export AerobrakingCampaignPropulsiveManeuverGuidanceModel
    include(joinpath(@__DIR__, "..", "..", "gnc", "guidance", "thruster_guidance", "thruster_guidance_models.jl"))
end
