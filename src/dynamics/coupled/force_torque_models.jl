"""
    Wrapper module for all dynamic effector models (all forces/torques)
"""
module DynamicEffectors
    function calcForceTorque end

    module GravityEffectors
        using ...ConfigTypes: ODEParams
        using ...AbstractTypes: AbstractForceTorqueModel
        import ..DynamicEffectors: calcForceTorque

        export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
        export aerobraking_gravity_force_ii

        include(joinpath(@__DIR__, "..", "..", "environment", "gravity", "gravity_models.jl"))
    end

    module AerodynamicEffectors
        using ...Structure
        using ...ConfigTypes: ODEParams, AeroScratchWorkspace
        using ...AbstractTypes: AbstractForceTorqueModel
        using ...ParallelPolicy
        using ...Kinematics
        using LinearAlgebra
        using StaticArrays
        import ..DynamicEffectors: calcForceTorque

        export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight

        include(joinpath(@__DIR__, "aerodynamic_wrench_models.jl"))
    end

    module PerturbationEffectors
        using ...ConfigTypes: ODEParams, NBodyScratchWorkspace, HarmonicsScratchWorkspace
        using ...AbstractTypes: AbstractPlanet, AbstractForceTorqueModel
        using ...Planets: Earth, Mars, Venus, Titan
        using ...ParallelPolicy
        using ...SimulationModel: SPICE_LOCK, SRPSunEphemerisCache, NBodyEphemerisCache, SpiceRhsMemo
        using StaticArrays
        import ..DynamicEffectors: calcForceTorque

        export NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
        export srp, srp_cannonball_accel

        include(joinpath(@__DIR__, "perturbations.jl"))
    end

    module ThrusterModels
        using ...AbstractTypes: AbstractThrusterModel

        export BaseThrusterModel

        include(joinpath(@__DIR__, "..", "..", "vehicle", "actuators", "thruster", "thruster_models.jl"))
    end

    module GuidanceModels
        using ...AbstractTypes: AbstractGuidanceModel

        export AerobrakingCampaignPropulsiveManeuverGuidanceModel

        include(joinpath(@__DIR__, "..", "..", "gnc", "guidance", "thruster_guidance", "thruster_guidance_models.jl"))
    end

    using .GravityEffectors: ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    using .GravityEffectors: aerobraking_gravity_force_ii, calcForceTorque!
    using .AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
    using .PerturbationEffectors: NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
    using .PerturbationEffectors: srp, srp_cannonball_accel
    using .ThrusterModels: BaseThrusterModel
    using .GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel

    # Public members to export
    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    export NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
    export srp, srp_cannonball_accel
    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
    export calcForceTorque
    export BaseThrusterModel
    export AerobrakingCampaignPropulsiveManeuverGuidanceModel
end
