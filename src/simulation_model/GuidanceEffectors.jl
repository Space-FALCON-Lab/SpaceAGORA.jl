module GuidanceEffectors
    using ..Analysis

    using ..ConfigTypes: ODEParams # Get the Planet struct
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractGuidanceModel
    using ..DynamicEffectors: AerobrakingCampaignPropulsiveManeuverGuidanceModel, BaseThrusterModel
    using ..LinearAlgebra       # Get deps from parent
    using ..StaticArrays        # Get deps from parent
    using ..Kinematics

    # Public members to export
    export calcGuidanceEffect!
    include("../guidance/Thruster_guidance_functions.jl") 
end