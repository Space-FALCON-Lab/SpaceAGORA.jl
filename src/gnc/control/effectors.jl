module ControlEffectors
    using ..Analysis

    using ..ConfigTypes: ODEParams # Get the Planet struct
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractThrusterModel
    using ..LinearAlgebra       # Get deps from parent
    using ..StaticArrays        # Get deps from parent
    using ..Kinematics
    using ..DynamicEffectors: BaseThrusterModel
    using ..GuidanceEffectors: AerobrakingCampaignPropulsiveManeuverGuidanceModel

    # Public members to export
    export calcControlForceTorque, calcControlEffect!, calcControlMassFlowRate

    include(joinpath(@__DIR__, "propulsive_maneuvers.jl"))
end
