module ControlHooks
    using ..Structure

    using ..ConfigTypes: ODEParams
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractThrusterModel
    using ..GuidanceHooks: AerobrakingCampaignPropulsiveManeuverGuidanceModel
    using ..GuidanceHooks: AerobrakingGuidanceInput, dispatch_aerobraking_guidance
    using ..AerobrakingPolicy: AerobrakingPolicyConfig, DefaultAerobrakingPolicySelector
    using ..ref_sys
    using ..DynamicEffectors: BaseThrusterModel
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics

    const config = Structure

    export calcControlForceTorque, calcControlEffect!, calcControlMassFlowRate
    export control_solarpanels_heatload, control_solarpanels_heatrate, control_solarpanels_openloop
    export asim_ctrl, asim_ctrl_plot

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    include(joinpath(@__DIR__, "propulsive_maneuvers.jl"))
    include(joinpath(@__DIR__, "aerobraking", "control_commands.jl"))
    include(joinpath(@__DIR__, "aerobraking", "constraint_tracking.jl"))
    include(joinpath(@__DIR__, "aerobraking", "tracking_executor.jl"))
end
