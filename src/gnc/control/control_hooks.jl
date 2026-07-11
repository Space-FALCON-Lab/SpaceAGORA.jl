module ControlHooks
    using ..Structure

    using ..ConfigTypes: ODEParams
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractThrusterModel
    using ..ThrusterModels: BaseThrusterModel
    using ..CommandTypes: PropulsiveManeuverCommand, PropulsiveBurnPlan
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..GuidanceHooks: AerobrakingGuidanceInput, dispatch_aerobraking_guidance
    using ..GuidanceHooks: AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
    using ..AerobrakingPolicy: AerobrakingPolicyConfig, DefaultAerobrakingPolicySelector
    using ..EnvironmentModels: getDensity
    using ..EphemeridesModels: ephemerides_requires_spice, planet_frame_lpi
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics

    const config = Structure

    export calcControlForceTorque, calcControlEffect!, calcControlMassFlowRate
    export AerobrakingEnergyDepletionControlModel, SolarPanelAngleOfAttackControlModel

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    include(joinpath(@__DIR__, "propulsive_maneuvers.jl"))
    include(joinpath(@__DIR__, "heat_rate_control.jl"))
    include(joinpath(@__DIR__, "heat_load_control.jl"))
    include(joinpath(@__DIR__, "struct_load_control.jl"))
    include(joinpath(@__DIR__, "targeting_control.jl"))
    include(joinpath(@__DIR__, "aerobraking", "control_commands.jl"))
    include(joinpath(@__DIR__, "aerobraking", "constraint_tracking.jl"))
    include(joinpath(@__DIR__, "aerobraking", "tracking_executor.jl"))
end
