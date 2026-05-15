module GuidanceHooks
    using ..Structure

    using ..ConfigTypes: ODEParams, Solution
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractGuidanceModel
    using ..GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel
    using ..GuidanceModels: RPOGuidanceModel, RPOPlan, RPOPlanBuffer, update_rpo_plan_buffer!
    using ..NavigationHooks: rpo_clearance_distance_to_station, rpo_clearance_to_station, rpo_path_clearance_stats
    using ..CommandTypes: PropulsiveManeuverCommand, AerobrakingControlCommand
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..AerobrakingPolicy: AbstractAerobrakingPolicySelector, AerobrakingPolicyConfig, E_EDG, T_EDG, select_strategy
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using Base.Threads: @threads, maxthreadid, threadid
    using Random

    const config = Structure
    const _PARENT = parentmodule(@__MODULE__)

    @inline function _control_module()
        return getfield(_PARENT, :ControlHooks)
    end
    @inline function _control_asim_ctrl(args...; kwargs...)
        return _control_module().asim_ctrl(args...; kwargs...)
    end
    @inline function _control_asim_ctrl_rf(args...; kwargs...)
        return _control_module().asim_ctrl_rf(args...; kwargs...)
    end
    @inline function _control_solarpanels_heatrate(args...; kwargs...)
        return _control_module().control_solarpanels_heatrate(args...; kwargs...)
    end

    export calcGuidanceEffect!
    export RPOPSOAdaptiveSettings, RPOPSOConfig, RPOPSOConfigurator, RPOPSOCullSettings
    export RPOPSOGIFSettings, RPOPSOObjectiveSettings
    export RPOPSOProbeSettings, RPOPSOReexploreSettings, RPOPSORefinementSettings
    export RPOPSORetimingSettings, RPOPSOScheduleSettings, RPOPSOSwarmSettings
    export RPOPSOStagnationSettings, rpo_pso_config
    export rpo_pso_plan_path, rpo_reference_from_path
    export AbstractAerobrakingStrategy, EEdgStrategy, TEdgStrategy
    export AerobrakingGuidanceInput, AerobrakingGuidanceOutput, AerobrakingControlCommand
    export compute_aerobraking_guidance, dispatch_aerobraking_guidance

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "hypr.jl"))
    include(joinpath(@__DIR__, "rpo", "rpo_reference_trajectory.jl"))
    include(joinpath(@__DIR__, "rpo", "rpo_guidance_hooks.jl"))
    include(joinpath(@__DIR__, "aerobraking", "interfaces.jl"))
    include(joinpath(@__DIR__, "aerobraking", "common", "closed_form_solution.jl"))
    include(joinpath(@__DIR__, "aerobraking", "common", "heat_rate_models.jl"))
    include(joinpath(@__DIR__, "aerobraking", "t_edg", "trajectory_predictor.jl"))
    include(joinpath(@__DIR__, "aerobraking", "t_edg", "eom_predictor.jl"))
    include(joinpath(@__DIR__, "aerobraking", "t_edg", "targeting_solver.jl"))
    include(joinpath(@__DIR__, "aerobraking", "e_edg", "switch_window_solver.jl"))
    include(joinpath(@__DIR__, "aerobraking", "e_edg", "second_switch_solver.jl"))
    include(joinpath(@__DIR__, "aerobraking", "e_edg", "energy_profile_solver.jl"))
    include(joinpath(@__DIR__, "aerobraking", "e_edg", "e_edg_strategy.jl"))
    include(joinpath(@__DIR__, "aerobraking", "t_edg", "t_edg_strategy.jl"))
    include(joinpath(@__DIR__, "aerobraking", "dispatcher.jl"))

    include(joinpath(@__DIR__, "thruster_guidance", "thruster_guidance_functions.jl"))
end
