module GuidanceHooks
    using ..Structure

    using ..ConfigTypes: ODEParams, Solution
    using ..AbstractTypes: AbstractPlanet, AbstractControlEffectorModel, AbstractGuidanceModel
    using ..GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel
    using ..GuidanceModels: ApoapsisTargetPeriapsisRaiseGuidanceModel
    using ..GuidanceModels: RPOGuidanceModel, RPOPlan, RPOPlanBuffer, update_rpo_plan_buffer!
    using ..NavigationHooks: RPOReferenceGeometry, RPOStationGeometry
    using ..NavigationHooks: rpo_clearance_distance_to_station, rpo_clearance_to_station, rpo_path_clearance_stats
    using ..CommandTypes: PropulsiveManeuverCommand, AerobrakingControlCommand
    using ..EphemeridesModels: planet_frame_lpi
    using ..GravityEffectors: aerobraking_gravity_force_ii
    using ..AerobrakingPolicy: AbstractAerobrakingPolicySelector, AerobrakingPolicyConfig, E_EDG, T_EDG, select_strategy
    using ..ReferenceSystems
    using ..LinearAlgebra
    using ..StaticArrays
    using ..Kinematics
    using ..HYPRUtils
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
    export RPOPSOEarlyStoppingSettings, RPOPSOObjectiveSettings
    export RPOPSOProbeSettings, RPOPSOReexploreSettings, RPOPSORefinementSettings
    export RPOPSORetimingSettings, RPOPSOScheduleSettings, RPOPSOSwarmSettings
    export RPOPSORRTConnectWarmstartSettings, RPOPSOStagnationSettings, rpo_pso_config
    export rpo_path_objective_components, rpo_pso_plan_path, rpo_reference_from_path
    export RPOReplanningConfig, RPOReplanningSphere
    export rpo_active_replanning_spheres, rpo_geometry_with_replanning_spheres
    export rpo_reference_tracking_error, rpo_remaining_reference_path
    export rpo_replanning_decision, rpo_replanning_sphere_center
    export RPOCHOMPSettings, RPORRTConnectSettings, RPORRTStarSettings, RPOSTOMPSettings, RPOTrajectoryOptimizerSettings
    export RPOLQMPCTrackingSettings, RPOPlannerComparisonCase, RPOPlannerComparisonConfig
    export normalize_rpo_comparison_planner_type, rpo_comparison_planner_label
    export rpo_740_mpc_final_pso_config, rpo_chomp_plan_path, rpo_rrt_connect_plan_path, rpo_rrt_connect_bezier_plan_path, rpo_rrt_star_plan_path, rpo_stomp_plan_path
    export rpo_plan_comparison_path, rpo_run_planner_comparison_batch
    export rpo_track_retimed_path_lqmpc, rpo_flatten_planner_results
    export rpo_comparison_cost_iteration_plot, rpo_write_planner_comparison_outputs
    export AbstractAerobrakingStrategy, EEdgStrategy, TEdgStrategy
    export AerobrakingGuidanceInput, AerobrakingGuidanceOutput, AerobrakingControlCommand
    export compute_aerobraking_guidance, dispatch_aerobraking_guidance
    export AerobrakingEnergyDepletionConfig, AerobrakingEnergyDepletionState
    export AerobrakingEnergyDepletionGuidanceModel

    include(joinpath(@__DIR__, "..", "internal", "bridge_helpers.jl"))
    include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "pso_parameters.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "path_retiming.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "path_sampling.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "path_costs.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "pso_adaptive_policy.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "pso_helpers.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "pso_refinement.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "rrt_connect.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr", "pso_path_planning.jl"))
    include(joinpath(@__DIR__, "rpo", "rpo_reference_trajectory.jl"))
    include(joinpath(@__DIR__, "rpo", "comparison_methods", "trajectory_optimizers.jl"))
    include(joinpath(@__DIR__, "rpo", "comparison_methods", "planner_comparison.jl"))
    include(joinpath(@__DIR__, "rpo", "hypr_planning", "replanning.jl"))
    include(joinpath(@__DIR__, "rpo", "rpo_guidance_hooks.jl"))
    include(joinpath(@__DIR__, "aerobraking", "interfaces.jl"))
    include(joinpath(@__DIR__, "target_energy_bracketing.jl"))
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
