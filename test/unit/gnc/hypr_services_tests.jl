module HYPRServiceContractTests
using Test, SpaceAGORA
const C=SpaceAGORA.HYPRServices
@testset "HYPR versioned service identity" begin
    @test C.CONTRACT_VERSION == v"1.0.0"
    @test C.SUPPORTED_HYPR_VERSIONS == (v"0.1.0",v"0.1.1")
    @test_throws C.CompatibilityError C.check_provider(:HYPR,v"0.1.2",v"1.0.0")
    @test C.require_version(v"1.0.0") === nothing
    @test_throws SpaceAGORA.SimulationModel.HYPRSupport.HYPRCompatibilityError C.require_version(v"2.0.0")
    @test_throws SpaceAGORA.SimulationModel.HYPRSupport.HYPRCompatibilityError C.check_provider(:HYPR,v"9.0.0",v"1.0.0")
    @test C.SwarmPolicy.hypr_iteration_weights === SpaceAGORA.SimulationModel.HYPRUtils.hypr_iteration_weights
    @test C.SwarmPolicy.hypr_material_improvement === SpaceAGORA.SimulationModel.HYPRUtils.hypr_material_improvement
    @test C.SwarmPolicy.hypr_protected_particle_mask === SpaceAGORA.SimulationModel.HYPRUtils.hypr_protected_particle_mask
    @test C.RPO.RPOAdaptiveSamplingSettings === SpaceAGORA.SimulationModel.GuidanceHooks.RPOAdaptiveSamplingSettings
    @test C.RPO.RPOPSOConfig === SpaceAGORA.SimulationModel.GuidanceHooks.RPOPSOConfig
    @test C.RPO.RPORRTConnectSettings === SpaceAGORA.SimulationModel.GuidanceHooks.RPORRTConnectSettings
    @test C.RPO.RPORRTStarSettings === SpaceAGORA.SimulationModel.GuidanceHooks.RPORRTStarSettings
    @test C.RPO.RPORetimeCurve === SpaceAGORA.SimulationModel.GuidanceHooks.RPORetimeCurve
    @test C.RPO._rpo_rrt_configured_result === SpaceAGORA.SimulationModel.GuidanceHooks._rpo_rrt_configured_result
    @test C.RPO._rpo_sampling_settings === SpaceAGORA.SimulationModel.GuidanceHooks._rpo_sampling_settings
    @test C.RPO.rpo_adaptive_pso_config === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_adaptive_pso_config
    @test C.RPO.rpo_adaptive_sampling_min_ds_m === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_adaptive_sampling_min_ds_m
    @test C.RPO.rpo_adaptive_segment_samples === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_adaptive_segment_samples
    @test C.RPO.rpo_arc_length_params === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_arc_length_params
    @test C.RPO.rpo_clearance_stats_from_samples === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_clearance_stats_from_samples
    @test C.RPO.rpo_estimate_geometry_complexity === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_estimate_geometry_complexity
    @test C.RPO.rpo_fit_bezier_fixed_endpoints === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_fit_bezier_fixed_endpoints
    @test C.RPO.rpo_fuel_proxy_dt_s === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_fuel_proxy_dt_s
    @test C.RPO.rpo_fuel_proxy_from_samples === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_fuel_proxy_from_samples
    @test C.RPO.rpo_hcw_fuel_proxy === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_hcw_fuel_proxy
    @test C.RPO.rpo_hypr_refinement_sampling_density_m === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_hypr_refinement_sampling_density_m
    @test C.RPO.rpo_hypr_sampling_density_m === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_hypr_sampling_density_m
    @test C.RPO.rpo_manuscript_adaptive_pso_config === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_manuscript_adaptive_pso_config
    @test C.RPO.rpo_manuscript_exploration_score === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_manuscript_exploration_score
    @test C.RPO.rpo_manuscript_path_cost_components === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_manuscript_path_cost_components
    @test C.RPO.rpo_normalized_path_cost_components === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_normalized_path_cost_components
    @test C.RPO.rpo_obstacle_sigmoid_penalty === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_obstacle_sigmoid_penalty
    @test C.RPO.rpo_obstacle_sigmoid_threshold === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_obstacle_sigmoid_threshold
    @test C.RPO.rpo_path_cost === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_path_cost
    @test C.RPO.rpo_path_cost_normalization_refs === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_path_cost_normalization_refs
    @test C.RPO.rpo_path_length === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_path_length
    @test C.RPO.rpo_position_to_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_position_to_path
    @test C.RPO.rpo_post_refine_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_post_refine_path
    @test C.RPO.rpo_probe_geometry_metrics === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_probe_geometry_metrics
    @test C.RPO.rpo_profile_hcw_fuel_proxy === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_profile_hcw_fuel_proxy
    @test C.RPO.rpo_pso_bounds === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_bounds
    @test C.RPO.rpo_pso_config === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_config
    @test C.RPO.rpo_pso_cull_swarm! === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_cull_swarm!
    @test C.RPO.rpo_pso_early_stopping_feasible === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_early_stopping_feasible
    @test C.RPO.rpo_pso_effective_safe_distance === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_effective_safe_distance
    @test C.RPO.rpo_pso_empty_warmstart_diagnostics === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_empty_warmstart_diagnostics
    @test C.RPO.rpo_pso_iteration_weights === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_iteration_weights
    @test C.RPO.rpo_pso_material_improvement === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_material_improvement
    @test C.RPO.rpo_pso_plan_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_plan_path
    @test C.RPO.rpo_pso_project_to_segment === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_project_to_segment
    @test C.RPO.rpo_pso_protected_particle_mask === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_protected_particle_mask
    @test C.RPO.rpo_pso_rrt_warmstart_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_rrt_warmstart_path
    @test C.RPO.rpo_pso_stagnation_count_after_learning === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_stagnation_count_after_learning
    @test C.RPO.rpo_pso_station_bounds === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_station_bounds
    @test C.RPO.rpo_pso_tapered_noise_scale === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_tapered_noise_scale
    @test C.RPO.rpo_pso_warmstart_bounds === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_pso_warmstart_bounds
    @test C.RPO.rpo_reference_from_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_reference_from_path
    @test C.RPO.rpo_refine_lower_degree === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refine_lower_degree
    @test C.RPO.rpo_refine_shortcut_refit === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refine_shortcut_refit
    @test C.RPO.rpo_refine_tighten_handles === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refine_tighten_handles
    @test C.RPO.rpo_refinement_bernstein === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_bernstein
    @test C.RPO.rpo_refinement_better === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_better
    @test C.RPO.rpo_refinement_clamp_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_clamp_path
    @test C.RPO.rpo_refinement_config === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_config
    @test C.RPO.rpo_refinement_project_to_segment === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_project_to_segment
    @test C.RPO.rpo_refinement_sample_params === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_sample_params
    @test C.RPO.rpo_refinement_segment_is_safe === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_segment_is_safe
    @test C.RPO.rpo_refinement_segment_samples === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_segment_samples
    @test C.RPO.rpo_refinement_shortcut_samples === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_refinement_shortcut_samples
    @test C.RPO.rpo_resample_polyline_points === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_resample_polyline_points
    @test C.RPO.rpo_retime_available_distance === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retime_available_distance
    @test C.RPO.rpo_retime_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retime_path
    @test C.RPO.rpo_retime_pointwise_speed === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retime_pointwise_speed
    @test C.RPO.rpo_retime_profile === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retime_profile
    @test C.RPO.rpo_retime_samples === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retime_samples
    @test C.RPO.rpo_retime_sampling_ds_m === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retime_sampling_ds_m
    @test C.RPO.rpo_retimed_reference === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retimed_reference
    @test C.RPO.rpo_retimed_reference_from_profile === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_retimed_reference_from_profile
    @test C.RPO.rpo_rrt_connect_bezier_plan_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_rrt_connect_bezier_plan_path
    @test C.RPO.rpo_rrt_connect_plan_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_rrt_connect_plan_path
    @test C.RPO.rpo_rrt_star_plan_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_rrt_star_plan_path
    @test C.RPO.rpo_sample_path === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_sample_path
    @test C.RPO.rpo_sample_path_bezier === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_sample_path_bezier
    @test C.RPO.rpo_sample_path_bezier_adaptive === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_sample_path_bezier_adaptive
    @test C.RPO.rpo_sample_path_bezier_adaptive_with_params === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_sample_path_bezier_adaptive_with_params
    @test C.RPO.rpo_sample_path_polyline === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_sample_path_polyline
    @test C.RPO.rpo_sample_path_polyline_adaptive === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_sample_path_polyline_adaptive
    @test C.RPO.rpo_sample_path_with_params === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_sample_path_with_params
    @test C.RPO.rpo_try_accept_refinement === SpaceAGORA.SimulationModel.GuidanceHooks.rpo_try_accept_refinement
    @test C.RPO.validate_rpo_pso_config === SpaceAGORA.SimulationModel.GuidanceHooks.validate_rpo_pso_config
    @test C.RPO.hypr_iteration_weights === SpaceAGORA.SimulationModel.HYPRUtils.hypr_iteration_weights
    @test C.RPO.hypr_material_improvement === SpaceAGORA.SimulationModel.HYPRUtils.hypr_material_improvement
    @test C.RPO.hypr_protected_particle_mask === SpaceAGORA.SimulationModel.HYPRUtils.hypr_protected_particle_mask
    @test C.RPO.rpo_clearance_distance_to_station === SpaceAGORA.SimulationModel.NavigationHooks.rpo_clearance_distance_to_station
    @test C.RPO.rpo_path_clearance_stats === SpaceAGORA.SimulationModel.NavigationHooks.rpo_path_clearance_stats
    @test C.RPO.RPOReferenceGeometry === SpaceAGORA.SimulationModel.NavigationHooks.RPOReferenceGeometry
    @test C.RobotArm.RobotArmHYPRConfig === SpaceAGORA.SimulationModel.RobotArmPlanning.RobotArmHYPRConfig
    @test C.RobotArm.RobotArmHYPRResult === SpaceAGORA.SimulationModel.RobotArmPlanning.RobotArmHYPRResult
    @test C.RobotArm.RobotArmPlan === SpaceAGORA.SimulationModel.RobotArmPlanning.RobotArmPlan
    @test C.RobotArm.RobotArmPlannerConfig === SpaceAGORA.SimulationModel.RobotArmPlanning.RobotArmPlannerConfig
    @test C.RobotArm.RobotArmRRTConnectTree === SpaceAGORA.SimulationModel.RobotArmPlanning.RobotArmRRTConnectTree
    @test C.RobotArm.RobotArmSphereObstacle === SpaceAGORA.SimulationModel.RobotArmPlanning.RobotArmSphereObstacle
    @test C.RobotArm._reference_times === SpaceAGORA.SimulationModel.RobotArmPlanning._reference_times
    @test C.RobotArm._robot_arm_control_points === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_control_points
    @test C.RobotArm._robot_arm_empty_rrt_warmstart_diagnostics === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_empty_rrt_warmstart_diagnostics
    @test C.RobotArm._robot_arm_flatten_internal_points === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_flatten_internal_points
    @test C.RobotArm._robot_arm_hypr_base_wrench_ratios === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_base_wrench_ratios
    @test C.RobotArm._robot_arm_hypr_cloth_base_wrench_ratios === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_cloth_base_wrench_ratios
    @test C.RobotArm._robot_arm_hypr_cloth_state_for_reaction === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_cloth_state_for_reaction
    @test C.RobotArm._robot_arm_hypr_cull_swarm! === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_cull_swarm!
    @test C.RobotArm._robot_arm_hypr_early_stopping_feasible === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_early_stopping_feasible
    @test C.RobotArm._robot_arm_hypr_iteration_weights === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_iteration_weights
    @test C.RobotArm._robot_arm_hypr_link_com_history === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_link_com_history
    @test C.RobotArm._robot_arm_hypr_material_improvement === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_material_improvement
    @test C.RobotArm._robot_arm_hypr_post_refine_points === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_post_refine_points
    @test C.RobotArm._robot_arm_hypr_reaction_scale === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_reaction_scale
    @test C.RobotArm._robot_arm_hypr_reference_times_from_scales === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_reference_times_from_scales
    @test C.RobotArm._robot_arm_hypr_refinement_better === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_refinement_better
    @test C.RobotArm._robot_arm_hypr_retime_reference === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_retime_reference
    @test C.RobotArm._robot_arm_hypr_rigid_base_wrench_ratios === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_hypr_rigid_base_wrench_ratios
    @test C.RobotArm._robot_arm_path_length === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_path_length
    @test C.RobotArm._robot_arm_path_smoothness === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_path_smoothness
    @test C.RobotArm._robot_arm_plan_from_q_reference === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_plan_from_q_reference
    @test C.RobotArm._robot_arm_resample_polyline_points === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_resample_polyline_points
    @test C.RobotArm._robot_arm_rrt_connect! === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_connect!
    @test C.RobotArm._robot_arm_rrt_connect_warmstart_path === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_connect_warmstart_path
    @test C.RobotArm._robot_arm_rrt_extend! === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_extend!
    @test C.RobotArm._robot_arm_rrt_join_paths === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_join_paths
    @test C.RobotArm._robot_arm_rrt_nearest_index === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_nearest_index
    @test C.RobotArm._robot_arm_rrt_path_score === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_path_score
    @test C.RobotArm._robot_arm_rrt_random_state === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_random_state
    @test C.RobotArm._robot_arm_rrt_segment_is_safe === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_segment_is_safe
    @test C.RobotArm._robot_arm_rrt_segment_samples === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_segment_samples
    @test C.RobotArm._robot_arm_rrt_shortcut_path === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_shortcut_path
    @test C.RobotArm._robot_arm_rrt_steer === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_steer
    @test C.RobotArm._robot_arm_rrt_warmstart_fields === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_rrt_warmstart_fields
    @test C.RobotArm._robot_arm_seed_control_points === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_seed_control_points
    @test C.RobotArm._robot_arm_segment_distance === SpaceAGORA.SimulationModel.RobotArmPlanning._robot_arm_segment_distance
    @test C.RobotArm._validate_robot_arm_hypr_config === SpaceAGORA.SimulationModel.RobotArmPlanning._validate_robot_arm_hypr_config
    @test C.RobotArm.plan_robot_arm_motion === SpaceAGORA.SimulationModel.RobotArmPlanning.plan_robot_arm_motion
    @test C.RobotArm.plan_robot_arm_motion_hypr === SpaceAGORA.SimulationModel.RobotArmPlanning.plan_robot_arm_motion_hypr
    @test C.RobotArm.robot_arm_clearance_stats_from_samples === SpaceAGORA.SimulationModel.RobotArmPlanning.robot_arm_clearance_stats_from_samples
    @test C.RobotArm.robot_arm_hypr_path_cost_components === SpaceAGORA.SimulationModel.RobotArmPlanning.robot_arm_hypr_path_cost_components
    @test C.RobotArm.robot_arm_sample_hypr_path === SpaceAGORA.SimulationModel.RobotArmPlanning.robot_arm_sample_hypr_path
    @test C.RobotArm.hypr_iteration_weights === SpaceAGORA.SimulationModel.HYPRUtils.hypr_iteration_weights
    @test C.RobotArm.hypr_material_improvement === SpaceAGORA.SimulationModel.HYPRUtils.hypr_material_improvement
    @test C.RobotArm.hypr_rrt_join_paths === SpaceAGORA.SimulationModel.HYPRUtils.hypr_rrt_join_paths
    @test C.RobotArm.hypr_rrt_nearest_index === SpaceAGORA.SimulationModel.HYPRUtils.hypr_rrt_nearest_index
    @test C.RobotArm.hypr_rrt_steer === SpaceAGORA.SimulationModel.HYPRUtils.hypr_rrt_steer
    @test C.RobotArm.ClothArmBasePose === SpaceAGORA.SimulationModel.Robotics.ClothArmBasePose
    @test C.RobotArm.ClothArmModel === SpaceAGORA.SimulationModel.Robotics.ClothArmModel
    @test C.RobotArm.cloth_fk === SpaceAGORA.SimulationModel.Robotics.cloth_fk
    @test C.RobotArm.cloth_ik === SpaceAGORA.SimulationModel.Robotics.cloth_ik
    @test C.Planner.RPOPlanningRequest === SpaceAGORA.RPOPlannerInterfaces.RPOPlanningRequest
    @test C.Planner.RPOPlanningResult === SpaceAGORA.RPOPlannerInterfaces.RPOPlanningResult
    @test C.Planner.RPOReference === SpaceAGORA.RPOPlannerInterfaces.RPOReference
    @test C.Planner.rpo_planning_budget === SpaceAGORA.RPOPlannerInterfaces.rpo_planning_budget
    @test C.Planner.validate_rpo_result === SpaceAGORA.RPOPlannerInterfaces.validate_rpo_result
    @test C.Planner.HYPRRPOPlanner === SpaceAGORA.HYPRRPOPlanning.HYPRRPOPlanner
    @test C.Planner._plan_hypr_rpo! === SpaceAGORA.HYPRRPOPlanning._plan_hypr_rpo!
    @test C.Cloth.ClothRobotArmDynamics.simulate_cloth_robot_arm_plan === SpaceAGORA.SimulationModel.ClothRobotArmDynamics.simulate_cloth_robot_arm_plan
    @test C.Cloth.ClothRobotArmDynamics.assign_coupled_cloth_robot_arm_rhs! === SpaceAGORA.SimulationModel.ClothRobotArmDynamics.assign_coupled_cloth_robot_arm_rhs!
    @test C.Cloth.ClothMultibody.compliant_state_parts === SpaceAGORA.SimulationModel.ClothMultibody.compliant_state_parts
    for (group, name) in C.IMPLEMENTED_FUNCTIONS
        @test getfield(getfield(C,group),name) isa Function
    end
    @test isempty(intersect(Set(C.IMPLEMENTED_FUNCTIONS),Set(C.CONSUMED_SERVICES)))
    for group in (:SwarmPolicy, :RPO, :RobotArm, :Planner)
        declared = Set(name for (g,name) in (C.IMPLEMENTED_FUNCTIONS..., C.CONSUMED_SERVICES...) if g == group)
        exposed = setdiff(Set(names(getfield(C, group))), Set([group]))
        @test declared == exposed
    end
    @test C.CLOTH_SERVICES == (
        (:ClothRobotArmDynamics, :simulate_cloth_robot_arm_plan),
        (:ClothRobotArmDynamics, :assign_coupled_cloth_robot_arm_rhs!),
        (:ClothMultibody, :compliant_state_parts),
    )
end
# Use an isolated support module so these refusal tests do not poison other suites.
const Fixture=Module(:HYPRProviderFixture)
Base.include(Fixture,joinpath(@__DIR__,"..","..","..","src","gnc","hypr","support.jl"))
const Support=Fixture.HYPRSupport
@testset "One runtime provider, fail closed after conflict" begin
    @test !Support.hypr_available()
    @test Support.activate!(:HYPR,v"0.1.0") === nothing
    @test Support.activate!(:HYPR,v"0.1.0") === nothing
    @test Support.hypr_available()
    @test_throws Support.HYPRCompatibilityError Support.activate!(:Other,v"0.1.0")
    @test !Support.hypr_available()
    @test_throws Support.HYPRCompatibilityError Support.require_hypr()
    @test_throws Support.HYPRCompatibilityError Support.activate!(:HYPR,v"0.1.0")
end
const OldFixture=Module(:OldHYPRProviderFixture)
Base.include(OldFixture,joinpath(@__DIR__,"..","..","..","src","gnc","hypr","support.jl"))
@testset "Old activation protocol is unsupported" begin
    @test_throws OldFixture.HYPRSupport.HYPRCompatibilityError OldFixture.HYPRSupport.activate!()
    @test !OldFixture.HYPRSupport.hypr_available()
end
end
