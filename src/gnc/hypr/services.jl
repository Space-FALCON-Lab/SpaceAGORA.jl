"""
    SpaceAGORA.HYPRServices

Versioned services and extension points for the optional HYPR package. Bindings
alias their existing owners; this facade does not duplicate implementation or
change public type identity. See the HYPR service-contract guide for the explicit
inventory and compatibility policy. Underscore names listed here are supported
only through this versioned boundary, not through arbitrary internal modules.
"""
module HYPRServices
import ..SimulationModel.HYPRSupport
const CONTRACT_VERSION = v"1.0.0"
const HYPR_VERSION = v"0.1.0"
"""Compatibility error type, preserving its existing owner."""
const CompatibilityError = HYPRSupport.HYPRCompatibilityError

"""Require exactly the interface version implemented by this SpaceAGORA build."""
function require_version(expected::VersionNumber)
    expected == CONTRACT_VERSION || throw(HYPRSupport.HYPRCompatibilityError(
        "HYPR requires SpaceAGORA service contract $expected; this build provides $CONTRACT_VERSION."))
    return nothing
end

"""Check compatibility before the optional extension defines planner methods."""
function check_provider(provider::Symbol, version::VersionNumber, contract::VersionNumber)
    require_version(contract)
    provider === :HYPR && version == HYPR_VERSION || throw(HYPRSupport.HYPRCompatibilityError(
        "Unsupported HYPR implementation/version; this contract supports HYPR $HYPR_VERSION."))
    HYPRSupport.check_provider(provider, version)
    return nothing
end

"""Activate one compatible implementation after its extension has initialized."""
function activate!(provider::Symbol, version::VersionNumber, contract::VersionNumber)
    check_provider(provider, version, contract)
    HYPRSupport.activate!(provider, version)
    return nothing
end

module SwarmPolicy
import ...SimulationModel.HYPRUtils:
    hypr_iteration_weights,
    hypr_material_improvement,
    hypr_protected_particle_mask
export hypr_iteration_weights, hypr_material_improvement, hypr_protected_particle_mask
end

module RPO
import ...SimulationModel.GuidanceHooks:
    RPOAdaptiveSamplingSettings,
    RPOPSOConfig,
    RPORRTConnectSettings,
    RPORRTStarSettings,
    RPORetimeCurve,
    _rpo_rrt_configured_result,
    _rpo_sampling_settings,
    rpo_adaptive_pso_config,
    rpo_adaptive_sampling_min_ds_m,
    rpo_adaptive_segment_samples,
    rpo_arc_length_params,
    rpo_clearance_stats_from_samples,
    rpo_estimate_geometry_complexity,
    rpo_fit_bezier_fixed_endpoints,
    rpo_fuel_proxy_dt_s,
    rpo_fuel_proxy_from_samples,
    rpo_hcw_fuel_proxy,
    rpo_hypr_refinement_sampling_density_m,
    rpo_hypr_sampling_density_m,
    rpo_manuscript_adaptive_pso_config,
    rpo_manuscript_exploration_score,
    rpo_manuscript_path_cost_components,
    rpo_normalized_path_cost_components,
    rpo_obstacle_sigmoid_penalty,
    rpo_obstacle_sigmoid_threshold,
    rpo_path_cost,
    rpo_path_cost_normalization_refs,
    rpo_path_length,
    rpo_position_to_path,
    rpo_post_refine_path,
    rpo_probe_geometry_metrics,
    rpo_profile_hcw_fuel_proxy,
    rpo_pso_bounds,
    rpo_pso_config,
    rpo_pso_cull_swarm!,
    rpo_pso_early_stopping_feasible,
    rpo_pso_effective_safe_distance,
    rpo_pso_empty_warmstart_diagnostics,
    rpo_pso_iteration_weights,
    rpo_pso_material_improvement,
    rpo_pso_plan_path,
    rpo_pso_project_to_segment,
    rpo_pso_protected_particle_mask,
    rpo_pso_rrt_warmstart_path,
    rpo_pso_stagnation_count_after_learning,
    rpo_pso_station_bounds,
    rpo_pso_tapered_noise_scale,
    rpo_pso_warmstart_bounds,
    rpo_reference_from_path,
    rpo_refine_lower_degree,
    rpo_refine_shortcut_refit,
    rpo_refine_tighten_handles,
    rpo_refinement_bernstein,
    rpo_refinement_better,
    rpo_refinement_clamp_path,
    rpo_refinement_config,
    rpo_refinement_project_to_segment,
    rpo_refinement_sample_params,
    rpo_refinement_segment_is_safe,
    rpo_refinement_segment_samples,
    rpo_refinement_shortcut_samples,
    rpo_resample_polyline_points,
    rpo_retime_available_distance,
    rpo_retime_path,
    rpo_retime_pointwise_speed,
    rpo_retime_profile,
    rpo_retime_samples,
    rpo_retime_sampling_ds_m,
    rpo_retimed_reference,
    rpo_retimed_reference_from_profile,
    rpo_rrt_connect_bezier_plan_path,
    rpo_rrt_connect_plan_path,
    rpo_rrt_star_plan_path,
    rpo_sample_path,
    rpo_sample_path_bezier,
    rpo_sample_path_bezier_adaptive,
    rpo_sample_path_bezier_adaptive_with_params,
    rpo_sample_path_polyline,
    rpo_sample_path_polyline_adaptive,
    rpo_sample_path_with_params,
    rpo_try_accept_refinement,
    validate_rpo_pso_config
import ...SimulationModel.HYPRUtils:
    hypr_iteration_weights,
    hypr_material_improvement,
    hypr_protected_particle_mask
import ...SimulationModel.NavigationHooks:
    rpo_clearance_distance_to_station,
    rpo_path_clearance_stats,
    RPOReferenceGeometry
export RPOAdaptiveSamplingSettings, RPOPSOConfig, RPORRTConnectSettings, RPORRTStarSettings, RPORetimeCurve, _rpo_rrt_configured_result, _rpo_sampling_settings, rpo_adaptive_pso_config, rpo_adaptive_sampling_min_ds_m, rpo_adaptive_segment_samples, rpo_arc_length_params, rpo_clearance_stats_from_samples, rpo_estimate_geometry_complexity, rpo_fit_bezier_fixed_endpoints, rpo_fuel_proxy_dt_s, rpo_fuel_proxy_from_samples, rpo_hcw_fuel_proxy, rpo_hypr_refinement_sampling_density_m, rpo_hypr_sampling_density_m, rpo_manuscript_adaptive_pso_config, rpo_manuscript_exploration_score, rpo_manuscript_path_cost_components, rpo_normalized_path_cost_components, rpo_obstacle_sigmoid_penalty, rpo_obstacle_sigmoid_threshold, rpo_path_cost, rpo_path_cost_normalization_refs, rpo_path_length, rpo_position_to_path, rpo_post_refine_path, rpo_probe_geometry_metrics, rpo_profile_hcw_fuel_proxy, rpo_pso_bounds, rpo_pso_config, rpo_pso_cull_swarm!, rpo_pso_early_stopping_feasible, rpo_pso_effective_safe_distance, rpo_pso_empty_warmstart_diagnostics, rpo_pso_iteration_weights, rpo_pso_material_improvement, rpo_pso_plan_path, rpo_pso_project_to_segment, rpo_pso_protected_particle_mask, rpo_pso_rrt_warmstart_path, rpo_pso_stagnation_count_after_learning, rpo_pso_station_bounds, rpo_pso_tapered_noise_scale, rpo_pso_warmstart_bounds, rpo_reference_from_path, rpo_refine_lower_degree, rpo_refine_shortcut_refit, rpo_refine_tighten_handles, rpo_refinement_bernstein, rpo_refinement_better, rpo_refinement_clamp_path, rpo_refinement_config, rpo_refinement_project_to_segment, rpo_refinement_sample_params, rpo_refinement_segment_is_safe, rpo_refinement_segment_samples, rpo_refinement_shortcut_samples, rpo_resample_polyline_points, rpo_retime_available_distance, rpo_retime_path, rpo_retime_pointwise_speed, rpo_retime_profile, rpo_retime_samples, rpo_retime_sampling_ds_m, rpo_retimed_reference, rpo_retimed_reference_from_profile, rpo_rrt_connect_bezier_plan_path, rpo_rrt_connect_plan_path, rpo_rrt_star_plan_path, rpo_sample_path, rpo_sample_path_bezier, rpo_sample_path_bezier_adaptive, rpo_sample_path_bezier_adaptive_with_params, rpo_sample_path_polyline, rpo_sample_path_polyline_adaptive, rpo_sample_path_with_params, rpo_try_accept_refinement, validate_rpo_pso_config, hypr_iteration_weights, hypr_material_improvement, hypr_protected_particle_mask, rpo_clearance_distance_to_station, rpo_path_clearance_stats, RPOReferenceGeometry
end

module RobotArm
import ...SimulationModel.RobotArmPlanning:
    RobotArmHYPRConfig,
    RobotArmHYPRResult,
    RobotArmPlan,
    RobotArmPlannerConfig,
    RobotArmRRTConnectTree,
    RobotArmSphereObstacle,
    _reference_times,
    _robot_arm_control_points,
    _robot_arm_empty_rrt_warmstart_diagnostics,
    _robot_arm_flatten_internal_points,
    _robot_arm_hypr_base_wrench_ratios,
    _robot_arm_hypr_cloth_base_wrench_ratios,
    _robot_arm_hypr_cloth_state_for_reaction,
    _robot_arm_hypr_cull_swarm!,
    _robot_arm_hypr_early_stopping_feasible,
    _robot_arm_hypr_iteration_weights,
    _robot_arm_hypr_link_com_history,
    _robot_arm_hypr_material_improvement,
    _robot_arm_hypr_post_refine_points,
    _robot_arm_hypr_reaction_scale,
    _robot_arm_hypr_reference_times_from_scales,
    _robot_arm_hypr_refinement_better,
    _robot_arm_hypr_retime_reference,
    _robot_arm_hypr_rigid_base_wrench_ratios,
    _robot_arm_path_length,
    _robot_arm_path_smoothness,
    _robot_arm_plan_from_q_reference,
    _robot_arm_resample_polyline_points,
    _robot_arm_rrt_connect!,
    _robot_arm_rrt_connect_warmstart_path,
    _robot_arm_rrt_extend!,
    _robot_arm_rrt_join_paths,
    _robot_arm_rrt_nearest_index,
    _robot_arm_rrt_path_score,
    _robot_arm_rrt_random_state,
    _robot_arm_rrt_segment_is_safe,
    _robot_arm_rrt_segment_samples,
    _robot_arm_rrt_shortcut_path,
    _robot_arm_rrt_steer,
    _robot_arm_rrt_warmstart_fields,
    _robot_arm_seed_control_points,
    _robot_arm_segment_distance,
    _validate_robot_arm_hypr_config,
    plan_robot_arm_motion,
    plan_robot_arm_motion_hypr,
    robot_arm_clearance_stats_from_samples,
    robot_arm_hypr_path_cost_components,
    robot_arm_sample_hypr_path
import ...SimulationModel.HYPRUtils:
    hypr_iteration_weights,
    hypr_material_improvement,
    hypr_rrt_join_paths,
    hypr_rrt_nearest_index,
    hypr_rrt_steer
import ...SimulationModel.Robotics:
    ClothArmBasePose,
    ClothArmModel,
    cloth_fk,
    cloth_ik
export RobotArmHYPRConfig, RobotArmHYPRResult, RobotArmPlan, RobotArmPlannerConfig, RobotArmRRTConnectTree, RobotArmSphereObstacle, _reference_times, _robot_arm_control_points, _robot_arm_empty_rrt_warmstart_diagnostics, _robot_arm_flatten_internal_points, _robot_arm_hypr_base_wrench_ratios, _robot_arm_hypr_cloth_base_wrench_ratios, _robot_arm_hypr_cloth_state_for_reaction, _robot_arm_hypr_cull_swarm!, _robot_arm_hypr_early_stopping_feasible, _robot_arm_hypr_iteration_weights, _robot_arm_hypr_link_com_history, _robot_arm_hypr_material_improvement, _robot_arm_hypr_post_refine_points, _robot_arm_hypr_reaction_scale, _robot_arm_hypr_reference_times_from_scales, _robot_arm_hypr_refinement_better, _robot_arm_hypr_retime_reference, _robot_arm_hypr_rigid_base_wrench_ratios, _robot_arm_path_length, _robot_arm_path_smoothness, _robot_arm_plan_from_q_reference, _robot_arm_resample_polyline_points, _robot_arm_rrt_connect!, _robot_arm_rrt_connect_warmstart_path, _robot_arm_rrt_extend!, _robot_arm_rrt_join_paths, _robot_arm_rrt_nearest_index, _robot_arm_rrt_path_score, _robot_arm_rrt_random_state, _robot_arm_rrt_segment_is_safe, _robot_arm_rrt_segment_samples, _robot_arm_rrt_shortcut_path, _robot_arm_rrt_steer, _robot_arm_rrt_warmstart_fields, _robot_arm_seed_control_points, _robot_arm_segment_distance, _validate_robot_arm_hypr_config, plan_robot_arm_motion, plan_robot_arm_motion_hypr, robot_arm_clearance_stats_from_samples, robot_arm_hypr_path_cost_components, robot_arm_sample_hypr_path, hypr_iteration_weights, hypr_material_improvement, hypr_rrt_join_paths, hypr_rrt_nearest_index, hypr_rrt_steer, ClothArmBasePose, ClothArmModel, cloth_fk, cloth_ik
end

module Planner
import ...RPOPlannerInterfaces:
    RPOPlanningRequest,
    RPOPlanningResult,
    RPOReference,
    rpo_planning_budget,
    validate_rpo_result
import ...HYPRRPOPlanning:
    HYPRRPOPlanner,
    _plan_hypr_rpo!
import Random: AbstractRNG
export RPOPlanningRequest, RPOPlanningResult, RPOReference, rpo_planning_budget, validate_rpo_result, HYPRRPOPlanner, _plan_hypr_rpo!, AbstractRNG
end

module Cloth
module ClothRobotArmDynamics
import ....SimulationModel.ClothRobotArmDynamics: simulate_cloth_robot_arm_plan, assign_coupled_cloth_robot_arm_rhs!
end
module ClothMultibody
import ....SimulationModel.ClothMultibody: compliant_state_parts
end
end

# Explicit versioned ownership inventory. Aliases retain original method/type owners.
const IMPLEMENTED_FUNCTIONS = (
    (:RobotArm, :_robot_arm_control_points),
    (:RobotArm, :_robot_arm_empty_rrt_warmstart_diagnostics),
    (:RobotArm, :_robot_arm_flatten_internal_points),
    (:RobotArm, :_robot_arm_hypr_base_wrench_ratios),
    (:RobotArm, :_robot_arm_hypr_cloth_base_wrench_ratios),
    (:RobotArm, :_robot_arm_hypr_cloth_state_for_reaction),
    (:RobotArm, :_robot_arm_hypr_cull_swarm!),
    (:RobotArm, :_robot_arm_hypr_early_stopping_feasible),
    (:RobotArm, :_robot_arm_hypr_iteration_weights),
    (:RobotArm, :_robot_arm_hypr_link_com_history),
    (:RobotArm, :_robot_arm_hypr_material_improvement),
    (:RobotArm, :_robot_arm_hypr_post_refine_points),
    (:RobotArm, :_robot_arm_hypr_reaction_scale),
    (:RobotArm, :_robot_arm_hypr_reference_times_from_scales),
    (:RobotArm, :_robot_arm_hypr_refinement_better),
    (:RobotArm, :_robot_arm_hypr_retime_reference),
    (:RobotArm, :_robot_arm_hypr_rigid_base_wrench_ratios),
    (:RobotArm, :_robot_arm_plan_from_q_reference),
    (:RobotArm, :_robot_arm_resample_polyline_points),
    (:RobotArm, :_robot_arm_rrt_connect!),
    (:RobotArm, :_robot_arm_rrt_connect_warmstart_path),
    (:RobotArm, :_robot_arm_rrt_extend!),
    (:RobotArm, :_robot_arm_rrt_join_paths),
    (:RobotArm, :_robot_arm_rrt_nearest_index),
    (:RobotArm, :_robot_arm_rrt_path_score),
    (:RobotArm, :_robot_arm_rrt_random_state),
    (:RobotArm, :_robot_arm_rrt_segment_is_safe),
    (:RobotArm, :_robot_arm_rrt_segment_samples),
    (:RobotArm, :_robot_arm_rrt_shortcut_path),
    (:RobotArm, :_robot_arm_rrt_steer),
    (:RobotArm, :_robot_arm_rrt_warmstart_fields),
    (:RobotArm, :_robot_arm_seed_control_points),
    (:RobotArm, :hypr_iteration_weights),
    (:RobotArm, :hypr_material_improvement),
    (:RobotArm, :plan_robot_arm_motion_hypr),
    (:RobotArm, :robot_arm_hypr_path_cost_components),
    (:SwarmPolicy, :hypr_iteration_weights),
    (:SwarmPolicy, :hypr_material_improvement),
    (:SwarmPolicy, :hypr_protected_particle_mask),
    (:Planner, :_plan_hypr_rpo!),
    (:RPO, :_rpo_rrt_configured_result),
    (:RPO, :_rpo_sampling_settings),
    (:RPO, :hypr_iteration_weights),
    (:RPO, :hypr_material_improvement),
    (:RPO, :hypr_protected_particle_mask),
    (:RPO, :rpo_adaptive_pso_config),
    (:RPO, :rpo_adaptive_sampling_min_ds_m),
    (:RPO, :rpo_estimate_geometry_complexity),
    (:RPO, :rpo_fit_bezier_fixed_endpoints),
    (:RPO, :rpo_fuel_proxy_dt_s),
    (:RPO, :rpo_fuel_proxy_from_samples),
    (:RPO, :rpo_manuscript_adaptive_pso_config),
    (:RPO, :rpo_manuscript_exploration_score),
    (:RPO, :rpo_manuscript_path_cost_components),
    (:RPO, :rpo_normalized_path_cost_components),
    (:RPO, :rpo_path_cost),
    (:RPO, :rpo_path_cost_normalization_refs),
    (:RPO, :rpo_position_to_path),
    (:RPO, :rpo_post_refine_path),
    (:RPO, :rpo_probe_geometry_metrics),
    (:RPO, :rpo_pso_bounds),
    (:RPO, :rpo_pso_cull_swarm!),
    (:RPO, :rpo_pso_early_stopping_feasible),
    (:RPO, :rpo_pso_effective_safe_distance),
    (:RPO, :rpo_pso_empty_warmstart_diagnostics),
    (:RPO, :rpo_pso_iteration_weights),
    (:RPO, :rpo_pso_material_improvement),
    (:RPO, :rpo_pso_plan_path),
    (:RPO, :rpo_pso_project_to_segment),
    (:RPO, :rpo_pso_protected_particle_mask),
    (:RPO, :rpo_pso_rrt_warmstart_path),
    (:RPO, :rpo_pso_stagnation_count_after_learning),
    (:RPO, :rpo_pso_station_bounds),
    (:RPO, :rpo_pso_tapered_noise_scale),
    (:RPO, :rpo_pso_warmstart_bounds),
    (:RPO, :rpo_reference_from_path),
    (:RPO, :rpo_refine_lower_degree),
    (:RPO, :rpo_refine_shortcut_refit),
    (:RPO, :rpo_refine_tighten_handles),
    (:RPO, :rpo_refinement_bernstein),
    (:RPO, :rpo_refinement_better),
    (:RPO, :rpo_refinement_clamp_path),
    (:RPO, :rpo_refinement_config),
    (:RPO, :rpo_refinement_project_to_segment),
    (:RPO, :rpo_refinement_sample_params),
    (:RPO, :rpo_refinement_segment_is_safe),
    (:RPO, :rpo_refinement_segment_samples),
    (:RPO, :rpo_refinement_shortcut_samples),
    (:RPO, :rpo_retime_available_distance),
    (:RPO, :rpo_retime_path),
    (:RPO, :rpo_retime_pointwise_speed),
    (:RPO, :rpo_retime_profile),
    (:RPO, :rpo_retime_sampling_ds_m),
    (:RPO, :rpo_retimed_reference),
    (:RPO, :rpo_rrt_connect_bezier_plan_path),
    (:RPO, :rpo_rrt_connect_plan_path),
    (:RPO, :rpo_rrt_star_plan_path),
    (:RPO, :rpo_sample_path),
    (:RPO, :rpo_sample_path_bezier_adaptive),
    (:RPO, :rpo_sample_path_bezier_adaptive_with_params),
    (:RPO, :rpo_sample_path_polyline_adaptive),
    (:RPO, :rpo_sample_path_with_params),
    (:RPO, :rpo_try_accept_refinement),
)
const CONSUMED_SERVICES = (
    (:RobotArm, :ClothArmBasePose),
    (:RobotArm, :ClothArmModel),
    (:RobotArm, :RobotArmHYPRConfig),
    (:RobotArm, :RobotArmHYPRResult),
    (:RobotArm, :RobotArmPlan),
    (:RobotArm, :RobotArmPlannerConfig),
    (:RobotArm, :RobotArmRRTConnectTree),
    (:RobotArm, :RobotArmSphereObstacle),
    (:RobotArm, :_reference_times),
    (:RobotArm, :_robot_arm_path_length),
    (:RobotArm, :_robot_arm_path_smoothness),
    (:RobotArm, :_robot_arm_segment_distance),
    (:RobotArm, :_validate_robot_arm_hypr_config),
    (:RobotArm, :cloth_fk),
    (:RobotArm, :cloth_ik),
    (:RobotArm, :hypr_rrt_join_paths),
    (:RobotArm, :hypr_rrt_nearest_index),
    (:RobotArm, :hypr_rrt_steer),
    (:RobotArm, :plan_robot_arm_motion),
    (:RobotArm, :robot_arm_clearance_stats_from_samples),
    (:RobotArm, :robot_arm_sample_hypr_path),
    (:Planner, :HYPRRPOPlanner),
    (:Planner, :RPOPlanningRequest),
    (:Planner, :RPOPlanningResult),
    (:Planner, :RPOReference),
    (:Planner, :rpo_planning_budget),
    (:Planner, :validate_rpo_result),
    (:Planner, :AbstractRNG),
    (:RPO, :RPOAdaptiveSamplingSettings),
    (:RPO, :RPOPSOConfig),
    (:RPO, :RPORRTConnectSettings),
    (:RPO, :RPORRTStarSettings),
    (:RPO, :RPOReferenceGeometry),
    (:RPO, :RPORetimeCurve),
    (:RPO, :rpo_adaptive_segment_samples),
    (:RPO, :rpo_arc_length_params),
    (:RPO, :rpo_clearance_distance_to_station),
    (:RPO, :rpo_clearance_stats_from_samples),
    (:RPO, :rpo_hcw_fuel_proxy),
    (:RPO, :rpo_hypr_refinement_sampling_density_m),
    (:RPO, :rpo_hypr_sampling_density_m),
    (:RPO, :rpo_obstacle_sigmoid_penalty),
    (:RPO, :rpo_obstacle_sigmoid_threshold),
    (:RPO, :rpo_path_clearance_stats),
    (:RPO, :rpo_path_length),
    (:RPO, :rpo_profile_hcw_fuel_proxy),
    (:RPO, :rpo_pso_config),
    (:RPO, :rpo_resample_polyline_points),
    (:RPO, :rpo_retime_samples),
    (:RPO, :rpo_retimed_reference_from_profile),
    (:RPO, :rpo_sample_path_bezier),
    (:RPO, :rpo_sample_path_polyline),
    (:RPO, :validate_rpo_pso_config),
)
# Entries here are (namespace within Cloth, service name).
const CLOTH_SERVICES = (
    (:ClothRobotArmDynamics, :simulate_cloth_robot_arm_plan),
    (:ClothRobotArmDynamics, :assign_coupled_cloth_robot_arm_rhs!),
    (:ClothMultibody, :compliant_state_parts),
)
end
