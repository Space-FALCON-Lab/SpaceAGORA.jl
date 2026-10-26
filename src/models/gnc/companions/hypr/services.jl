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
const CONTRACT_VERSION = v"3.0.0"
const HYPR_VERSION = v"0.5.0"
const SUPPORTED_HYPR_VERSIONS = (HYPR_VERSION,)
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
    provider === :HYPR && version in SUPPORTED_HYPR_VERSIONS || throw(HYPRSupport.HYPRCompatibilityError(
        "Unsupported HYPR implementation/version; this contract supports HYPR versions $SUPPORTED_HYPR_VERSIONS."))
    HYPRSupport.check_provider(provider, version)
    return nothing
end

"""Activate one compatible implementation after its extension has initialized."""
function activate!(provider::Symbol, version::VersionNumber, contract::VersionNumber)
    check_provider(provider, version, contract)
    HYPRSupport.activate!(provider, version)
    return nothing
end

module RPO
import ...SimulationModel.GuidanceHooks:
    RPOAdaptiveSamplingSettings,
    RPOPSOConfig,
    RPORRTConnectSettings,
    RPORRTStarSettings,
    RPORetimeCurve,
    rpo_adaptive_sampling_min_ds_m,
    rpo_adaptive_segment_samples,
    rpo_arc_length_params,
    rpo_clearance_stats_from_samples,
    rpo_fuel_proxy_from_samples,
    rpo_hcw_fuel_proxy,
    rpo_hypr_refinement_sampling_density_m,
    rpo_hypr_sampling_density_m,
    rpo_normalized_path_cost_components,
    rpo_obstacle_sigmoid_penalty,
    rpo_obstacle_sigmoid_threshold,
    rpo_path_cost_normalization_refs,
    rpo_path_length,
    rpo_post_refine_path,
    rpo_profile_hcw_fuel_proxy,
    rpo_pso_config,
    rpo_pso_plan_path,
    rpo_reference_from_path,
    rpo_resample_polyline_points,
    rpo_retime_profile,
    rpo_retime_samples,
    rpo_retime_sampling_ds_m,
    rpo_retimed_reference_from_profile,
    rpo_rrt_connect_bezier_plan_path,
    rpo_rrt_connect_plan_path,
    rpo_rrt_star_plan_path,
    rpo_sample_path,
    rpo_sample_path_bezier_adaptive,
    rpo_sample_path_bezier_adaptive_with_params,
    rpo_sample_path_polyline,
    rpo_sample_path_polyline_adaptive,
    rpo_sample_path_with_params,
    validate_rpo_pso_config
import ...SimulationModel.NavigationHooks:
    rpo_clearance_distance_to_station,
    rpo_path_clearance_stats,
    RPOReferenceGeometry
export RPOAdaptiveSamplingSettings, RPOPSOConfig, RPORRTConnectSettings, RPORRTStarSettings, RPORetimeCurve, rpo_adaptive_sampling_min_ds_m, rpo_adaptive_segment_samples, rpo_arc_length_params, rpo_clearance_stats_from_samples, rpo_fuel_proxy_from_samples, rpo_hcw_fuel_proxy, rpo_hypr_refinement_sampling_density_m, rpo_hypr_sampling_density_m, rpo_normalized_path_cost_components, rpo_obstacle_sigmoid_penalty, rpo_obstacle_sigmoid_threshold, rpo_path_cost_normalization_refs, rpo_path_length, rpo_post_refine_path, rpo_profile_hcw_fuel_proxy, rpo_pso_config, rpo_pso_plan_path, rpo_reference_from_path, rpo_resample_polyline_points, rpo_retime_profile, rpo_retime_samples, rpo_retime_sampling_ds_m, rpo_retimed_reference_from_profile, rpo_rrt_connect_bezier_plan_path, rpo_rrt_connect_plan_path, rpo_rrt_star_plan_path, rpo_sample_path, rpo_sample_path_bezier_adaptive, rpo_sample_path_bezier_adaptive_with_params, rpo_sample_path_polyline, rpo_sample_path_polyline_adaptive, rpo_sample_path_with_params, validate_rpo_pso_config, rpo_clearance_distance_to_station, rpo_path_clearance_stats, RPOReferenceGeometry
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
    _robot_arm_path_length,
    _robot_arm_path_smoothness,
    _validate_robot_arm_hypr_config,
    plan_robot_arm_motion_hypr,
    robot_arm_clearance_stats_from_samples,
    robot_arm_hypr_path_cost_components,
    robot_arm_sample_hypr_path
import ...SimulationModel.HYPRUtils:
    hypr_rrt_join_paths,
    hypr_rrt_nearest_index,
    hypr_rrt_steer
import ...SimulationModel.Robotics:
    ClothArmBasePose,
    ClothArmModel,
    cloth_fk,
    cloth_ik
export RobotArmHYPRConfig, RobotArmHYPRResult, RobotArmPlan, RobotArmPlannerConfig, RobotArmRRTConnectTree, RobotArmSphereObstacle, _reference_times, _robot_arm_path_length, _robot_arm_path_smoothness, _validate_robot_arm_hypr_config, plan_robot_arm_motion_hypr, robot_arm_clearance_stats_from_samples, robot_arm_hypr_path_cost_components, robot_arm_sample_hypr_path, hypr_rrt_join_paths, hypr_rrt_nearest_index, hypr_rrt_steer, ClothArmBasePose, ClothArmModel, cloth_fk, cloth_ik
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
    (:RobotArm, :plan_robot_arm_motion_hypr),
    (:RobotArm, :robot_arm_hypr_path_cost_components),
    (:Planner, :_plan_hypr_rpo!),
    (:RPO, :rpo_adaptive_sampling_min_ds_m),
    (:RPO, :rpo_fuel_proxy_from_samples),
    (:RPO, :rpo_normalized_path_cost_components),
    (:RPO, :rpo_path_cost_normalization_refs),
    (:RPO, :rpo_post_refine_path),
    (:RPO, :rpo_pso_plan_path),
    (:RPO, :rpo_reference_from_path),
    (:RPO, :rpo_retime_profile),
    (:RPO, :rpo_retime_sampling_ds_m),
    (:RPO, :rpo_rrt_connect_bezier_plan_path),
    (:RPO, :rpo_rrt_connect_plan_path),
    (:RPO, :rpo_rrt_star_plan_path),
    (:RPO, :rpo_sample_path),
    (:RPO, :rpo_sample_path_bezier_adaptive),
    (:RPO, :rpo_sample_path_bezier_adaptive_with_params),
    (:RPO, :rpo_sample_path_polyline_adaptive),
    (:RPO, :rpo_sample_path_with_params),
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
    (:RobotArm, :_validate_robot_arm_hypr_config),
    (:RobotArm, :cloth_fk),
    (:RobotArm, :cloth_ik),
    (:RobotArm, :hypr_rrt_join_paths),
    (:RobotArm, :hypr_rrt_nearest_index),
    (:RobotArm, :hypr_rrt_steer),
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
