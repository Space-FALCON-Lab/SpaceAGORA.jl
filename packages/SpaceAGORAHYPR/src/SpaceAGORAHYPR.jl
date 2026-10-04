"""Optional HYPR search and configured planning for SpaceAGORA."""
module SpaceAGORAHYPR
import SpaceAGORA
module SwarmPolicy
import SpaceAGORA.SimulationModel.HYPRUtils: hypr_iteration_weights, hypr_material_improvement, hypr_protected_particle_mask
include("swarm_policy.jl")
end
module RPO
import SpaceAGORA
using LinearAlgebra, Random, StaticArrays
using SpaceAGORA.SimulationModel.HYPRUtils
using Logging: Logging
using Base.Threads: @threads, maxthreadid, threadid
using SpaceAGORA.SimulationModel.NavigationHooks
import SpaceAGORA.SimulationModel.GuidanceHooks:
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
include("rpo/path_retiming.jl")
include("rpo/sampling_adapters.jl")
include("rpo/metric_adapters.jl")
include("rpo/path_costs.jl")
include("rpo/pso_adaptive_policy.jl")
include("rpo/pso_helpers.jl")
include("rpo/pso_refinement.jl")
include("rpo/rrt_adapters.jl")
include("rpo/pso_path_planning.jl")
include("rpo/rpo_reference_trajectory.jl")
end
module RobotArm
import SpaceAGORA
using LinearAlgebra, Random, StaticArrays
using SpaceAGORA.SimulationModel.HYPRUtils
using SpaceAGORA.SimulationModel.Robotics
import SpaceAGORA.SimulationModel.RobotArmPlanning:
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
include("robot/swarm_and_retiming.jl")
include("robot/rrt_warmstart.jl")
include("robot/planner_core.jl")
end
module PlannerAdapter
using LinearAlgebra: norm
import SpaceAGORA
const P = SpaceAGORA.RPOPlannerInterfaces
const S = SpaceAGORA.SimulationModel
const G = S.GuidanceHooks
import SpaceAGORA.HYPRRPOPlanning: HYPRRPOPlanner, _plan_hypr_rpo!
include("rpo_planner_execution.jl")
end
function __init__()
    SpaceAGORA.SimulationModel.HYPRSupport.activate!()
end
end
