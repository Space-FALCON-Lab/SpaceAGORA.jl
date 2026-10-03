# Compatibility function identities only. The optional companion owns execution.
# Private names remain available qualified for existing research clients.
using ..HYPRSupport
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_control_points end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_empty_rrt_warmstart_diagnostics end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_flatten_internal_points end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_base_wrench_ratios end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_cloth_base_wrench_ratios end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_cloth_state_for_reaction end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_cull_swarm! end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_early_stopping_feasible end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_iteration_weights end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_link_com_history end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_material_improvement end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_post_refine_points end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_reaction_scale end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_reference_times_from_scales end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_refinement_better end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_retime_reference end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_hypr_rigid_base_wrench_ratios end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_plan_from_q_reference end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_resample_polyline_points end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_connect! end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_connect_warmstart_path end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_extend! end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_join_paths end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_nearest_index end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_path_score end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_random_state end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_segment_is_safe end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_segment_samples end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_shortcut_path end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_steer end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_rrt_warmstart_fields end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _robot_arm_seed_control_points end
"""
    plan_robot_arm_motion_hypr(model, base_pose, q_start, target;
        planner_config=RobotArmPlannerConfig(), hypr_config=RobotArmHYPRConfig(),
        obstacles=RobotArmSphereObstacle[], rng=Random.default_rng())

Plan a robot-arm trajectory with HYPR sampling, PSO search, optional RRT warm
start, and retiming. `q_start` is the initial joint configuration; `target` is the
end-effector target passed to inverse kinematics. Load `SpaceAGORAHYPR` before
calling this function. The returned `RobotArmHYPRResult` holds the joint reference
and end-effector trajectory in its `plan` field, together with planning metrics.
"""
function plan_robot_arm_motion_hypr end
plan_robot_arm_motion_hypr(args...; kwargs...) = HYPRSupport.unavailable(plan_robot_arm_motion_hypr, args)
"""
    robot_arm_hypr_path_cost_components(points, model, base_pose, obstacles, cfg;
        cost_cutoff=Inf)

Evaluate HYPR objective components for a candidate robot-arm joint path. Columns
of `points` are joint-space control points; `cfg` is a `RobotArmHYPRConfig`.
Returns the sampled path's cost and feasibility components. Load
`SpaceAGORAHYPR` before calling this function.
"""
function robot_arm_hypr_path_cost_components end
robot_arm_hypr_path_cost_components(args...; kwargs...) = HYPRSupport.unavailable(robot_arm_hypr_path_cost_components, args)
