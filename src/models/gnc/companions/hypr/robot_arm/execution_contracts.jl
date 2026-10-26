# Compatibility function identities only. The optional companion owns execution.
# HYPR's own helpers are defined in its extension, not declared here.
using ..HYPRSupport
"""
    plan_robot_arm_motion_hypr(model, base_pose, q_start, target;
        planner_config=RobotArmPlannerConfig(), hypr_config=RobotArmHYPRConfig(),
        obstacles=RobotArmSphereObstacle[], rng)

Plan a robot-arm trajectory with HYPR sampling, PSO search, optional RRT warm
start, and retiming. `q_start` is the initial joint configuration; `target` is the
end-effector target passed to inverse kinematics. `rng` (an `AbstractRNG`) is
required and owns every random draw; there is no default global stream. Load `SpaceAGORA` and `HYPR` before
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
`SpaceAGORA` and `HYPR` before calling this function.
"""
function robot_arm_hypr_path_cost_components end
robot_arm_hypr_path_cost_components(args...; kwargs...) = HYPRSupport.unavailable(robot_arm_hypr_path_cost_components, args)
