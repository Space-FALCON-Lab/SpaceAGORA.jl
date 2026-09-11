"""Robot-arm planning module for quintic references and HYPR-based joint-space plans."""
module RobotArmPlanning

using LinearAlgebra
using Random
using StaticArrays
using ..Robotics
using ..HYPRUtils

export RobotArmPlannerConfig, RobotArmPlan, plan_robot_arm_motion, robot_arm_plan_sample
export RobotArmSphereObstacle, RobotArmHYPRConfig, RobotArmHYPRResult
export plan_robot_arm_motion_hypr, robot_arm_sample_hypr_path
export robot_arm_clearance_stats_from_samples, robot_arm_hypr_path_cost_components

"""Configuration for quintic joint-space robot-arm planning."""
Base.@kwdef struct RobotArmPlannerConfig
    dt_s::Float64 = 0.1
    duration_s::Float64 = 12.0
    ik_tol_m::Float64 = 1.0e-4
    ik_max_iters::Int = 100
    ik_damping::Float64 = 1.0e-3
end

"""Reference trajectory and model metadata for a planned robot-arm motion."""
struct RobotArmPlan
    model::ClothArmModel
    base_pose::ClothArmBasePose
    t_ref_s::Vector{Float64}
    q_ref::Matrix{Float64}
    dq_ref::Matrix{Float64}
    ddq_ref::Matrix{Float64}
    ee_ref::Matrix{Float64}
    q_start::Vector{Float64}
    q_goal::Vector{Float64}
    target::SVector{3, Float64}
    final_error_m::Float64
    planner::Symbol
end

"""Evaluate the smooth quintic blend position at a normalized time."""
function _quintic_scalar(s::Float64)
    σ = clamp(s, 0.0, 1.0)
    return 10σ^3 - 15σ^4 + 6σ^5
end

"""Evaluate the time derivative of the quintic blend."""
function _quintic_scalar_dot(s::Float64, Tf::Float64)
    σ = clamp(s, 0.0, 1.0)
    return (30σ^2 - 60σ^3 + 30σ^4) / Tf
end

"""Evaluate the second time derivative of the quintic blend."""
function _quintic_scalar_ddot(s::Float64, Tf::Float64)
    σ = clamp(s, 0.0, 1.0)
    return (60σ - 180σ^2 + 120σ^3) / Tf^2
end

"""Create the monotonic sample-time grid for a robot-arm reference trajectory."""
function _reference_times(dt_s::Float64, duration_s::Float64)
    dt_s > 0.0 || throw(ArgumentError("dt_s must be positive."))
    duration_s > 0.0 || throw(ArgumentError("duration_s must be positive."))
    t = collect(0.0:dt_s:duration_s)
    if isempty(t) || t[end] < duration_s
        push!(t, duration_s)
    end
    return t
end

include(joinpath(@__DIR__, "robot_arm_hypr.jl"))

"""Plan a quintic joint-space motion from an initial state to a target surface point."""
function plan_robot_arm_motion(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_start,
    target;
    config::RobotArmPlannerConfig=RobotArmPlannerConfig(),
    planner::Symbol=:cloth_quintic,
    hypr_config::Union{RobotArmHYPRConfig, Nothing}=nothing,
    obstacles::AbstractVector{RobotArmSphereObstacle}=RobotArmSphereObstacle[],
    rng=Random.default_rng(),
)
    if planner == :hypr
        hypr = plan_robot_arm_motion_hypr(
            model,
            base_pose,
            q_start,
            target;
            planner_config=config,
            hypr_config=something(hypr_config, RobotArmHYPRConfig()),
            obstacles=obstacles,
            rng=rng,
        )
        return hypr.plan
    elseif planner != :cloth_quintic
        throw(ArgumentError("Unsupported robot arm planner: $planner. Use :cloth_quintic or :hypr."))
    end
    q0 = Float64.(collect(q_start))
    q_goal = cloth_ik(
        model,
        base_pose,
        SVector{3, Float64}(target);
        q_seed=q0,
        max_iters=config.ik_max_iters,
        position_tol_m=config.ik_tol_m,
        damping=config.ik_damping,
    )
    t_ref = _reference_times(config.dt_s, config.duration_s)
    n = length(q0)
    nt = length(t_ref)
    q_ref = zeros(n, nt)
    dq_ref = zeros(n, nt)
    ddq_ref = zeros(n, nt)
    ee_ref = zeros(3, nt)
    Δq = q_goal .- q0
    for (k, t) in pairs(t_ref)
        s = t / config.duration_s
        a = _quintic_scalar(s)
        adot = _quintic_scalar_dot(s, config.duration_s)
        addot = _quintic_scalar_ddot(s, config.duration_s)
        q_ref[:, k] .= q0 .+ a .* Δq
        dq_ref[:, k] .= adot .* Δq
        ddq_ref[:, k] .= addot .* Δq
        ee_ref[:, k] .= cloth_fk(model, base_pose, q_ref[:, k]).end_effector_position
    end
    target_v = SVector{3, Float64}(target)
    final_error = norm(SVector{3, Float64}(ee_ref[:, end]) - target_v)
    return RobotArmPlan(
        model,
        base_pose,
        t_ref,
        q_ref,
        dq_ref,
        ddq_ref,
        ee_ref,
        q0,
        q_goal,
        target_v,
        final_error,
        planner,
    )
end

"""Sample a robot-arm plan at an arbitrary time."""
function robot_arm_plan_sample(plan::RobotArmPlan, t_s::Real)
    t = Float64(t_s)
    if t <= plan.t_ref_s[1]
        idx = 1
        return (
            q=copy(plan.q_ref[:, idx]),
            dq=copy(plan.dq_ref[:, idx]),
            ddq=copy(plan.ddq_ref[:, idx]),
            ee=SVector{3, Float64}(plan.ee_ref[:, idx]),
        )
    elseif t >= plan.t_ref_s[end]
        idx = length(plan.t_ref_s)
        return (
            q=copy(plan.q_ref[:, idx]),
            dq=copy(plan.dq_ref[:, idx]),
            ddq=copy(plan.ddq_ref[:, idx]),
            ee=SVector{3, Float64}(plan.ee_ref[:, idx]),
        )
    end
    hi = searchsortedfirst(plan.t_ref_s, t)
    lo = hi - 1
    α = (t - plan.t_ref_s[lo]) / (plan.t_ref_s[hi] - plan.t_ref_s[lo])
    q = (1 - α) .* plan.q_ref[:, lo] .+ α .* plan.q_ref[:, hi]
    dq = (1 - α) .* plan.dq_ref[:, lo] .+ α .* plan.dq_ref[:, hi]
    ddq = (1 - α) .* plan.ddq_ref[:, lo] .+ α .* plan.ddq_ref[:, hi]
    ee = (1 - α) .* plan.ee_ref[:, lo] .+ α .* plan.ee_ref[:, hi]
    return (q=q, dq=dq, ddq=ddq, ee=SVector{3, Float64}(ee))
end

end # module RobotArmPlanning
