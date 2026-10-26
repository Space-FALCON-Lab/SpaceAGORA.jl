"""Robot-arm planning module for quintic references and HYPR-based joint-space plans."""
module RobotArmPlanning

using LinearAlgebra
using Random
using StaticArrays
using ..Robotics
using ..HYPRUtils
import ...SimulationLifecycle

export RobotArmPlannerConfig, RobotArmPlan, RobotArmPlanBuffer, plan_robot_arm_motion, robot_arm_plan_sample
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

"""
    RobotArmPlanBuffer(; plan=nothing, updated_at_s=0.0)
    RobotArmPlanBuffer(model, base_pose, q_start)

The link between a robot arm's planner and its controller: the `RobotArmPlan`
the arm tracks and the time (s into the run) the plan's time axis starts from.
A `RobotArmControlEffector` configured with `plan_buffer=` reads it at each tick,
and the run couples the arm's state to it as to a controller's own `plan`: the
arm's state is sized and started from the plan the buffer holds when the run
starts, so it must hold one. The second form holds the arm of `model` at its
`base_pose` still at the joint angles `q_start` (rad).

An `AlgorithmGuidance` of the joint-reference family (`JointReferenceFamily`)
configured with `link=` the same buffer installs each accepted plan into it, at
the time of its request. At each run's start, a plan and time a planner installed
in the previous run return to the buffer's plan and time at that run's start,
while a plan or time a user set since is kept as the new configuration.
"""
Base.@kwdef mutable struct RobotArmPlanBuffer
    plan::Union{Nothing, RobotArmPlan} = nothing
    updated_at_s::Float64 = 0.0
    # The plan and time the previous run started from, and a copy of the last
    # install in that run (`_buffer_run_start`).
    run_start::Union{Nothing, NamedTuple{(:plan, :updated_at_s), Tuple{Union{Nothing, RobotArmPlan}, Float64}}} = nothing
    run_written::Union{Nothing, NamedTuple{(:plan, :updated_at_s), Tuple{RobotArmPlan, Float64}}} = nothing
end

# The plan and time the next run starts from, without changing the buffer: an entry
# still equal to the previous run's last install returns to its value at that
# run's start, an entry a user changed since is kept (`_run_start_entries`).
_buffer_run_start(buffer::RobotArmPlanBuffer) = SimulationLifecycle._run_start_entries(
    buffer.run_start, buffer.run_written, (plan=buffer.plan, updated_at_s=buffer.updated_at_s))

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

include(joinpath(@__DIR__, "..", "..", "companions", "hypr", "robot_arm.jl"))
include(joinpath(@__DIR__, "path_geometry.jl"))

"""
    plan_robot_arm_motion(model, base_pose, q_start, target; config, planner=:cloth_quintic,
                          hypr_config=nothing, obstacles=RobotArmSphereObstacle[], rng=nothing)

Plan a joint-space motion from an initial state to a target surface point. The
default `:cloth_quintic` planner is deterministic and draws no random numbers.
`planner=:hypr` samples, so it requires an explicit `rng` (an `AbstractRNG`); there
is no default global stream.
"""
function plan_robot_arm_motion(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_start,
    target;
    config::RobotArmPlannerConfig=RobotArmPlannerConfig(),
    planner::Symbol=:cloth_quintic,
    hypr_config::Union{RobotArmHYPRConfig, Nothing}=nothing,
    obstacles::AbstractVector{RobotArmSphereObstacle}=RobotArmSphereObstacle[],
    rng::Union{Nothing, Random.AbstractRNG}=nothing,
)
    if planner == :hypr
        rng === nothing && throw(ArgumentError(
            "plan_robot_arm_motion(...; planner=:hypr) samples and needs an explicit rng, for example rng=MersenneTwister(seed)."))
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

# A plan that keeps the arm at `q_start` (its goal is where those joints put the
# end effector, so the blend does not move).
function RobotArmPlanBuffer(model::ClothArmModel, base_pose::ClothArmBasePose, q_start)
    q = Float64.(collect(q_start))
    target = cloth_fk(model, base_pose, q).end_effector_position
    return RobotArmPlanBuffer(plan=plan_robot_arm_motion(model, base_pose, q, target), updated_at_s=0.0)
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
