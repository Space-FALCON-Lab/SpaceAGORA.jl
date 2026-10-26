# The joint-reference family of the GNC algorithm contract (gnc_algorithms.jl): an
# algorithm proposes a joint-space reference (a `RobotArmPlan`) for one
# spacecraft's robot arm, trusted code checks its joint limits, its time grid, its
# goal and its sampled obstacle clearance, and an accepted reference is installed
# in the `RobotArmPlanBuffer` that the arm's controller reads. A refused proposal
# fails the request. `AlgorithmGuidance` runs any planner of the family.
module JointReferencePlanning

using LinearAlgebra: norm
using StaticArrays: SVector
using ..SimulationModel
import ..GNCAlgorithms as GA
import ..GNCRuntime as GR
import ..AlgorithmLifecycle as AL
using ..SimulationModel.RobotArmPlanning: RobotArmPlan, RobotArmPlanBuffer, RobotArmPlannerConfig,
    RobotArmSphereObstacle, plan_robot_arm_motion, robot_arm_clearance_stats_from_samples
using ..SimulationModel.Robotics: ClothArmModel, ClothArmBasePose, cloth_fk
const H = SimulationModel.ControlHooks

export AbstractJointReferencePlanner, JointReferenceFamily, JointReferenceRequest, JointReferenceSettings,
    QuinticJointPlanner, validate_joint_reference, joint_reference_is_current

"""
    JointReferenceFamily()

The joint-reference family of the GNC algorithm contract. Its request is a
`JointReferenceRequest`, a snapshot of one spacecraft's arm and of the motion
asked for; its product a `RobotArmPlan`, proposed in an `AlgorithmResult`; its
trusted validation `validate_joint_reference` (joint limits, a uniform time grid,
the goal and the obstacle clearance at the reference's samples); its install
target the `RobotArmPlanBuffer` a `RobotArmControlEffector` reads (its
`plan_buffer`); its currency rule `joint_reference_is_current`; and its rejection
policy `:fail`: a refused proposal stops the run with an `AlgorithmGuidanceError`
whose reason is the validation's.

An `AlgorithmGuidance` runs any planner of the family (an
`AbstractJointReferencePlanner`, such as `QuinticJointPlanner`), with
`JointReferenceSettings` as its `settings` and the controller's plan buffer as its
`link`. The buffer is a GNC link: the binding produces it, the controller consumes
it, and before the run the engine checks that one controller reads it, for the
binding's spacecraft and arm. The buffer holds the plan the arm starts from, such
as a hold at the start joints; each run starts from it again. The controller and
the arm's coupled state track an installed plan from the time it is installed at.
"""
struct JointReferenceFamily <: GA.AbstractAlgorithmFamily end

"""
    AbstractJointReferencePlanner <: AbstractGNCAlgorithm

An algorithm of the joint-reference family (`JointReferenceFamily`, inherited). A
planner implements `propose!(state, planner, request::JointReferenceRequest, rng)`
and returns an `AlgorithmResult` whose product, for a `:candidate`, is the
`RobotArmPlan` it proposes. `QuinticJointPlanner` is a subtype.
"""
abstract type AbstractJointReferencePlanner <: GA.AbstractGNCAlgorithm end
GA.algorithm_family(::AbstractJointReferencePlanner) = JointReferenceFamily()

_vector3(x, name) = length(x) == 3 ? SVector{3, Float64}(x[1], x[2], x[3]) :
    throw(ArgumentError("$name must contain 3 values."))

function _joint_vector(model::ClothArmModel, q, name)
    length(q) == length(model.joints) ||
        throw(ArgumentError("$name must have one entry per joint of the arm ($(length(model.joints)))."))
    v = Float64.(collect(q))
    all(isfinite, v) || throw(ArgumentError("$name must be finite."))
    return v
end

function _positive(x, name)
    (isfinite(x) && x > 0) || throw(ArgumentError("$name must be finite and positive."))
    return Float64(x)
end

function _nonnegative(x, name)
    (isfinite(x) && x >= 0) || throw(ArgumentError("$name must be finite and nonnegative."))
    return Float64(x)
end

"""
    JointReferenceSettings(; model, base_pose, target, obstacles=RobotArmSphereObstacle[], dt_s=0.1,
                           duration_s=12.0, clearance_m=0.0, goal_tolerance_m=5.0e-3)

The joint-reference family's settings for an `AlgorithmGuidance`: the arm, as a
`RobotArmModel` and its `RobotArmBasePose`; the end-effector `target` position (m,
in the frame the base pose is given in); the sphere `obstacles` the arm must clear
by `clearance_m` (m); the reference's sample step `dt_s` and duration
`duration_s` (s); and the distance `goal_tolerance_m` (m) within which the
reference's last sample must reach the target. Each request starts from the arm's
joints as the spacecraft's state measures them.
"""
struct JointReferenceSettings
    model::ClothArmModel
    base_pose::ClothArmBasePose
    target::SVector{3, Float64}
    obstacles::Vector{RobotArmSphereObstacle}
    dt_s::Float64
    duration_s::Float64
    clearance_m::Float64
    goal_tolerance_m::Float64
    function JointReferenceSettings(; model::ClothArmModel, base_pose::ClothArmBasePose, target,
            obstacles::AbstractVector{RobotArmSphereObstacle}=RobotArmSphereObstacle[], dt_s::Real=0.1,
            duration_s::Real=12.0, clearance_m::Real=0.0, goal_tolerance_m::Real=5.0e-3)
        new(model, base_pose, _vector3(target, "target"), collect(obstacles), _positive(dt_s, "dt_s"),
            _positive(duration_s, "duration_s"), _nonnegative(clearance_m, "clearance_m"),
            _positive(goal_tolerance_m, "goal_tolerance_m"))
    end
end

"""
    JointReferenceRequest(; request_id, reason, spacecraft_index, time_s, model, base_pose, q_start, target,
                          obstacles=RobotArmSphereObstacle[], dt_s, duration_s, clearance_m=0.0,
                          goal_tolerance_m=5.0e-3)

The joint-reference family's owned request snapshot, in SI units, for the arm of
the spacecraft at `spacecraft_index` (its position in the run, from 1) at elapsed
time `time_s`: the arm (`model` and `base_pose`), its start joints `q_start`
(rad), the end-effector `target` (m), the sphere `obstacles` and the clearance
`clearance_m` (m) the reference must keep from them, the sample step `dt_s` and
duration `duration_s` (s) of the reference, and the goal tolerance
`goal_tolerance_m` (m). Trusted family code builds it (`family_request`), with the
arm's measured joints as the start; the planner receives its own copy. The
reference's time axis starts at `time_s`.
"""
struct JointReferenceRequest <: GA.AbstractAlgorithmRequest
    request_id::Int
    reason::Symbol
    spacecraft_index::Int
    time_s::Float64
    model::ClothArmModel
    base_pose::ClothArmBasePose
    q_start::Vector{Float64}
    target::SVector{3, Float64}
    obstacles::Vector{RobotArmSphereObstacle}
    dt_s::Float64
    duration_s::Float64
    clearance_m::Float64
    goal_tolerance_m::Float64
    function JointReferenceRequest(; request_id::Integer, reason::Symbol, spacecraft_index::Integer, time_s::Real,
            model::ClothArmModel, base_pose::ClothArmBasePose, q_start, target,
            obstacles::AbstractVector{RobotArmSphereObstacle}=RobotArmSphereObstacle[], dt_s::Real, duration_s::Real,
            clearance_m::Real=0.0, goal_tolerance_m::Real=5.0e-3)
        request_id > 0 || throw(ArgumentError("request_id must be positive."))
        spacecraft_index > 0 || throw(ArgumentError("spacecraft_index must be positive."))
        isfinite(time_s) || throw(ArgumentError("time_s must be finite."))
        new(Int(request_id), reason, Int(spacecraft_index), Float64(time_s), model, base_pose,
            _joint_vector(model, q_start, "q_start"), _vector3(target, "target"), collect(obstacles),
            _positive(dt_s, "dt_s"), _positive(duration_s, "duration_s"), _nonnegative(clearance_m, "clearance_m"),
            _positive(goal_tolerance_m, "goal_tolerance_m"))
    end
end

# ---- Trusted validation ------------------------------------------------------

# The same arm: equal links, joints and mount offset (the vectors are compared by
# value, so a copied request matches).
_same_arm(a::ClothArmModel, b::ClothArmModel) =
    a.links == b.links && a.joints == b.joints && a.mount_offset_body == b.mount_offset_body

# A uniform grid from 0 to `duration_s`: every step `dt_s` but the last, which is
# no longer than `dt_s` and ends at the duration (the grid `plan_robot_arm_motion`
# samples), to a tolerance of 1e-9 of the duration.
function _uniform_grid(t::AbstractVector{Float64}, dt_s::Float64, duration_s::Float64)
    tol = 1e-9 * max(duration_s, 1.0)
    length(t) >= 2 || return false
    t[1] == 0.0 || return false
    for k in 2:(length(t) - 1)
        abs(t[k] - (k - 1) * dt_s) <= tol || return false
    end
    last_step = t[end] - t[end - 1]
    return abs(t[end] - duration_s) <= tol && 0.0 < last_step <= dt_s + tol
end

"""
    validate_joint_reference(request::JointReferenceRequest, result) -> ValidationResult

The joint-reference family's trusted validation of a proposal `result` against the
lifecycle's own copy of `request`. A result for another request (a different
`request_id`) is refused with `:request_mismatch`, whatever its status. A
`:candidate` for this request is accepted, with reason `:accepted` and metrics
`min_clearance_m`, `final_error_m` and `samples`, when its product is a
`RobotArmPlan` that:

1. is for the request's arm and base pose (`:model_mismatch`,
   `:base_pose_mismatch`);
2. has one row per joint and one column per sample in its joint and end-effector
   references, at least two samples, and finite times and references
   (`:invalid_shape`, `:nonfinite_reference`);
3. is sampled on the request's grid: from 0 every `dt_s`, with a last step of at
   most `dt_s`, to `duration_s` (`:invalid_timing`);
4. starts at the request's start joints, to 1e-9 rad (`:start_mismatch`);
5. keeps every joint within its limits at every sample (`:joint_limit`);
6. reaches the target within `goal_tolerance_m` at its last sample, by the arm's
   forward kinematics (`:goal_not_reached`);
7. clears every obstacle by `clearance_m` at every sample, link by link
   (`:insufficient_clearance`).

Any other result is refused with its status as the reason (`:infeasible`,
`:unsupported`, `:failed`), and a candidate whose product is not a `RobotArmPlan`
with `:invalid_product`. Validation installs nothing and draws from no random
stream.
"""
function validate_joint_reference(request::JointReferenceRequest, result)
    refuse(reason, metrics=NamedTuple()) = GA.ValidationResult(false, reason, metrics)
    result.request_id == request.request_id ||
        return refuse(:request_mismatch, (result_request_id=result.request_id, request_id=request.request_id))
    result.status === :candidate || return refuse(result.status, (termination=result.termination,))
    plan = result.product
    plan isa RobotArmPlan || return refuse(:invalid_product)
    model, base = request.model, request.base_pose
    _same_arm(plan.model, model) || return refuse(:model_mismatch)
    plan.base_pose == base || return refuse(:base_pose_mismatch)
    n, m = length(model.joints), length(plan.t_ref_s)
    (m >= 2 && length(plan.q_start) == n && all(A -> size(A) == (n, m), (plan.q_ref, plan.dq_ref, plan.ddq_ref)) &&
        size(plan.ee_ref) == (3, m)) || return refuse(:invalid_shape, (joints=n, samples=m))
    all(A -> all(isfinite, A), (plan.t_ref_s, plan.q_ref, plan.dq_ref, plan.ddq_ref, plan.ee_ref)) ||
        return refuse(:nonfinite_reference)
    _uniform_grid(plan.t_ref_s, request.dt_s, request.duration_s) ||
        return refuse(:invalid_timing, (dt_s=request.dt_s, duration_s=request.duration_s, last_time_s=plan.t_ref_s[end]))
    start_error = maximum(abs, view(plan.q_ref, :, 1) .- request.q_start)
    start_error <= 1e-9 || return refuse(:start_mismatch, (start_error_rad=start_error,))
    for (j, joint) in enumerate(model.joints), k in 1:m
        q = plan.q_ref[j, k]
        joint.lower_rad <= q <= joint.upper_rad || return refuse(:joint_limit,
            (joint=j, time_s=plan.t_ref_s[k], angle_rad=q, lower_rad=joint.lower_rad, upper_rad=joint.upper_rad))
    end
    final_error = norm(cloth_fk(model, base, view(plan.q_ref, :, m)).end_effector_position - request.target)
    final_error <= request.goal_tolerance_m || return refuse(:goal_not_reached, (final_error_m=final_error,))
    clearance = robot_arm_clearance_stats_from_samples(model, base, plan.q_ref, request.obstacles, request.clearance_m)
    clearance.violation_count == 0 || return refuse(:insufficient_clearance,
        (min_clearance_m=clearance.min_clearance, clearance_m=request.clearance_m))
    return GA.ValidationResult(true, :accepted, (min_clearance_m=clearance.min_clearance, final_error_m=final_error,
        samples=m))
end

"""
    joint_reference_is_current(request::JointReferenceRequest, plan::RobotArmPlan, time_s) -> Bool

Whether an installed `plan` may be used at `time_s`: from its request's time on.
After its last sample the controller holds the final joints.
"""
joint_reference_is_current(request::JointReferenceRequest, plan::RobotArmPlan, time_s) =
    isfinite(time_s) && Float64(time_s) >= request.time_s

# ---- The family's methods of the general contract -----------------------------

function GA.preflight_family(::JointReferenceFamily, planner, settings, args, spacecraft_index; isolate_state::Bool=true)
    settings isa JointReferenceSettings || throw(ArgumentError(
        "The joint-reference family needs JointReferenceSettings as the binding's settings, not $(repr(settings))."))
    return nothing
end

function GA.family_link(::JointReferenceFamily, link, settings, p, spacecraft_index)
    link isa RobotArmPlanBuffer || throw(ArgumentError(
        "The joint-reference family installs into a RobotArmPlanBuffer: configure the binding's `link` as the " *
        "plan_buffer of the arm's RobotArmControlEffector."))
    return link
end

# The start joints are the arm's measured joints. The arm controller that reads the
# plan buffer couples the arm's state to the spacecraft's (the buffer holds the
# plan it starts from), so the state carries them.
function GA.family_request(::JointReferenceFamily, settings::JointReferenceSettings, u, p, t, spacecraft_index, reason,
        request_id)
    i = Int(spacecraft_index)
    model = settings.model
    n = length(model.joints)
    sc_view = hasproperty(u, :sc) && length(u.sc) >= i ? u.sc[i] : nothing
    measured = sc_view === nothing ? nothing : H._robot_arm_measured_joints(model, settings.base_pose.quaternion, n, sc_view)
    measured === nothing && throw(ArgumentError(
        "Spacecraft $(i)'s state carries no arm joint state: its RobotArmControlEffector must read the binding's " *
        "plan buffer, which couples the arm's state."))
    return JointReferenceRequest(; request_id=request_id, reason=reason, spacecraft_index=i, time_s=t, model=model,
        base_pose=settings.base_pose, q_start=measured[1:n], target=settings.target, obstacles=settings.obstacles,
        dt_s=settings.dt_s, duration_s=settings.duration_s, clearance_m=settings.clearance_m,
        goal_tolerance_m=settings.goal_tolerance_m)
end

GA.validate_proposal(::JointReferenceFamily, request::JointReferenceRequest, result, trusted) =
    validate_joint_reference(request, result)

# The installed plan and its time origin replace the buffer's in one step.
function GA.install!(buffer::RobotArmPlanBuffer, ::JointReferenceFamily, plan::RobotArmPlan,
        request::JointReferenceRequest, t; validation, diagnostics)
    buffer.plan, buffer.updated_at_s = plan, Float64(t)
    buffer.run_written = (plan=deepcopy(plan), updated_at_s=Float64(t))
    return nothing
end

GA.is_current(::JointReferenceFamily, request::JointReferenceRequest, plan::RobotArmPlan, time_s) =
    joint_reference_is_current(request, plan, time_s)
GA.on_rejection(::JointReferenceFamily) = :fail
GA.family_contract_version(::JointReferenceFamily) = v"1.0.0"
GA.algorithm_result_type(::JointReferenceFamily) = GA.AlgorithmResult

# ---- The plan buffer as a GNC link --------------------------------------------

# One controller reads a plan buffer; at most one model writes it.
GR.link_cardinality(::RobotArmPlanBuffer) = (producers=0:1, consumers=1:1)

# The controller's MPC fits the arm whose plan the buffer starts with, whether or
# not a model writes the buffer; a planner acts on the controller's spacecraft, for
# that arm.
function GR.check_link(buffer::RobotArmPlanBuffer, producers, consumers, args)
    controller = only(consumers)
    controller isa H.RobotArmControlEffector || throw(ArgumentError(
        "A RobotArmPlanBuffer is read by a RobotArmControlEffector, not a $(nameof(typeof(controller)))."))
    start, mpc = buffer.plan, controller.controller
    if start !== nothing && mpc isa H.RobotArmJointMPCController
        joints = length(start.q_start)
        length(mpc.joint_inertia_kg_m2) == joints || throw(ArgumentError(
            "The RobotArmControlEffector's MPC controller is built for $(length(mpc.joint_inertia_kg_m2)) joints, " *
            "the arm of its plan buffer's plan has $(joints)."))
    end
    for producer in producers
        producer isa AL.AlgorithmGuidance || continue
        producer.spacecraft_index == controller.spacecraft_index || throw(ArgumentError(
            "The joint-reference planner plans for spacecraft $(producer.spacecraft_index), but the " *
            "RobotArmControlEffector reading its plan buffer acts on spacecraft $(controller.spacecraft_index)."))
        settings = producer.settings
        settings isa JointReferenceSettings || continue
        start === nothing || (_same_arm(start.model, settings.model) && start.base_pose == settings.base_pose) ||
            throw(ArgumentError("The plan buffer starts with a plan for another arm or base pose than the " *
                "joint-reference planner's settings; the arm's state is sized from that plan."))
    end
    return nothing
end

# Each run starts from the buffer's configured plan and time: a plan the previous
# run installed returns to the one that run started from, a plan a user set since
# is kept (`_buffer_run_start`). The engine builds the arm's initial state before
# this runs, by the same rule, so state and buffer start alike.
function GR.link_initialize!(buffer::RobotArmPlanBuffer, producers, consumers, u, p, t)
    start = SimulationModel.RobotArmPlanning._buffer_run_start(buffer)
    buffer.plan, buffer.updated_at_s = start.plan, start.updated_at_s
    buffer.run_start, buffer.run_written = deepcopy(start), nothing
    return nothing
end

# ---- A planner -----------------------------------------------------------------

"""
    QuinticJointPlanner(; ik_tol_m=1.0e-4, ik_max_iters=100, ik_damping=1.0e-3)

A planner of the joint-reference family (`AbstractJointReferencePlanner`): the
quintic joint-space blend of `plan_robot_arm_motion`, from the request's start
joints to the inverse-kinematics solution for its target (solved to `ik_tol_m`,
in at most `ik_max_iters` damped steps of damping `ik_damping`), sampled every
`dt_s` over `duration_s`. It does not avoid obstacles, so the family's validation
refuses a reference that does not clear them. When the inverse kinematics does not
reach the target, its result is `:infeasible`, with termination
`:target_unreachable`. It draws no random numbers (`consumes_rng` is `false`).
"""
struct QuinticJointPlanner <: AbstractJointReferencePlanner
    ik_tol_m::Float64
    ik_max_iters::Int
    ik_damping::Float64
    function QuinticJointPlanner(; ik_tol_m::Real=1.0e-4, ik_max_iters::Integer=100, ik_damping::Real=1.0e-3)
        ik_max_iters > 0 || throw(ArgumentError("ik_max_iters must be positive."))
        new(_positive(ik_tol_m, "ik_tol_m"), Int(ik_max_iters), _nonnegative(ik_damping, "ik_damping"))
    end
end

GA.consumes_rng(::QuinticJointPlanner) = false

function GA.propose!(state, planner::QuinticJointPlanner, request::JointReferenceRequest, rng)
    config = RobotArmPlannerConfig(dt_s=request.dt_s, duration_s=request.duration_s, ik_tol_m=planner.ik_tol_m,
        ik_max_iters=planner.ik_max_iters, ik_damping=planner.ik_damping)
    plan = try
        plan_robot_arm_motion(request.model, request.base_pose, request.q_start, request.target; config=config)
    catch err
        # The inverse kinematics reports a target it cannot reach this way.
        (err isa ErrorException && startswith(err.msg, "cloth_ik failed")) || rethrow()
        return GA.AlgorithmResult(request_id=request.request_id, status=:infeasible, termination=:target_unreachable)
    end
    return GA.AlgorithmResult(request_id=request.request_id, status=:candidate, termination=:planned, product=plan,
        diagnostics=(final_error_m=plan.final_error_m, q_goal=copy(plan.q_goal)))
end

end # module JointReferencePlanning
