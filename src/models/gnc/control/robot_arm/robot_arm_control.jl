"""Cached robot-arm control command held between controller updates."""
Base.@kwdef mutable struct RobotArmHeldActuation
    joint_torque_nm::Vector{Float64} = Float64[]
    base_force_ii::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    base_torque_body::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
end

"""Precomputed joint-space MPC controller for robot-arm tracking."""
mutable struct RobotArmJointMPCController
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    horizon::Int
    H::Matrix{Float64}
    E::Matrix{Float64}
    F::Matrix{Float64}
    U_prev::Vector{Float64}
    joint_inertia_kg_m2::Vector{Float64}
    dt_s::Float64                 # the step the joint model was discretized with
end

"""
Control effector that tracks a robot-arm plan and emits joint torque commands for
the spacecraft at run position `spacecraft_index`.

The plan is `plan`, with its time axis starting at `updated_at_s` (s into the
run), or, when `plan_buffer` is set, the plan of that `RobotArmPlanBuffer` from its
own `updated_at_s`: the link a joint-reference planner (an `AlgorithmGuidance` of
`JointReferenceFamily`) installs into. Configure one of the two; a buffer must
hold the plan the arm starts from. The run couples the arm's state to the plan
the controller tracks. At each run's start (`initialize_control!`) the MPC warm
start returns to zero, and a held command the controller wrote in the previous run
returns to its value at that run's start, while one a user set since is kept as the
new configuration: a held command configured with no plan is the command for the
whole run.
"""
mutable struct RobotArmControlEffector <: AbstractControlEffectorModel
    spacecraft_index::Int
    plan::Union{Nothing, RobotArmPlan}
    controller::Any
    updated_at_s::Float64
    joint_kp::Float64
    joint_kd::Float64
    control_dt_s::Float64
    k_translation_n_m::Float64
    c_translation_n_s_m::Float64
    k_rotation_n_m_rad::Float64
    c_rotation_n_m_s_rad::Float64
    held::RobotArmHeldActuation
    plan_buffer::Union{Nothing, RobotArmPlanBuffer}
    run_start::Union{Nothing, RobotArmHeldActuation}     # the held command the previous run started from
    run_written::Union{Nothing, RobotArmHeldActuation}   # a copy of the controller's last write of it in that run
end
# `run_start` and `run_written` start empty; they are not keywords.
function RobotArmControlEffector(; spacecraft_index=1, plan=nothing, controller=nothing,
        updated_at_s=0.0, joint_kp=0.0, joint_kd=0.0, control_dt_s=0.1, k_translation_n_m=5.0e3,
        c_translation_n_s_m=30.0, k_rotation_n_m_rad=15.0, c_rotation_n_m_s_rad=0.5, held=RobotArmHeldActuation(),
        plan_buffer=nothing)
    return RobotArmControlEffector(spacecraft_index, plan, controller, updated_at_s, joint_kp, joint_kd,
        control_dt_s, k_translation_n_m, c_translation_n_s_m, k_rotation_n_m_rad, c_rotation_n_m_s_rad, held,
        plan_buffer, nothing, nothing)
end

"""Estimate default joint inertias from link masses and lengths."""
_robot_arm_default_joint_inertia(plan::RobotArmPlan) = _robot_arm_default_joint_inertia(plan.model)
function _robot_arm_default_joint_inertia(model::ClothArmModel)
    links = model.links
    n = length(links)
    inertia = fill(1.0e-4, n)
    @inbounds for i in 1:n
        distance_to_joint = 0.0
        total = 0.0
        for j in i:n
            link = links[j]
            com_distance = distance_to_joint + norm(link.com_offset_parent)
            total += link.mass_kg * max(com_distance, link.radius_m)^2
            distance_to_joint += norm(link.vector_parent)
        end
        inertia[i] = max(total, 1.0e-6)
    end
    return inertia
end

"""
    init_robot_arm_joint_mpc(plan::RobotArmPlan; dt_s=0.1, horizon=8, q_weight=40.0, dq_weight=6.0,
                             torque_weight=0.1, terminal_weight_scale=4.0, joint_inertia_kg_m2=nothing)
    init_robot_arm_joint_mpc(model::RobotArmModel; kwargs...)

Initialize the robot-arm joint-space MPC gain and limits, for the arm of `plan` or
for `model` itself (one joint per joint of the arm), as a controller of a
`RobotArmControlEffector` whose plans arrive during the run through its
`plan_buffer`. `dt_s` is the control step (s); the default joint inertias come
from the arm's link masses and lengths.
"""
init_robot_arm_joint_mpc(plan::RobotArmPlan; kwargs...) =
    _init_robot_arm_joint_mpc(length(plan.q_start), plan.model; kwargs...)
init_robot_arm_joint_mpc(model::ClothArmModel; kwargs...) = _init_robot_arm_joint_mpc(length(model.joints), model; kwargs...)

function _init_robot_arm_joint_mpc(
    n::Int,
    model::ClothArmModel;
    dt_s::Real=0.1,
    horizon::Integer=8,
    q_weight::Real=40.0,
    dq_weight::Real=6.0,
    torque_weight::Real=0.1,
    terminal_weight_scale::Real=4.0,
    joint_inertia_kg_m2=nothing,
)
    horizon_i = Int(horizon)
    horizon_i > 0 || throw(ArgumentError("horizon must be positive."))
    dt = Float64(dt_s)
    dt > 0.0 || throw(ArgumentError("dt_s must be positive."))
    inertia = joint_inertia_kg_m2 === nothing ?
        _robot_arm_default_joint_inertia(model) :
        Float64.(collect(joint_inertia_kg_m2))
    length(inertia) == n || throw(ArgumentError("joint_inertia_kg_m2 must have one entry per joint."))
    all(>(0.0), inertia) || throw(ArgumentError("joint_inertia_kg_m2 entries must be positive."))

    Ad = [
        Matrix{Float64}(I, n, n) dt .* Matrix{Float64}(I, n, n)
        zeros(n, n) Matrix{Float64}(I, n, n)
    ]
    invJ = Diagonal(1.0 ./ inertia)
    Bd = [
        0.5 * dt^2 .* Matrix{Float64}(invJ)
        dt .* Matrix{Float64}(invJ)
    ]
    nx = 2n
    nu = n
    Abar, Bbar = rpo_prediction_mats(Ad, Bd, horizon_i)
    Q = Diagonal(vcat(fill(Float64(q_weight), n), fill(Float64(dq_weight), n)))
    Qf = Float64(terminal_weight_scale) .* Q
    R = Diagonal(fill(Float64(torque_weight), nu))
    Qbar = rpo_block_diag([fill(Matrix{Float64}(Q), max(horizon_i - 1, 0)); [Matrix{Float64}(Qf)]])
    Rbar = rpo_block_diag(fill(Matrix{Float64}(R), horizon_i))
    H = Bbar' * Qbar * Bbar + Rbar
    H = 0.5 .* (H .+ H')
    E = Bbar' * Qbar * Abar
    F = Bbar' * Qbar
    return RobotArmJointMPCController(Ad, Bd, horizon_i, H, E, F, zeros(nu * horizon_i), inertia, dt)
end

"""Build a finite-horizon robot-arm joint reference preview."""
function robot_arm_joint_mpc_reference_preview(plan::RobotArmPlan, t_elapsed_s::Real, dt_s::Real, horizon::Integer)
    n = length(plan.q_start)
    out = zeros(2n, Int(horizon) + 1)
    @inbounds for j in 0:Int(horizon)
        sample = robot_arm_plan_sample(plan, Float64(t_elapsed_s) + j * Float64(dt_s))
        out[1:n, j + 1] .= sample.q
        out[(n + 1):(2n), j + 1] .= sample.dq
    end
    return out
end

"""Compute the first joint acceleration command from the robot-arm MPC law."""
function robot_arm_joint_mpc_control(ctrl::RobotArmJointMPCController, x, x_ref)
    nx = size(ctrl.Ad, 1)
    nu = size(ctrl.Bd, 2)
    ref = Matrix{Float64}(x_ref)
    size(ref, 1) == nx || throw(ArgumentError("arm MPC reference has $(size(ref, 1)) rows, expected $nx."))
    if size(ref, 2) < ctrl.horizon + 1
        padded = repeat(ref[:, end], 1, ctrl.horizon + 1)
        padded[:, 1:size(ref, 2)] .= ref
        ref = padded
    end
    ref_stack = vec(ref[:, 2:(ctrl.horizon + 1)])
    rhs = ctrl.F * ref_stack - ctrl.E * Vector{Float64}(x)
    U = ctrl.H \ rhs
    ctrl.U_prev[1:(end - nu)] .= U[(nu + 1):end]
    ctrl.U_prev[(end - nu + 1):end] .= U[(end - nu + 1):end]
    return Vector{Float64}(U[1:nu])
end

"""Return a quaternion conjugate for robot-arm control calculations."""
@inline function _robot_arm_control_quat_conj(q)
    qv = project_unit_quaternion(q)
    return SVector{4, Float64}(-qv[1], -qv[2], -qv[3], qv[4])
end

"""Project a relative quaternion rotation onto a joint axis."""
@inline function _robot_arm_control_axis_angle_about(q_rel, axis)
    q = project_unit_quaternion(q_rel)
    q = q[4] < 0.0 ? -q : q
    a = SVector{3, Float64}(axis)
    return 2.0 * atan(dot(SVector{3, Float64}(q[1], q[2], q[3]), a), q[4])
end

"""Estimate current robot-arm joint positions and velocities from spacecraft view state."""
robot_arm_measured_joint_state(plan::RobotArmPlan, sc_view) =
    _robot_arm_measured_joints(plan.model, plan.base_pose.quaternion, length(plan.q_start), sc_view)

# The joint angles and rates of an arm of `n` joints, measured from the spacecraft's
# state view (`nothing` when it carries no arm state), with `base_quaternion` the
# arm base's attitude when the view has none.
function _robot_arm_measured_joints(model::ClothArmModel, base_quaternion, n::Int, sc_view)
    hasproperty(sc_view, :arm_q) || return nothing
    hasproperty(sc_view, :arm_ω) || return nothing
    q_meas = zeros(n)
    dq_meas = zeros(n)
    parent_q = hasproperty(sc_view, :q) ? project_unit_quaternion(sc_view.q) : base_quaternion
    parent_ω_world = hasproperty(sc_view, :ω) ?
        rot(parent_q) * SVector{3, Float64}(sc_view.ω) :
        SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for i in 1:n
        child_q = project_unit_quaternion(sc_view.arm_q[:, i])
        rel_q = quat_mult(_robot_arm_control_quat_conj(parent_q), child_q)
        axis_parent = model.joints[i].axis_parent
        axis_world = rot(parent_q) * axis_parent
        child_ω_world = rot(child_q) * SVector{3, Float64}(sc_view.arm_ω[:, i])
        q_meas[i] = _robot_arm_control_axis_angle_about(rel_q, axis_parent)
        dq_meas[i] = dot(child_ω_world - parent_ω_world, axis_world)
        parent_q = child_q
        parent_ω_world = child_ω_world
    end
    return vcat(q_meas, dq_meas)
end

"""Sample the robot-arm control reference at elapsed plan time."""
function _robot_arm_control_reference_state(plan::RobotArmPlan, t_elapsed_s::Real)
    sample = robot_arm_plan_sample(plan, t_elapsed_s)
    return vcat(sample.q, sample.dq)
end

"""Extract the spacecraft view state needed by robot-arm control."""
function _robot_arm_control_spacecraft_state(u, idx::Int)
    return hasproperty(u, :sc) && length(u.sc) >= idx ? u.sc[idx] : nothing
end

# The plan the controller tracks and the time its axis starts from: the plan
# buffer's when one is configured, otherwise the model's own.
_robot_arm_reference(model::RobotArmControlEffector) = model.plan_buffer === nothing ?
    (model.plan, model.updated_at_s) : (model.plan_buffer.plan, model.plan_buffer.updated_at_s)

"""Apply the robot-arm controller at the current simulation time and cache the joint actuation command."""
function update_control!(model::RobotArmControlEffector, u, p::ODEParams, t::Float64, sat_idx::Int)
    sat_idx == model.spacecraft_index || return nothing
    plan, origin_s = _robot_arm_reference(model)
    plan === nothing && return nothing
    plan_t = max(0.0, t - origin_s)
    sample = robot_arm_plan_sample(plan, plan_t)
    joint_torque = zeros(length(sample.q))
    if model.controller isa RobotArmJointMPCController
        sc_view = _robot_arm_control_spacecraft_state(u, model.spacecraft_index)
        x_meas = sc_view === nothing ? nothing : robot_arm_measured_joint_state(plan, sc_view)
        x = x_meas === nothing ? _robot_arm_control_reference_state(plan, plan_t) : x_meas
        x_ref = robot_arm_joint_mpc_reference_preview(plan, plan_t, model.control_dt_s, model.controller.horizon)
        joint_torque .= robot_arm_joint_mpc_control(model.controller, x, x_ref)
    end
    model.held = RobotArmHeldActuation(
        joint_torque_nm=joint_torque,
        base_force_ii=SVector{3, Float64}(0.0, 0.0, 0.0),
        base_torque_body=SVector{3, Float64}(0.0, 0.0, 0.0),
    )
    model.run_written = deepcopy(model.held)
    return nothing
end

# A controller tracks its own plan or its plan buffer's, not both; a buffer holds
# the plan the arm starts from, which sizes the arm's state.
function SimulationLifecycle.preflight_control(model::RobotArmControlEffector, args; isolate_state)
    model.plan_buffer === nothing && return nothing
    model.plan === nothing || throw(ArgumentError(
        "RobotArmControlEffector for spacecraft $(model.spacecraft_index) has both a plan and a plan_buffer; " *
        "configure one: plan= for a fixed plan, plan_buffer= for the plans a joint-reference planner installs."))
    model.plan_buffer.plan === nothing && throw(ArgumentError(
        "RobotArmControlEffector for spacecraft $(model.spacecraft_index) reads an empty plan_buffer: give it the " *
        "plan the arm starts from, from which the run sizes the arm's state; RobotArmPlanBuffer(model, base_pose, " *
        "q_start) holds the arm still at q_start."))
    return nothing
end

# Each run starts the held command from the configuration (G-8): one the previous
# run wrote returns to its value at that run's start, one a user set since is kept
# (a held command configured with no plan is the command for the whole run). The
# MPC warm start returns to zero.
function SimulationLifecycle.initialize_control!(model::RobotArmControlEffector, u, p, t)
    written = model.run_written === nothing ? nothing : (held=model.run_written,)
    start = model.run_start === nothing ? nothing : (held=model.run_start,)
    model.held = SimulationLifecycle._run_start_entries(start, written, (held=model.held,)).held
    model.run_start, model.run_written = deepcopy(model.held), nothing
    model.controller isa RobotArmJointMPCController && fill!(model.controller.U_prev, 0.0)
    return nothing
end

# The held command with the record of the run's last write of it (G-8) and the MPC
# warm start go into checkpoint files; a resumed run puts them back after its
# run-start reset.
function SimulationLifecycle.checkpoint_state(model::RobotArmControlEffector)
    ctrl = model.controller
    return (held=deepcopy(model.held), run_written=deepcopy(model.run_written),
        U_prev=ctrl isa RobotArmJointMPCController ? copy(ctrl.U_prev) : nothing)
end
function SimulationLifecycle.restore_checkpoint_state!(model::RobotArmControlEffector, data)
    ctrl = model.controller
    (ctrl isa RobotArmJointMPCController) == (data.U_prev !== nothing) || throw(ArgumentError(
        "The checkpoint's robot-arm controller state does not match this model's controller; resume with the same controller."))
    model.held, model.run_written = data.held, data.run_written
    ctrl isa RobotArmJointMPCController && copyto!(ctrl.U_prev, data.U_prev)
    return nothing
end

"""Return the held force and torque command for a control effector."""
function control_wrench(model::RobotArmControlEffector, x::StateSample, p::ODEParams, t::Float64)
    x.spacecraft_index == model.spacecraft_index || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    return model.held.base_force_ii, model.held.base_torque_body
end

# No control_propellant_rate method: the arm uses no propellant, and the
# dynamics then take the rate as zero.
