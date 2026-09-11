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
end

"""Control effector that tracks a robot-arm plan and emits joint torque commands."""
Base.@kwdef mutable struct RobotArmControlEffector <: AbstractControlEffectorModel
    spacecraft_idx::Int = 1
    plan::Union{Nothing, RobotArmPlan} = nothing
    controller::Any = nothing
    updated_at_s::Float64 = 0.0
    joint_kp::Float64 = 0.0
    joint_kd::Float64 = 0.0
    control_dt_s::Float64 = 0.1
    k_translation_n_m::Float64 = 5.0e3
    c_translation_n_s_m::Float64 = 30.0
    k_rotation_n_m_rad::Float64 = 15.0
    c_rotation_n_m_s_rad::Float64 = 0.5
    held::RobotArmHeldActuation = RobotArmHeldActuation()
end

"""Estimate default joint inertias from link masses and lengths."""
function _robot_arm_default_joint_inertia(plan::RobotArmPlan)
    links = plan.model.links
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

"""Initialize the robot-arm joint-space MPC gain and limits."""
function init_robot_arm_joint_mpc(
    plan::RobotArmPlan;
    dt_s::Real=0.1,
    horizon::Integer=8,
    q_weight::Real=40.0,
    dq_weight::Real=6.0,
    torque_weight::Real=0.1,
    terminal_weight_scale::Real=4.0,
    joint_inertia_kg_m2=nothing,
)
    n = length(plan.q_start)
    horizon_i = Int(horizon)
    horizon_i > 0 || throw(ArgumentError("horizon must be positive."))
    dt = Float64(dt_s)
    dt > 0.0 || throw(ArgumentError("dt_s must be positive."))
    inertia = joint_inertia_kg_m2 === nothing ?
        _robot_arm_default_joint_inertia(plan) :
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
    return RobotArmJointMPCController(Ad, Bd, horizon_i, H, E, F, zeros(nu * horizon_i), inertia)
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
function robot_arm_measured_joint_state(plan::RobotArmPlan, sc_view)
    hasproperty(sc_view, :arm_q) || return nothing
    hasproperty(sc_view, :arm_ω) || return nothing
    n = length(plan.q_start)
    q_meas = zeros(n)
    dq_meas = zeros(n)
    parent_q = hasproperty(sc_view, :q) ? project_unit_quaternion(sc_view.q) : plan.base_pose.quaternion
    parent_ω_world = hasproperty(sc_view, :ω) ?
        rot(parent_q) * SVector{3, Float64}(sc_view.ω) :
        SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for i in 1:n
        child_q = project_unit_quaternion(sc_view.arm_q[:, i])
        rel_q = quat_mult(_robot_arm_control_quat_conj(parent_q), child_q)
        axis_parent = plan.model.joints[i].axis_parent
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

"""Apply the robot-arm controller at the current simulation time and cache the joint actuation command."""
function calcControlEffect!(model::RobotArmControlEffector, u, p::ODEParams, t::Float64, sat_idx::Int)
    sat_idx == model.spacecraft_idx || return nothing
    model.plan === nothing && return nothing
    plan_t = max(0.0, t - model.updated_at_s)
    sample = robot_arm_plan_sample(model.plan, plan_t)
    joint_torque = zeros(length(sample.q))
    if model.controller isa RobotArmJointMPCController
        sc_view = _robot_arm_control_spacecraft_state(u, model.spacecraft_idx)
        x_meas = sc_view === nothing ? nothing : robot_arm_measured_joint_state(model.plan, sc_view)
        x = x_meas === nothing ? _robot_arm_control_reference_state(model.plan, plan_t) : x_meas
        x_ref = robot_arm_joint_mpc_reference_preview(model.plan, plan_t, model.control_dt_s, model.controller.horizon)
        joint_torque .= robot_arm_joint_mpc_control(model.controller, x, x_ref)
    end
    model.held = RobotArmHeldActuation(
        joint_torque_nm=joint_torque,
        base_force_ii=SVector{3, Float64}(0.0, 0.0, 0.0),
        base_torque_body=SVector{3, Float64}(0.0, 0.0, 0.0),
    )
    return nothing
end

"""Return the held force and torque command for a control effector."""
function calcControlForceTorque(model::RobotArmControlEffector, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)
    i == model.spacecraft_idx || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    return model.held.base_force_ii, model.held.base_torque_body
end

"""Return propellant mass flow for a control effector."""
function calcControlMassFlowRate(model::RobotArmControlEffector, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64
    return 0.0
end
