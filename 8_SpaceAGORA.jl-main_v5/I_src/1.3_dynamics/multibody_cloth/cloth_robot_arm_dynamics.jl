"""Coupled cloth robot-arm dynamics module for compliant arm realization and simulation."""
module ClothRobotArmDynamics

using LinearAlgebra
using StaticArrays
using ..ClothMultibody
using ..RobotArmPlanning
using ..Robotics

export ClothRobotArmDynamicsMode, ClothRobotArmReferenceState, cloth_reference_state
export ClothRobotArmSimulation, cloth_robot_arm_multibody, cloth_robot_arm_initial_state
export cloth_robot_arm_rest_quaternions, cloth_robot_arm_end_effector
export cloth_robot_arm_actuators
export coupled_cloth_robot_arm_state_shape, initialize_coupled_cloth_robot_arm_state!
export assign_coupled_cloth_robot_arm_rhs!
export simulate_cloth_robot_arm_plan

const ClothRobotArmDynamicsMode = Symbol
const _Q_IDENTITY = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

"""Reference joint, velocity, and acceleration sample for cloth-arm coupling."""
struct ClothRobotArmReferenceState
    t_s::Float64
    state::ClothArmState
    end_effector_position::SVector{3, Float64}
end

"""Sample a robot-arm plan into the reference state used by cloth coupling."""
function cloth_reference_state(plan::RobotArmPlan, t_s::Real)::ClothRobotArmReferenceState
    sample = robot_arm_plan_sample(plan, t_s)
    state = cloth_fk_state(plan.model, plan.base_pose, sample.q; dq=sample.dq)
    return ClothRobotArmReferenceState(Float64(t_s), state, sample.ee)
end

"""Trajectory output from simulating a cloth-coupled robot arm."""
struct ClothRobotArmSimulation
    multibody_model::CompliantMultibodyModel
    trajectory::CompliantMultibodyTrajectory
    end_effector_positions::Matrix{Float64}
    reference_end_effector_positions::Matrix{Float64}
    tracking_error_m::Vector{Float64}
    joint_compliance_torques_body::Array{Float64, 3}
    joint_actuator_torques_body::Array{Float64, 3}
    joint_total_torques_body::Array{Float64, 3}
end

"""Normalize a quaternion and return it as a static vector."""
@inline function _unit_quat(q)::SVector{4, Float64}
    qv = SVector{4, Float64}(q)
    nq = norm(qv)
    return isfinite(nq) && nq > eps(Float64) ? qv / nq : SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
end

"""Return the conjugate of a quaternion."""
@inline function _quat_conj(q)::SVector{4, Float64}
    qv = _unit_quat(q)
    return SVector{4, Float64}(-qv[1], -qv[2], -qv[3], qv[4])
end

"""Multiply two quaternions using the project convention."""
@inline function _quat_mul(a, b)::SVector{4, Float64}
    aq = _unit_quat(a)
    bq = _unit_quat(b)
    ax, ay, az, aw = aq
    bx, by, bz, bw = bq
    return _unit_quat(SVector{4, Float64}(
        aw * bx + bw * ax + ay * bz - az * by,
        aw * by + bw * ay + az * bx - ax * bz,
        aw * bz + bw * az + ax * by - ay * bx,
        aw * bw - ax * bx - ay * by - az * bz,
    ))
end

"""Convert a quaternion to a world-frame rotation matrix."""
@inline function _rot(q_raw)::SMatrix{3, 3, Float64}
    q = _unit_quat(q_raw)
    x, y, z, w = q
    return SMatrix{3, 3, Float64}(
        1 - 2y^2 - 2z^2, 2x * y + 2z * w, 2x * z - 2y * w,
        2x * y - 2z * w, 1 - 2x^2 - 2z^2, 2y * z + 2x * w,
        2x * z + 2y * w, 2y * z - 2x * w, 1 - 2x^2 - 2y^2,
    )
end

"""Return a link inertia matrix from explicit inertia or cylindrical approximation."""
function _link_inertia(link)
    len = norm(link.vector_parent)
    r = link.radius_m
    m = link.mass_kg
    return SMatrix{3, 3, Float64}(
        0.5 * m * r^2, 0.0, 0.0,
        0.0, (1 / 12) * m * (3r^2 + len^2), 0.0,
        0.0, 0.0, (1 / 12) * m * (3r^2 + len^2),
    )
end

"""Create a 3-by-3 diagonal static matrix from a scalar."""
@inline _diag3(v::Real) = SMatrix{3, 3, Float64}(Float64(v) * I)

"""Normalize scalar, vector, tuple, or matrix compliance input into a 3-by-3 matrix."""
function _compliance_matrix(value, name::Symbol)::SMatrix{3, 3, Float64}
    if value isa Real
        return _diag3(value)
    elseif value isa AbstractMatrix
        size(value) == (3, 3) || throw(ArgumentError("$(name) matrix entries must be 3x3."))
        return SMatrix{3, 3, Float64}(value)
    elseif value isa AbstractVector
        length(value) == 3 || throw(ArgumentError("$(name) vector entries must have length 3 for diagonal matrices."))
        return SMatrix{3, 3, Float64}(Diagonal(Float64.(value)))
    elseif value isa Tuple
        length(value) == 3 || throw(ArgumentError("$(name) tuple entries must have length 3 for diagonal matrices."))
        return SMatrix{3, 3, Float64}(Diagonal(collect(Float64.(value))))
    else
        throw(ArgumentError("$(name) must be a scalar, a 3x3 matrix, or a per-joint vector of those values."))
    end
end

"""Expand joint compliance values into one matrix per robot-arm joint."""
function _joint_compliance_matrices(value, n::Int, name::Symbol)::Vector{SMatrix{3, 3, Float64}}
    if value isa Real || value isa AbstractMatrix || value isa Tuple
        return fill(_compliance_matrix(value, name), n)
    elseif value isa AbstractVector && length(value) == n
        return [_compliance_matrix(value[i], name) for i in 1:n]
    else
        throw(ArgumentError("$(name) must be a scalar, a 3x3 matrix, or a vector with one entry per joint."))
    end
end

"""Compute world velocity of a body-fixed offset point for coupled arm dynamics."""
@inline function _body_offset_velocity_world(v_world, q, ω_body, p_body)
    return SVector{3, Float64}(v_world) + _rot(q) * cross(SVector{3, Float64}(ω_body), SVector{3, Float64}(p_body))
end

"""Convert a quaternion error into an axis-angle error vector."""
@inline function _axis_angle_error(q_err_raw)::SVector{3, Float64}
    q_err = _unit_quat(q_err_raw)
    q_err = q_err[4] < 0.0 ? -q_err : q_err
    v = SVector{3, Float64}(q_err[1], q_err[2], q_err[3])
    nv = norm(v)
    if nv <= 1.0e-12
        return 2.0 * v
    end
    θ = 2.0 * atan(nv, q_err[4])
    return (θ / nv) * v
end

"""Multiply two quaternions without normalizing the result."""
@inline function _quat_raw_mul(a::SVector{4, Float64}, b::SVector{4, Float64})::SVector{4, Float64}
    ax, ay, az, aw = a
    bx, by, bz, bw = b
    return SVector{4, Float64}(
        aw * bx + bw * ax + ay * bz - az * by,
        aw * by + bw * ay + az * bx - ax * bz,
        aw * bz + bw * az + ax * by - ay * bx,
        aw * bw - ax * bx - ay * by - az * bz,
    )
end

"""Compute the rest child-to-parent orientation between two body quaternions."""
function _rest_child_parent_quat(parent_q, child_q)
    return _quat_mul(_quat_conj(parent_q), child_q)
end

"""Compute rest child-parent orientations for the robot-arm links at a reference time."""
function cloth_robot_arm_rest_quaternions(plan::RobotArmPlan, t_s::Real)
    sample = robot_arm_plan_sample(plan, t_s)
    pose = cloth_fk(plan.model, plan.base_pose, sample.q)
    rests = Vector{SVector{4, Float64}}(undef, length(plan.model.links))
    parent_q = plan.base_pose.quaternion
    for i in eachindex(plan.model.links)
        child_q = pose.link_quaternions[i]
        rests[i] = _rest_child_parent_quat(parent_q, child_q)
        parent_q = child_q
    end
    return rests
end

"""Pack the initial cloth-coupled robot-arm multibody state from a plan."""
function cloth_robot_arm_initial_state(plan::RobotArmPlan; t_s::Real=0.0)
    sample = robot_arm_plan_sample(plan, t_s)
    pose = cloth_fk(plan.model, plan.base_pose, sample.q)
    return compliant_state_vector(
        pose.link_com_positions,
        pose.link_quaternions;
        velocities=fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(pose.link_com_positions)),
        angular_velocities=fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(pose.link_com_positions)),
    )
end

"""Return body and state dimensions for a coupled cloth robot-arm plan."""
function coupled_cloth_robot_arm_state_shape(plan::RobotArmPlan)
    n = length(plan.model.links)
    return (
        arm_r=zeros(3, n),
        arm_q=zeros(4, n),
        arm_v=zeros(3, n),
        arm_ω=zeros(3, n),
    )
end

"""Write the cloth-coupled robot-arm initial state into a spacecraft view."""
function initialize_coupled_cloth_robot_arm_state!(sc_view, plan::RobotArmPlan; t_s::Real=0.0)
    hasproperty(sc_view, :arm_r) || return nothing
    sample = robot_arm_plan_sample(plan, t_s)
    base_q = hasproperty(sc_view, :q) ? _unit_quat(sc_view.q) : plan.base_pose.quaternion
    base_r = SVector{3, Float64}(sc_view.pos)
    base_v = hasproperty(sc_view, :vel) ? SVector{3, Float64}(sc_view.vel) : SVector{3, Float64}(0.0, 0.0, 0.0)
    base_ω = hasproperty(sc_view, :ω) ? SVector{3, Float64}(sc_view.ω) : SVector{3, Float64}(0.0, 0.0, 0.0)
    pose = cloth_fk(plan.model, ClothArmBasePose(base_r, base_q), sample.q)
    R_base = _rot(base_q)
    ω_world = R_base * base_ω
    for i in eachindex(plan.model.links)
        r = pose.link_com_positions[i]
        q = pose.link_quaternions[i]
        R_link = _rot(q)
        sc_view.arm_r[:, i] .= r
        sc_view.arm_q[:, i] .= q
        sc_view.arm_v[:, i] .= base_v + cross(ω_world, r - base_r)
        sc_view.arm_ω[:, i] .= R_link' * ω_world
    end
    return nothing
end

"""Build the compliant multibody model corresponding to a robot-arm plan."""
function cloth_robot_arm_multibody(
    plan::RobotArmPlan;
    k_translation_n_m=5.0e3,
    c_translation_n_s_m=30.0,
    k_rotation_n_m_rad=15.0,
    c_rotation_n_m_s_rad=0.5,
)
    arm = plan.model
    pose0 = cloth_fk(arm, plan.base_pose, plan.q_ref[:, 1])
    bodies = CompliantBody[
        CompliantBody(link.name, link.mass_kg, _link_inertia(link)) for link in arm.links
    ]
    n = length(arm.links)
    Kx = _joint_compliance_matrices(k_translation_n_m, n, :k_translation_n_m)
    Cx = _joint_compliance_matrices(c_translation_n_s_m, n, :c_translation_n_s_m)
    Kr = _joint_compliance_matrices(k_rotation_n_m_rad, n, :k_rotation_n_m_rad)
    Cr = _joint_compliance_matrices(c_rotation_n_m_s_rad, n, :c_rotation_n_m_s_rad)
    rests0 = cloth_robot_arm_rest_quaternions(plan, 0.0)
    joints = CompliantJoint[]
    for i in eachindex(arm.links)
        link = arm.links[i]
        if i == 1
            parent = 0
            parent_point = arm.mount_offset_body
            parent_q = plan.base_pose.quaternion
        else
            parent = i - 1
            parent_link = arm.links[i - 1]
            parent_point = parent_link.vector_parent - parent_link.com_offset_parent
            parent_q = pose0.link_quaternions[i - 1]
        end
        child_point = -link.com_offset_parent
        push!(joints, CompliantJoint(
            Symbol("cloth_arm_joint_$i"),
            parent,
            i,
            parent_point,
            child_point,
            Kx[i],
            Cx[i],
            Kr[i],
            Cr[i],
            rests0[i],
        ))
    end
    return CompliantMultibodyModel(bodies, joints, plan.base_pose.position, plan.base_pose.quaternion)
end

"""Build compliant joint actuators that track a robot-arm reference state."""
function cloth_robot_arm_actuators(
    plan::RobotArmPlan;
    torque_limit_n_m=Inf,
    kp_n_m_rad=0.0,
    kd_n_m_s_rad=0.0,
    feedforward_torque_child_body=(0.0, 0.0, 0.0),
    efficiency::Real=1.0,
)
    return CompliantJointActuator[
        CompliantJointActuator(
            Symbol("cloth_arm_actuator_$i"),
            i;
            torque_limit_n_m=torque_limit_n_m,
            kp_n_m_rad=kp_n_m_rad,
            kd_n_m_s_rad=kd_n_m_s_rad,
            feedforward_torque_child_body=feedforward_torque_child_body,
            efficiency=efficiency,
        ) for i in eachindex(plan.model.links)
    ]
end

"""Read one coupled robot-arm body state from a spacecraft view."""
@inline function _coupled_body_state(sc_view, i::Int)
    return (
        r=SVector{3, Float64}(sc_view.arm_r[:, i]),
        q=_unit_quat(sc_view.arm_q[:, i]),
        v=SVector{3, Float64}(sc_view.arm_v[:, i]),
        ω=SVector{3, Float64}(sc_view.arm_ω[:, i]),
    )
end

"""Return parent attachment kinematics for spacecraft-coupled arm dynamics."""
function _coupled_parent_kinematics(sc_view, parent::Int, p_body)
    if parent == 0
        q = hasproperty(sc_view, :q) ? _unit_quat(sc_view.q) : _Q_IDENTITY
        r = SVector{3, Float64}(sc_view.pos)
        v = SVector{3, Float64}(sc_view.vel)
        ω = hasproperty(sc_view, :ω) ? SVector{3, Float64}(sc_view.ω) : SVector{3, Float64}(0.0, 0.0, 0.0)
        R = _rot(q)
        return (
            r=r,
            q=q,
            R=R,
            ω=ω,
            ω_world=R * ω,
            point=r + R * p_body,
            point_velocity=_body_offset_velocity_world(v, q, ω, p_body),
        )
    end
    s = _coupled_body_state(sc_view, parent)
    R = _rot(s.q)
    return (
        r=s.r,
        q=s.q,
        R=R,
        ω=s.ω,
        ω_world=R * s.ω,
        point=s.r + R * p_body,
        point_velocity=_body_offset_velocity_world(s.v, s.q, s.ω, p_body),
    )
end

"""Write coupled cloth robot-arm derivatives into the simulation RHS vector."""
function assign_coupled_cloth_robot_arm_rhs!(
    du_view,
    sc_view,
    plan::RobotArmPlan,
    t_s::Real,
    base_forces_world,
    base_torques_body;
    k_translation_n_m=5.0e3,
    c_translation_n_s_m=30.0,
    k_rotation_n_m_rad=15.0,
    c_rotation_n_m_s_rad=0.5,
    joint_actuators::AbstractVector{CompliantJointActuator}=CompliantJointActuator[],
)
    hasproperty(sc_view, :arm_r) || return nothing
    n = length(plan.model.links)
    n == 0 && return nothing
    Kx = _joint_compliance_matrices(k_translation_n_m, n, :k_translation_n_m)
    Cx = _joint_compliance_matrices(c_translation_n_s_m, n, :c_translation_n_s_m)
    Kr = _joint_compliance_matrices(k_rotation_n_m_rad, n, :k_rotation_n_m_rad)
    Cr = _joint_compliance_matrices(c_rotation_n_m_s_rad, n, :c_rotation_n_m_s_rad)
    rests = cloth_robot_arm_rest_quaternions(plan, t_s)
    forces = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n]
    torques_body = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n]
    base_q = hasproperty(sc_view, :q) ? _unit_quat(sc_view.q) : _Q_IDENTITY
    R_base = _rot(base_q)

    for i in eachindex(plan.model.links)
        link = plan.model.links[i]
        if i == 1
            parent_idx = 0
            parent_point = plan.model.mount_offset_body
        else
            parent_idx = i - 1
            parent_link = plan.model.links[i - 1]
            parent_point = parent_link.vector_parent - parent_link.com_offset_parent
        end
        child_point = -link.com_offset_parent
        parent = _coupled_parent_kinematics(sc_view, parent_idx, parent_point)
        child = _coupled_parent_kinematics(sc_view, i, child_point)

        Δr = parent.point - child.point
        Δv = parent.point_velocity - child.point_velocity
        force_parent = -Kx[i] * Δr - Cx[i] * Δv
        force_child = -force_parent

        if parent_idx == 0
            base_forces_world .+= force_parent
            base_torques_body .+= cross(parent_point, R_base' * force_parent)
        else
            forces[parent_idx] += force_parent
            torques_body[parent_idx] += cross(parent_point, parent.R' * force_parent)
        end
        forces[i] += force_child
        torques_body[i] += cross(child_point, child.R' * force_child)

        desired_child_q = _quat_mul(parent.q, rests[i])
        q_err = _quat_mul(desired_child_q, _quat_conj(child.q))
        ϕ_world = _axis_angle_error(q_err)
        Δω_world = child.ω_world - parent.ω_world
        torque_child_world = Kr[i] * ϕ_world - Cr[i] * Δω_world
        for actuator in joint_actuators
            actuator.joint == i || continue
            pd_world = actuator.kp_n_m_rad * ϕ_world - actuator.kd_n_m_s_rad * Δω_world
            raw_child_body = child.R' * pd_world + actuator.feedforward_torque_child_body
            limited_child_body = clamp.(raw_child_body, -actuator.torque_limit_n_m, actuator.torque_limit_n_m)
            torque_child_world += child.R * (actuator.efficiency * limited_child_body)
        end
        torque_parent_world = -torque_child_world
        if parent_idx == 0
            base_torques_body .+= R_base' * torque_parent_world
        else
            torques_body[parent_idx] += parent.R' * torque_parent_world
        end
        torques_body[i] += child.R' * torque_child_world
    end

    for i in eachindex(plan.model.links)
        link = plan.model.links[i]
        state = _coupled_body_state(sc_view, i)
        J = _link_inertia(link)
        qdot = 0.5 .* _quat_raw_mul(state.q, SVector{4, Float64}(state.ω[1], state.ω[2], state.ω[3], 0.0))
        du_view.arm_r[:, i] .= state.v
        du_view.arm_q[:, i] .= qdot
        du_view.arm_v[:, i] .= forces[i] / link.mass_kg
        du_view.arm_ω[:, i] .= J \ (torques_body[i] - cross(state.ω, J * state.ω))
    end
    return nothing
end

"""Return the end-effector position and orientation from a coupled arm state."""
function cloth_robot_arm_end_effector(plan::RobotArmPlan, x::AbstractVector)
    n = length(plan.model.links)
    n > 0 || return plan.base_pose.position
    state = compliant_state_parts(x, n)
    link = plan.model.links[end]
    tip_body = link.vector_parent - link.com_offset_parent
    return state.r + _rot(state.q) * tip_body
end

"""Simulate a robot-arm plan with compliant cloth-style multibody dynamics."""
function simulate_cloth_robot_arm_plan(
    plan::RobotArmPlan;
    dt_s::Real=(length(plan.t_ref_s) >= 2 ? plan.t_ref_s[2] - plan.t_ref_s[1] : 0.1),
    duration_s::Real=plan.t_ref_s[end],
    integrator::Symbol=:implicit_midpoint,
    k_translation_n_m=5.0e3,
    c_translation_n_s_m=30.0,
    k_rotation_n_m_rad=15.0,
    c_rotation_n_m_s_rad=0.5,
    joint_actuators::Union{Nothing, AbstractVector{CompliantJointActuator}}=nothing,
    actuator_torque_limit_n_m=Inf,
    actuator_kp_n_m_rad=0.0,
    actuator_kd_n_m_s_rad=0.0,
)
    model = cloth_robot_arm_multibody(
        plan;
        k_translation_n_m=k_translation_n_m,
        c_translation_n_s_m=c_translation_n_s_m,
        k_rotation_n_m_rad=k_rotation_n_m_rad,
        c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    )
    actuators = joint_actuators === nothing ?
        cloth_robot_arm_actuators(
            plan;
            torque_limit_n_m=actuator_torque_limit_n_m,
            kp_n_m_rad=actuator_kp_n_m_rad,
            kd_n_m_s_rad=actuator_kd_n_m_s_rad,
        ) :
        collect(joint_actuators)
    dt = Float64(dt_s)
    duration = Float64(duration_s)
    times = collect(0.0:dt:duration)
    if isempty(times) || times[end] < duration
        push!(times, duration)
    end
    states = Vector{Vector{Float64}}(undef, length(times))
    states[1] = cloth_robot_arm_initial_state(plan; t_s=0.0)
    for k in 2:length(times)
        t_prev = times[k - 1]
        rests = cloth_robot_arm_rest_quaternions(plan, t_prev)
        stepper = if integrator === :implicit_midpoint
            step_compliant_multibody_implicit_midpoint
        elseif integrator === :rk4
            step_compliant_multibody_rk4
        else
            throw(ArgumentError("Unsupported Cloth robot arm integrator $(repr(integrator)); expected :implicit_midpoint or :rk4."))
        end
        states[k] = stepper(
            model,
            states[k - 1],
            t_prev,
            times[k] - t_prev;
            joint_rest_quaternions=rests,
            joint_actuators=actuators,
        )
    end
    ee = zeros(3, length(times))
    ee_ref = zeros(3, length(times))
    err = zeros(length(times))
    nj = length(model.joints)
    compliance_torque = zeros(3, nj, length(times))
    actuator_torque = zeros(3, nj, length(times))
    total_torque = zeros(3, nj, length(times))
    for (k, t) in pairs(times)
        ee[:, k] .= cloth_robot_arm_end_effector(plan, states[k])
        ee_ref[:, k] .= robot_arm_plan_sample(plan, t).ee
        err[k] = norm(ee[:, k] - ee_ref[:, k])
        loads = compliant_joint_loads(
            model,
            states[k];
            joint_rest_quaternions=cloth_robot_arm_rest_quaternions(plan, t),
            joint_actuators=actuators,
        )
        for (j, load) in pairs(loads)
            compliance_torque[:, j, k] .= load.compliance_torque_child_body
            actuator_torque[:, j, k] .= load.actuator_torque_child_body
            total_torque[:, j, k] .= load.total_torque_child_body
        end
    end
    return ClothRobotArmSimulation(
        model,
        CompliantMultibodyTrajectory(times, states),
        ee,
        ee_ref,
        err,
        compliance_torque,
        actuator_torque,
        total_torque,
    )
end

end # module ClothRobotArmDynamics
