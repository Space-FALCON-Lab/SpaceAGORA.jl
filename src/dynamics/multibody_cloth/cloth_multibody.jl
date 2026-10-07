"""Compliant multibody dynamics module for cloth-like panels, topologies, joint loads, and integration."""
module ClothMultibody

using LinearAlgebra
using StaticArrays

export CompliantBody, CompliantJoint, CompliantMultibodyModel, CompliantMultibodyTrajectory
export CompliantTopologyNode, CompliantTopologyEdge, CompliantTopologyBuild
export CompliantJointActuator, CompliantJointLoad
export rectangular_prism_inertia, thin_panel_inertia
export build_compliant_topology, build_rectangular_compliant_grid
export compliant_state_vector, compliant_state_parts, compliant_multibody_dynamics
export compliant_joint_loads
export CompliantMountKinematics, compliant_mount_kinematics, compliant_joint_loads_in_place!
export compliant_body_derivatives!, compliant_rest_state, compliant_state_in_mount_frame
export compliant_state_to_inertial!, compliant_mount_wrench, compliant_world_to_body
export CompiledCompliantModel, compile_compliant_model
export step_compliant_multibody_rk4, step_compliant_multibody_implicit_midpoint
export simulate_compliant_multibody

const _Q_IDENTITY = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

"""Rigid body properties used by the compliant multibody cloth model."""
struct CompliantBody
    name::Symbol
    mass_kg::Float64
    inertia_body_kg_m2::SMatrix{3, 3, Float64}
end

"""Compliant connection between a parent body and child body."""
struct CompliantJoint
    name::Symbol
    parent::Int
    child::Int
    parent_point_body::SVector{3, Float64}
    child_point_body::SVector{3, Float64}
    k_translation_n_m::SMatrix{3, 3, Float64}
    c_translation_n_s_m::SMatrix{3, 3, Float64}
    k_rotation_n_m_rad::SMatrix{3, 3, Float64}
    c_rotation_n_m_s_rad::SMatrix{3, 3, Float64}
    rest_child_parent_quat::SVector{4, Float64}
end

"""Body and joint collection defining a compliant multibody system."""
struct CompliantMultibodyModel
    bodies::Vector{CompliantBody}
    joints::Vector{CompliantJoint}
    base_position::SVector{3, Float64}
    base_quaternion::SVector{4, Float64}
end

"""Time history returned by compliant multibody simulation."""
struct CompliantMultibodyTrajectory
    t_s::Vector{Float64}
    states::Vector{Vector{Float64}}
end

"""Node specification used to build compliant multibody topologies."""
struct CompliantTopologyNode
    name::Symbol
    mass_kg::Float64
    inertia_body_kg_m2::SMatrix{3, 3, Float64}
    position::SVector{3, Float64}
    quaternion::SVector{4, Float64}
    velocity::SVector{3, Float64}
    angular_velocity::SVector{3, Float64}
end

"""Edge specification used to connect compliant topology nodes."""
struct CompliantTopologyEdge
    name::Symbol
    parent::Int
    child::Int
    parent_point_body::SVector{3, Float64}
    child_point_body::SVector{3, Float64}
    k_translation_n_m::SMatrix{3, 3, Float64}
    c_translation_n_s_m::SMatrix{3, 3, Float64}
    k_rotation_n_m_rad::SMatrix{3, 3, Float64}
    c_rotation_n_m_s_rad::SMatrix{3, 3, Float64}
    rest_child_parent_quat::Union{Nothing, SVector{4, Float64}}
end

"""Built compliant topology with model, initial state, and name-to-index lookup."""
struct CompliantTopologyBuild
    model::CompliantMultibodyModel
    initial_state::Vector{Float64}
    node_names::Vector{Symbol}
    joint_names::Vector{Symbol}
end

"""Commanded joint actuator torque applied in the compliant multibody model."""
struct CompliantJointActuator
    name::Symbol
    joint::Int
    torque_limit_n_m::SVector{3, Float64}
    kp_n_m_rad::SMatrix{3, 3, Float64}
    kd_n_m_s_rad::SMatrix{3, 3, Float64}
    feedforward_torque_child_body::SVector{3, Float64}
    efficiency::Float64
end

"""Force and torque diagnostics for one compliant joint."""
struct CompliantJointLoad
    name::Symbol
    parent::Int
    child::Int
    translation_force_parent_world::SVector{3, Float64}
    translation_force_child_world::SVector{3, Float64}
    compliance_torque_parent_world::SVector{3, Float64}
    compliance_torque_child_world::SVector{3, Float64}
    actuator_torque_parent_world::SVector{3, Float64}
    actuator_torque_child_world::SVector{3, Float64}
    compliance_torque_child_body::SVector{3, Float64}
    actuator_torque_child_body::SVector{3, Float64}
    total_torque_child_body::SVector{3, Float64}
end

"""Normalize a quaternion and return it as a static vector."""
@inline function _unit_quat(q)::SVector{4, Float64}
    qv = SVector{4, Float64}(q)
    nq = norm(qv)
    return isfinite(nq) && nq > eps(Float64) ? qv / nq : _Q_IDENTITY
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

"""Convert a quaternion error into an axis-angle error vector."""
function _axis_angle_error(q_err_raw)::SVector{3, Float64}
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

"""Compute world velocity of a body-fixed offset point."""
@inline function _body_offset_velocity(v_world, q, ω_body, p_body)
    return SVector{3, Float64}(v_world) + _rot(q) * cross(SVector{3, Float64}(ω_body), p_body)
end

@inline _vec3(v)::SVector{3, Float64} = SVector{3, Float64}(v)
@inline _diag3(v::Real)::SMatrix{3, 3, Float64} = SMatrix{3, 3, Float64}(Float64(v) * I)
@inline _mat3(m)::SMatrix{3, 3, Float64} = m isa Real ? _diag3(m) : SMatrix{3, 3, Float64}(m)

"""Return the body-frame inertia matrix for a rectangular prism."""
function rectangular_prism_inertia(mass_kg::Real, dimensions_m)::SMatrix{3, 3, Float64}
    d = _vec3(dimensions_m)
    m = Float64(mass_kg)
    all(d .> 0.0) || throw(ArgumentError("rectangular prism dimensions must be positive."))
    return SMatrix{3, 3, Float64}(
        (m / 12.0) * (d[2]^2 + d[3]^2), 0.0, 0.0,
        0.0, (m / 12.0) * (d[1]^2 + d[3]^2), 0.0,
        0.0, 0.0, (m / 12.0) * (d[1]^2 + d[2]^2),
    )
end

thin_panel_inertia(mass_kg::Real, width_m::Real, height_m::Real; thickness_m::Real=1.0e-3) =
    rectangular_prism_inertia(mass_kg, (width_m, height_m, thickness_m))

function CompliantTopologyNode(
    name::Symbol;
    mass_kg::Real,
    inertia_body_kg_m2,
    position,
    quaternion=_Q_IDENTITY,
    velocity=(0.0, 0.0, 0.0),
    angular_velocity=(0.0, 0.0, 0.0),
)
    mass = Float64(mass_kg)
    mass > 0.0 || throw(ArgumentError("node mass_kg must be positive."))
    return CompliantTopologyNode(
        name,
        mass,
        _mat3(inertia_body_kg_m2),
        _vec3(position),
        _unit_quat(quaternion),
        _vec3(velocity),
        _vec3(angular_velocity),
    )
end

function CompliantTopologyEdge(
    name::Symbol,
    parent::Integer,
    child::Integer;
    parent_point_body,
    child_point_body,
    k_translation_n_m=0.0,
    c_translation_n_s_m=0.0,
    k_rotation_n_m_rad=0.0,
    c_rotation_n_m_s_rad=0.0,
    rest_child_parent_quat=nothing,
)
    return CompliantTopologyEdge(
        name,
        Int(parent),
        Int(child),
        _vec3(parent_point_body),
        _vec3(child_point_body),
        _mat3(k_translation_n_m),
        _mat3(c_translation_n_s_m),
        _mat3(k_rotation_n_m_rad),
        _mat3(c_rotation_n_m_s_rad),
        rest_child_parent_quat === nothing ? nothing : _unit_quat(rest_child_parent_quat),
    )
end

function CompliantJointActuator(
    name::Symbol,
    joint::Integer;
    torque_limit_n_m=Inf,
    kp_n_m_rad=0.0,
    kd_n_m_s_rad=0.0,
    feedforward_torque_child_body=(0.0, 0.0, 0.0),
    efficiency::Real=1.0,
)
    limit = torque_limit_n_m isa Real ?
        SVector{3, Float64}(fill(Float64(torque_limit_n_m), 3)) :
        abs.(_vec3(torque_limit_n_m))
    all(limit .>= 0.0) || throw(ArgumentError("actuator torque limits must be nonnegative."))
    eff = Float64(efficiency)
    eff >= 0.0 || throw(ArgumentError("actuator efficiency must be nonnegative."))
    return CompliantJointActuator(
        name,
        Int(joint),
        limit,
        _mat3(kp_n_m_rad),
        _mat3(kd_n_m_s_rad),
        _vec3(feedforward_torque_child_body),
        eff,
    )
end

"""Compute the rest child-to-parent orientation between two body quaternions."""
function _rest_child_parent_quat(parent_q, child_q)
    return _quat_mul(_quat_conj(parent_q), child_q)
end

"""Build a compliant multibody model and initial state from topology nodes and edges."""
function build_compliant_topology(
    nodes::AbstractVector{CompliantTopologyNode},
    edges::AbstractVector{CompliantTopologyEdge};
    base_position=(0.0, 0.0, 0.0),
    base_quaternion=_Q_IDENTITY,
)::CompliantTopologyBuild
    n = length(nodes)
    n > 0 || throw(ArgumentError("at least one topology node is required."))
    base_q = _unit_quat(base_quaternion)
    bodies = CompliantBody[
        CompliantBody(node.name, node.mass_kg, node.inertia_body_kg_m2) for node in nodes
    ]
    joints = CompliantJoint[]
    for edge in edges
        1 <= edge.child <= n || throw(ArgumentError("edge $(edge.name) child index $(edge.child) is outside 1:$n."))
        0 <= edge.parent <= n || throw(ArgumentError("edge $(edge.name) parent index $(edge.parent) is outside 0:$n."))
        edge.parent != edge.child || throw(ArgumentError("edge $(edge.name) cannot connect a node to itself."))
        parent_q = edge.parent == 0 ? base_q : nodes[edge.parent].quaternion
        child_q = nodes[edge.child].quaternion
        rest = edge.rest_child_parent_quat === nothing ?
            _rest_child_parent_quat(parent_q, child_q) :
            edge.rest_child_parent_quat
        push!(joints, CompliantJoint(
            edge.name,
            edge.parent,
            edge.child,
            edge.parent_point_body,
            edge.child_point_body,
            edge.k_translation_n_m,
            edge.c_translation_n_s_m,
            edge.k_rotation_n_m_rad,
            edge.c_rotation_n_m_s_rad,
            rest,
        ))
    end
    x0 = compliant_state_vector(
        [node.position for node in nodes],
        [node.quaternion for node in nodes];
        velocities=[node.velocity for node in nodes],
        angular_velocities=[node.angular_velocity for node in nodes],
    )
    return CompliantTopologyBuild(
        CompliantMultibodyModel(bodies, joints, _vec3(base_position), base_q),
        x0,
        [node.name for node in nodes],
        [edge.name for edge in edges],
    )
end

"""Return either a scalar grid parameter or its per-index value."""
@inline function _grid_value(v, idx::Int)
    if v isa AbstractVector
        return v[idx]
    elseif v isa AbstractMatrix
        return v[idx]
    else
        return v
    end
end

"""Build a rectangular compliant grid of panel bodies and spring-damper joints."""
function build_rectangular_compliant_grid(
    rows::Integer,
    cols::Integer;
    spacing_m=(1.0, 1.0),
    tile_size_m=spacing_m,
    thickness_m::Real=1.0e-3,
    mass_kg=1.0,
    k_translation_n_m=5.0e3,
    c_translation_n_s_m=30.0,
    k_rotation_n_m_rad=15.0,
    c_rotation_n_m_s_rad=0.5,
    anchor_index::Union{Nothing, Integer}=nothing,
    anchor_k_translation_n_m=k_translation_n_m,
    anchor_c_translation_n_s_m=c_translation_n_s_m,
    anchor_k_rotation_n_m_rad=k_rotation_n_m_rad,
    anchor_c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
    base_position=(0.0, 0.0, 0.0),
    base_quaternion=_Q_IDENTITY,
)::CompliantTopologyBuild
    rcount = Int(rows)
    ccount = Int(cols)
    rcount > 0 && ccount > 0 || throw(ArgumentError("grid rows and cols must be positive."))
    spacing = SVector{2, Float64}(spacing_m)
    tile = SVector{2, Float64}(tile_size_m)
    all(spacing .> 0.0) || throw(ArgumentError("grid spacing must be positive."))
    all(tile .> 0.0) || throw(ArgumentError("tile size must be positive."))

    nodes = CompliantTopologyNode[]
    for r in 1:rcount, c in 1:ccount
        idx = (r - 1) * ccount + c
        x = (c - 1 - 0.5 * (ccount - 1)) * spacing[1]
        y = (r - 1 - 0.5 * (rcount - 1)) * spacing[2]
        m = Float64(_grid_value(mass_kg, idx))
        push!(nodes, CompliantTopologyNode(
            Symbol("panel_$(r)_$(c)");
            mass_kg=m,
            inertia_body_kg_m2=thin_panel_inertia(m, tile[1], tile[2]; thickness_m=thickness_m),
            position=SVector{3, Float64}(x, y, 0.0),
        ))
    end

    edges = CompliantTopologyEdge[]
    for r in 1:rcount, c in 1:ccount
        idx = (r - 1) * ccount + c
        if c < ccount
            right = idx + 1
            push!(edges, CompliantTopologyEdge(
                Symbol("panel_edge_$(r)_$(c)_right"),
                idx,
                right;
                parent_point_body=SVector{3, Float64}(0.5 * spacing[1], 0.0, 0.0),
                child_point_body=SVector{3, Float64}(-0.5 * spacing[1], 0.0, 0.0),
                k_translation_n_m=k_translation_n_m,
                c_translation_n_s_m=c_translation_n_s_m,
                k_rotation_n_m_rad=k_rotation_n_m_rad,
                c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
            ))
        end
        if r < rcount
            below = idx + ccount
            push!(edges, CompliantTopologyEdge(
                Symbol("panel_edge_$(r)_$(c)_down"),
                idx,
                below;
                parent_point_body=SVector{3, Float64}(0.0, 0.5 * spacing[2], 0.0),
                child_point_body=SVector{3, Float64}(0.0, -0.5 * spacing[2], 0.0),
                k_translation_n_m=k_translation_n_m,
                c_translation_n_s_m=c_translation_n_s_m,
                k_rotation_n_m_rad=k_rotation_n_m_rad,
                c_rotation_n_m_s_rad=c_rotation_n_m_s_rad,
            ))
        end
    end

    if anchor_index !== nothing
        aidx = Int(anchor_index)
        1 <= aidx <= length(nodes) || throw(ArgumentError("anchor_index must be inside the grid node range."))
        push!(edges, CompliantTopologyEdge(
            :base_anchor,
            0,
            aidx;
            parent_point_body=nodes[aidx].position,
            child_point_body=SVector{3, Float64}(0.0, 0.0, 0.0),
            k_translation_n_m=anchor_k_translation_n_m,
            c_translation_n_s_m=anchor_c_translation_n_s_m,
            k_rotation_n_m_rad=anchor_k_rotation_n_m_rad,
            c_rotation_n_m_s_rad=anchor_c_rotation_n_m_s_rad,
        ))
    end

    return build_compliant_topology(nodes, edges; base_position=base_position, base_quaternion=base_quaternion)
end

"""Pack compliant body positions, quaternions, and velocities into a flat state vector."""
function compliant_state_vector(
    positions,
    quaternions;
    velocities=fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(positions)),
    angular_velocities=fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(positions)),
)
    n = length(positions)
    length(quaternions) == n || throw(ArgumentError("quaternions length must match positions."))
    length(velocities) == n || throw(ArgumentError("velocities length must match positions."))
    length(angular_velocities) == n || throw(ArgumentError("angular_velocities length must match positions."))
    x = zeros(13n)
    for i in 1:n
        b = 13 * (i - 1)
        x[(b + 1):(b + 3)] .= SVector{3, Float64}(positions[i])
        x[(b + 4):(b + 7)] .= _unit_quat(quaternions[i])
        x[(b + 8):(b + 10)] .= SVector{3, Float64}(velocities[i])
        x[(b + 11):(b + 13)] .= SVector{3, Float64}(angular_velocities[i])
    end
    return x
end

"""Return the state-vector slices for one compliant body."""
function compliant_state_parts(x::AbstractVector, i::Int)
    b = 13 * (i - 1)
    return (
        r=SVector{3, Float64}(x[(b + 1):(b + 3)]),
        q=_unit_quat(x[(b + 4):(b + 7)]),
        v=SVector{3, Float64}(x[(b + 8):(b + 10)]),
        ω=SVector{3, Float64}(x[(b + 11):(b + 13)]),
    )
end

"""Normalize all body quaternions in a compliant state vector."""
function _normalize_state_quaternions!(x::AbstractVector)
    n = length(x) ÷ 13
    for i in 1:n
        b = 13 * (i - 1)
        x[(b + 4):(b + 7)] .= _unit_quat(x[(b + 4):(b + 7)])
    end
    return x
end

"""Kinematics of a body-fixed attachment point from a body state (position, attitude, velocity, body rate)."""
@inline function _body_kinematics(r, q, v, ω, p_body)
    R = _rot(q)
    return (
        r=r,
        q=q,
        R=R,
        ω=ω,
        ω_world=R * ω,
        point=r + R * p_body,
        point_velocity=_body_offset_velocity(v, q, ω, p_body),
    )
end

"""Return parent attachment kinematics for a compliant joint."""
function _parent_kinematics(model::CompliantMultibodyModel, x::AbstractVector, parent::Int, p_body)
    if parent == 0
        q = model.base_quaternion
        r = model.base_position
        ω = SVector{3, Float64}(0.0, 0.0, 0.0)
        v = SVector{3, Float64}(0.0, 0.0, 0.0)
        R = _rot(q)
        return (
            r=r,
            q=q,
            R=R,
            ω=ω,
            ω_world=R * ω,
            point=r + R * p_body,
            point_velocity=v,
        )
    end
    s = compliant_state_parts(x, parent)
    return _body_kinematics(s.r, s.q, s.v, s.ω, p_body)
end

"""Return the rest orientation for a compliant joint instance."""
@inline function _joint_rest(joint::CompliantJoint, jidx::Int, joint_rest_quaternions)
    return joint_rest_quaternions === nothing ? joint.rest_child_parent_quat : _unit_quat(joint_rest_quaternions[jidx])
end

"""Transform actuator torque from child joint coordinates into world coordinates."""
function _actuator_torque_child_world(
    actuator,
    child,
    ϕ_world::SVector{3, Float64},
    Δω_world::SVector{3, Float64},
)::SVector{3, Float64}
    pd_world = actuator.kp_n_m_rad * ϕ_world - actuator.kd_n_m_s_rad * Δω_world
    raw_child_body = child.R' * pd_world + actuator.feedforward_torque_child_body
    limited_child_body = clamp.(raw_child_body, -actuator.torque_limit_n_m, actuator.torque_limit_n_m)
    return child.R * (actuator.efficiency * limited_child_body)
end

"""
Evaluate one compliant joint from parent and child attachment kinematics. Shared by
`compliant_joint_loads` and `compliant_joint_loads_in_place!`; allocation-free.
"""
@inline function _joint_load_kernel(joint, jidx::Int, parent, child, rest, joint_actuators)
    Δr = parent.point - child.point
    Δv = parent.point_velocity - child.point_velocity
    force_parent = -joint.k_translation_n_m * Δr - joint.c_translation_n_s_m * Δv
    force_child = -force_parent

    desired_child_q = _quat_mul(parent.q, rest)
    q_err = _quat_mul(desired_child_q, _quat_conj(child.q))
    ϕ_world = _axis_angle_error(q_err)
    Δω_world = child.ω_world - parent.ω_world
    compliance_child_world = joint.k_rotation_n_m_rad * ϕ_world - joint.c_rotation_n_m_s_rad * Δω_world
    actuator_child_world = SVector{3, Float64}(0.0, 0.0, 0.0)
    for actuator in joint_actuators
        actuator.joint == jidx || continue
        actuator_child_world += _actuator_torque_child_world(actuator, child, ϕ_world, Δω_world)
    end
    return (
        force_parent=force_parent,
        force_child=force_child,
        compliance_child_world=compliance_child_world,
        actuator_child_world=actuator_child_world,
    )
end

"""Compute compliant joint forces, torques, and diagnostics for a state."""
function compliant_joint_loads(
    model::CompliantMultibodyModel,
    x::AbstractVector;
    joint_rest_quaternions::Union{Nothing, AbstractVector}=nothing,
    joint_actuators::AbstractVector{CompliantJointActuator}=CompliantJointActuator[],
)::Vector{CompliantJointLoad}
    n = length(model.bodies)
    length(x) == 13n || throw(ArgumentError("state length must be 13 * number of bodies."))
    loads = CompliantJointLoad[]
    for (jidx, joint) in pairs(model.joints)
        parent = _parent_kinematics(model, x, joint.parent, joint.parent_point_body)
        child = _parent_kinematics(model, x, joint.child, joint.child_point_body)
        rest = _joint_rest(joint, jidx, joint_rest_quaternions)
        k = _joint_load_kernel(joint, jidx, parent, child, rest, joint_actuators)

        compliance_child_world = k.compliance_child_world
        actuator_child_world = k.actuator_child_world
        compliance_parent_world = -compliance_child_world
        actuator_parent_world = -actuator_child_world
        compliance_child_body = child.R' * compliance_child_world
        actuator_child_body = child.R' * actuator_child_world
        push!(loads, CompliantJointLoad(
            joint.name,
            joint.parent,
            joint.child,
            k.force_parent,
            k.force_child,
            compliance_parent_world,
            compliance_child_world,
            actuator_parent_world,
            actuator_child_world,
            compliance_child_body,
            actuator_child_body,
            compliance_child_body + actuator_child_body,
        ))
    end
    return loads
end

"""Evaluate the compliant multibody state derivative."""
function compliant_multibody_dynamics(
    model::CompliantMultibodyModel,
    x::AbstractVector,
    t::Real=0.0;
    joint_rest_quaternions::Union{Nothing, AbstractVector}=nothing,
    joint_actuators::AbstractVector{CompliantJointActuator}=CompliantJointActuator[],
    external_forces_world=fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(model.bodies)),
    external_torques_body=fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(model.bodies)),
)
    n = length(model.bodies)
    length(x) == 13n || throw(ArgumentError("state length must be 13 * number of bodies."))
    forces = [SVector{3, Float64}(external_forces_world[i]) for i in 1:n]
    torques_body = [SVector{3, Float64}(external_torques_body[i]) for i in 1:n]

    loads = compliant_joint_loads(model, x; joint_rest_quaternions=joint_rest_quaternions, joint_actuators=joint_actuators)
    for (jidx, joint) in pairs(model.joints)
        load = loads[jidx]
        parent = _parent_kinematics(model, x, joint.parent, joint.parent_point_body)
        child = _parent_kinematics(model, x, joint.child, joint.child_point_body)
        if joint.parent != 0
            forces[joint.parent] += load.translation_force_parent_world
            torques_body[joint.parent] +=
                cross(joint.parent_point_body, parent.R' * load.translation_force_parent_world) +
                parent.R' * (load.compliance_torque_parent_world + load.actuator_torque_parent_world)
        end
        forces[joint.child] += load.translation_force_child_world
        torques_body[joint.child] +=
            cross(joint.child_point_body, child.R' * load.translation_force_child_world) +
            child.R' * (load.compliance_torque_child_world + load.actuator_torque_child_world)
    end

    dx = zeros(Float64, 13n)
    for i in 1:n
        body = model.bodies[i]
        s = compliant_state_parts(x, i)
        b = 13 * (i - 1)
        dx[(b + 1):(b + 3)] .= s.v
        dx[(b + 4):(b + 7)] .= 0.5 .* _quat_raw_mul(s.q, SVector{4, Float64}(s.ω[1], s.ω[2], s.ω[3], 0.0))
        dx[(b + 8):(b + 10)] .= forces[i] / body.mass_kg
        J = body.inertia_body_kg_m2
        dx[(b + 11):(b + 13)] .= J \ (torques_body[i] - cross(s.ω, J * s.ω))
    end
    return dx
end

"""Advance compliant multibody dynamics with one RK4 step."""
function step_compliant_multibody_rk4(
    model::CompliantMultibodyModel,
    x::AbstractVector,
    t::Real,
    dt::Real;
    dynamics_kwargs...,
)
    h = Float64(dt)
    x0 = Vector{Float64}(x)
    f = (xx, tt) -> compliant_multibody_dynamics(model, xx, tt; dynamics_kwargs...)
    k1 = f(x0, t)
    k2 = f(x0 .+ 0.5h .* k1, t + 0.5h)
    k3 = f(x0 .+ 0.5h .* k2, t + 0.5h)
    k4 = f(x0 .+ h .* k3, t + h)
    xn = x0 .+ (h / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    return _normalize_state_quaternions!(xn)
end

"""Advance compliant multibody dynamics with one implicit-midpoint step."""
function step_compliant_multibody_implicit_midpoint(
    model::CompliantMultibodyModel,
    x::AbstractVector,
    t::Real,
    dt::Real;
    tol::Float64=1.0e-9,
    max_iters::Int=12,
    dynamics_kwargs...,
)
    h = Float64(dt)
    x0 = Vector{Float64}(x)
    xn = x0 .+ h .* compliant_multibody_dynamics(model, x0, t; dynamics_kwargs...)

    # Return the implicit-midpoint nonlinear residual for a candidate next state.
    function residual(z)
        mid = 0.5 .* (x0 .+ z)
        return z .- x0 .- h .* compliant_multibody_dynamics(model, mid, t + 0.5h; dynamics_kwargs...)
    end

    r = residual(xn)
    for _ in 1:max_iters
        maximum(abs.(r)) <= tol && return _normalize_state_quaternions!(xn)
        m = length(xn)
        J = Matrix{Float64}(I, m, m)
        for j in 1:m
            step = sqrt(eps(Float64)) * max(1.0, abs(xn[j]))
            z = copy(xn)
            z[j] += step
            J[:, j] .= (residual(z) .- r) ./ step
        end
        Δ = -(J \ r)
        α = 1.0
        accepted = false
        rnorm = maximum(abs.(r))
        while α >= 1.0e-3
            trial = copy(xn)
            trial .+= α .* Δ
            _normalize_state_quaternions!(trial)
            r_trial = residual(trial)
            if maximum(abs.(r_trial)) < rnorm
                xn = trial
                r = r_trial
                accepted = true
                break
            end
            α *= 0.5
        end
        accepted || break
    end
    maximum(abs.(r)) <= 100tol ||
        throw(ErrorException("implicit midpoint failed to converge; residual=$(maximum(abs.(r)))."))
    return _normalize_state_quaternions!(xn)
end

"""Simulate compliant multibody dynamics over a fixed time step sequence."""
function simulate_compliant_multibody(
    model::CompliantMultibodyModel,
    x0::AbstractVector;
    dt_s::Real,
    duration_s::Real,
    integrator::Symbol=:implicit_midpoint,
    dynamics_kwargs...,
)::CompliantMultibodyTrajectory
    dt = Float64(dt_s)
    duration = Float64(duration_s)
    dt > 0.0 || throw(ArgumentError("dt_s must be positive."))
    duration >= 0.0 || throw(ArgumentError("duration_s must be nonnegative."))
    times = collect(0.0:dt:duration)
    if isempty(times) || times[end] < duration
        push!(times, duration)
    end
    states = Vector{Vector{Float64}}(undef, length(times))
    states[1] = _normalize_state_quaternions!(Vector{Float64}(x0))
    for k in 2:length(times)
        stepper = if integrator === :implicit_midpoint
            step_compliant_multibody_implicit_midpoint
        elseif integrator === :rk4
            step_compliant_multibody_rk4
        else
            throw(ArgumentError("Unsupported Cloth multibody integrator $(repr(integrator)); expected :implicit_midpoint or :rk4."))
        end
        states[k] = stepper(model, states[k - 1], times[k - 1], times[k] - times[k - 1]; dynamics_kwargs...)
    end
    return CompliantMultibodyTrajectory(times, states)
end

# ---------------------------------------------------------------------------
# Mounted use: a compliant model attached to a moving frame (spacecraft link)
# ---------------------------------------------------------------------------
#
# When a model is mounted on a moving link, joints with `parent == 0` attach to the MOUNT frame
# instead of the fixed base: the model's `base_position`/`base_quaternion` are not used. The body
# states live in column-per-body matrices (absolute inertial position and velocity, scalar-last
# attitude quaternion, body-frame angular velocity), the layout of the engine state.

"""Inertial kinematics of a mount frame: origin `r`, attitude `q`, origin velocity `v`, angular velocity `ω` in mount-frame components."""
struct CompliantMountKinematics
    r::SVector{3, Float64}
    q::SVector{4, Float64}
    v::SVector{3, Float64}
    ω::SVector{3, Float64}
end

"""
    compliant_mount_kinematics(body_pos, body_vel, body_q, body_ω, link_offset, link_q, mount_point, mount_q)

Mount-frame kinematics for a mount at `mount_point` (link frame) with orientation `mount_q` relative to
the link frame, on a link whose frame origin sits at `link_offset` (body frame) with attitude `link_q`
relative to a dynamic body with inertial pose `body_pos`/`body_q`, velocity `body_vel` and body-frame
rate `body_ω`. Quaternions are scalar-last, active body-to-inertial.
"""
@inline function compliant_mount_kinematics(body_pos, body_vel, body_q, body_ω, link_offset, link_q, mount_point, mount_q)
    qb = _unit_quat(body_q)
    Rb = _rot(qb)
    ω_world = Rb * SVector{3, Float64}(body_ω)
    e_link = Rb * SVector{3, Float64}(link_offset)
    q_link = _quat_mul(qb, link_q)
    e_mount = _rot(q_link) * SVector{3, Float64}(mount_point)
    q_mount = _quat_mul(q_link, mount_q)
    e = e_link + e_mount
    return CompliantMountKinematics(
        SVector{3, Float64}(body_pos) + e,
        q_mount,
        SVector{3, Float64}(body_vel) + cross(ω_world, e),
        _rot(q_mount)' * ω_world,
    )
end

"""
    compliant_mount_wrench(mount, body_ref_pos, force_world, torque_world) -> (force_world, torque_world)

Move a reaction (force at the mount origin plus a torque, both inertial) to the reference point
`body_ref_pos` of the body the mount is fixed to: the force is unchanged and the torque gains
`(mount.r - body_ref_pos) × force`.
"""
@inline function compliant_mount_wrench(mount::CompliantMountKinematics, body_ref_pos, force_world, torque_world)
    f = SVector{3, Float64}(force_world)
    return f, cross(mount.r - SVector{3, Float64}(body_ref_pos), f) + SVector{3, Float64}(torque_world)
end

"""Components in the body frame of attitude `q` (active, scalar-last) of an inertial vector."""
@inline compliant_world_to_body(q, v_world) = _rot(q)' * SVector{3, Float64}(v_world)

@inline function _mount_parent_kinematics(mount::CompliantMountKinematics, p_body)
    return _body_kinematics(mount.r, mount.q, mount.v, mount.ω, p_body)
end

@inline function _column_kinematics(pos, quat, vel, angvel, c::Int, p_body)
    return _body_kinematics(
        SVector{3, Float64}(pos[1, c], pos[2, c], pos[3, c]),
        _unit_quat(SVector{4, Float64}(quat[1, c], quat[2, c], quat[3, c], quat[4, c])),
        SVector{3, Float64}(vel[1, c], vel[2, c], vel[3, c]),
        SVector{3, Float64}(angvel[1, c], angvel[2, c], angvel[3, c]),
        p_body,
    )
end

@inline function _add_column!(A::AbstractMatrix, c::Int, v)
    A[1, c] += v[1]
    A[2, c] += v[2]
    A[3, c] += v[3]
    return nothing
end

struct CompiledJoint
    parent::Int
    child::Int
    parent_point_body::SVector{3, Float64}
    child_point_body::SVector{3, Float64}
    k_translation_n_m::SMatrix{3, 3, Float64, 9}
    c_translation_n_s_m::SMatrix{3, 3, Float64, 9}
    k_rotation_n_m_rad::SMatrix{3, 3, Float64, 9}
    c_rotation_n_m_s_rad::SMatrix{3, 3, Float64, 9}
end

struct CompiledActuator
    joint::Int
    torque_limit_n_m::SVector{3, Float64}
    kp_n_m_rad::SMatrix{3, 3, Float64, 9}
    kd_n_m_s_rad::SMatrix{3, 3, Float64, 9}
    feedforward_torque_child_body::SVector{3, Float64}
    efficiency::Float64
end

"""
    CompiledCompliantModel

Concretely typed copy of a model's bodies, joints and joint actuators for the allocation-free in-place
functions. The public types keep abstractly typed matrix fields (their numerics are part of the
existing behavior), which would force dynamic dispatch in a hot loop; this mirror has the same
field names and is built once per run by [`compile_compliant_model`](@ref).
"""
struct CompiledCompliantModel
    bodies::Vector{NamedTuple{(:mass_kg, :inertia_body_kg_m2), Tuple{Float64, SMatrix{3, 3, Float64, 9}}}}
    joints::Vector{CompiledJoint}
    actuators::Vector{CompiledActuator}
end

"""Build the concretely typed [`CompiledCompliantModel`](@ref) of a model and its joint actuators."""
function compile_compliant_model(model::CompliantMultibodyModel, actuators::AbstractVector{CompliantJointActuator}=CompliantJointActuator[])
    bodies = [(mass_kg=b.mass_kg, inertia_body_kg_m2=SMatrix{3, 3, Float64, 9}(b.inertia_body_kg_m2)) for b in model.bodies]
    joints = CompiledJoint[CompiledJoint(
        j.parent, j.child, j.parent_point_body, j.child_point_body,
        SMatrix{3, 3, Float64, 9}(j.k_translation_n_m), SMatrix{3, 3, Float64, 9}(j.c_translation_n_s_m),
        SMatrix{3, 3, Float64, 9}(j.k_rotation_n_m_rad), SMatrix{3, 3, Float64, 9}(j.c_rotation_n_m_s_rad)) for j in model.joints]
    acts = CompiledActuator[CompiledActuator(
        a.joint, a.torque_limit_n_m, SMatrix{3, 3, Float64, 9}(a.kp_n_m_rad), SMatrix{3, 3, Float64, 9}(a.kd_n_m_s_rad),
        a.feedforward_torque_child_body, a.efficiency) for a in actuators]
    return CompiledCompliantModel(bodies, joints, acts)
end

"""
    compliant_joint_loads_in_place!(forces, torques_body, model, col0, pos, quat, vel, angvel, mount, rest)
        -> (force_on_mount_world, torque_on_mount_world)

Allocation-free counterpart of [`compliant_joint_loads`](@ref) with the same joint math, on a
[`CompiledCompliantModel`](@ref) (which carries the joint actuators). Body `i` of
`model` is column `col0 + i` of the 3-by-N `forces` (world), `torques_body` (body frame), `pos`, `vel`,
`angvel` and the 4-by-N `quat`. The joint forces and torques are written to those columns (they are
zeroed first); `rest[j]` is the rest orientation of joint `j`. Joints with parent 0 attach to `mount`;
their reaction is returned as the total force on the mount frame (world) and the torque about the mount
origin (world), to be applied to whatever the mount is fixed to.
"""
function compliant_joint_loads_in_place!(
    forces::AbstractMatrix{Float64},
    torques_body::AbstractMatrix{Float64},
    model::CompiledCompliantModel,
    col0::Int,
    pos::AbstractMatrix{Float64},
    quat::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    angvel::AbstractMatrix{Float64},
    mount::CompliantMountKinematics,
    rest::AbstractVector{SVector{4, Float64}},
)
    @inbounds for i in eachindex(model.bodies), k in 1:3
        forces[k, col0 + i] = 0.0
        torques_body[k, col0 + i] = 0.0
    end
    f_mount = SVector{3, Float64}(0.0, 0.0, 0.0)
    t_mount = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for jidx in eachindex(model.joints)
        joint = model.joints[jidx]
        parent = joint.parent == 0 ?
            _mount_parent_kinematics(mount, joint.parent_point_body) :
            _column_kinematics(pos, quat, vel, angvel, col0 + joint.parent, joint.parent_point_body)
        cc = col0 + joint.child
        child = _column_kinematics(pos, quat, vel, angvel, cc, joint.child_point_body)
        k = _joint_load_kernel(joint, jidx, parent, child, rest[jidx], model.actuators)

        compliance_parent_world = -k.compliance_child_world
        actuator_parent_world = -k.actuator_child_world
        if joint.parent == 0
            f_mount += k.force_parent
            t_mount += cross(parent.point - mount.r, k.force_parent) + compliance_parent_world + actuator_parent_world
        else
            pc = col0 + joint.parent
            _add_column!(forces, pc, k.force_parent)
            _add_column!(torques_body, pc,
                cross(joint.parent_point_body, parent.R' * k.force_parent) +
                parent.R' * (compliance_parent_world + actuator_parent_world))
        end
        _add_column!(forces, cc, k.force_child)
        _add_column!(torques_body, cc,
            cross(joint.child_point_body, child.R' * k.force_child) +
            child.R' * (k.compliance_child_world + k.actuator_child_world))
    end
    return f_mount, t_mount
end

"""
    compliant_body_derivatives!(dpos, dquat, dvel, dangvel, model, col0, forces, torques_body, quat, vel, angvel)

Rigid-body derivative of every body of `model` (columns `col0 + 1:col0 + n`) from the total world-frame
`forces` and body-frame `torques_body` (joint loads plus anything the caller added), with the same
equations as [`compliant_multibody_dynamics`](@ref). Allocation-free.
"""
function compliant_body_derivatives!(
    dpos::AbstractMatrix{Float64}, dquat::AbstractMatrix{Float64}, dvel::AbstractMatrix{Float64}, dangvel::AbstractMatrix{Float64},
    model::CompiledCompliantModel, col0::Int,
    forces::AbstractMatrix{Float64}, torques_body::AbstractMatrix{Float64},
    quat::AbstractMatrix{Float64}, vel::AbstractMatrix{Float64}, angvel::AbstractMatrix{Float64},
)
    @inbounds for i in eachindex(model.bodies)
        c = col0 + i
        body = model.bodies[i]
        q = _unit_quat(SVector{4, Float64}(quat[1, c], quat[2, c], quat[3, c], quat[4, c]))
        ω = SVector{3, Float64}(angvel[1, c], angvel[2, c], angvel[3, c])
        qdot = 0.5 * _quat_raw_mul(q, SVector{4, Float64}(ω[1], ω[2], ω[3], 0.0))
        J = body.inertia_body_kg_m2
        α = J \ (SVector{3, Float64}(torques_body[1, c], torques_body[2, c], torques_body[3, c]) - cross(ω, J * ω))
        for k in 1:3
            dpos[k, c] = vel[k, c]
            dvel[k, c] = forces[k, c] / body.mass_kg
            dangvel[k, c] = α[k]
        end
        for k in 1:4
            dquat[k, c] = qdot[k]
        end
    end
    return nothing
end

"""
    compliant_state_in_mount_frame(model, x) -> Vector{Float64}

Express a state of `model` given in the frame where its base sits at `base_position`/`base_quaternion`
(as produced by `build_compliant_topology`) relative to that base, i.e. as mount-frame coordinates.
The identity when the base is at the origin with identity attitude.
"""
function compliant_state_in_mount_frame(model::CompliantMultibodyModel, x::AbstractVector)
    n = length(model.bodies)
    length(x) == 13n || throw(ArgumentError("state length must be 13 * number of bodies."))
    out = Vector{Float64}(x)
    qb = _unit_quat(model.base_quaternion)
    if model.base_position == SVector{3, Float64}(0.0, 0.0, 0.0) && qb == _Q_IDENTITY
        return out
    end
    Rt = _rot(qb)'
    qc = _quat_conj(qb)
    for i in 1:n
        b = 13 * (i - 1)
        s = compliant_state_parts(x, i)
        out[(b + 1):(b + 3)] .= Rt * (s.r - model.base_position)
        out[(b + 4):(b + 7)] .= _quat_mul(qc, s.q)
        out[(b + 8):(b + 10)] .= Rt * s.v
    end
    return out
end

"""
    compliant_rest_state(model) -> Vector{Float64}

Zero-velocity state at the model's rest geometry in the mount frame, built by walking a spanning tree
from the joints attached to the mount (`parent == 0`): each child takes its joint's rest orientation and
the position that makes the joint points coincide. Models with closed loops (meshes) have one geometry
per spanning tree; pass `initial_state` (for example a topology build's) for those.
"""
function compliant_rest_state(model::CompliantMultibodyModel)
    n = length(model.bodies)
    pos = Vector{SVector{3, Float64}}(undef, n)
    quat = Vector{SVector{4, Float64}}(undef, n)
    done = falses(n)
    progress = true
    while progress && !all(done)
        progress = false
        for joint in model.joints
            done[joint.child] && continue
            if joint.parent == 0
                q_parent = _Q_IDENTITY
                point = joint.parent_point_body
            elseif done[joint.parent]
                q_parent = quat[joint.parent]
                point = pos[joint.parent] + _rot(q_parent) * joint.parent_point_body
            else
                continue
            end
            q = _quat_mul(q_parent, joint.rest_child_parent_quat)
            quat[joint.child] = q
            pos[joint.child] = point - _rot(q) * joint.child_point_body
            done[joint.child] = true
            progress = true
        end
    end
    all(done) || throw(ArgumentError("compliant model bodies $(findall(!, done)) are not reachable from the mount frame through its joints; pass an explicit initial_state."))
    return compliant_state_vector(pos, quat)
end

"""
    compliant_state_to_inertial!(att_r, att_q, att_v, att_ω, col0, x_mount, mount)

Write the mount-frame state `x_mount` (13 numbers per body, relative to a fixed mount) into columns
`col0 + 1:col0 + n` of the inertial engine arrays for a mount that moves with `mount` kinematics: the
body gets the rigid-body velocity of its mount-frame position plus its relative velocity, and its
body-frame rate adds the mount's rotation.
"""
function compliant_state_to_inertial!(att_r, att_q, att_v, att_ω, col0::Int, x_mount::AbstractVector, mount::CompliantMountKinematics)
    n = length(x_mount) ÷ 13
    Rm = _rot(mount.q)
    ω_world = Rm * mount.ω
    for i in 1:n
        s = compliant_state_parts(x_mount, i)
        c = col0 + i
        rel = Rm * s.r
        q = _quat_mul(mount.q, s.q)
        r = mount.r + rel
        v = mount.v + cross(ω_world, rel) + Rm * s.v
        ω = s.ω + _rot(q)' * ω_world
        for k in 1:3
            att_r[k, c] = r[k]
            att_v[k, c] = v[k]
            att_ω[k, c] = ω[k]
        end
        for k in 1:4
            att_q[k, c] = q[k]
        end
    end
    return nothing
end

end # module ClothMultibody
