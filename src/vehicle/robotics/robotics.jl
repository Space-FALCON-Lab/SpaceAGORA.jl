"""Robot and cloth-arm kinematics module for serial-arm models, FK, IK, and surface targeting."""
module Robotics

using LinearAlgebra
using StaticArrays

export ClothArmBasePose, ClothArmLink, ClothArmJoint, ClothArmModel, ClothArmPose, ClothArmState
export default_cloth_arm_model, cloth_fk, cloth_fk_state, cloth_end_effector_pose
export cloth_ik, cloth_total_reach, closest_surface_target

const _Q_IDENTITY = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

"""Base pose for a cloth-arm kinematic chain in world coordinates."""
struct ClothArmBasePose
    position::SVector{3, Float64}
    quaternion::SVector{4, Float64}
end

"""Base pose for a cloth-arm kinematic chain in world coordinates."""
ClothArmBasePose(position::Union{AbstractVector{<:Real}, NTuple{3, <:Real}}) =
    ClothArmBasePose(SVector{3, Float64}(position), _Q_IDENTITY)

ClothArmBasePose(
    position::Union{AbstractVector{<:Real}, NTuple{3, <:Real}},
    quaternion::Union{AbstractVector{<:Real}, NTuple{4, <:Real}},
) =
    ClothArmBasePose(SVector{3, Float64}(position), _unit_quat(quaternion))

"""Rigid link geometry and mass properties for a cloth-arm model."""
struct ClothArmLink
    name::Symbol
    vector_parent::SVector{3, Float64}
    com_offset_parent::SVector{3, Float64}
    radius_m::Float64
    mass_kg::Float64
end

"""Single revolute joint axis and position limits for a cloth-arm model."""
struct ClothArmJoint
    name::Symbol
    axis_parent::SVector{3, Float64}
    lower_rad::Float64
    upper_rad::Float64
end

"""Serial cloth-arm model composed of links and joints."""
struct ClothArmModel
    links::Vector{ClothArmLink}
    joints::Vector{ClothArmJoint}
    mount_offset_body::SVector{3, Float64}
end

"""Forward-kinematics output for cloth-arm joint and link frames."""
struct ClothArmPose
    base_position::SVector{3, Float64}
    base_quaternion::SVector{4, Float64}
    joint_origins::Vector{SVector{3, Float64}}
    joint_axes_world::Vector{SVector{3, Float64}}
    link_com_positions::Vector{SVector{3, Float64}}
    link_tip_positions::Vector{SVector{3, Float64}}
    link_quaternions::Vector{SVector{4, Float64}}
    end_effector_position::SVector{3, Float64}
    end_effector_quaternion::SVector{4, Float64}
end

"""Cloth-arm pose plus joint velocity and acceleration state."""
struct ClothArmState
    pose::ClothArmPose
    q::Vector{Float64}
    dq::Vector{Float64}
    link_linear_velocities::Vector{SVector{3, Float64}}
    link_angular_velocities::Vector{SVector{3, Float64}}
end

"""Normalize a quaternion and return it as a static vector."""
@inline function _unit_quat(q)::SVector{4, Float64}
    qv = SVector{4, Float64}(q)
    nq = norm(qv)
    return isfinite(nq) && nq > eps(Float64) ? qv / nq : _Q_IDENTITY
end

"""Multiply two quaternions using the project convention."""
@inline function _quat_mul(a::SVector{4, Float64}, b::SVector{4, Float64})::SVector{4, Float64}
    ax, ay, az, aw = a
    bx, by, bz, bw = b
    return SVector{4, Float64}(
        aw * bx + bw * ax + ay * bz - az * by,
        aw * by + bw * ay + az * bx - ax * bz,
        aw * bz + bw * az + ax * by - ay * bx,
        aw * bw - ax * bx - ay * by - az * bz,
    )
end

"""Create a unit quaternion from an axis-angle rotation."""
@inline function _quat_from_axis_angle(axis::SVector{3, Float64}, θ::Float64)::SVector{4, Float64}
    na = norm(axis)
    if !(isfinite(na) && na > eps(Float64))
        return _Q_IDENTITY
    end
    u = axis / na
    s = sin(0.5 * θ)
    return SVector{4, Float64}(s * u[1], s * u[2], s * u[3], cos(0.5 * θ))
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

"""Construct the default serial cloth-arm model used by examples and tests."""
function default_cloth_arm_model(;
    link_lengths_m=(0.18, 0.16, 0.12),
    link_radii_m=(0.018, 0.016, 0.014),
    link_masses_kg=(0.15, 0.10, 0.05),
    joint_axes=((0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0)),
    joint_limit_rad=deg2rad(175.0),
    mount_offset_body=(0.0, 0.0, 0.0),
)
    n = length(link_lengths_m)
    length(link_radii_m) == n || throw(ArgumentError("link_radii_m must match link_lengths_m."))
    length(link_masses_kg) == n || throw(ArgumentError("link_masses_kg must match link_lengths_m."))
    length(joint_axes) == n || throw(ArgumentError("joint_axes must match link_lengths_m."))
    links = ClothArmLink[
        ClothArmLink(
            Symbol("link_$i"),
            SVector{3, Float64}(Float64(link_lengths_m[i]), 0.0, 0.0),
            SVector{3, Float64}(0.5 * Float64(link_lengths_m[i]), 0.0, 0.0),
            Float64(link_radii_m[i]),
            Float64(link_masses_kg[i]),
        ) for i in 1:n
    ]
    joints = ClothArmJoint[
        ClothArmJoint(
            Symbol("joint_$i"),
            _normalize_axis(SVector{3, Float64}(joint_axes[i])),
            -Float64(joint_limit_rad),
            Float64(joint_limit_rad),
        ) for i in 1:n
    ]
    return ClothArmModel(links, joints, SVector{3, Float64}(mount_offset_body))
end

"""Normalize a joint axis and reject zero-length inputs."""
@inline function _normalize_axis(axis::SVector{3, Float64})::SVector{3, Float64}
    na = norm(axis)
    na > eps(Float64) || throw(ArgumentError("joint axis must be nonzero."))
    return axis / na
end

"""Return the total reach of all links in a cloth-arm model."""
function cloth_total_reach(model::ClothArmModel)::Float64
    return sum(norm(link.vector_parent) for link in model.links)
end

"""Validate that a joint vector length matches the arm model."""
function _validate_joint_vector(model::ClothArmModel, q)
    length(q) == length(model.joints) ||
        throw(ArgumentError("expected $(length(model.joints)) joint coordinates, got $(length(q))."))
    return Float64.(collect(q))
end

"""Evaluate cloth-arm forward kinematics for a base pose and joint vector."""
function cloth_fk(model::ClothArmModel, base_pose::ClothArmBasePose, q)::ClothArmPose
    qv = _validate_joint_vector(model, q)
    R_base = _rot(base_pose.quaternion)
    r_curr = base_pose.position + R_base * model.mount_offset_body
    Q_curr = base_pose.quaternion

    joint_origins = SVector{3, Float64}[]
    joint_axes_world = SVector{3, Float64}[]
    link_com_positions = SVector{3, Float64}[]
    link_tip_positions = SVector{3, Float64}[]
    link_quaternions = SVector{4, Float64}[]

    for i in eachindex(model.links)
        joint = model.joints[i]
        link = model.links[i]
        R_parent = _rot(Q_curr)
        axis_world = R_parent * joint.axis_parent
        Q_joint = _quat_from_axis_angle(joint.axis_parent, qv[i])
        Q_link = _unit_quat(_quat_mul(Q_curr, Q_joint))
        R_link = _rot(Q_link)

        push!(joint_origins, r_curr)
        push!(joint_axes_world, axis_world)
        push!(link_com_positions, r_curr + R_link * link.com_offset_parent)
        r_curr = r_curr + R_link * link.vector_parent
        push!(link_tip_positions, r_curr)
        push!(link_quaternions, Q_link)
        Q_curr = Q_link
    end

    return ClothArmPose(
        base_pose.position,
        base_pose.quaternion,
        joint_origins,
        joint_axes_world,
        link_com_positions,
        link_tip_positions,
        link_quaternions,
        r_curr,
        Q_curr,
    )
end

"""Evaluate cloth-arm forward kinematics and attach joint velocity and acceleration."""
function cloth_fk_state(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q;
    dq=zeros(length(model.joints)),
)::ClothArmState
    qv = _validate_joint_vector(model, q)
    dqv = _validate_joint_vector(model, dq)
    pose = cloth_fk(model, base_pose, qv)
    return ClothArmState(
        pose,
        qv,
        dqv,
        fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(model.links)),
        fill(SVector{3, Float64}(0.0, 0.0, 0.0), length(model.links)),
    )
end

"""Return the end-effector position and orientation from a cloth-arm pose or state."""
cloth_end_effector_pose(state::ClothArmState) =
    state.pose.end_effector_position, state.pose.end_effector_quaternion

"""Return the end-effector position and orientation from a cloth-arm pose or state."""
cloth_end_effector_pose(pose::ClothArmPose) =
    pose.end_effector_position, pose.end_effector_quaternion

"""Compute the geometric end-effector position Jacobian for the cloth arm."""
function _ee_position_jacobian(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q::Vector{Float64},
    ee::SVector{3, Float64};
    h::Float64=1.0e-6,
)
    n = length(q)
    J = zeros(3, n)
    q_trial = copy(q)
    for i in 1:n
        step = max(h, h * abs(q[i]))
        q_trial[i] = q[i] + step
        hi = cloth_fk(model, base_pose, q_trial).end_effector_position
        q_trial[i] = q[i] - step
        lo = cloth_fk(model, base_pose, q_trial).end_effector_position
        q_trial[i] = q[i]
        J[:, i] .= (hi - lo) / (2step)
    end
    return J
end

"""Solve cloth-arm inverse kinematics with damped least-squares iterations."""
function cloth_ik(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    target;
    q_seed=zeros(length(model.joints)),
    max_iters::Int=100,
    position_tol_m::Float64=1.0e-4,
    damping::Float64=1.0e-3,
    step_limit_rad::Float64=0.25,
)
    q = _validate_joint_vector(model, q_seed)
    target_v = SVector{3, Float64}(target)
    for _ in 1:max_iters
        pose = cloth_fk(model, base_pose, q)
        err = target_v - pose.end_effector_position
        norm(err) <= position_tol_m && return q
        J = _ee_position_jacobian(model, base_pose, q, pose.end_effector_position)
        Δq = (J' * J + damping^2 * I) \ (J' * Vector(err))
        Δnorm = norm(Δq)
        if Δnorm > step_limit_rad
            Δq .*= step_limit_rad / Δnorm
        end
        q .+= Δq
        for i in eachindex(q)
            q[i] = clamp(q[i], model.joints[i].lower_rad, model.joints[i].upper_rad)
        end
    end
    pose = cloth_fk(model, base_pose, q)
    err = norm(target_v - pose.end_effector_position)
    err <= 10position_tol_m || throw(ErrorException("cloth_ik failed to reach target; final error=$(err) m."))
    return q
end

"""Pick the surface point nearest a reference and offset it by an optional standoff."""
function closest_surface_target(points, reference_position; standoff_m::Real=0.0)
    pts = Matrix{Float64}(points)
    size(pts, 1) == 3 || throw(ArgumentError("surface points must be a 3 x N matrix."))
    size(pts, 2) > 0 || throw(ArgumentError("surface points must contain at least one point."))
    ref = SVector{3, Float64}(reference_position)
    best_idx = 1
    best_d2 = Inf
    for i in axes(pts, 2)
        p = SVector{3, Float64}(pts[:, i])
        d2 = sum(abs2, p - ref)
        if d2 < best_d2
            best_d2 = d2
            best_idx = i
        end
    end
    surface_point = SVector{3, Float64}(pts[:, best_idx])
    normal = ref - surface_point
    nrm = norm(normal)
    normal = nrm > eps(Float64) ? normal / nrm : SVector{3, Float64}(1.0, 0.0, 0.0)
    return (
        target=surface_point + Float64(standoff_m) * normal,
        surface_point=surface_point,
        surface_normal=normal,
        index=best_idx,
    )
end

end # module Robotics
