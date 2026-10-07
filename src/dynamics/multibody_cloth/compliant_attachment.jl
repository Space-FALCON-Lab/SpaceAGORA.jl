# Included inside `SpacecraftModels` (src/vehicle/spacecraft/model.jl) after `Link`, so that
# `SpacecraftModel` can hold attachments. The dynamics live in
# src/dynamics/multibody_cloth/compliant_attachment_dynamics.jl.

using ..ClothMultibody: CompliantMultibodyModel, CompliantTopologyBuild, CompliantJointActuator
using ..ClothMultibody: compliant_rest_state, compliant_state_in_mount_frame

"""
    CompliantAttachment(; model, link, mount_point=(0,0,0), mount_quaternion=(0,0,0,1),
                        initial_state=nothing, joint_actuators=[], rest_schedule=nothing)

A [`CompliantMultibodyModel`](@ref) mounted on a spacecraft `Link`. The attachment bodies (a cloth
panel mesh, a flexible appendage) move under their own compliant joint springs, dampers and actuators
and exchange forces and torques with the link they are mounted on, in both directions. `link` may be the
root, a link merged through a `:fixed` joint, or a link of a moving (articulated) body; the reaction
reaches the dynamics through that link's body.

- `model`: a `CompliantMultibodyModel`, or a `CompliantTopologyBuild` (its `initial_state` is then the
  default `initial_state`). At least one joint must have `parent == 0`: those joints attach to the mount
  frame. The model's `base_position` and `base_quaternion` are IGNORED: the mount frame replaces the
  fixed base. (A topology build's own state is converted from the base pose to the mount frame.)
- `link`: the spacecraft `Link` it is mounted on. It must be one of the spacecraft's links (or its root).
- `mount_point`: mount frame origin in the link frame (m). The link frame origin is the link center of
  mass (the spacecraft position for the root link) and its axes are the link body axes.
- `mount_quaternion`: mount frame orientation relative to the link frame, scalar-last, active.
- `initial_state`: 13 numbers per body (position, scalar-last attitude quaternion, velocity, body-frame
  angular velocity) in the MOUNT frame, relative to a mount that is at rest. The run adds the mount's
  inertial pose and rigid-body velocity. Default: the build's state, or zero velocities at the model's
  rest geometry (`compliant_rest_state`).
- `joint_actuators`: [`CompliantJointActuator`](@ref)s of the model's joints.
- `rest_schedule`: `nothing` (each joint keeps its own `rest_child_parent_quat`) or a callable
  `rest_schedule(out::Vector{SVector{4,Float64}}, t) -> nothing` that writes the rest quaternion of every
  joint for time `t` (seconds since the start) into `out`. It is called on every RHS evaluation, so keep
  it allocation-free for allocation-free runs.

Attachment bodies are NOT links of the spacecraft: they are not in its `links`, `dry_mass` or `mass`,
and carry their own masses (the system total is the spacecraft plus the attachments). In this version
the only loads on attachment bodies are gravity at each body's own position (position-only gravity
effectors) and the compliant joint loads: no aerodynamic, SRP, thermal or other effector acts on them.
"""
struct CompliantAttachment
    model::CompliantMultibodyModel
    link::Link
    mount_point::SVector{3, Float64}
    mount_quaternion::SVector{4, Float64}
    initial_state::Vector{Float64}
    joint_actuators::Vector{CompliantJointActuator}
    rest_schedule::Any
end

function CompliantAttachment(;
    model::Union{CompliantMultibodyModel, CompliantTopologyBuild},
    link::Link,
    mount_point=SVector{3, Float64}(0.0, 0.0, 0.0),
    mount_quaternion=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
    initial_state::Union{Nothing, AbstractVector{<:Real}}=nothing,
    joint_actuators::AbstractVector{CompliantJointActuator}=CompliantJointActuator[],
    rest_schedule=nothing,
)
    m = model isa CompliantTopologyBuild ? model.model : model
    nb = length(m.bodies)
    nb > 0 || throw(ArgumentError("CompliantAttachment model has no bodies."))
    any(j -> j.parent == 0, m.joints) ||
        throw(ArgumentError("CompliantAttachment model has no joint with parent == 0: nothing attaches it to the mount, so it would float free of the spacecraft."))
    mp = SVector{3, Float64}(mount_point)
    all(isfinite, mp) || throw(ArgumentError("CompliantAttachment mount_point must be finite."))
    mq = SVector{4, Float64}(mount_quaternion)
    nq = norm(mq)
    (isfinite(nq) && nq > 0.0) || throw(ArgumentError("CompliantAttachment mount_quaternion must be finite and nonzero."))
    mq = mq / nq
    for act in joint_actuators
        1 <= act.joint <= length(m.joints) ||
            throw(ArgumentError("CompliantAttachment actuator $(act.name) refers to joint $(act.joint), outside 1:$(length(m.joints))."))
    end
    x0 = if initial_state !== nothing
        Vector{Float64}(initial_state)
    elseif model isa CompliantTopologyBuild
        compliant_state_in_mount_frame(m, model.initial_state)
    else
        compliant_rest_state(m)
    end
    length(x0) == 13nb ||
        throw(ArgumentError("CompliantAttachment initial_state has $(length(x0)) entries; expected 13 per body = $(13nb)."))
    all(isfinite, x0) || throw(ArgumentError("CompliantAttachment initial_state must be finite."))
    rest_schedule === nothing || applicable(rest_schedule, Vector{SVector{4, Float64}}(undef, length(m.joints)), 0.0) ||
        throw(ArgumentError("CompliantAttachment rest_schedule must be callable as rest_schedule(out::Vector{SVector{4,Float64}}, t)."))
    return CompliantAttachment(m, link, mp, mq, x0, Vector{CompliantJointActuator}(joint_actuators), rest_schedule)
end

"""Number of bodies of a compliant attachment."""
attachment_body_count(a::CompliantAttachment)::Int = length(a.model.bodies)
