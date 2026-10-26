"""
Run-time support for `CompliantAttachment`s: per-spacecraft runtime, state layout, initialization and
the loads and derivatives of the attachment bodies.

Frames and conventions
- State per spacecraft with attachments: `att_r` (3-by-N), `att_q` (4-by-N), `att_v` (3-by-N),
  `att_ω` (3-by-N), N the body count over all attachments in attachment order. `att_r` and `att_v` are
  position and velocity RELATIVE to the spacecraft's translational state (`pos`/`vel`: the bus center of
  mass for a rigid spacecraft, the root composite center of mass for an articulated one), in INERTIAL
  axes. `att_q` is the absolute scalar-last active body-to-inertial quaternion and `att_ω` the absolute
  body-frame angular velocity. Relative positions keep the spring forces free of the roundoff of
  orbital-radius positions (k * ulp(7e6 m)).
- The mount frame is `mount_point`/`mount_quaternion` on the mounted link; its inertial pose and
  rigid-body velocity follow the link, which follows its dynamic body (the bus for a rigid spacecraft,
  the body of the link in an articulated one). Joints with `parent == 0` attach to the mount frame.
- Reaction on the mount: a force at the mount origin and a torque, applied to the link's body.
"""
module CompliantAttachmentDynamics

using LinearAlgebra
using StaticArrays
using ..ClothMultibody
using ..SpacecraftModels: SpacecraftModel, CompliantAttachment, Link
using ..ArticulatedBody: ArticulatedTree, ArticulatedWorkspace, articulated_kinematics

export AttachmentRuntime, build_attachment_runtime, attachment_state_shape, initialize_attachment_state!
export spacecraft_has_attachments, attachment_total_body_count, attachment_total_mass
export relative_gravity
export attachment_loads!, finish_attachments!, apply_attachments_rigid!, apply_attachments_articulated!

"""Evaluates a rest schedule at the runtime clock without boxing: `thunk()` writes into the joint rest buffer."""
struct ScheduleThunk{F}
    f::F
    buf::Vector{SVector{4, Float64}}
    time::Base.RefValue{Float64}
end

function (s::ScheduleThunk)()
    s.f(s.buf, s.time[])
    return nothing
end

"""
Per-spacecraft attachment runtime built once per run: the attachments, their column offsets in the
`att_*` state, rest-orientation buffers, rest-schedule thunks, the mounted link's frame in its dynamic
body, and force/torque scratch. One thread touches a satellite's runtime at a time.
"""
struct AttachmentRuntime
    n_att::Int
    n_bodies::Int
    attachments::Vector{CompliantAttachment}
    compiled::Vector{ClothMultibody.CompiledCompliantModel}
    col0::Vector{Int}
    rest::Vector{Vector{SVector{4, Float64}}}
    thunks::Vector{Any}
    time::Base.RefValue{Float64}
    link_body::Vector{Int}
    link_offset::Vector{SVector{3, Float64}}
    link_q::Vector{SVector{4, Float64}}
    forces::Matrix{Float64}
    torques::Matrix{Float64}
    base_gravity::Vector{SVector{3, Float64}}   # gravity at the base position for the current RHS call
end

spacecraft_has_attachments(sc::SpacecraftModel)::Bool = !isempty(sc.attachments)
attachment_total_body_count(sc::SpacecraftModel)::Int = sum((length(a.model.bodies) for a in sc.attachments); init=0)

"""Total mass of all attachment bodies (kg); not part of the spacecraft's `dry_mass` or `mass`."""
attachment_total_mass(sc::SpacecraftModel)::Float64 =
    sum((b.mass_kg for a in sc.attachments for b in a.model.bodies); init=0.0)

# Frame of `link` in the dynamic body that carries it: (body index, offset in body frame, attitude in body frame).
function _link_frame(sc::SpacecraftModel, tree::Union{Nothing, ArticulatedTree}, link::Link)
    if tree === nothing
        link === sc.root && return 1, SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        q = SVector{4, Float64}(link.q)
        return 1, SVector{3, Float64}(link.r), q / norm(q)
    end
    links = Link[l for l in sc.links]
    any(l -> l === sc.root, links) || pushfirst!(links, sc.root)
    idx = findfirst(l -> l === link, links)
    idx === nothing && throw(ArgumentError("Attachment link is not one of the spacecraft's links."))
    return tree.body_of_link[idx], tree.link_com_in_body[idx], tree.link_q_in_body[idx]
end

"""
    build_attachment_runtime(sc, tree=nothing) -> Union{Nothing, AttachmentRuntime}

`nothing` for a spacecraft without attachments. `tree` is the articulated tree of the spacecraft, or
`nothing` for a rigid one.
"""
function build_attachment_runtime(sc::SpacecraftModel, tree::Union{Nothing, ArticulatedTree}=nothing)
    isempty(sc.attachments) && return nothing
    n = length(sc.attachments)
    time = Ref(0.0)
    col0 = zeros(Int, n)
    rest = Vector{Vector{SVector{4, Float64}}}(undef, n)
    thunks = Vector{Any}(undef, n)
    lb = zeros(Int, n)
    lo = Vector{SVector{3, Float64}}(undef, n)
    lq = Vector{SVector{4, Float64}}(undef, n)
    total = 0
    for (k, a) in pairs(sc.attachments)
        col0[k] = total
        total += length(a.model.bodies)
        rest[k] = SVector{4, Float64}[j.rest_child_parent_quat for j in a.model.joints]
        thunks[k] = a.rest_schedule === nothing ? nothing : ScheduleThunk(a.rest_schedule, rest[k], time)
        lb[k], lo[k], lq[k] = _link_frame(sc, tree, a.link)
    end
    compiled = [ClothMultibody.compile_compliant_model(a.model, a.joint_actuators) for a in sc.attachments]
    return AttachmentRuntime(n, total, Vector{CompliantAttachment}(sc.attachments), compiled, col0, rest, thunks, time, lb, lo, lq,
        zeros(3, total), zeros(3, total), [SVector{3, Float64}(0.0, 0.0, 0.0)])
end

"""State arrays of a spacecraft with attachments (for the ComponentVector layout)."""
function attachment_state_shape(sc::SpacecraftModel)
    n = attachment_total_body_count(sc)
    return (att_r=zeros(3, n), att_q=zeros(4, n), att_v=zeros(3, n), att_ω=zeros(3, n))
end

"""
    initialize_attachment_state!(sc_view, sc, tree=nothing, rt=build_attachment_runtime(sc, tree))

Write the initial attachment state into `sc_view` from each attachment's mount-frame `initial_state`,
the mount's inertial pose and its rigid-body velocity. `sc_view` already holds the spacecraft's
initial `pos`, `vel`, `q`, `ω` (and `joint_q`/`joint_qd` for an articulated spacecraft).
"""
function initialize_attachment_state!(sc_view, sc::SpacecraftModel, tree::Union{Nothing, ArticulatedTree}=nothing,
        rt::Union{Nothing, AttachmentRuntime}=build_attachment_runtime(sc, tree))
    rt === nothing && return nothing
    q = SVector{4, Float64}(sc_view.q); ω = SVector{3, Float64}(sc_view.ω)
    # Relative state: the mount kinematics are taken about the base, so position and velocity are zero there.
    z = SVector{3, Float64}(0.0, 0.0, 0.0)
    kin = nothing
    if tree !== nothing
        base = (pos=z, vel=z, q=q, ω=ω)
        kin = articulated_kinematics(tree, base, Float64.(sc_view.joint_q), Float64.(sc_view.joint_qd))
    end
    for k in 1:rt.n_att
        a = rt.attachments[k]
        b = rt.link_body[k]
        mount = if kin === nothing
            ClothMultibody.compliant_mount_kinematics(z, z, q, ω, rt.link_offset[k], rt.link_q[k], a.mount_point, a.mount_quaternion)
        else
            ClothMultibody.compliant_mount_kinematics(kin.pos[b], kin.vel[b], kin.quat[b], kin.ω[b], rt.link_offset[k], rt.link_q[k], a.mount_point, a.mount_quaternion)
        end
        ClothMultibody.compliant_state_to_inertial!(sc_view.att_r, sc_view.att_q, sc_view.att_v, sc_view.att_ω, rt.col0[k], a.initial_state, mount)
    end
    return nothing
end

"""
    relative_gravity(gravity, r_base, ρ, g_base) -> g(r_base + ρ) - g(r_base)

Gravity difference between a body at `r_base + ρ` and the base, `g_base = gravity(r_base)`. The fallback
subtracts two evaluations. The engine adds a method for its gravity callable that differences the
point-mass part analytically (Encke), so the difference does not cancel two 8 m/s^2 terms. This module is
included before the gravity code, hence the hook.
"""
@inline relative_gravity(gravity, r_base, ρ, g_base) = gravity(r_base + ρ) - g_base

"""
    attachment_loads!(du_view, sc_view, rt, k, mount, gravity, base_pos, base_gravity) -> (force_on_mount, torque_on_mount)

Evaluate attachment `k` in coordinates relative to the base (`mount` is relative too): its rest schedule
(at the runtime clock `rt.time`), the compliant joint loads, and the body derivatives written to
`du_view.att_*`. Returns the reaction on the mount (force at the mount origin and torque about it,
inertial). `att_v'` is written WITHOUT the base acceleration: [`finish_attachments!`](@ref) subtracts it
once the base acceleration, which includes these reactions, is known.

Only the gravity DIFFERENCE `g(r_base + r_rel) - g(r_base)` enters here (about 1e-6 m/s^2 for a 1 m
offset in LEO, against 8 m/s^2 for each term). [`relative_gravity`](@ref) forms it, by Encke
differencing for the engine's gravity models, so the large common part never reaches the relative
acceleration. The remaining `g(r_base) - a_base` is the negative of the
non-gravitational base acceleration and is formed in `finish_attachments!`.
"""
function attachment_loads!(du_view, sc_view, rt::AttachmentRuntime, k::Int, mount::ClothMultibody.CompliantMountKinematics, gravity,
        base_pos::SVector{3, Float64}, base_gravity::SVector{3, Float64})
    model = rt.compiled[k]
    col0 = rt.col0[k]
    thunk = rt.thunks[k]
    thunk === nothing || thunk()
    att_r = sc_view.att_r; att_q = sc_view.att_q; att_v = sc_view.att_v; att_ω = sc_view.att_ω
    f_mount, t_mount = ClothMultibody.compliant_joint_loads_in_place!(
        rt.forces, rt.torques, model, col0, att_r, att_q, att_v, att_ω, mount, rt.rest[k])
    @inbounds for i in eachindex(model.bodies)
        c = col0 + i
        g = relative_gravity(gravity, base_pos, SVector{3, Float64}(att_r[1, c], att_r[2, c], att_r[3, c]), base_gravity)
        m = model.bodies[i].mass_kg
        rt.forces[1, c] += m * g[1]
        rt.forces[2, c] += m * g[2]
        rt.forces[3, c] += m * g[3]
    end
    ClothMultibody.compliant_body_derivatives!(
        du_view.att_r, du_view.att_q, du_view.att_v, du_view.att_ω, model, col0,
        rt.forces, rt.torques, att_q, att_v, att_ω)
    return f_mount, t_mount
end

"""
    finish_attachments!(du_view, rt, base_acceleration)

Second stage of the attachment derivatives, once the base acceleration (inertial, including the
attachment reactions and gravity at the base) is known: `att_v' += g(r_base) - a_base`, the part of the
relative acceleration that is common to all bodies. Allocation-free.
"""
function finish_attachments!(du_view, rt::AttachmentRuntime, base_acceleration)
    d = rt.base_gravity[1] - SVector{3, Float64}(base_acceleration[1], base_acceleration[2], base_acceleration[3])
    dv = du_view.att_v
    @inbounds for c in 1:rt.n_bodies
        dv[1, c] += d[1]
        dv[2, c] += d[2]
        dv[3, c] += d[3]
    end
    return nothing
end

"""
    apply_attachments_rigid!(du_view, sc_view, rt, t, forces, torques, gravity)

Rigid spacecraft, stage 1: evaluate every attachment (relative to the bus) and add the reactions to the
bus `forces` (inertial) and `torques` (body frame, about the bus position) before the bus equations of
motion. Call [`finish_attachments!`](@ref) with the bus acceleration afterwards.
"""
function apply_attachments_rigid!(du_view, sc_view, rt::AttachmentRuntime, t::Float64, forces, torques, gravity)
    rt.time[] = t
    z = SVector{3, Float64}(0.0, 0.0, 0.0)
    pos = SVector{3, Float64}(sc_view.pos[1], sc_view.pos[2], sc_view.pos[3])
    q = SVector{4, Float64}(sc_view.q[1], sc_view.q[2], sc_view.q[3], sc_view.q[4])
    ω = SVector{3, Float64}(sc_view.ω[1], sc_view.ω[2], sc_view.ω[3])
    g_base = SVector{3, Float64}(gravity(pos))
    rt.base_gravity[1] = g_base
    for k in 1:rt.n_att
        a = rt.attachments[k]
        mount = ClothMultibody.compliant_mount_kinematics(z, z, q, ω, rt.link_offset[k], rt.link_q[k], a.mount_point, a.mount_quaternion)
        f, τ = attachment_loads!(du_view, sc_view, rt, k, mount, gravity, pos, g_base)
        F, τ_ref = ClothMultibody.compliant_mount_wrench(mount, z, f, τ)
        τb = ClothMultibody.compliant_world_to_body(q, τ_ref)
        forces[1] += F[1]; forces[2] += F[2]; forces[3] += F[3]
        torques[1] += τb[1]; torques[2] += τb[2]; torques[3] += τb[3]
    end
    return nothing
end

"""
    apply_attachments_articulated!(du_view, sc_view, rt, ws, t, gravity)

Articulated spacecraft, stage 1: with the body kinematics of `ws` computed RELATIVE to the root COM
(`articulated_kinematics!` with zero base position and velocity), evaluate every attachment and write the
reaction on each mounted body into `ws.ext_force` and `ws.ext_torque` (inertial, about the body COM),
which `articulated_dynamics!` takes as `body_force_world`/`body_torque_world`. Call
[`finish_attachments!`](@ref) with the returned base acceleration afterwards.
"""
function apply_attachments_articulated!(du_view, sc_view, rt::AttachmentRuntime, ws::ArticulatedWorkspace{Float64}, t::Float64, gravity)
    rt.time[] = t
    z = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for b in eachindex(ws.ext_force)
        ws.ext_force[b] = z
        ws.ext_torque[b] = z
    end
    pos = SVector{3, Float64}(sc_view.pos[1], sc_view.pos[2], sc_view.pos[3])
    g_base = SVector{3, Float64}(gravity(pos))
    rt.base_gravity[1] = g_base
    for k in 1:rt.n_att
        a = rt.attachments[k]
        b = rt.link_body[k]
        mount = ClothMultibody.compliant_mount_kinematics(
            ws.kpos[b], ws.kvel[b], ws.kquat[b], ws.kwb[b], rt.link_offset[k], rt.link_q[k], a.mount_point, a.mount_quaternion)
        f, τ = attachment_loads!(du_view, sc_view, rt, k, mount, gravity, pos, g_base)
        F, τ_ref = ClothMultibody.compliant_mount_wrench(mount, ws.kpos[b], f, τ)
        ws.ext_force[b] += F
        ws.ext_torque[b] += τ_ref
    end
    return nothing
end

end # module CompliantAttachmentDynamics
