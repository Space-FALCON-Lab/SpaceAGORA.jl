"""
Articulated-body dynamics for a free-floating root body with joint-space children.

Stage A (standalone library, no engine integration): `build_articulated_tree` turns a
`SpacecraftModel` whose `Joint`s carry `joint_type` into an immutable `ArticulatedTree`;
`articulated_dynamics!` returns the root accelerations and the joint accelerations.

Conventions
- Quaternions are scalar-last `[x, y, z, w]`. A state quaternion `q` is the ACTIVE body-to-inertial
  rotation: `v_inertial = R_a(q) * v_body`, composition `q_child_inertial = q_parent ⊗ q_rel`, and the
  engine's `rot(q)` (`QuaternionMath.rot`) is the transpose `R_a(q)'` (inertial-to-body). The
  standalone cloth model's `_rot(q)` is `R_a(q)`.
- Root state: `pos`/`vel` are the inertial position and velocity of the ROOT BODY COM, `q` the root
  attitude, `ω` the root body-frame angular velocity. The root body is the root link plus any links
  attached through `:fixed` joints plus the propellant (a point mass AT the composite COM of the
  dry group: it adds mass and no inertia and never moves the COM), so its COM is the composite
  COM of the dry group. The system COM differs from the root COM whenever
  moving bodies are present.
- Joint coordinates: `:hinge`/`:slide` one scalar each (rad, m) about/along `axis` (parent frame),
  `:ball` a scalar-last quaternion (4 entries in `joint_q`, 3 angular-velocity entries in
  `joint_qd`). Joint rates and accelerations of a ball joint are the relative angular velocity
  and acceleration expressed in the PARENT body frame. Coordinate 0 is the configured geometry.
- Algorithm: the equations of motion are assembled from per-body Jacobians (`M = Σ Jᵀ 𝕄 J`,
  bias from zero-acceleration kinematics) and solved by a hand-written Cholesky factorization.
  For the target size (<= ~20 dof) this is O(n^3) with tiny n, much simpler than Featherstone's ABA,
  is allocation-free with a preallocated workspace, and is generic over the element type.
"""
module ArticulatedBody

using LinearAlgebra
using StaticArrays
using ..SpacecraftModels: SpacecraftModel, Link, Joint

export ArticulatedTree, ArticulatedWorkspace, ArticulatedBaseState
export build_articulated_tree, articulated_dynamics!, articulated_forward_kinematics
export articulated_kinematics, articulated_kinematics!, articulated_link_poses, articulated_joint_qdot, articulated_joint_qdot!
export articulated_potential_energy, articulated_moving_mass, articulated_has_moving_joints
export ArticulatedRuntime

const _Q_ID = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
const _I3 = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

# ---------------------------------------------------------------------------
# Generic quaternion / rotation helpers (no Float64 annotations: duals and BigFloat flow through)
# ---------------------------------------------------------------------------

@inline function _qmul(a, b)
    return SVector(
        a[4] * b[1] + b[4] * a[1] + a[2] * b[3] - a[3] * b[2],
        a[4] * b[2] + b[4] * a[2] + a[3] * b[1] - a[1] * b[3],
        a[4] * b[3] + b[4] * a[3] + a[1] * b[2] - a[2] * b[1],
        a[4] * b[4] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
    )
end

@inline _qconj(q) = SVector(-q[1], -q[2], -q[3], q[4])

@inline function _qnormalize(q)
    n = sqrt(q[1]^2 + q[2]^2 + q[3]^2 + q[4]^2)
    return q / n
end

# Active (body-to-inertial) rotation matrix of a (possibly unnormalized) scalar-last quaternion.
@inline function _rotmat(q)
    x, y, z, w = q[1], q[2], q[3], q[4]
    s = 2 / (x * x + y * y + z * z + w * w)
    xx = s * x * x; yy = s * y * y; zz = s * z * z
    xy = s * x * y; xz = s * x * z; yz = s * y * z
    wx = s * w * x; wy = s * w * y; wz = s * w * z
    return SMatrix{3, 3}(
        1 - yy - zz, xy + wz, xz - wy,
        xy - wz, 1 - xx - zz, yz + wx,
        xz + wy, yz - wx, 1 - xx - yy,
    )
end

# Rodrigues rotation about a unit axis.
@inline function _axis_rotmat(a, θ)
    c = cos(θ); s = sin(θ); t = 1 - c
    return SMatrix{3, 3}(
        c + t * a[1] * a[1], t * a[1] * a[2] + s * a[3], t * a[1] * a[3] - s * a[2],
        t * a[1] * a[2] - s * a[3], c + t * a[2] * a[2], t * a[2] * a[3] + s * a[1],
        t * a[1] * a[3] + s * a[2], t * a[2] * a[3] - s * a[1], c + t * a[3] * a[3],
    )
end

@inline function _axis_quat(a, θ)
    h = θ / 2
    s = sin(h)
    return SVector(a[1] * s, a[2] * s, a[3] * s, cos(h))
end

# Axis-angle (rotation vector) of a unit quaternion, shortest rotation.
@inline function _axis_angle(q)
    sgn = q[4] < 0 ? -one(q[4]) : one(q[4])
    x, y, z, w = sgn * q[1], sgn * q[2], sgn * q[3], sgn * q[4]
    nv2 = x * x + y * y + z * z
    if nv2 < 1.0e-18
        f = 2 / w
    else
        nv = sqrt(nv2)
        f = 2 * atan(nv, w) / nv
    end
    return SVector(f * x, f * y, f * z)
end

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    ArticulatedBaseState(pos, vel, q, ω)

Root-body state: inertial COM position and velocity, scalar-last attitude quaternion (active
body-to-inertial) and body-frame angular velocity.
"""
struct ArticulatedBaseState{T <: Real}
    pos::SVector{3, T}
    vel::SVector{3, T}
    q::SVector{4, T}
    ω::SVector{3, T}
end

function ArticulatedBaseState(pos, vel, q, ω)
    T = promote_type(eltype(pos), eltype(vel), eltype(q), eltype(ω))
    return ArticulatedBaseState{T}(SVector{3, T}(pos), SVector{3, T}(vel), SVector{4, T}(q), SVector{3, T}(ω))
end

"""
    ArticulatedTree

Immutable description of the dynamic tree built by [`build_articulated_tree`](@ref). Dynamic bodies
are numbered `1:nb`, body 1 is the root, and every parent index is smaller than its child's. Fields
indexed per body (`parent`, `jtype`, `axis`, `d1`, `d2`, `C0`, `qC0`, `mass`, `inertia`, `k`, `c`,
`rest`, `kmat`, `cmat`, `restq`, `qoff`, `voff`) hold the root at index 1 with placeholder values.

- `parent[b]`: parent body (0 for the root); `jtype[b]`: `:root`, `:hinge`, `:slide` or `:ball`.
- `axis[b]`: hinge/slide axis, unit, in the parent body frame.
- `d1[b]`: parent COM to joint point, parent frame, at coordinate 0. `d2[b]`: child COM to joint
  point, child frame. `C0[b]`/`qC0[b]`: child-to-parent rotation at coordinate 0.
- `mass`, `inertia`: composite mass and inertia about the composite COM in the body frame (the frame
  of the body's head link).
- `qoff[b]`/`voff[b]`: 1-based start of the body's entries in `joint_q`/`joint_qd` (0 for the root).
- `nq`, `nv`: joint coordinate and joint rate counts.
- Link bookkeeping for output: `body_of_link`, `link_com_in_body`, `link_q_in_body`, `root_com_bus`.
- `joint_q`/`joint_qd` list the NON-fixed joints in the order of `spacecraft.joints`;
  `joint_of_body[b]` is the index into `spacecraft.joints` of the joint of body `b` (0 for the
  root). `q0`/`qd0` are the initial joint coordinates and rates from `Joint(initial_q, initial_qd)`.
- `mass[1]` includes the propellant mass given at build time as a point mass at the root COM;
  the engine overrides the root mass at runtime (`root_mass` keyword of
  [`articulated_dynamics!`](@ref)).
"""
struct ArticulatedTree
    nb::Int
    nq::Int
    nv::Int
    parent::Vector{Int}
    jtype::Vector{Symbol}
    axis::Vector{SVector{3, Float64}}
    d1::Vector{SVector{3, Float64}}
    d2::Vector{SVector{3, Float64}}
    C0::Vector{SMatrix{3, 3, Float64, 9}}
    qC0::Vector{SVector{4, Float64}}
    mass::Vector{Float64}
    inertia::Vector{SMatrix{3, 3, Float64, 9}}
    k::Vector{Float64}
    c::Vector{Float64}
    rest::Vector{Float64}
    kmat::Vector{SMatrix{3, 3, Float64, 9}}
    cmat::Vector{SMatrix{3, 3, Float64, 9}}
    restq::Vector{SVector{4, Float64}}
    qoff::Vector{Int}
    voff::Vector{Int}
    body_of_link::Vector{Int}
    link_com_in_body::Vector{SVector{3, Float64}}
    link_q_in_body::Vector{SVector{4, Float64}}
    root_com_bus::SVector{3, Float64}
    q0::Vector{Float64}
    qd0::Vector{Float64}
    joint_of_body::Vector{Int}
end

"""
    ArticulatedWorkspace(tree, T=Float64)

Preallocated scratch space for [`articulated_dynamics!`](@ref), built once per tree. Use
`T = Float64` on the hot path; build one with another element type (BigFloat, ForwardDiff duals) to
differentiate or extend precision (that path may allocate).
"""
struct ArticulatedWorkspace{T}
    R::Vector{SMatrix{3, 3, T, 9}}
    x::Vector{SVector{3, T}}
    v::Vector{SVector{3, T}}
    w::Vector{SVector{3, T}}
    ab::Vector{SVector{3, T}}
    alb::Vector{SVector{3, T}}
    Jv::Array{T, 3}
    Jw::Array{T, 3}
    A::Matrix{T}
    M::Matrix{T}
    rhs::Vector{T}
    u::Vector{T}
    qdd::Vector{T}
    # Kinematics output buffers (see `articulated_kinematics!`): per dynamic body inertial COM
    # position, attitude quaternion, COM velocity, body-frame and world-frame angular velocity.
    kpos::Vector{SVector{3, T}}
    kquat::Vector{SVector{4, T}}
    kvel::Vector{SVector{3, T}}
    kwb::Vector{SVector{3, T}}
    kww::Vector{SVector{3, T}}
    # External wrench per dynamic body about its COM, world frame (see `articulated_dynamics!`);
    # written by the caller (the engine) and zeroed by it.
    ext_force::Vector{SVector{3, T}}
    ext_torque::Vector{SVector{3, T}}
end

function ArticulatedWorkspace(tree::ArticulatedTree, ::Type{T}=Float64) where {T}
    nb = tree.nb
    n = 6 + tree.nv
    z3 = zero(SVector{3, T})
    return ArticulatedWorkspace{T}(
        fill(zero(SMatrix{3, 3, T, 9}), nb), fill(z3, nb), fill(z3, nb), fill(z3, nb), fill(z3, nb), fill(z3, nb),
        zeros(T, 3, n, nb), zeros(T, 3, n, nb), zeros(T, 3, n), zeros(T, n, n), zeros(T, n), zeros(T, n),
        zeros(T, tree.nv),
        fill(z3, nb), fill(zero(SVector{4, T}), nb), fill(z3, nb), fill(z3, nb), fill(z3, nb),
        fill(z3, nb), fill(z3, nb),
    )
end

# ---------------------------------------------------------------------------
# Tree construction
# ---------------------------------------------------------------------------

function _link_config(l::Link, is_root::Bool)
    # Root link r/q are inertial-frame values, not a bus-frame offset: the bus frame is the root frame.
    is_root && return (zero(SVector{3, Float64}), _Q_ID)
    q = SVector{4, Float64}(l.q)
    nq = norm(q)
    (isfinite(nq) && nq > 0) || throw(ArgumentError("Link quaternion must be finite and nonzero."))
    return (SVector{3, Float64}(l.r), q / nq)
end

"""
    build_articulated_tree(sc::SpacecraftModel; prop_mass=sc.prop_mass) -> ArticulatedTree

Build the dynamic tree from `sc.joints`: `link1` is the parent and `link2` the child of every joint.
Validates a proper tree (one parent per non-root link, no cycles, every non-root link reachable, the
root never a child) and that `p1ᵇ` (parent frame) and `p2ᵇ` (child frame) map to the same point in
the configured geometry within 1e-9 m (non-fixed joints only; `:fixed` joints merge at the
configured geometry as given). Links joined by `:fixed` joints are merged into their parent body
(composite mass, COM and inertia); `prop_mass` adds root mass only (a point mass at the root
composite COM: no inertia, no COM shift). Link `r`/`q` of
non-root links are the configured offset and attitude relative to the bus (root) frame; the root
link's own `r`/`q` are ignored (the bus frame origin is the root link COM).
"""
function build_articulated_tree(sc::SpacecraftModel; prop_mass::Real=sc.prop_mass)
    links = Link[l for l in sc.links]
    any(l -> l === sc.root, links) || pushfirst!(links, sc.root)
    nl = length(links)
    index = IdDict{Link, Int}()
    for (i, l) in pairs(links)
        haskey(index, l) && throw(ArgumentError("Link $i appears twice in the spacecraft link list."))
        index[l] = i
    end
    iroot = index[sc.root]
    prop_mass >= 0 || throw(ArgumentError("prop_mass must be >= 0, got $prop_mass."))

    # Topology.
    parent_link = zeros(Int, nl)
    joint_of = zeros(Int, nl)
    for (jk, jt) in pairs(sc.joints)
        haskey(index, jt.link1) || throw(ArgumentError("Joint $jk: link1 is not part of the spacecraft links."))
        haskey(index, jt.link2) || throw(ArgumentError("Joint $jk: link2 is not part of the spacecraft links."))
        p = index[jt.link1]; c = index[jt.link2]
        p == c && throw(ArgumentError("Joint $jk connects link $p to itself."))
        c == iroot && throw(ArgumentError("Joint $jk: the root link ($c) cannot be the child of a joint."))
        parent_link[c] == 0 ||
            throw(ArgumentError("Joint $jk: link $c already has a parent via joint $(joint_of[c]); the structure must be a tree."))
        parent_link[c] = p
        joint_of[c] = jk
    end
    for l in 1:nl
        l == iroot && continue
        parent_link[l] != 0 ||
            throw(ArgumentError("Link $l is not reachable from the root: no joint has it as link2."))
        a = l
        for _ in 1:nl
            a = parent_link[a]
            (a == iroot || a == 0) && break
        end
        a == iroot || throw(ArgumentError("Link $l is part of a cycle and is not reachable from the root."))
    end

    # Breadth-first order from the root.
    order = Int[iroot]
    for l in order, c in 1:nl
        parent_link[c] == l && push!(order, c)
    end
    length(order) == nl || throw(ArgumentError("Not every link is reachable from the root."))

    cfg = [_link_config(links[l], l == iroot) for l in 1:nl]
    rbus = [cfg[l][1] for l in 1:nl]
    qbus = [cfg[l][2] for l in 1:nl]
    Rbus = [_rotmat(qbus[l]) for l in 1:nl]

    # Attachment-point consistency in the configured geometry.
    for l in 1:nl
        l == iroot && continue
        jt = sc.joints[joint_of[l]]
        jt.joint_type === :fixed && continue
        p = parent_link[l]
        P1 = rbus[p] + Rbus[p] * SVector{3, Float64}(jt.p1ᵇ)
        P2 = rbus[l] + Rbus[l] * SVector{3, Float64}(jt.p2ᵇ)
        err = norm(P1 - P2)
        err <= 1.0e-9 ||
            throw(ArgumentError("Joint $(joint_of[l]) (link $p -> link $l): p1ᵇ and p2ᵇ map to different points in the configured geometry (mismatch $(err) m > 1e-9 m)."))
    end

    # Group links into dynamic bodies. Head of a body: the root, or the child of a moving joint.
    body_of_link = zeros(Int, nl)
    heads = Int[]
    for l in order
        if l == iroot
            push!(heads, l); body_of_link[l] = 1
        elseif sc.joints[joint_of[l]].joint_type === :fixed
            body_of_link[l] = body_of_link[parent_link[l]]
        else
            push!(heads, l); body_of_link[l] = length(heads)
        end
    end
    nb = length(heads)

    mass = zeros(nb); comb = fill(zero(SVector{3, Float64}), nb)
    for l in 1:nl
        b = body_of_link[l]
        mass[b] += links[l].m
        comb[b] += links[l].m * rbus[l]
    end
    mass[1] += prop_mass # point mass at the root composite COM: mass only
    for b in 1:nb
        mass[b] > 0 || throw(ArgumentError("Dynamic body $b has non-positive mass."))
    end
    for b in 1:nb
        comb[b] /= (b == 1 ? mass[b] - prop_mass : mass[b])   # COM of the dry group
    end
    inertia = fill(zero(SMatrix{3, 3, Float64, 9}), nb)
    link_com = Vector{SVector{3, Float64}}(undef, nl)
    link_q = Vector{SVector{4, Float64}}(undef, nl)
    for l in 1:nl
        b = body_of_link[l]; h = heads[b]
        RB = Rbus[h]'
        d = RB * (rbus[l] - comb[b])
        Rl = RB * Rbus[l]
        link_com[l] = d
        link_q[l] = _qmul(_qconj(qbus[h]), qbus[l])
        inertia[b] += Rl * SMatrix{3, 3, Float64, 9}(links[l].inertia) * Rl' + links[l].m * (dot(d, d) * _I3 - d * d')
    end

    parent = zeros(Int, nb)
    jtype = fill(:root, nb)
    axis = fill(zero(SVector{3, Float64}), nb)
    d1 = fill(zero(SVector{3, Float64}), nb); d2 = fill(zero(SVector{3, Float64}), nb)
    C0 = fill(_I3, nb); qC0 = fill(_Q_ID, nb)
    k = zeros(nb); c = zeros(nb); rest = zeros(nb)
    kmat = fill(zero(SMatrix{3, 3, Float64, 9}), nb); cmat = fill(zero(SMatrix{3, 3, Float64, 9}), nb)
    restq = fill(_Q_ID, nb)
    qoff = zeros(Int, nb); voff = zeros(Int, nb)
    joint_of_body = zeros(Int, nb)
    nq = 0; nv = 0
    q0 = Float64[]; qd0 = Float64[]
    # Degrees of freedom follow the order of `sc.joints`.
    for (jk, jt) in pairs(sc.joints)
        jt.joint_type === :fixed && continue
        b = body_of_link[index[jt.link2]]
        joint_of_body[b] = jk
        if jt.joint_type === :ball
            qoff[b] = nq + 1; nq += 4
            voff[b] = nv + 1; nv += 3
            append!(q0, SVector{4, Float64}(jt.initial_q)); append!(qd0, SVector{3, Float64}(jt.initial_qd))
        else
            qoff[b] = nq + 1; nq += 1
            voff[b] = nv + 1; nv += 1
            push!(q0, jt.initial_q); push!(qd0, jt.initial_qd)
        end
    end
    for b in 2:nb
        h = heads[b]
        jt = sc.joints[joint_of[h]]
        pl = parent_link[h]
        A = body_of_link[pl]
        hA = heads[A]
        parent[b] = A
        jtype[b] = jt.joint_type
        P = rbus[pl] + Rbus[pl] * SVector{3, Float64}(jt.p1ᵇ)
        d1[b] = Rbus[hA]' * (P - comb[A])
        d2[b] = Rbus[h]' * (P - comb[b])
        axis[b] = (Rbus[hA]' * Rbus[pl]) * SVector{3, Float64}(jt.axis)
        C0[b] = Rbus[hA]' * Rbus[h]
        qC0[b] = _qmul(_qconj(qbus[hA]), qbus[h])
        if jt.joint_type === :ball
            kmat[b] = SMatrix{3, 3, Float64, 9}(jt.stiffness)
            cmat[b] = SMatrix{3, 3, Float64, 9}(jt.damping)
            restq[b] = SVector{4, Float64}(jt.rest)
        else
            k[b] = jt.stiffness; c[b] = jt.damping; rest[b] = jt.rest
        end
    end
    return ArticulatedTree(nb, nq, nv, parent, jtype, axis, d1, d2, C0, qC0, mass, inertia, k, c, rest, kmat, cmat,
        restq, qoff, voff, body_of_link, link_com, link_q, comb[1], q0, qd0, joint_of_body)
end

# ---------------------------------------------------------------------------
# Joint-space helpers
# ---------------------------------------------------------------------------

# Relative pose of body b w.r.t. its parent at joint coordinates `jq`:
# (Rrel child-to-parent, qrel, joint-point offset from parent COM in the parent frame).
@inline function _joint_pose(tree::ArticulatedTree, b::Int, jq)
    jt = tree.jtype[b]
    qo = tree.qoff[b]
    if jt === :hinge
        a = tree.axis[b]
        θ = jq[qo]
        return _axis_rotmat(a, θ) * tree.C0[b], _qmul(_axis_quat(a, θ), tree.qC0[b]), tree.d1[b]
    elseif jt === :slide
        s = jq[qo]
        return tree.C0[b], tree.qC0[b], tree.d1[b] + s * tree.axis[b]
    else # :ball
        qb = _qnormalize(SVector(jq[qo], jq[qo + 1], jq[qo + 2], jq[qo + 3]))
        return _rotmat(qb) * tree.C0[b], _qmul(qb, tree.qC0[b]), tree.d1[b]
    end
end

"""
    articulated_joint_qdot(tree, joint_q, joint_qd) -> Vector

Time derivative of `joint_q` for integration: the rate itself for hinge/slide, and
`q̇ = ½ [ω_rel, 0] ⊗ q` for ball joints (`ω_rel` in the parent frame, scalar-last).
"""
function articulated_joint_qdot(tree::ArticulatedTree, joint_q::AbstractVector, joint_qd::AbstractVector)
    T = promote_type(eltype(joint_q), eltype(joint_qd))
    return articulated_joint_qdot!(zeros(T, tree.nq), tree, joint_q, joint_qd)
end

"""
    articulated_joint_qdot!(out, tree, joint_q, joint_qd) -> out

In-place, allocation-free form of [`articulated_joint_qdot`](@ref).
"""
function articulated_joint_qdot!(out::AbstractVector, tree::ArticulatedTree, joint_q::AbstractVector, joint_qd::AbstractVector)
    T = promote_type(eltype(joint_q), eltype(joint_qd))
    @inbounds for b in 2:tree.nb
        qo = tree.qoff[b]; vo = tree.voff[b]
        if tree.jtype[b] === :ball
            q = SVector(joint_q[qo], joint_q[qo + 1], joint_q[qo + 2], joint_q[qo + 3])
            ω = SVector(joint_qd[vo], joint_qd[vo + 1], joint_qd[vo + 2], zero(T))
            dq = _qmul(ω, q) / 2
            for i in 1:4
                out[qo + i - 1] = dq[i]
            end
        else
            out[qo] = joint_qd[vo]
        end
    end
    return out
end

"""
    articulated_potential_energy(tree, joint_q) -> Real

Joint spring energy: `½ k (q - rest)²` for hinge/slide and `½ ϕᵀ K ϕ` for ball joints, with
`ϕ` the axis-angle of `rest⁻¹ ⊗ q`. Gravity is not included.
"""
function articulated_potential_energy(tree::ArticulatedTree, joint_q::AbstractVector)
    U = zero(eltype(joint_q))
    for b in 2:tree.nb
        qo = tree.qoff[b]
        if tree.jtype[b] === :ball
            q = _qnormalize(SVector(joint_q[qo], joint_q[qo + 1], joint_q[qo + 2], joint_q[qo + 3]))
            ϕ = _axis_angle(_qmul(_qconj(tree.restq[b]), q))
            U += dot(ϕ, tree.kmat[b] * ϕ) / 2
        else
            U += tree.k[b] * (joint_q[qo] - tree.rest[b])^2 / 2
        end
    end
    return U
end

"""
    articulated_moving_mass(tree) -> Float64

Total mass of the non-root dynamic bodies (constant in a run). The engine sets the root mass to the
spacecraft's total mass minus this value.
"""
articulated_moving_mass(tree::ArticulatedTree) = sum(@view tree.mass[2:end]; init=0.0)

"""
    articulated_has_moving_joints(sc) -> Bool

Whether any joint of `sc` is not `:fixed`. Only such spacecraft are articulated.
"""
articulated_has_moving_joints(sc::SpacecraftModel) = any(j -> j.joint_type !== :fixed, sc.joints)

"""
    ArticulatedRuntime

Per-spacecraft run data held by the engine: the immutable `tree`, a Float64 workspace used by one
thread at a time, and the constant `moving_mass`.
"""
struct ArticulatedRuntime
    tree::ArticulatedTree
    ws::ArticulatedWorkspace{Float64}
    moving_mass::Float64
end

function ArticulatedRuntime(tree::ArticulatedTree)
    return ArticulatedRuntime(tree, ArticulatedWorkspace(tree), articulated_moving_mass(tree))
end

# ---------------------------------------------------------------------------
# Kinematics outputs
# ---------------------------------------------------------------------------

"""
    articulated_forward_kinematics(tree, base_pos, base_q, joint_q) -> (pos, quat)

Inertial COM position and scalar-last attitude quaternion (active body-to-inertial, same convention
as the root state `q`) of every dynamic body (index 1 is the root). The body frame is the frame of
the body's head link. Use [`articulated_link_poses`](@ref) for every individual link.
"""
function articulated_forward_kinematics(tree::ArticulatedTree, base_pos, base_q, joint_q::AbstractVector)
    T = promote_type(eltype(base_pos), eltype(base_q), eltype(joint_q))
    pos = Vector{SVector{3, T}}(undef, tree.nb)
    quat = Vector{SVector{4, T}}(undef, tree.nb)
    return _forward_kinematics!(pos, quat, tree, base_pos, base_q, joint_q)
end

function _forward_kinematics!(pos, quat, tree::ArticulatedTree, base_pos, base_q, joint_q::AbstractVector)
    T = eltype(eltype(pos))
    pos[1] = SVector{3, T}(base_pos)
    quat[1] = _qnormalize(SVector{4, T}(base_q))
    for b in 2:tree.nb
        p = tree.parent[b]
        Rrel, qrel, off = _joint_pose(tree, b, joint_q)
        Rp = _rotmat(quat[p])
        q = _qnormalize(_qmul(quat[p], qrel))
        quat[b] = q
        pos[b] = pos[p] + Rp * off - _rotmat(q) * tree.d2[b]
    end
    return pos, quat
end

"""
    articulated_kinematics(tree, base_state, joint_q, joint_qd) -> NamedTuple

`(pos, quat, vel, ω)` per dynamic body: inertial COM position and attitude quaternion, inertial COM
velocity, and body-frame angular velocity.
"""
function articulated_kinematics(tree::ArticulatedTree, base_state, joint_q::AbstractVector, joint_qd::AbstractVector)
    T = promote_type(eltype(base_state.pos), eltype(base_state.vel), eltype(base_state.q), eltype(base_state.ω),
        eltype(joint_q), eltype(joint_qd))
    nb = tree.nb
    k = (pos=Vector{SVector{3, T}}(undef, nb), quat=Vector{SVector{4, T}}(undef, nb), vel=Vector{SVector{3, T}}(undef, nb),
        ω=Vector{SVector{3, T}}(undef, nb), ω_world=Vector{SVector{3, T}}(undef, nb))
    articulated_kinematics!(k, tree, base_state, joint_q, joint_qd)
    return (pos=k.pos, quat=k.quat, vel=k.vel, ω=k.ω)
end

"""
    articulated_kinematics!(buf, tree, base_state, joint_q, joint_qd) -> buf

In-place, allocation-free [`articulated_kinematics`](@ref): `buf` is a NamedTuple of per-body vectors
`(pos, quat, vel, ω, ω_world)` (the world-frame angular velocity is the extra scratch output).
"""
function articulated_kinematics!(buf, tree::ArticulatedTree, base_state, joint_q::AbstractVector, joint_qd::AbstractVector)
    pos = buf.pos; quat = buf.quat; vel = buf.vel; wb = buf.ω; ww = buf.ω_world
    T = eltype(eltype(pos))
    _forward_kinematics!(pos, quat, tree, base_state.pos, base_state.q, joint_q)
    nb = tree.nb
    vel[1] = SVector{3, T}(base_state.vel)
    wb[1] = SVector{3, T}(base_state.ω)
    ww[1] = _rotmat(quat[1]) * wb[1]
    for b in 2:nb
        p = tree.parent[b]
        Rp = _rotmat(quat[p]); R = _rotmat(quat[b])
        Rrel, _, off = _joint_pose(tree, b, joint_q)
        vo = tree.voff[b]
        wrel = zero(SVector{3, T}); sd = zero(T)
        if tree.jtype[b] === :hinge
            wrel = tree.axis[b] * joint_qd[vo]
        elseif tree.jtype[b] === :ball
            wrel = SVector(joint_qd[vo], joint_qd[vo + 1], joint_qd[vo + 2])
        else
            sd = joint_qd[vo]
        end
        ww[b] = ww[p] + Rp * wrel
        wb[b] = R' * ww[b]
        e = Rp * off
        f = R * tree.d2[b]
        vel[b] = vel[p] + cross(ww[p], e) + sd * (Rp * tree.axis[b]) - cross(ww[b], f)
    end
    return buf
end

"""
    articulated_link_poses(tree, pos, quat) -> (link_pos, link_quat)

Inertial COM position and attitude of every spacecraft link (including links merged through
`:fixed` joints), in the order of the spacecraft's link list used by `build_articulated_tree`
(the root link is added at the front if the list lacks it), from per-body `pos`/`quat` as
returned by [`articulated_forward_kinematics`](@ref).
"""
function articulated_link_poses(tree::ArticulatedTree, pos, quat)
    T = promote_type(eltype(eltype(pos)), eltype(eltype(quat)))
    nl = length(tree.body_of_link)
    lp = Vector{SVector{3, T}}(undef, nl)
    lq = Vector{SVector{4, T}}(undef, nl)
    for l in 1:nl
        b = tree.body_of_link[l]
        lp[l] = pos[b] + _rotmat(quat[b]) * tree.link_com_in_body[l]
        lq[l] = _qnormalize(_qmul(quat[b], tree.link_q_in_body[l]))
    end
    return lp, lq
end

# ---------------------------------------------------------------------------
# Dynamics
# ---------------------------------------------------------------------------

@inline _col(A::Array{T, 3}, j, b) where {T} = SVector{3, T}(A[1, j, b], A[2, j, b], A[3, j, b])
@inline function _setcol!(A::Array{T, 3}, j, b, v) where {T}
    A[1, j, b] = v[1]; A[2, j, b] = v[2]; A[3, j, b] = v[3]
    return nothing
end

# In-place Cholesky solve of the SPD system M x = b (lower triangle of M is overwritten by L).
function _chol_solve!(x::AbstractVector, M::AbstractMatrix, b::AbstractVector, n::Int)
    @inbounds for j in 1:n
        s = M[j, j]
        for k in 1:(j - 1)
            s -= M[j, k] * M[j, k]
        end
        s > 0 || error("Articulated mass matrix is not positive definite (pivot $j = $s).")
        ljj = sqrt(s)
        M[j, j] = ljj
        for i in (j + 1):n
            t = M[i, j]
            for k in 1:(j - 1)
                t -= M[i, k] * M[j, k]
            end
            M[i, j] = t / ljj
        end
    end
    @inbounds for i in 1:n
        t = b[i]
        for k in 1:(i - 1)
            t -= M[i, k] * x[k]
        end
        x[i] = t / M[i, i]
    end
    @inbounds for i in n:-1:1
        t = x[i]
        for k in (i + 1):n
            t -= M[k, i] * x[k]
        end
        x[i] = t / M[i, i]
    end
    return x
end

"""
    articulated_dynamics!(ws, tree, base_state, joint_q, joint_qd,
                          base_force_world, base_torque_body, body_gravity)
        -> (base_accel_world, base_angular_accel_body, joint_qdd)

Forward dynamics of the articulated tree. `base_state` provides `pos`, `vel`, `q`, `ω` of the root
body COM (see [`ArticulatedBaseState`](@ref)); `base_force_world` and `base_torque_body` are the
external NON-gravity loads on the root (about its COM). `body_gravity(r_inertial) -> g_inertial` is
evaluated at every body's own COM, the root included, so gravity-gradient effects emerge from the
body offsets. Joint springs and dampers act in joint space (`τ = -k(q - rest) - c q̇`; for ball
joints `τ = -R_rest K ϕ - C ω_rel` with `ϕ` the axis-angle of `rest⁻¹ ⊗ q`).

`body_force_world` and `body_torque_world` (default `nothing`) are optional per-dynamic-body external
wrenches: vectors of 3-vectors indexed by dynamic body (index 1 is the root), a force through the
body COM and a torque about it, both in the INERTIAL frame (the engine hands the attachment reactions
over this way). They enter as generalized forces `Jᵀ [F; τ]` through each body's Jacobian columns,
in addition to `base_force_world`/`base_torque_body` (which stay root-only).

`relative_gravity` (default `nothing`) switches to the Encke form: a callable `ρ -> g(x_root + ρ) -
g(x_root)` replaces `body_gravity` for the bodies other than the root, the system is solved with those
relative accelerations only, and the root gravity (`root_gravity`, default `body_gravity(base pos)`) is
added to the returned root acceleration. This is exact: a uniform field `g` gives `u = [g; 0]` because
the root translation Jacobian columns are the identity and have no rotational part. The common ~8 m/s^2
therefore never has to cancel inside the Cholesky solve.

`root_mass` (default `tree.mass[1]`) overrides the root body mass at runtime so mass flow needs no
tree rebuild; the root inertia stays the configured composite (propellant carries no inertia).

Returns the root COM acceleration (inertial), the root body-rate derivative, and the joint
accelerations (a vector owned by `ws`, overwritten on the next call). For `T = Float64` and a
workspace built once per tree the call performs no allocation; the math is generic over `eltype(ws)`.
"""
function articulated_dynamics!(ws::ArticulatedWorkspace{T}, tree::ArticulatedTree, base_state,
        joint_q::AbstractVector, joint_qd::AbstractVector, base_force_world, base_torque_body, body_gravity;
        root_mass=tree.mass[1], body_force_world=nothing, body_torque_world=nothing,
        relative_gravity=nothing, root_gravity=nothing) where {T}
    nb = tree.nb
    nv = tree.nv
    n = 6 + nv
    R = ws.R; x = ws.x; v = ws.v; w = ws.w; ab = ws.ab; alb = ws.alb
    Jv = ws.Jv; Jw = ws.Jw

    # Root.
    R1 = _rotmat(SVector{4, T}(base_state.q))
    ω1 = SVector{3, T}(base_state.ω)
    R[1] = R1
    # With `relative_gravity` the positions are kept ROOT-RELATIVE (the root at the origin): they are used
    # only for gravity, and ρ_b never forms an orbital-magnitude position.
    x[1] = relative_gravity === nothing ? SVector{3, T}(base_state.pos) : zero(SVector{3, T})
    v[1] = SVector{3, T}(base_state.vel)
    w[1] = R1 * ω1
    ab[1] = zero(SVector{3, T})
    alb[1] = zero(SVector{3, T})
    @inbounds for j in 1:n, i in 1:3
        Jv[i, j, 1] = zero(T)
        Jw[i, j, 1] = zero(T)
    end
    @inbounds for k in 1:3
        Jv[k, k, 1] = one(T)
        for i in 1:3
            Jw[i, 3 + k, 1] = R1[i, k]
        end
    end

    # Forward pass: poses, velocities, zero-acceleration (bias) accelerations, Jacobians.
    @inbounds for b in 2:nb
        p = tree.parent[b]
        Rp = R[p]; wp = w[p]
        Rrel, _, off = _joint_pose(tree, b, joint_q)
        Rb = Rp * Rrel
        R[b] = Rb
        jt = tree.jtype[b]
        vo = tree.voff[b]
        col = 6 + vo
        ax = Rp * tree.axis[b]                      # world axis (hinge/slide)
        wrel = zero(SVector{3, T}); sd = zero(T)
        if jt === :hinge
            wrel = ax * joint_qd[vo]
        elseif jt === :ball
            wrel = Rp * SVector{3, T}(joint_qd[vo], joint_qd[vo + 1], joint_qd[vo + 2])
        else
            sd = joint_qd[vo]
        end
        wb = wp + wrel
        w[b] = wb
        alb[b] = alb[p] + cross(wp, wrel)
        e = Rp * off
        f = Rb * tree.d2[b]
        edot = cross(wp, e) + sd * ax
        x[b] = x[p] + e - f
        v[b] = v[p] + edot - cross(wb, f)
        ebias = cross(alb[p], e) + cross(wp, edot) + sd * cross(wp, ax)
        fdot = cross(wb, f)
        fbias = cross(alb[b], f) + cross(wb, fdot)
        ab[b] = ab[p] + ebias - fbias
        for j in 1:n
            jwp = _col(Jw, j, p)
            dw = zero(SVector{3, T})
            if jt === :hinge
                j == col && (dw = ax)
            elseif jt === :ball
                if col <= j <= col + 2
                    kk = j - col + 1
                    dw = SVector{3, T}(Rp[1, kk], Rp[2, kk], Rp[3, kk])
                end
            end
            jwb = jwp + dw
            jvb = _col(Jv, j, p) - cross(e, jwp) + cross(f, jwb)
            if jt === :slide && j == col
                jvb += ax
            end
            _setcol!(Jw, j, b, jwb)
            _setcol!(Jv, j, b, jvb)
        end
    end

    # Assemble M u̇ = rhs.
    M = ws.M; rhs = ws.rhs; A = ws.A
    @inbounds for j in 1:n
        rhs[j] = zero(T)
        for i in 1:n
            M[i, j] = zero(T)
        end
    end
    @inbounds for b in 1:nb
        Rb = R[b]
        Iw = Rb * SMatrix{3, 3, T, 9}(tree.inertia[b]) * Rb'
        mb = b == 1 ? T(root_mass) : T(tree.mass[b])
        g = relative_gravity === nothing ? SVector{3, T}(body_gravity(x[b])) :
            (b == 1 ? zero(SVector{3, T}) : SVector{3, T}(relative_gravity(x[b])))
        Fb = mb * (ab[b] - g)
        Tb = Iw * alb[b] + cross(w[b], Iw * w[b])
        for j in 1:n
            jw = _col(Jw, j, b)
            jv = _col(Jv, j, b)
            rhs[j] -= dot(jv, Fb) + dot(jw, Tb)
            if body_force_world !== nothing
                # External wrench on body b about its COM: generalized force Jᵀ [F; τ].
                rhs[j] += dot(jv, SVector{3, T}(body_force_world[b])) + dot(jw, SVector{3, T}(body_torque_world[b]))
            end
            Aj = Iw * jw
            A[1, j] = Aj[1]; A[2, j] = Aj[2]; A[3, j] = Aj[3]
        end
        for j in 1:n
            jv = _col(Jv, j, b)
            for k in 1:j
                kv = _col(Jv, k, b)
                M[j, k] += mb * dot(jv, kv) + Jw[1, j, b] * A[1, k] + Jw[2, j, b] * A[2, k] + Jw[3, j, b] * A[3, k]
            end
        end
    end
    @inbounds for k in 1:3
        rhs[k] += T(base_force_world[k])
        rhs[3 + k] += T(base_torque_body[k])
    end
    # Joint-space spring/damper generalized forces.
    @inbounds for b in 2:nb
        jt = tree.jtype[b]
        qo = tree.qoff[b]; vo = tree.voff[b]
        if jt === :ball
            qb = _qnormalize(SVector{4, T}(joint_q[qo], joint_q[qo + 1], joint_q[qo + 2], joint_q[qo + 3]))
            ϕ = _axis_angle(_qmul(_qconj(tree.restq[b]), qb))
            Rr = _rotmat(tree.restq[b])
            ωr = SVector{3, T}(joint_qd[vo], joint_qd[vo + 1], joint_qd[vo + 2])
            Q = -(Rr * (tree.kmat[b] * ϕ)) - tree.cmat[b] * ωr
            for i in 1:3
                rhs[6 + vo + i - 1] += Q[i]
            end
        else
            rhs[6 + vo] += -tree.k[b] * (joint_q[qo] - tree.rest[b]) - tree.c[b] * joint_qd[vo]
        end
    end
    @inbounds for j in 1:n, i in 1:(j - 1)
        M[i, j] = M[j, i]
    end
    u = _chol_solve!(ws.u, M, rhs, n)
    qdd = ws.qdd
    @inbounds for i in 1:nv
        qdd[i] = u[6 + i]
    end
    a_root = SVector{3, T}(u[1], u[2], u[3])
    if relative_gravity !== nothing
        a_root += root_gravity === nothing ? SVector{3, T}(body_gravity(SVector{3, T}(base_state.pos))) : SVector{3, T}(root_gravity)
    end
    return a_root, SVector{3, T}(u[4], u[5], u[6]), qdd
end

end # module ArticulatedBody
