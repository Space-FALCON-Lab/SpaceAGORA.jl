module ArticulatedBodyTests

using Test
using LinearAlgebra
using StaticArrays
using OrdinaryDiffEq
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const AB = SM.ArticulatedBody
const CM = SM.ClothMultibody

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

rotq(axis, ang) = (a = normalize(collect(Float64, axis)); SVector{4, Float64}(a[1] * sin(ang / 2), a[2] * sin(ang / 2), a[3] * sin(ang / 2), cos(ang / 2)))
Rmat(q) = AB._rotmat(SVector{4, Float64}(q))   # active body-to-inertial

function mklink(; root=false, m, dims, r=(0.0, 0.0, 0.0), q=(0.0, 0.0, 0.0, 1.0))
    return SM.Link(root=root, m=m, dims=MVector{3, Float64}(dims...), r=MVector{3, Float64}(r...),
        q=MVector{4, Float64}(q...))
end

# Joint at point P (bus frame) between a parent and child link at their configured geometry.
function mkjoint(par::SM.Link, child::SM.Link, P; kwargs...)
    rp = par.root ? zeros(3) : collect(par.r)
    Rp = par.root ? Matrix(1.0I, 3, 3) : Matrix(Rmat(par.q))
    Rc = Matrix(Rmat(child.q))
    p1 = SVector{3, Float64}(Rp' * (collect(P) - rp))
    p2 = SVector{3, Float64}(Rc' * (collect(P) - collect(child.r)))
    return SM.Joint(par, p1, child, p2; kwargs...)
end

# State layout: pos(3) vel(3) q(4) ω(3) joint_q(nq) joint_qd(nv)
function pack(tree, pos, vel, q, ω, jq, jqd)
    return vcat(collect(pos), collect(vel), collect(q), collect(ω), collect(jq), collect(jqd))
end

function unpack(tree, u)
    nq, nv = tree.nq, tree.nv
    base = AB.ArticulatedBaseState(SVector{3}(u[1:3]), SVector{3}(u[4:6]), SVector{4}(u[7:10]), SVector{3}(u[11:13]))
    return base, u[14:(13 + nq)], u[(14 + nq):(13 + nq + nv)]
end

function make_rhs(tree; force=SVector(0.0, 0.0, 0.0), torque=SVector(0.0, 0.0, 0.0), gravity=(r -> SVector(0.0, 0.0, 0.0)))
    ws = AB.ArticulatedWorkspace(tree)
    return function (du, u, p, t)
        base, jq, jqd = unpack(tree, u)
        a, α, qdd = AB.articulated_dynamics!(ws, tree, base, jq, jqd, force, torque, gravity)
        nq, nv = tree.nq, tree.nv
        du[1:3] .= base.vel
        du[4:6] .= a
        qn = SVector{4, Float64}(u[7:10])
        du[7:10] .= AB._qmul(qn, SVector(base.ω[1], base.ω[2], base.ω[3], 0.0)) / 2
        du[11:13] .= α
        du[14:(13 + nq)] .= AB.articulated_joint_qdot(tree, jq, jqd)
        du[(14 + nq):(13 + nq + nv)] .= qdd
        return nothing
    end
end

# ---------------------------------------------------------------------------
# Convention
# ---------------------------------------------------------------------------

@testset "quaternion convention" begin
    q = rotq([0.3, -0.5, 0.8], 1.1)
    # Engine rot(q) is inertial-to-body; the cloth _rot(q) and the articulated _rotmat(q) are body-to-inertial.
    @test SM.QuaternionMath.rot(q) ≈ Rmat(q)' atol = 1e-14
    @test CM._rot(q) ≈ Rmat(q) atol = 1e-14
    # Engine kinematics q̇ = ½ q ⊗ [ω_body, 0] matches a body-to-inertial quaternion.
    ω = SVector(0.1, -0.2, 0.3)
    @test SM.DynamicsRotational.quaternion_derivative(ω, q) ≈ AB._qmul(q, SVector(ω..., 0.0)) / 2 atol = 1e-14
end

# ---------------------------------------------------------------------------
# (f) Joint constructor and tree validation
# ---------------------------------------------------------------------------

@testset "Joint constructors" begin
    a = mklink(m = 1.0, dims = (1, 1, 1)); b = mklink(m = 1.0, dims = (1, 1, 1), r = (1, 0, 0))
    j = SM.Joint(a, b)
    @test j.joint_type === :fixed && j.stiffness == 0.0 && j.damping == 0.0 && j.rest == 0.0
    @test SM.Joint(a, SVector(0.5, 0, 0.0), b, SVector(-0.5, 0, 0.0)).joint_type === :fixed
    @test SM.Joint(; link1 = a, link2 = b).joint_type === :fixed
    @test SM.Joint(j).joint_type === :fixed
    h = SM.Joint(a, b; joint_type = :hinge, axis = [0, 0, 2.0], stiffness = 3, damping = 0.5, rest = 0.1)
    @test h.axis == SVector(0.0, 0.0, 1.0) && h.stiffness == 3.0 && h.damping == 0.5 && h.rest == 0.1
    hc = SM.Joint(h)
    @test hc.joint_type === :hinge && hc.axis == h.axis && hc.stiffness == 3.0 && hc.rest == 0.1
    bl = SM.Joint(a, b; joint_type = :ball, stiffness = 2.0)
    @test bl.stiffness == SMatrix{3, 3, Float64, 9}(2.0I) && bl.damping == zero(SMatrix{3, 3, Float64, 9})
    @test bl.rest == SVector(0.0, 0.0, 0.0, 1.0)
    @test SM.Joint(a, b; joint_type = :ball, stiffness = SMatrix{3, 3}(1.0, 0, 0, 0, 2.0, 0, 0, 0, 3.0)).joint_type === :ball
    @test SM.Joint(; link1 = a, link2 = b, joint_type = :slide, axis = (1, 0, 0)).joint_type === :slide

    @test_throws ArgumentError SM.Joint(a, b; joint_type = :hinge)                          # axis required
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :slide)
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :hinge, axis = [0, 0, 0])        # zero axis
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :hinge, axis = [0, 0, 1], stiffness = -1)
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :hinge, axis = [0, 0, 1], damping = -1e-3)
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :hinge, axis = [0, 0, 1], stiffness = SMatrix{3, 3}(1.0I))
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :ball, damping = SMatrix{3, 3}(-1.0I))
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :screw)
    @test_throws ArgumentError SM.Joint(a, b; joint_type = :ball, rest = [0, 0, 0, 0])
end

@testset "Tree validation" begin
    bus = mklink(root = true, m = 5.0, dims = (1, 1, 1))
    A = mklink(m = 1.0, dims = (1, 1, 1), r = (0, 1, 0))
    B = mklink(m = 1.0, dims = (1, 1, 1), r = (0, 2, 0))
    P = (0.0, 0.5, 0.0)
    ok = SpaceAGORA.SpacecraftModel(; joints = [mkjoint(bus, A, P; joint_type = :hinge, axis = [0, 0, 1])], links = [bus, A], root = bus)
    tree = AB.build_articulated_tree(ok)
    @test tree.nb == 2 && tree.nq == 1 && tree.nv == 1 && tree.parent == [0, 1]

    # Mismatched attachment point, named in the error.
    bad = SM.Joint(bus, SVector(0.0, 0.5, 0.0), A, SVector(0.0, 0.5 + 1e-3, 0.0); joint_type = :hinge, axis = [0, 0, 1])
    sc = SpaceAGORA.SpacecraftModel(; joints = [bad], links = [bus, A], root = bus)
    err = try AB.build_articulated_tree(sc); nothing catch e; e end
    @test err isa ArgumentError && occursin("Joint 1", err.msg) && occursin("p1", err.msg)

    # Root as a child.
    sc = SpaceAGORA.SpacecraftModel(; joints = [SM.Joint(A, SVector(0.0, 0, 0), bus, SVector(0.0, 0, 0))], links = [bus, A], root = bus)
    @test_throws ArgumentError AB.build_articulated_tree(sc)

    # Cycle among non-root links (neither reachable from the root).
    sc = SpaceAGORA.SpacecraftModel(; joints = [mkjoint(A, B, (0.0, 1.5, 0.0)), mkjoint(B, A, (0.0, 1.5, 0.0))], links = [bus, A, B], root = bus)
    err = try AB.build_articulated_tree(sc); nothing catch e; e end
    @test err isa ArgumentError && occursin("cycle", err.msg)

    # Two parents for one link.
    sc = SpaceAGORA.SpacecraftModel(; joints = [mkjoint(bus, B, (0.0, 1.0, 0.0)), mkjoint(A, B, (0.0, 1.5, 0.0))], links = [bus, A, B], root = bus)
    @test_throws ArgumentError AB.build_articulated_tree(sc)

    # Unreachable link.
    sc = SpaceAGORA.SpacecraftModel(; joints = [mkjoint(bus, A, P)], links = [bus, A, B], root = bus)
    @test_throws ArgumentError AB.build_articulated_tree(sc)
end

# ---------------------------------------------------------------------------
# (a) All-fixed tree equals one rigid body
# ---------------------------------------------------------------------------

@testset "all-fixed tree is a single rigid body" begin
    bus = mklink(root = true, m = 12.0, dims = (1.0, 0.8, 0.6))
    A = mklink(m = 2.0, dims = (0.1, 1.2, 0.5), r = (0.2, 1.0, -0.1), q = rotq([1, 2, 3], 0.7))
    B = mklink(m = 1.5, dims = (0.3, 0.4, 0.5), r = (-0.5, -0.9, 0.3), q = rotq([0, 1, 1], -0.4))
    C = mklink(m = 0.7, dims = (0.2, 0.2, 0.9), r = (0.6, 1.4, 0.1), q = rotq([1, 0, 0], 0.5))
    prop = 3.0
    joints = [mkjoint(bus, A, (0.0, 0.5, 0.0)), mkjoint(bus, B, (-0.4, -0.5, 0.2)), mkjoint(A, C, (0.5, 1.2, 0.0))]
    sc = SpaceAGORA.SpacecraftModel(; joints = joints, links = [bus, A, B, C], root = bus, prop_mass = prop)
    tree = AB.build_articulated_tree(sc; prop_mass = prop)
    @test tree.nb == 1 && tree.nq == 0 && tree.nv == 0

    # Independent composite properties about the composite COM, in the bus frame.
    items = [(m = bus.m, r = zeros(3), R = Matrix(1.0I, 3, 3), I = Matrix(bus.inertia))]
    for l in (A, B, C)
        push!(items, (m = l.m, r = collect(l.r), R = Matrix(SM.QuaternionMath.rot(l.q)'), I = Matrix(l.inertia)))
    end
    mdry = sum(i.m for i in items)
    com = sum(i.m * i.r for i in items) / mdry
    Itot = sum(i.R * i.I * i.R' + i.m * (dot(i.r - com, i.r - com) * I - (i.r - com) * (i.r - com)') for i in items)
    # Propellant is a point mass at the dry composite COM: mass only, no inertia, no COM shift.
    mt = mdry + prop
    @test tree.mass[1] ≈ mt rtol = 1e-14
    @test Matrix(tree.inertia[1]) ≈ Itot rtol = 1e-13
    @test collect(tree.root_com_bus) ≈ com atol = 1e-14

    F = SVector(1.3, -0.7, 2.1); T = SVector(0.4, 0.2, -0.9)
    ω = SVector(0.3, -0.2, 0.5)
    base = AB.ArticulatedBaseState(SVector(1e6, 0, 0), SVector(0, 7e3, 0), rotq([1, 1, 0], 0.4), ω)
    ws = AB.ArticulatedWorkspace(tree)
    a, α, qdd = AB.articulated_dynamics!(ws, tree, base, Float64[], Float64[], F, T, r -> SVector(0.0, 0.0, 0.0))
    α_ref = Itot \ (collect(T) - cross(collect(ω), Itot * collect(ω)))
    @test isempty(qdd)
    @test norm(a - F / mt) / norm(F / mt) < 1e-12
    @test norm(collect(α) - α_ref) / norm(α_ref) < 1e-12
    @info "(a) all-fixed relative errors" acc = norm(a - F / mt) / norm(F / mt) ang = norm(collect(α) - α_ref) / norm(α_ref)
end

# ---------------------------------------------------------------------------
# (b) two-body hinge, free flight
# ---------------------------------------------------------------------------
#
# Planar motion about the hinge axis z. Bus COM to joint: a (along y), joint to panel COM: b.
# Absolute angles: bus φ, panel φ + θ. Relative COM vector r = a ĥ(φ) + b ĥ(φ+θ), so for small angles
#   ṙ_x = -(a+b) φ̇ - b θ̇,   T = ½ Ib φ̇² + ½ Ip (φ̇+θ̇)² + ½ μ ((a+b) φ̇ + b θ̇)²,  μ = M m / (M + m).
# Mass matrix in (φ̇, θ̇):
#   M11 = Ib + Ip + μ (a+b)²,  M12 = Ip + μ b (a+b),  M22 = Ip + μ b².
# With zero total angular momentum L = M11 φ̇ + M12 θ̇ = 0 the base rotation follows the hinge,
# φ̇ = -(M12/M11) θ̇, and the hinge sees the free-free effective inertia
#   I_eff = M22 - M12² / M11,   ω0² = k / I_eff,   damped decay rate λ = c / (2 I_eff).

function hinge_setup(; k, c, θ0)
    bus = mklink(root = true, m = 10.0, dims = (1.0, 1.0, 1.0))
    a, b = 0.5, 0.6
    panel = mklink(m = 2.0, dims = (0.05, 1.0, 0.5), r = (0.0, a + b, 0.0))
    joint = mkjoint(bus, panel, (0.0, a, 0.0); joint_type = :hinge, axis = [0, 0, 1], stiffness = k, damping = c)
    sc = SpaceAGORA.SpacecraftModel(; joints = [joint], links = [bus, panel], root = bus)
    tree = AB.build_articulated_tree(sc)
    Ib = bus.inertia[3, 3]; Ip = panel.inertia[3, 3]
    μ = bus.m * panel.m / (bus.m + panel.m)
    M11 = Ib + Ip + μ * (a + b)^2; M12 = Ip + μ * b * (a + b); M22 = Ip + μ * b^2
    Ieff = M22 - M12^2 / M11
    u0 = pack(tree, zeros(3), zeros(3), [0, 0, 0, 1.0], zeros(3), [θ0], [0.0])
    return tree, u0, Ieff
end

function crossings(f, ts)   # roots of f(t) between samples, by bisection
    out = Float64[]
    fp = f(ts[1])
    for i in 2:length(ts)
        fi = f(ts[i])
        if fp * fi < 0
            lo, hi = ts[i - 1], ts[i]; flo = fp
            for _ in 1:60
                mid = (lo + hi) / 2; fm = f(mid)
                if flo * fm <= 0; hi = mid else lo = mid; flo = fm end
            end
            push!(out, (lo + hi) / 2)
        end
        fp = fi
    end
    return out
end

@testset "hinge free-free frequency and decay" begin
    k = 4.0
    tree, u0, Ieff = hinge_setup(k = k, c = 0.0, θ0 = 0.01)
    ω0 = sqrt(k / Ieff); Tper = 2π / ω0
    tend = 12Tper
    prob = ODEProblem(make_rhs(tree), u0, (0.0, tend))
    sol = solve(prob, Vern9(); reltol = 1e-13, abstol = 1e-15)
    ts = range(0, tend, length = 2400)
    zc = crossings(t -> sol(t)[14], ts)
    n = length(zc)
    @test n >= 20
    Tmeas = 2 * (zc[end] - zc[1]) / (n - 1)
    ferr = abs(Tmeas - Tper) / Tper
    @info "(b) hinge frequency" Tpred = Tper Tmeas ferr
    @test ferr < 0.005

    # Damped: decay rate of the extrema.
    c = 0.1 * 2 * sqrt(k * Ieff)         # damping ratio 0.1
    tree, u0, Ieff = hinge_setup(k = k, c = c, θ0 = 0.01)
    λ = c / (2Ieff)
    tend = 12Tper
    sold = solve(ODEProblem(make_rhs(tree), u0, (0.0, tend)), Vern9(); reltol = 1e-13, abstol = 1e-15)
    te = crossings(t -> sold(t)[15], range(0, tend, length = 2400))   # θ̇ = 0 at extrema
    amp = [abs(sold(t)[14]) for t in te]
    A = hcat(ones(length(te)), te)
    slope = (A \ log.(amp))[2]
    derr = abs(-slope - λ) / λ
    @info "(b) hinge decay" λpred = λ λmeas = -slope derr
    @test derr < 0.02
end

# ---------------------------------------------------------------------------
# (c) conservation, mixed chain
# ---------------------------------------------------------------------------

function chain_spacecraft(; prop = 0.5)
    bus = mklink(root = true, m = 20.0, dims = (1.0, 0.8, 0.6))
    A = mklink(m = 3.0, dims = (0.1, 1.2, 0.6), r = (0.0, 1.0, 0.2), q = rotq([1, 0, 0], 0.3))
    B = mklink(m = 2.0, dims = (0.2, 0.9, 0.4), r = (0.3, 2.2, 0.3), q = rotq([0, 1, 0], -0.5))
    C = mklink(m = 1.5, dims = (0.3, 0.3, 0.6), r = (0.4, 3.0, 0.1), q = rotq([1, 1, 0], 0.8))
    D = mklink(m = 2.0, dims = (0.4, 0.4, 0.4), r = (-0.8, 0.0, 0.1), q = rotq([0, 0, 1], 0.2))
    E = mklink(m = 0.8, dims = (0.1, 0.5, 0.3), r = (0.5, 1.2, -0.3), q = rotq([1, 1, 1], -0.6))
    joints = [
        mkjoint(bus, A, (0.0, 0.6, 0.1); joint_type = :hinge, axis = [0.3, 0.2, 0.9], stiffness = 6.0),
        mkjoint(A, B, (0.1, 1.7, 0.25); joint_type = :slide, axis = [1.0, 0.5, 0.0], stiffness = 15.0),
        mkjoint(B, C, (0.35, 2.7, 0.2); joint_type = :ball, stiffness = 3.0, rest = rotq([0, 0, 1], 0.2)),
        mkjoint(bus, D, (-0.5, 0.0, 0.05)),
        mkjoint(A, E, (0.3, 1.0, -0.1)),
    ]
    sc = SpaceAGORA.SpacecraftModel(; joints = joints, links = [bus, A, B, C, D, E], root = bus, prop_mass = prop)
    return sc
end

function invariants(tree, u)
    base, jq, jqd = unpack(tree, u)
    kin = AB.articulated_kinematics(tree, base, jq, jqd)
    mt = sum(tree.mass)
    X = sum(tree.mass[b] * kin.pos[b] for b in 1:tree.nb) / mt
    V = sum(tree.mass[b] * kin.vel[b] for b in 1:tree.nb) / mt
    P = mt * V
    L = zero(SVector{3, Float64}); KE = 0.0
    for b in 1:tree.nb
        R = Rmat(kin.quat[b])
        ωb = kin.ω[b]
        L += tree.mass[b] * cross(kin.pos[b] - X, kin.vel[b] - V) + R * (tree.inertia[b] * ωb)
        KE += tree.mass[b] * dot(kin.vel[b], kin.vel[b]) / 2 + dot(ωb, tree.inertia[b] * ωb) / 2
    end
    return P, L, KE + AB.articulated_potential_energy(tree, jq)
end

@testset "conservation, hinge+slide+ball chain" begin
    sc = chain_spacecraft()
    tree = AB.build_articulated_tree(sc)
    @test tree.nb == 4 && tree.nq == 6 && tree.nv == 5        # D merges into the root, E into body A
    @test tree.mass[1] ≈ 20.0 + 0.5 + 2.0
    @test SM.articulated_moving_mass(tree) ≈ sum(tree.mass[2:end])
    qb = rotq([1, 2, 1], 0.5)
    u0 = pack(tree, [1e3, -2e3, 5e2], [0.3, -0.2, 0.1], rotq([1, -1, 2], 0.9), [0.05, 0.1, -0.04],
        [0.4, 0.15, qb...], [0.5, 0.3, 0.2, -0.1, 0.15])
    P0, L0, E0 = invariants(tree, u0)
    prob = ODEProblem(make_rhs(tree), u0, (0.0, 100.0))
    sol = solve(prob, Vern9(); reltol = 1e-13, abstol = 1e-14)
    dP = 0.0; dL = 0.0; dE = 0.0
    for t in range(0, 100, length = 201)
        P, L, E = invariants(tree, sol(t))
        dP = max(dP, norm(P - P0) / norm(P0)); dL = max(dL, norm(L - L0) / norm(L0)); dE = max(dE, abs(E - E0) / abs(E0))
    end
    @info "(c) conservation drifts over 100 s" linear_momentum = dP angular_momentum = dL energy = dE
    @test dP < 1e-9
    @test dL < 1e-9
    @test dE < 1e-9
    # Energy exchange actually happened (the test is not vacuous).
    _, jq, _ = unpack(tree, sol(50.0))
    @test abs(jq[1] - 0.4) > 1e-2 || abs(jq[2] - 0.15) > 1e-2
end

# ---------------------------------------------------------------------------
# (d) cross-check against the compliant multibody model
# ---------------------------------------------------------------------------

@testset "hinge vs compliant multibody" begin
    k = 4.0; θ0 = 0.05
    tree, u0, Ieff = hinge_setup(k = k, c = 0.0, θ0 = θ0)
    ω0 = sqrt(k / Ieff); Tper = 2π / ω0; tend = 4Tper
    sola = solve(ODEProblem(make_rhs(tree), u0, (0.0, tend)), Vern9(); reltol = 1e-12, abstol = 1e-14)

    bus_m, bus_I = 10.0, mklink(root = true, m = 10.0, dims = (1.0, 1.0, 1.0)).inertia
    pan = mklink(m = 2.0, dims = (0.05, 1.0, 0.5))
    a, b = 0.5, 0.6
    Rz(θ) = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
    P = SVector(0.0, a, 0.0)
    pan_pos = P + SVector{3}(Rz(θ0) * [0, b, 0.0])
    nodes = [
        CM.CompliantTopologyNode(:bus, bus_m, bus_I, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)),
        CM.CompliantTopologyNode(:panel, pan.m, pan.inertia, pan_pos, rotq([0, 0, 1], θ0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)),
    ]
    kt = 1.0e6
    ct = 2 * sqrt(kt * 10.0 * 2.0 / 12.0) * 0.5
    krot = SMatrix{3, 3, Float64, 9}(1e4, 0, 0, 0, 1e4, 0, 0, 0, k)
    edge = CM.CompliantTopologyEdge(:hinge, 1, 2, P, SVector(0.0, -b, 0.0),
        SMatrix{3, 3, Float64, 9}(kt * I), SMatrix{3, 3, Float64, 9}(ct * I), krot, zero(SMatrix{3, 3, Float64, 9}),
        SVector(0.0, 0.0, 0.0, 1.0))
    build = CM.build_compliant_topology(nodes, [edge])
    f = (du, x, p, t) -> (du .= CM.compliant_multibody_dynamics(build.model, x, t); nothing)
    solc = solve(ODEProblem(f, build.initial_state, (0.0, tend)), Vern9(); reltol = 1e-11, abstol = 1e-13)
    ang(x, i) = (s = CM.compliant_state_parts(x, i); 2 * atan(s.q[3], s.q[4]))
    maxdiff = 0.0
    for t in range(0, tend, length = 400)
        θc = ang(solc(t), 2) - ang(solc(t), 1)
        maxdiff = max(maxdiff, abs(θc - sola(t)[14]))
    end
    rel = maxdiff / θ0
    @info "(d) hinge angle, articulated vs compliant (max abs diff / amplitude)" rel
    @test rel < 0.01
end

# ---------------------------------------------------------------------------
# (e) generic element type
# ---------------------------------------------------------------------------

@testset "BigFloat dynamics call" begin
    sc = chain_spacecraft()
    tree = AB.build_articulated_tree(sc)
    qb = rotq([1, 2, 1], 0.5)
    u = pack(tree, [1e3, -2e3, 5e2], [0.3, -0.2, 0.1], rotq([1, -1, 2], 0.9), [0.05, 0.1, -0.04],
        [0.4, 0.15, qb...], [0.5, 0.3, 0.2, -0.1, 0.15])
    grav(r) = -3.986e14 * r / norm(r)^3
    F = SVector(0.2, 0.1, -0.3); T = SVector(0.01, 0.02, 0.03)
    base, jq, jqd = unpack(tree, u)
    a, α, qdd = AB.articulated_dynamics!(AB.ArticulatedWorkspace(tree), tree, base, jq, jqd, F, T, grav)
    a0, α0, qdd0 = collect(a), collect(α), copy(qdd)
    bf(x) = BigFloat.(x)
    wsb = AB.ArticulatedWorkspace(tree, BigFloat)
    baseb = AB.ArticulatedBaseState(bf(base.pos), bf(base.vel), bf(base.q), bf(base.ω))
    gravb(r) = -BigFloat(3.986e14) * r / norm(r)^3
    ab, αb, qddb = AB.articulated_dynamics!(wsb, tree, baseb, bf(jq), bf(jqd), bf(F), bf(T), gravb)
    @test eltype(ab) == BigFloat && eltype(qddb) == BigFloat
    @test norm(Float64.(ab) - a0) / norm(a0) < 1e-10
    @test norm(Float64.(αb) - α0) / norm(α0) < 1e-10
    @test norm(Float64.(qddb) - qdd0) / norm(qdd0) < 1e-10
    # ForwardDiff is not a dependency of this project, so BigFloat stands in for duals.
end

# ---------------------------------------------------------------------------
# (g) allocations
# ---------------------------------------------------------------------------

function alloc_probe(ws, tree, base, jq, jqd, F, T, grav)
    return @allocated AB.articulated_dynamics!(ws, tree, base, jq, jqd, F, T, grav)
end

@testset "zero allocations, warmed Float64 call" begin
    sc = chain_spacecraft()
    tree = AB.build_articulated_tree(sc)
    ws = AB.ArticulatedWorkspace(tree)
    qb = rotq([1, 2, 1], 0.5)
    u = pack(tree, [1e3, -2e3, 5e2], [0.3, -0.2, 0.1], rotq([1, -1, 2], 0.9), [0.05, 0.1, -0.04],
        [0.4, 0.15, qb...], [0.5, 0.3, 0.2, -0.1, 0.15])
    base, jq, jqd = unpack(tree, u)
    grav = r -> SVector(0.0, 0.0, -9.8)
    F = SVector(0.2, 0.1, -0.3); T = SVector(0.01, 0.02, 0.03)
    AB.articulated_dynamics!(ws, tree, base, jq, jqd, F, T, grav)
    alloc = alloc_probe(ws, tree, base, jq, jqd, F, T, grav)
    @info "(g) allocations of a warmed call" alloc
    @test alloc == 0
end

end # module
