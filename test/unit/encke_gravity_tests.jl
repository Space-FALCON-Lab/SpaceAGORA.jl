module EnckeGravityTests

# Encke differencing of point-mass gravity (`encke_point_mass_difference`) and its use in
# `articulated_dynamics!` (root-relative gravity form).

using Test
using LinearAlgebra
using StaticArrays
using Random
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const AB = SM.ArticulatedBody
const encke = SM.DynamicEffectors.GravityEffectors.encke_point_mass_difference

naive(μ, r, ρ) = (g(x) = -μ * x / norm(x)^3; g(r + ρ) - g(r))

# BigFloat reference from the SAME rounded inputs.
function bigref(μ, r, ρ)
    setprecision(BigFloat, 256) do
        return naive(BigFloat(μ), BigFloat.(r), BigFloat.(ρ))
    end
end

relerr(x, ref) = setprecision(BigFloat, 256) do
    Float64(norm(BigFloat.(x) - ref) / norm(ref))
end

function directions(T, r)
    rh = r / norm(r)
    perp = normalize(cross(rh, T.([0.3, -0.5, 0.8])))
    rng = MersenneTwister(20261006)
    rand_dirs = [normalize(T.(randn(rng, 3))) for _ in 1:3]
    return [rh, -rh, perp, rand_dirs...]
end

function sweep(T, scales)
    μ = T(3.986004418e14)
    r = SVector{3, T}(7.0e6, 1.2e5, -3.4e5)
    worst = 0.0; worst_naive = 0.0
    for s in scales, d in directions(T, r)
        ρ = SVector{3, T}(T(s) * norm(r) * d)
        ref = bigref(μ, r, ρ)
        worst = max(worst, relerr(encke(μ, r, ρ), ref))
        worst_naive = max(worst_naive, relerr(naive(μ, r, ρ), ref))
    end
    return worst, worst_naive
end

@testset "Encke point-mass difference" begin
    scales = 10.0 .^ (-12:-1)
    w64, n64 = sweep(Float64, scales)
    w32, n32 = sweep(Float32, scales)
    @info "Encke max relative error (units of eps) and naive max relative error" Float64 = w64 / eps(Float64) Float32 = w32 / eps(Float32) naive_Float64 = n64 naive_Float32 = n32
    @test w64 <= 16 * eps(Float64)
    @test w32 <= 16 * eps(Float32)
    # The point of the change: plain differencing loses the difference at small ρ/|r|.
    @test n32 > 1e3 * w32
    @test n64 > 1e3 * w64
    w32b, n32b = sweep(Float32, [1e-5])
    @info "Float32, rho/r = 1e-5" encke = w32b naive = n32b
    @test n32b > 1e3 * w32b   # observed ~1e-2 against ~1e-6: eps*r/rho is ~1e-2 here, not 1e1
    @test w32b <= 16 * eps(Float32)
    @test encke(3.986e14, SVector(7.0e6, 0.0, 0.0), SVector(0.0, 0.0, 0.0)) == SVector(0.0, 0.0, 0.0)
    @test eltype(encke(BigFloat(1), SVector{3, BigFloat}(7, 0, 0), SVector{3, BigFloat}(1, 0, 0))) == BigFloat
end

@testset "Encke allocations" begin
    μ = 3.986004418e14
    r = SVector(7.0e6, 1.2e5, -3.4e5); ρ = SVector(1.0, -2.0, 0.5)
    encke(μ, r, ρ)
    @test (@allocated encke(μ, r, ρ)) == 0 skip=Base.JLOptions().code_coverage != 0  # coverage instrumentation adds allocations
end

# ---------------------------------------------------------------------------
# articulated_dynamics!: absolute versus root-relative (Encke) form
# ---------------------------------------------------------------------------

rotq(axis, ang) = (a = normalize(collect(Float64, axis)); SVector{4, Float64}(a[1] * sin(ang / 2), a[2] * sin(ang / 2), a[3] * sin(ang / 2), cos(ang / 2)))
Rmat(q) = AB._rotmat(SVector{4, Float64}(q))
mklink(; root=false, m, dims, r=(0.0, 0.0, 0.0), q=(0.0, 0.0, 0.0, 1.0)) =
    SM.Link(root=root, m=m, dims=MVector{3, Float64}(dims...), r=MVector{3, Float64}(r...), q=MVector{4, Float64}(q...))
function mkjoint(par, child, P; kwargs...)
    rp = par.root ? zeros(3) : collect(par.r)
    Rp = par.root ? Matrix(1.0I, 3, 3) : Matrix(Rmat(par.q))
    Rc = Matrix(Rmat(child.q))
    return SM.Joint(par, SVector{3, Float64}(Rp' * (collect(P) - rp)), child, SVector{3, Float64}(Rc' * (collect(P) - collect(child.r))); kwargs...)
end

function test_tree(k = 1.0)
    bus = mklink(root = true, m = 20.0, dims = (1.0, 0.8, 0.6))
    A = mklink(m = 3.0, dims = (0.1, 1.2, 0.6), r = (0.0, 1.0, 0.2), q = rotq([1, 0, 0], 0.3))
    B = mklink(m = 2.0, dims = (0.2, 0.9, 0.4), r = (0.3, 2.2, 0.3), q = rotq([0, 1, 0], -0.5))
    joints = [
        mkjoint(bus, A, (0.0, 0.6, 0.1); joint_type = :hinge, axis = [0.3, 0.2, 0.9], stiffness = 6.0k),
        mkjoint(A, B, (0.1, 1.7, 0.25); joint_type = :ball, stiffness = 3.0k),
    ]
    return AB.build_articulated_tree(SpaceAGORA.SpacecraftModel(; joints = joints, links = [bus, A, B], root = bus, prop_mass = 0.0))
end

@testset "articulated_dynamics!: absolute vs Encke form, springs x$k" for k in (1.0, 0.0)
    tree = test_tree(k)
    μ = 3.986004418e14
    pos = [6.9e6, 1.2e6, -2.0e6]
    q0 = rotq([1, -1, 2], 0.9)
    qb = rotq([1, 2, 1], 0.5)
    jq = [0.4, qb...]
    jqd = k == 0.0 ? zeros(4) : [0.5, 0.3, 0.2, -0.1]
    F = k == 0.0 ? zero(SVector{3, Float64}) : SVector(0.2, 0.1, -0.3); T = k == 0.0 ? zero(SVector{3, Float64}) : SVector(0.01, 0.02, 0.03)
    mkbase(::Type{S}) where {S} = AB.ArticulatedBaseState(S.(pos), S.([0.3, -0.2, 7.5e3]), S.(q0), S.([0.05, 0.1, -0.04]))

    function run(S, relative)
        ws = AB.ArticulatedWorkspace(tree, S)
        base = mkbase(S)
        μS = S(μ)
        g(r) = -μS * r / norm(r)^3
        r0 = SVector{3, S}(base.pos)
        kw = relative ? (; relative_gravity = ρ -> encke(μS, r0, ρ), root_gravity = g(r0)) : (;)
        a, α, qdd = AB.articulated_dynamics!(ws, tree, base, S.(jq), S.(jqd), S.(F), S.(T), g; kw...)
        return collect(a), collect(α), collect(qdd)
    end

    setprecision(BigFloat, 256) do
        refa, refα, refq = run(BigFloat, false)
        enb = run(BigFloat, true)
        for (x, y) in zip(enb, (refa, refα, refq))
            @test maximum(abs.(x - y) ./ max.(abs.(y), 1e-30)) < 1e-60
        end
        e(x, ref) = Float64(norm(BigFloat.(x) - ref) / norm(ref))
        abs64 = run(Float64, false)
        enc64 = run(Float64, true)
        err_abs = max(e(abs64[3], refq), e(abs64[2], refα))
        err_enc = max(e(enc64[3], refq), e(enc64[2], refα))
        @info "articulated_dynamics! Float64 (k=$k) relative error vs BigFloat (qdd, alpha)" absolute_form = err_abs encke_form = err_enc agree_root = norm(enc64[1] - abs64[1]) / norm(abs64[1]) agree_qdd = norm(enc64[3] - abs64[3]) / norm(abs64[3])
        @test norm(enc64[1] - abs64[1]) / norm(abs64[1]) < 1e-13
        k == 1.0 && @test norm(enc64[3] - abs64[3]) / norm(abs64[3]) < 1e-13
        @test err_enc < 1e-13
        # Gravity-gradient-driven joint accelerations (free joints): observed 3e-15 (Encke) against 2e-12.
        k == 0.0 && @test err_enc < err_abs / 100
        k == 1.0 && @test err_abs < 1e-13
    end
end

end # module
