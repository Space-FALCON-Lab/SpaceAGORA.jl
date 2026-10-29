using Test
using StaticArrays
using Random
using SpaceAGORA
const SM_HBP = SpaceAGORA.SimulationModel
const PERT_HBP = SM_HBP.DynamicEffectors.PerturbationEffectors

# `_harmonics_flat_batch_kernel!` must reproduce `_harmonics_scalar_force_ii`
# bit for bit — not approximately. The batch kernel carries no @turbo/@fastmath
# and nests degree, then order, then batch precisely so that a given satellite's
# accumulation runs in the scalar kernel's sequence, and the flat constellation
# route's determinism contract rests on that.
#
# The batch kernel holds a rolling three-row window over the associated-Legendre
# triangle rather than the whole triangle, so a row slot is reused every three
# degrees and again on the next call, which maps different rows onto the same
# slots. Exactly one column per row is read but never written and so has to be
# zero; which column that is depends on the recursion's bound,
# `min(max(M, 1) + 1, l)`. The `max(M, 1)` makes M = 0 a genuinely distinct case
# — the zonal-only model writes one column further than `M + 1` — and a
# regression there escaped both the square-degree benchmarks and the
# serial-vs-flat parity probe, because every case they cover has M >= 1.
#
# Hence the (L, M) matrix below, and hence `===` rather than `≈`.

const HBP_GRAVITY_FILE = joinpath(
    normpath(joinpath(@__DIR__, "..", "..", "..")),
    "data", "planets", "gravity", "EarthGGM05C.csv",
)

const HBP_LPI = SMatrix{3, 3, Float64, 9}(
    0.9362934, -0.3513767, 0.0,
    0.3513767,  0.9362934, 0.0,
    0.0,        0.0,       1.0,
)

function hbp_states(n::Int, seed::Int)
    rng = MersenneTwister(seed)
    states = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        r = 6.771e6 + 2.0e6 * rand(rng)
        u = 2 * rand(rng) - 1
        lon = 2pi * rand(rng)
        s = sqrt(1 - u^2)
        v = zeros(Float64, 8)
        v[1] = r * s * cos(lon)
        v[2] = r * s * sin(lon)
        v[3] = r * u
        v[5] = 7.5e3
        v[7] = 500.0 + 100 * rand(rng)   # mass
        states[i] = v
    end
    return states
end

# Run the batch kernel over `states` and return its per-satellite inertial force.
function hbp_batch_forces(model, states::Vector{Vector{Float64}}, calls::Int)
    B = length(states)
    work_items = collect(1:B)
    slots = zeros(Float64, 6, 1, B)
    pool = PERT_HBP._get_harmonics_batch_pool(model, 1, B)
    # More than one call on the same workspace: the second and later calls are
    # the ones that can observe a slot left dirty by the previous call.
    for _ in 1:calls
        PERT_HBP._harmonics_flat_batch_kernel!(
            slots, 1, model, states, work_items, 1, B, HBP_LPI, pool[1],
        )
    end
    return slots
end

function hbp_scalar_force(model, state::Vector{Float64})
    workspace = PERT_HBP._make_harmonics_scratch_workspace(model)
    pos = SVector{3, Float64}(state[1], state[2], state[3])
    force, _ = PERT_HBP._harmonics_scalar_force_ii(model, workspace, pos, state[7], HBP_LPI)
    return force
end

@testset "Harmonics batch kernel matches the scalar kernel bit for bit" begin
    @test isfile(HBP_GRAVITY_FILE)
    planet = SM_HBP.Earth()

    # M = 0 is zonal-only; M < L is an order-truncated field; M == L is the
    # square case the benchmarks and the parity probe already exercise.
    for (L, M) in ((20, 0), (12, 3), (20, 5), (30, 7), (50, 2), (50, 10), (50, 49), (20, 20), (50, 50))
        model = SM_HBP.GravitationalHarmonicsModel(L, M, HBP_GRAVITY_FILE, planet)
        states = hbp_states(48, 7 * L + M + 1)
        slots = hbp_batch_forces(model, states, 2)
        mismatches = 0
        for i in eachindex(states)
            expected = hbp_scalar_force(model, states[i])
            for k in 1:3
                slots[k, 1, i] === expected[k] || (mismatches += 1)
            end
        end
        @test mismatches == 0
    end
end

@testset "Harmonics batch kernel is exact across batch sizes on a reused workspace" begin
    planet = SM_HBP.Earth()
    # The workspace pool is keyed by the model and sized by the largest batch it
    # has seen, so a later, smaller batch runs against a workspace whose slots
    # still hold the wider run's values.
    for (L, M) in ((20, 0), (40, 8), (30, 30))
        model = SM_HBP.GravitationalHarmonicsModel(L, M, HBP_GRAVITY_FILE, planet)
        for B in (64, 7, 33, 1)
            states = hbp_states(B, 31 * L + 5 * M + B)
            slots = hbp_batch_forces(model, states, 2)
            mismatches = 0
            for i in eachindex(states)
                expected = hbp_scalar_force(model, states[i])
                for k in 1:3
                    slots[k, 1, i] === expected[k] || (mismatches += 1)
                end
            end
            @test mismatches == 0
        end
    end
end


@testset "Harmonics batch pool snapshot/restore" begin
    planet = SM_HBP.Earth()
    model_a = SM_HBP.GravitationalHarmonicsModel(4, 4, HBP_GRAVITY_FILE, planet)
    model_b = SM_HBP.GravitationalHarmonicsModel(4, 4, HBP_GRAVITY_FILE, planet)
    snap0 = PERT_HBP.harmonics_batch_pool_snapshot()
    try
        pool_a = PERT_HBP._get_harmonics_batch_pool(model_a, 2, 16)
        snap = PERT_HBP.harmonics_batch_pool_snapshot()
        # Visit other shapes, as a sweep would.
        PERT_HBP._get_harmonics_batch_pool(model_a, 8, 512)
        PERT_HBP._get_harmonics_batch_pool(model_b, 4, 64)
        PERT_HBP.harmonics_batch_pool_restore!(snap)
        now = PERT_HBP.harmonics_batch_pool_snapshot()
        @test now[objectid(model_a)] == [pool_a]
        @test length(pool_a) == 2 && size(pool_a[1].A, 1) == 16
        # (Compared against the snapshot rather than for absence: objectid keys can
        # be reused by a later model after GC, so a stale entry from an earlier
        # test may legitimately share model_b's key.)
        @test get(now, objectid(model_b), nothing) === get(snap, objectid(model_b), nothing)
        @test keys(now) == keys(snap)
    finally
        PERT_HBP.harmonics_batch_pool_restore!(snap0)
    end
end

@testset "Harmonics batch pool is right-sized per shape" begin
    planet = SM_HBP.Earth()
    model = SM_HBP.GravitationalHarmonicsModel(4, 4, HBP_GRAVITY_FILE, planet)
    ratio = PERT_HBP.HARMONICS_BATCH_POOL_MAX_ROW_RATIO
    maxshapes = PERT_HBP.HARMONICS_BATCH_POOL_MAX_SHAPES
    snap0 = PERT_HBP.harmonics_batch_pool_snapshot()
    try
        big = PERT_HBP._get_harmonics_batch_pool(model, 4, 1024)
        @test size(big[1].A, 1) == 1024
        # A much smaller request is NOT served from the oversized pool...
        small = PERT_HBP._get_harmonics_batch_pool(model, 16, 256)
        @test small !== big
        @test length(small) == 16 && size(small[1].A, 1) == 256
        # ...and the earlier shape is still there: alternating does not rebuild.
        @test PERT_HBP._get_harmonics_batch_pool(model, 4, 1024) === big
        @test PERT_HBP._get_harmonics_batch_pool(model, 16, 256) === small
        # A request within the ratio, with enough workers, reuses an existing pool.
        @test PERT_HBP._get_harmonics_batch_pool(model, 8, 256) === small
        @test PERT_HBP._get_harmonics_batch_pool(model, 4, 512) === big   # 1024 rows == ratio * 512
        # Too few workers in a fitting pool: rebuilt with enough.
        wider = PERT_HBP._get_harmonics_batch_pool(model, 32, 256)
        @test length(wider) >= 32
        # Bounded: never more than MAX_SHAPES pools per model.
        for b in (10, 40, 160, 640, 2560, 10240)
            PERT_HBP._get_harmonics_batch_pool(model, 2, b)
        end
        @test length(PERT_HBP.harmonics_batch_pool_snapshot()[objectid(model)]) <= maxshapes
        # Per-solve cached reference honors the same fit rule.
        ref = Ref{Any}(nothing)
        p1 = PERT_HBP._get_harmonics_batch_pool_cached!(ref, model, 2, 100)
        @test PERT_HBP._get_harmonics_batch_pool_cached!(ref, model, 2, 100) === p1
        p2 = PERT_HBP._get_harmonics_batch_pool_cached!(ref, model, 2, 10)
        @test size(p2[1].A, 1) <= ratio * 10
    finally
        PERT_HBP.harmonics_batch_pool_restore!(snap0)
    end
end

@testset "Harmonics batch kernel result is independent of the pool's row count" begin
    planet = SM_HBP.Earth()
    model = SM_HBP.GravitationalHarmonicsModel(8, 8, HBP_GRAVITY_FILE, planet)
    states = hbp_states(33, 77)
    work_items = collect(1:33)
    out = Dict{Int, Array{Float64, 3}}()
    for rows in (33, 64, 4096)
        ws = PERT_HBP._make_harmonics_batch_workspace(model, rows)
        slots = zeros(Float64, 6, 1, 33)
        for _ in 1:2
            PERT_HBP._harmonics_flat_batch_kernel!(slots, 1, model, states, work_items, 1, 33, HBP_LPI, ws)
        end
        out[rows] = slots
    end
    @test out[33] == out[64]
    @test out[33] == out[4096]
end
