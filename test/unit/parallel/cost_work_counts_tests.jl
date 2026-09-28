using Test
using SpaceAGORA

const PC = SpaceAGORA.SimulationModel.ParallelCost
const SM = SpaceAGORA.SimulationModel

# Independently re-derive the ALF recurrence count by simulating the kernel's
# `max_recur_row` state machine, rather than by restating the closed form.
#
# This is the non-circular half of the test: `_alf_recurrence_iterations`
# collapses a lazily-filled, two-entry-point loop into sum_{r=3}^{L+2}
# min(M'+1, r-2), and the collapse is where an off-by-one would hide. The
# simulation below mirrors `_harmonics_flat_batch_kernel!` phase 4 line for
# line and counts what it would actually execute.
function simulate_alf_iterations(L::Int, M::Int)::Int
    max_recur_row = 2
    total = 0
    for l in 1:L
        row = l + 1
        if row > max_recur_row
            jmax = min(max(M, 1) + 1, l - 1)
            total += max(0, jmax)
            max_recur_row = row
        end
        next_row = row + 1
        if next_row > max_recur_row
            jmax_next = min(max(M, 1) + 1, l)
            total += max(0, jmax_next)
            max_recur_row = next_row
        end
    end
    return total
end

@testset "ALF recurrence count matches the kernel state machine" begin
    for L in 1:60, M in (0, 1, 2, 4, 8, 16, L)
        M > L && continue
        @test PC._alf_recurrence_iterations(L, M) == simulate_alf_iterations(L, M)
    end
end

@testset "Harmonics counts follow the active-order structure" begin
    # A zonal-only field has no active tesseral/sectoral orders, so its cost is
    # linear in L. A full field is quadratic. Conflating the two -- which the
    # triangular (L+1)(L+2)/2 count does -- overstates a zonal L=50 model by
    # better than an order of magnitude, and that is precisely the shape of the
    # gravity_*_l50_vacuum family the router mis-routes today.
    zonal(L) = [Int[] for _ in 1:(L + 2)]
    full(L) = [collect(1:min(row - 1, L)) for row in 1:(L + 2)]

    zonal_terms = Int[]
    full_terms = Int[]
    for L in (10, 20, 40, 80)
        z = PC._alf_recurrence_iterations(L, 0) + (L + 2) + 2 + 3L
        push!(zonal_terms, z)
        f = PC._alf_recurrence_iterations(L, L) + (L + 2) + (L + 2) + 3L +
            sum(length, full(L)[2:(L + 1)])
        push!(full_terms, f)
    end

    # Zonal: doubling L roughly doubles the count (linear).
    for i in 2:length(zonal_terms)
        ratio = zonal_terms[i] / zonal_terms[i - 1]
        @test 1.6 < ratio < 2.4
    end
    # Full field: doubling L roughly quadruples it (quadratic).
    for i in 2:length(full_terms)
        ratio = full_terms[i] / full_terms[i - 1]
        @test 3.0 < ratio < 4.6
    end
    # And the two must diverge, not merely differ by a constant.
    @test full_terms[end] / zonal_terms[end] > 4 * (full_terms[1] / zonal_terms[1])
end

@testset "WorkCounts algebra" begin
    a = PC.WorkCounts(simd_terms = 3.0, coeff_touches = 5.0, unknown_effectors = 1)
    b = PC.WorkCounts(simd_terms = 4.0, scalar_items = 2.0, in_domain = false)
    c = a + b
    @test c.simd_terms == 7.0
    @test c.scalar_items == 2.0
    @test c.coeff_touches == 5.0
    @test c.unknown_effectors == 1
    # in_domain is conjunctive: one unrepresentable term disqualifies the sum.
    @test c.in_domain == false
    @test (a + zero(PC.WorkCounts)).simd_terms == a.simd_terms
end

@testset "Unknown effectors fall through to the probe, never to a guess" begin
    struct _UnmodelledEffector end
    @test PC.effector_cost_terms(_UnmodelledEffector()) === nothing
    # The default must not fabricate zero-cost counts: a caller distinguishing
    # "free" from "unmeasured" depends on getting `nothing` here.
    @test PC.effector_cost_terms(_UnmodelledEffector()) !== PC.WorkCounts()
end

@testset "Closed-form gravity terms" begin
    @test PC.effector_cost_terms(SM.InverseSquaredGravityModel()).coeff_touches == 0.0
    @test PC.effector_cost_terms(SM.InverseSquaredGravityModel()).scalar_items > 0.0
end

@testset "Flat-queue predicates agree with the engine's" begin
    # ParallelCost duplicates _batchable_effector / _harmonics_prepass_effector
    # because it is included below SimulationEngine and cannot call upward.
    # This test is what makes that duplication safe: if the engine's predicates
    # change and the mirrored copies do not, queue_nodes silently drifts from
    # the node count the RHS actually builds, and the cost model mis-ranks the
    # routing candidates with no other symptom.
    SE = SpaceAGORA.SimulationEngine
    E = SM.Earth()
    harm_file = joinpath(@__DIR__, "..", "..", "..", "data",
                         "Gravity_harmonics_data", "EarthGGM05C.csv")

    effectors = Any[
        SM.InverseSquaredGravityModel(),
        SM.InverseSquaredJ2GravityModel(),
        SM.NBodyGravityModel(["Sun"], "Earth"),
        SM.SolarRadiationPressureModel(1.8, 10.0),
        # gravity_gradient flips _batchable_effector for the central terms, so
        # both variants must be checked -- a type-only mirror would pass the
        # default case and silently disagree here.
        SM.InverseSquaredGravityModel(gravity_gradient=true),
        SM.InverseSquaredJ2GravityModel(gravity_gradient=true),
    ]
    if isfile(harm_file)
        push!(effectors, SM.GravitationalHarmonicsModel(4, 4, harm_file, E))
    end

    for e in effectors
        @test PC._cost_batchable_effector(e) == SE._batchable_effector(e)
        @test PC._cost_harmonics_prepass_effector(e) == SE._harmonics_prepass_effector(e)
    end

    # The mirrored node predicate must reproduce the engine's queue-only count.
    tup = Tuple(effectors)
    expected = SE._count_flat_queue_only_effectors(tup)
    got = count(e -> PC.flat_queue_node_effector(e), effectors)
    @test got == expected

    # And a solver partition puts every effector back on the queue.
    @test all(e -> PC.flat_queue_node_effector(e; partition_active=true), effectors)
end
