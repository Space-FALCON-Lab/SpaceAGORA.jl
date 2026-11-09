module ChangingTreeInterlinkTests

using Test
using DataFrames
using LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "..", "examples", "changing_tree_interlinks.jl"))
const C = ChangingTreeInterlinks

@testset "Changing tree configuration and event extraction" begin
    short = C.K.configuration()
    long = C.K.configuration(duration=C.DURATION_S)
    @test short.mission_configuration.mission_time == 60.0
    @test long.mission_configuration.mission_time == 2400.0
    @test long.interlink_model.record_available === true
    @test keys(long.interlink_model.graph.edge_index_lookup) == keys(short.interlink_model.graph.edge_index_lookup)
    for (index, (a, b)) in enumerate(zip(short.dynamics_model.spacecraft, long.dynamics_model.spacecraft))
        @test a.initial_condition.a == b.initial_condition.a
        @test a.initial_condition.ν == b.initial_condition.ν
        @test a.n_terminal == b.n_terminal == 1
        @test a.initial_condition.i == b.initial_condition.i == 0.0
        @test a.initial_condition.e == b.initial_condition.e == 0.0
        @test b.initial_condition.a - long.environment_model.planet.Rp_e == C.K.ALTITUDES_M[index]
    end
    for duration in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError C.K.configuration(; duration)
    end
    comparison = DataFrame(time_s=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        available=["1-3;2-3;2-4", "1-3;2-3;2-4", "1-3;2-3;2-4", "1-3;2-3;2-4", "2-3", "empty"],
        active=["1-3;2-4", "1-3;2-4", "2-3", "empty", "empty", "empty"])
    changes = C.status_changes(comparison)
    @test changes.time_s == [0.0, 2.0, 3.0, 4.0, 5.0]
    @test changes.active == ["1-3;2-4", "2-3", "empty", "empty", "empty"]
    @test nrow(comparison) == 6
    @test C.status_changes(comparison[1:1, :]) == comparison[1:1, :]
    @test_throws ArgumentError C.status_changes(comparison[1:0, :])
end

@testset "Changing-partner endpoint velocity checks" begin
    positions = [[0.0, 0, 0], [0.0, 1000, 0], [1000.0, 0, 0], [1000.0, 1000, 0]]
    pairs = [collect(C.K.PREDICTED_PAIRS), [(2, 3)], Tuple{Int, Int}[], Tuple{Int, Int}[]]
    history = [(; time=Float64(index - 1), active=C.K.link_key.(active))
        for (index, active) in enumerate(pairs)]
    # Fixed geometry makes the impulses exact: outer pairs, middle pair, then idle.
    outer = C.K.FORCE_N * [-1.0, 0, 0]
    middle = C.K.FORCE_N / sqrt(2.0) * [-1.0, 1, 0]
    impulses = [[zeros(3) for _ in 1:4], [outer, outer, -outer, -outer],
        [outer, outer + middle, -outer - middle, -outer],
        [outer, outer + middle, -outer - middle, -outer]]
    initial_dv = [0.01, -0.02, 0.03]
    states = [(; sc=[(; pos=positions[index], mass=C.K.MASSES_KG[index],
        interlink_dv=initial_dv + impulse[index] / C.K.MASSES_KG[index]) for index in 1:4])
        for impulse in impulses]
    dv = C.K.check_endpoint_dv(history, states; fixed_prediction=false)
    @test dv.active_duration_s == [1.0, 2.0, 2.0, 1.0]
    @test all(dv.error_norm_m_per_s .<= dv.tolerance_m_per_s)
    @test isapprox(dv.recorded_dvx_m_per_s,
        [last(impulses)[index][1] / C.K.MASSES_KG[index] for index in 1:4]; atol=1e-15, rtol=1e-12)
    @test isapprox(dv.recorded_dvy_m_per_s,
        [last(impulses)[index][2] / C.K.MASSES_KG[index] for index in 1:4]; atol=1e-15, rtol=1e-12)
    @test dv.recorded_dvz_m_per_s == zeros(4)
    for (i, j) in C.K.PREDICTED_PAIRS
        @test norm(last(impulses)[i] + last(impulses)[j]) > C.K.FORCE_N / 2
    end
end

end # module ChangingTreeInterlinkTests
