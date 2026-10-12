using Test
include(joinpath(@__DIR__, "..", "II_examples", "walker_interlinks.jl"))

@testset "Seeded Walker geometry and scheduling" begin
    first_case = walker_constellation(satellites=24, planes=6, range_m=20e6)
    repeated = walker_constellation(satellites=24, planes=6, range_m=20e6)
    different = walker_constellation(satellites=24, planes=6, seed=9)
    @test first_case.metadata == repeated.metadata
    @test first_case.metadata.raan_offset != different.metadata.raan_offset
    @test first_case.metadata.phasing in 0:5
    @test length(first_case.spacecraft) == 24
    @test length(first_case.model.linkgraph) == 23
    @test length(unique(vehicle.id for vehicle in first_case.spacecraft)) == 24
    @test all(vehicle.n_terminal == 1 for vehicle in first_case.spacecraft)
    orbits = [vehicle.initial_condition for vehicle in first_case.spacecraft]
    @test all(ic.e == 0.0 && ic.i ≈ deg2rad(53.0) for ic in orbits)
    for plane in 0:5
        base = orbits[4plane + 1]
        @test mod(base.Ω - orbits[1].Ω, 2pi) ≈ 2pi * plane / 6
        for slot in 1:3
            @test mod(orbits[4plane + slot + 1].ν - base.ν, 2pi) ≈ 2pi * slot / 4
        end
    end
    initial = walker_snapshot(first_case)
    state = deepcopy(initial)
    motion = sqrt(first_case.planet.μ / orbits[1].a^3)
    advance_walker_snapshot!(state, initial, motion, 2pi / motion)
    @test all(isapprox(a.pos, b.pos; atol=1e-7) for (a, b) in zip(initial.sc, state.sc))
    @test all(isapprox(a.vel, b.vel; atol=1e-7) for (a, b) in zip(initial.sc, state.sc))
    benchmark = benchmark_walker(first_case; duration=120.0, interval=60.0)
    @test size(benchmark.metrics, 1) == 3
    @test all(benchmark.metrics.selected .== 1)
    @test size(schedule_table(benchmark.history), 1) == 3
    @test_throws ArgumentError walker_constellation(satellites=1)
    @test_throws ArgumentError walker_constellation(satellites=25, planes=6)
    @test_throws ArgumentError walker_constellation(altitude=-1.0)
    @test_throws ArgumentError walker_constellation(target_idx=4001)
    @test_throws ArgumentError benchmark_walker(first_case; interval=0.0)
    @test_throws ArgumentError walker_configuration(first_case; duration=0.0)
end
