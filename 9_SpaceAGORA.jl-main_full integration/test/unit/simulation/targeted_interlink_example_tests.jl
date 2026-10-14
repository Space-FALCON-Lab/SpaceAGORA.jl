module TargetedInterlinkExampleTests
using Test
include(joinpath(@__DIR__, "..", "..", "..", "examples", "targeted_40sat_sma_interlinks.jl"))

@testset "30 helper / 10 target constellation" begin
    args = targeted_40sat_configuration(save_results=false)
    sc = args.dynamics_model.spacecraft
    @test length(sc) == 40
    @test args.scheduling_policy_model.target_idxs == collect(31:40)
    @test length(args.interlink_model.graph.edges) == 300
    for (ring, altitude, spacing, offset) in ((1:30, 1000e3, 12.0, 0.0), (31:40, 1100e3, 36.0, 1.0))
        phases = Float64[]
        for (k, i) in enumerate(ring)
            ic = sc[i].initial_condition
            @test ic.a == args.environment_model.planet.Rp_e + altitude
            @test ic.e == ic.i == ic.Ω == 0.0
            @test rad2deg(ic.ν) ≈ offset + spacing * (k - 1) atol=1e-12
            @test sc[i].n_terminal == 1
            push!(phases, rad2deg(ic.ν))
        end
        @test all(isapprox.(diff(vcat(phases, first(phases) + 360)), spacing; atol=1e-12))
    end
    @test all(key -> key[1][1] in 1:30 && key[2][1] in 31:40,
        keys(args.interlink_model.graph.edge_index_lookup))
    @test_throws ArgumentError targeted_40sat_configuration(orbits=0.0)
    solution = run_targeted_40sat_example(save_results=false, orbits=0.01, animate=false)
    @test all(isfinite, solution.u[end])
    @test solution.t[end] ≈ args.mission_configuration.mission_time / 2000
    history = solution.prob.p.args.interlink_model.history
    @test length(history) == solution.stats.naccept + 1
    @test any(sample -> !isempty(sample.active), history)
    @test all(sample -> length(sample.active) <= 10, history)
    frames = interlink_animation_frames(solution; frame_count=3)
    @test length(frames.endpoints) == 300
    @test all(status -> length(status) == 300 && all(in(0:2), status), frames.status)
    mktempdir() do dir
        path = save_interlink_animation(solution; output=dir, frame_count=2)
        @test isfile(path)
        @test String(read(path)[1:6]) in ("GIF87a", "GIF89a")
    end
end
end
