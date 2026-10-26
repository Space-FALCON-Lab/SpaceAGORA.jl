using Test
using SpaceAGORA

const PU = SpaceAGORA.SimulationModel.ParallelCost

@testset "universal scalability law" begin
    @testset "shape" begin
        # beta = 0 is Amdahl: saturating, never turning over.
        @test PU.usl_speedup(0.0, 0.0, 8) ≈ 8.0
        @test PU.usl_speedup(0.1, 0.0, 12) < 12.0
        @test PU.usl_peak_workers(0.1, 0.0) == Inf
        # beta > 0 turns the curve back down, which is the behaviour a
        # saturating model cannot express and the measured RHS curves require.
        peak = PU.usl_peak_workers(0.05, 0.01)
        @test 9.0 < peak < 10.0
        @test PU.usl_speedup(0.05, 0.01, 30) < PU.usl_speedup(0.05, 0.01, Int(round(peak)))
        # Never below one: declining to split is always available.
        @test PU.usl_speedup(0.9, 0.9, 64) >= 1.0
        @test PU.usl_speedup(0.5, 0.5, 1) == 1.0
    end

    @testset "fit inverts the model" begin
        for (a, b) in ((0.05, 0.01), (0.2, 0.0), (0.0, 0.005))
            ks = [1, 2, 4, 8, 12]
            sp = [PU.usl_speedup(a, b, k) for k in ks]
            fa, fb = PU._fit_usl(ks, sp)
            @test fa ≈ a atol = 1e-6
            @test fb ≈ b atol = 1e-6
        end
        # Degenerate ladders must not throw or produce a negative beta, which
        # would predict speedup rising without bound.
        @test PU._fit_usl(Int[], Float64[]) == (0.0, 0.0)
        @test PU._fit_usl([1], [1.0]) == (0.0, 0.0)
        a, b = PU._fit_usl([2], [1.5])
        @test b == 0.0 && a >= 0.0
        # A ladder that scales better than linearly (noise) clamps rather than
        # inverting the sign of either parameter.
        a2, b2 = PU._fit_usl([2, 4, 8], [2.2, 4.5, 9.0])
        @test a2 >= 0.0 && b2 >= 0.0
    end
end
