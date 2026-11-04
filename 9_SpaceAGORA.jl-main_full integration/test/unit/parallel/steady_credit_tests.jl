using Test
using SpaceAGORA
const SCamp = SpaceAGORA.SimulationCampaigns

_s(i, fin) = SCamp.MonteCarloSampleResult(index = i, seed = i, success = true, elapsed_s = 0.1, finished_ns = fin)

@testset "steady_per_sample_s reads the second half of the completions" begin
    # Four slow (cold) completions, then four fast ones: the steady figure is
    # the fast tail, the mean is dominated by the cold half.
    fins = [1.0e9, 2.0e9, 3.0e9, 4.0e9, 4.1e9, 4.2e9, 4.3e9, 4.4e9]
    r = SCamp.MonteCarloResult([_s(i, f) for (i, f) in enumerate(fins)], 4.4, 4)
    @test SCamp.steady_per_sample_s(r) ≈ 0.1 atol = 1e-9
    @test r.elapsed_s / 8 ≈ 0.55
    # Fewer than four stamps: the mean.
    r3 = SCamp.MonteCarloResult([_s(i, f) for (i, f) in enumerate(fins[1:3])], 3.0, 3)
    @test SCamp.steady_per_sample_s(r3) ≈ 1.0
    # Unstamped samples (worker-side constructions never collected): the mean.
    ru = SCamp.MonteCarloResult([_s(i, NaN) for i in 1:8], 4.4, 4)
    @test SCamp.steady_per_sample_s(ru) ≈ 0.55
    # Order of the vector does not matter; the stamps are sorted.
    rr = SCamp.MonteCarloResult([_s(i, f) for (i, f) in enumerate(reverse(fins))], 4.4, 4)
    @test SCamp.steady_per_sample_s(rr) ≈ 0.1 atol = 1e-9
end

@testset "every dispatcher stamps the samples it collects" begin
    seeds = collect(1:8)
    spec = SCamp.MonteCarloSpec(seeds = seeds, threads = 1)
    serial = SCamp._run_monte_carlo_serial(x -> x, seeds, spec)
    @test all(s -> isfinite(s.finished_ns), serial)
    @test issorted(s.finished_ns for s in serial)
    if Base.Threads.nthreads() >= 2
        spec2 = SCamp.MonteCarloSpec(seeds = seeds, threads = 2)
        thr = SCamp._run_monte_carlo_threaded(x -> x, seeds, spec2, 2)
        @test length(thr) == 8 && all(s -> isfinite(s.finished_ns), thr)
        mixed = SCamp._run_monte_carlo_mixed(x -> x, seeds, spec2, Int[], 2)   # local slots only
        @test length(mixed) == 8 && all(s -> isfinite(s.finished_ns), mixed)
        @test SCamp.steady_per_sample_s(SCamp.MonteCarloResult(mixed, 0.01, 2)) >= 0.0
    end
    # A re-indexed sample keeps its stamp.
    s = SCamp._reindexed_sample(serial[1], 42)
    @test s.index == 42 && s.finished_ns == serial[1].finished_ns
end
