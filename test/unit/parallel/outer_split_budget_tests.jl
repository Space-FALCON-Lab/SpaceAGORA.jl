using Test
using SpaceAGORA

# A threaded outer split must tell each concurrent sample its share of the
# thread pool. Left unset, a sample resolves its inner budget to the WHOLE pool
# and, under the V2 static width rule, threads its callbacks at that width
# beside every sibling doing the same. Measured on B15 mcgrid_16sat_8mc at
# 12 threads, 8 concurrent samples: 10.33 s/sample with the pool advertised,
# 3.07 s with the share advertised. These tests pin the advertisement.

const PPol = SpaceAGORA.SimulationModel.ParallelPolicy
const SCamp = SpaceAGORA.SimulationCampaigns

@testset "outer_split_env_pairs advertises the per-sample share" begin
    n = Base.Threads.nthreads()
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => nothing,
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        pairs = Dict(SCamp.outer_split_env_pairs(8))
        @test pairs["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "1"
        @test pairs["SPACEAGORA_INNER_THREAD_BUDGET"] == string(max(1, fld(n, 8)))
        # The share is never zero, however wide the split.
        @test parse(Int, Dict(SCamp.outer_split_env_pairs(10_000))["SPACEAGORA_INNER_THREAD_BUDGET"]) == 1
    end
    # An explicit user budget always wins: the split declares itself active
    # but does not overwrite what the caller asked for.
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "3") do
        pairs = Dict(SCamp.outer_split_env_pairs(8))
        @test pairs["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "1"
        @test !haskey(pairs, "SPACEAGORA_INNER_THREAD_BUDGET")
    end
end

@testset "run_monte_carlo(threads=N) samples see fld(nthreads, N)" begin
    n = Base.Threads.nthreads()
    probe = seed -> (budget = PPol.effective_inner_thread_budget(),
                     outer  = PPol.outer_parallel_active())
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => nothing,
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        # Serial: no split, no advertisement.
        r1 = SCamp.run_monte_carlo(probe, [1]; threads = 1)
        @test r1.samples[1].value.outer == false
        @test r1.samples[1].value.budget == n
        if n >= 2
            r2 = SCamp.run_monte_carlo(probe, 1:4; threads = 2)
            for s in r2.samples
                @test s.value.outer == true
                @test s.value.budget == max(1, fld(n, 2))
            end
            # The advertisement is scoped to the dispatch, not leaked.
            @test !haskey(ENV, "SPACEAGORA_INNER_THREAD_BUDGET")
            @test PPol.outer_parallel_active() == false
        else
            @test_skip "needs julia --threads>=2 to exercise the threaded split"
        end
    end
end
