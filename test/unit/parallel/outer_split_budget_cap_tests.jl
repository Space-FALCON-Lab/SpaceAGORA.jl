using Test
using SpaceAGORA

# An outer split caps an inherited inner thread budget; it does not obey it.
#
# `SPACEAGORA_INNER_THREAD_BUDGET` used to outrank the split's own share
# outright, on the reasoning that an explicit user budget always wins. That is
# the wrong way round when the inherited value is the WIDER of the two: it was
# written by something that did not know how many samples would run beside
# each other, and the split is the only thing that does. Honoring a whole-pool
# budget under a W-wide split hands each of W concurrent samples the whole
# pool -- the overstatement `outer_split_env_pairs` is documented to prevent --
# and nothing downstream can see it, because the route and the width are both
# still right.
#
# Measured on this repo's 24-logical-core workstation, mcgrid_8sat_16mc at
# 32 threads with 16 samples (R6, warm store, idle box): 1.57 s per campaign
# and 43.9 GB summed sample allocation with the share declared, against 2.50 s
# and 69.9 GB with an inherited budget of 32 honored -- 1.59x on both, with
# policy_threads_enabled_total 0 against 874.
#
# A NARROWER inherited budget still wins: lowering concurrency is always safe.

const SCampCap = SpaceAGORA.SimulationCampaigns
const PPolCap = SpaceAGORA.SimulationModel.ParallelPolicy

@testset "capped_inner_thread_budget caps a wider inherited budget" begin
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
        # Nothing inherited: the split declares its own share.
        @test SCampCap.capped_inner_thread_budget(4) == 4
        @test SCampCap.capped_inner_thread_budget(0) == 1
    end
    # Wider than the share: capped to the share. This is the regression.
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "32") do
        @test SCampCap.capped_inner_thread_budget(2) == 2
        @test SCampCap.capped_inner_thread_budget(31) == 31
    end
    # At or below the share: left alone, so an explicit narrower budget stands.
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "2") do
        @test SCampCap.capped_inner_thread_budget(8) === nothing
        @test SCampCap.capped_inner_thread_budget(2) === nothing
    end
    # Non-positive is how the policy layer spells "the whole pool", so it is
    # the widest value there is and must be capped like any other.
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "0") do
        @test SCampCap.capped_inner_thread_budget(3) == 3
    end
    # Unparseable: the split declines to rewrite a setting it cannot read,
    # rather than silently discarding what the caller asked for.
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "wide") do
        @test SCampCap.capped_inner_thread_budget(3) === nothing
    end
end

@testset "outer_split_env_pairs caps a wider inherited budget" begin
    n = Base.Threads.nthreads()
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(max(2, 4 * n))) do
        pairs = Dict(SCampCap.outer_split_env_pairs(8))
        @test pairs["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "1"
        # Pre-fix this key was absent, leaving every concurrent sample on the
        # inherited whole-pool budget.
        @test haskey(pairs, "SPACEAGORA_INNER_THREAD_BUDGET")
        @test pairs["SPACEAGORA_INNER_THREAD_BUDGET"] == string(max(1, fld(n, 8)))
    end
end

@testset "a threaded split's samples never see more than their share" begin
    n = Base.Threads.nthreads()
    if n < 2
        @test_skip "needs julia --threads>=2 to exercise the threaded split"
    else
        probe = seed -> PPolCap.effective_inner_thread_budget()
        share = max(1, fld(n, 2))
        # An inherited budget wider than the split can afford. Pre-fix every
        # sample resolved to `n`; post-fix none exceeds its share.
        withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(4 * n),
                "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
            r = SCampCap.run_monte_carlo(probe, 1:4; threads = 2)
            @test all(s -> s.value == share, r.samples)
            # Scoped to the dispatch: what the caller set is restored.
            @test ENV["SPACEAGORA_INNER_THREAD_BUDGET"] == string(4 * n)
        end
        # A narrower inherited budget is still honored.
        withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "1",
                "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
            r = SCampCap.run_monte_carlo(probe, 1:4; threads = 2)
            @test all(s -> s.value == 1, r.samples)
        end
    end
end

@testset "the adaptive route's samples never see more than their share" begin
    n = Base.Threads.nthreads()
    if n < 2
        @test_skip "needs julia --threads>=2 to exercise the adaptive split"
    else
        PPro = SpaceAGORA.ParallelProfiles
        features = PPro.OuterRouteFeatures(
            category = "montecarlo", n_sats = 8, mission_time_s = 3600.0,
            harmonics_degree = 50, density_family = "exponential",
            dynamic_effector_count = 2, montecarlo_samples = 8,
        )
        state = PPro.OuterRouteState()
        probe = seed -> (budget = PPolCap.effective_inner_thread_budget(),
                         outer = PPolCap.outer_parallel_active())
        withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1",
                # One process worker withdraws the process route, so this
                # exercises the threads route the collapse was observed on.
                "SPACEAGORA_PERF_PROCS" => "1",
                "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
                "SPACEAGORA_INNER_THREAD_BUDGET" => string(4 * n),
                "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => tempname()) do
            r = SCampCap.run_monte_carlo(probe, collect(1:8); threads = :auto,
                                         route_features = features, route_state = state)
            @test all(s -> s.value.outer, r.samples)
            # Whatever width the router picked, no sample may claim more than
            # the pool divided by the number of samples running beside it.
            width = max(1, r.route === :process ? max(1, r.local_slots) : r.threads)
            @test all(s -> s.value.budget <= max(1, fld(n, width)), r.samples)
            @test all(s -> s.value.budget < 4 * n, r.samples)
        end
    end
end
