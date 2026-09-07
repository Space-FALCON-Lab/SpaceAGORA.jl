using Test
using SpaceAGORA
const PPr = SpaceAGORA.ParallelProfiles
const SCamp = SpaceAGORA.SimulationCampaigns
const SEng = SpaceAGORA.SimulationEngine

_feat(n) = SCamp.campaign_route_features(samples = n, n_sats = 16, density_family = "exponential", mission_time_s = 3600.0)
# Two affordable workers; the thread budget is set below the coordinator's
# thread count so the rounds comparison has both outcomes in reach.
_tuning(; kw...) = PPr.OuterRouteTuning(;
    mixed_dispatch = true, memory_aware = false, mc_route_by_core_budget = true,
    process_max_workers = 2, outer_thread_budget = max(2, Base.Threads.nthreads() ÷ 2),
    explore_route_ties = true, split_race = true, kw...)
_kw(t) = (tuning = t, machine_class = :small, threads_available = true, parallel_enabled = true)

@testset "V2 Monte Carlo rule: fewer rounds wins, a tie defaults to the pool" begin
    t = _tuning()
    for n in (2, 4, 8, 16, 64, 256)
        r = PPr.mc_route_rounds(_feat(n), t)
        @test r !== nothing
        expect = r.process <= r.threads ? :process : :threads
        @test PPr.default_outer_route(_feat(n); tuning = t, machine_class = :small, threads_available = true) === expect
        @test PPr.mc_route_tie(_feat(n), t) == (r.process == r.threads)
    end
    # Mixed dispatch off: no rounds, the raw capacity comparison as before.
    @test PPr.mc_route_rounds(_feat(64), _tuning(mixed_dispatch = false)) === nothing
    @test !PPr.mc_route_tie(_feat(64), _tuning(mixed_dispatch = false))
    # A 1-thread coordinator: the rule is moot, process stays process.
    @test PPr.default_outer_route(_feat(64); tuning = t, machine_class = :small, threads_available = false) === :process
end

@testset "a rounds tie is measured: one campaign on the other arm, then the bandit" begin
    if Base.Threads.nthreads() < 4
        @test_skip "needs julia --threads>=4"
    else
        t = _tuning()
        n_tie = findfirst(n -> PPr.mc_route_tie(_feat(n), t), 2:(8 * Base.Threads.nthreads()))
        if n_tie === nothing
            @test_skip "no rounds tie reachable on this machine's core budget"
        else
            f = _feat(n_tie + 1)
            n = f.montecarlo_samples
            rec(st, route, per) = PPr.record_outer_route_feedback!(st, f; route = route, successes = n,
                failures = 0, elapsed_success_s = per * n, tuning = t)
            st = PPr.OuterRouteState()
            # Each parallel arm gets tie_explore_min_campaigns (2) campaigns: the
            # first is the cold one and is evicted by the second.
            @test PPr.select_outer_route!(st, f; _kw(t)...) === :process        # cold: the pool
            rec(st, :process, 3.0)                                              # cold reading
            @test PPr.select_outer_route!(st, f; _kw(t)...) === :process        # the pool again, warm
            rec(st, :process, 1.0)
            @test PPr.select_outer_route!(st, f; _kw(t)...) === :threads        # the other arm
            rec(st, :threads, 0.5)
            @test PPr.select_outer_route!(st, f; _kw(t)...) === :threads
            rec(st, :threads, 0.5)
            @test PPr.select_outer_route!(st, f; _kw(t)...) === :threads        # measured faster
            # The other way round: a pool whose warm campaign is the faster wins
            # even though its cold campaign was the slowest reading of all.
            st2 = PPr.OuterRouteState()
            rec(st2, :process, 3.0); rec(st2, :process, 0.3)
            rec(st2, :threads, 0.6); rec(st2, :threads, 0.6)
            @test PPr.select_outer_route!(st2, f; _kw(t)...) === :process
            # Off, the tie is never explored: the pool answers every time.
            st3 = PPr.OuterRouteState()
            t_off = _tuning(explore_route_ties = false)
            @test PPr.select_outer_route!(st3, f; _kw(t_off)...) === :process
            PPr.record_outer_route_feedback!(st3, f; route = :process, successes = f.montecarlo_samples,
                failures = 0, elapsed_success_s = 1.0 * f.montecarlo_samples, tuning = t_off)
            @test PPr.select_outer_route!(st3, f; _kw(t_off)...) === :process
        end
    end
end

@testset "local slots are bounded by the threads route's measured width" begin
    if Base.Threads.nthreads() < 4
        @test_skip "needs julia --threads>=4"
    else
        withenv("SPACEAGORA_OUTER_ROUTE_STATE_PATH" => joinpath(mktempdir(), "no_such_state.toml")) do
            t = _tuning()
            n = 64
            f = _feat(n)
            rf = SCamp._campaign_features_for_routing(f, n)
            st = PPr.OuterRouteState()
            cold = SCamp._campaign_route_plan(rf, n; state = st, tuning = t)
            if cold.route === :process && cold.local_slots > 1
                cands = PPr.outer_split_candidates(:threads; budget = Base.Threads.nthreads(), n_units = n, tuning = t)
                w0 = first(cands)
                # One measured threads width -- the only arm, hence the best.
                PPr.record_outer_split_feedback!(st, f; route = :threads, workers = w0, successes = n,
                    failures = 0, elapsed_success_s = Float64(n), tuning = t, weight = 2)
                warm = SCamp._campaign_route_plan(rf, n; state = st, tuning = t)
                if warm.route === :process
                    @test warm.local_slots <= w0
                    @test warm.local_slots <= cold.local_slots
                    @test warm.local_slots_at(warm.threads) == warm.local_slots
                else
                    @test warm.route === :threads
                end
            else
                @test cold.route in (:threads, :process, :none)
            end
        end
    end
end

@testset "split race is gated to the threads route" begin
    withenv("SPACEAGORA_OUTER_ROUTE_STATE_PATH" => joinpath(mktempdir(), "no_such_state.toml")) do
        st = PPr.OuterRouteState()
        t = _tuning(process_max_workers = 4)
        plan = SCamp._campaign_route_plan(SCamp._campaign_features_for_routing(_feat(256), 256), 256; state = st, tuning = t)
        if plan.route === :process
            @test plan.split_race == false
        else
            @test plan.route in (:threads, :none)
        end
    end
end

@testset "GC debt is set by a threaded dispatch and cleared by the collection" begin
    SCamp._GC_DEBT[] = false
    if Base.Threads.nthreads() >= 2
        SCamp._run_campaign_with_route_env(identity, SCamp.MonteCarloSpec(seeds = 1:4, threads = 2),
            (route = :threads, threads = 2, inner_thread_budget = 1, gc_first = true, local_slots = 0))
        @test SCamp._GC_DEBT[] == true
    else
        @test_skip "needs julia --threads>=2"
    end
end

@testset "satellite_batch calibration rungs carry a width" begin
    @test SEng._make_calib_satellite_batch_plan().allotment == 1
    @test SEng._make_calib_satellite_batch_plan(6).allotment == 6
    @test SEng._make_calib_satellite_batch_plan(0).allotment == 1
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "12", "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
        @test SEng._rhs_plan_width(SEng._make_calib_satellite_batch_plan()) == 12
        @test SEng._rhs_plan_width(SEng._make_calib_satellite_batch_plan(6)) == 6
        @test SEng._rhs_plan_width(SEng._make_calib_satellite_batch_plan(64)) == 12   # never above the budget
    end
end
