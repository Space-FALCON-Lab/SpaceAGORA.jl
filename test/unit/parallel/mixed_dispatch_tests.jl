using Test
using SpaceAGORA

# Mixed dispatch (V2): the process route fills the coordinator's spare threads
# with samples from the same queue the pool workers consume. Pool workers are
# --threads=1, so W workers on a T-thread coordinator otherwise leave T cores
# idle under the process route and W idle under the threads route.

const PPr = SpaceAGORA.ParallelProfiles
const SCamp = SpaceAGORA.SimulationCampaigns

_mc_feat(; samples = 64, n_sats = 16) = SCamp.campaign_route_features(
    samples = samples, n_sats = n_sats, density_family = "exponential", mission_time_s = 3600.0)

@testset "mixed_local_slots fills usable cores and keeps thread 1 for the feeders" begin
    n = Base.Threads.nthreads()
    usable = PPr.usable_core_budget()
    f = _mc_feat()
    on  = PPr.OuterRouteTuning(mixed_dispatch = true,  memory_aware = false)
    off = PPr.OuterRouteTuning(mixed_dispatch = false, memory_aware = false)
    @test PPr.mixed_local_slots(f, off, 2) == 0          # switch off: process-only
    @test PPr.mixed_local_slots(f, on, 0) == 0           # no pool: nothing to mix with
    if n > 1
        for w in 1:usable
            s = PPr.mixed_local_slots(f, on, w)
            @test s == max(0, min(n - 1, usable - w))
            @test w + s <= usable                        # never past the core budget
            @test s <= n - 1                             # thread 1 stays free
        end
    else
        @test PPr.mixed_local_slots(f, on, 1) == 0       # a 1-thread coordinator only feeds
    end
    # Capacity is workers plus slots; equals the worker count with the switch off.
    @test PPr.mixed_capacity(f, off) == PPr.effective_process_workers(f, off)
    @test PPr.mixed_capacity(f, on) == PPr.effective_process_workers(f, on) +
        PPr.mixed_local_slots(f, on, PPr.effective_process_workers(f, on))
end

@testset "the V2 Monte Carlo rule compares mixed capacity, not worker count" begin
    n = Base.Threads.nthreads()
    usable = PPr.usable_core_budget()
    f = _mc_feat()
    # Two affordable workers against an n-thread budget. Process-only loses the
    # core comparison whenever n > 2; mixed wins it whenever 2 + min(n-1, usable-2) >= n.
    base = (memory_aware = false, mc_route_by_core_budget = true, process_max_workers = 2,
            outer_thread_budget = n)
    t_off = PPr.OuterRouteTuning(; base..., mixed_dispatch = false)
    t_on  = PPr.OuterRouteTuning(; base..., mixed_dispatch = true)
    r_off = PPr.default_outer_route(f; tuning = t_off, machine_class = :small, threads_available = true)
    r_on  = PPr.default_outer_route(f; tuning = t_on,  machine_class = :small, threads_available = true)
    if n > 2
        @test r_off === :threads
    end
    # With mixed dispatch on the rule counts rounds (mc_route_rounds): the pool
    # wins outright when it saves a round and is the cold answer at a tie
    # (which is then measured, see mc_route_tie_tests).
    rounds = PPr.mc_route_rounds(f, t_on)
    if n > 1 && rounds !== nothing
        @test r_on === (rounds.process <= rounds.threads ? :process : :threads)
    end
end

@testset "_run_monte_carlo_mixed: every sample once, in order, on local slots" begin
    seeds = collect(101:112)
    spec = SCamp.MonteCarloSpec(seeds = seeds, threads = 1)
    out = SCamp._run_monte_carlo_mixed(x -> x * 2, seeds, spec, Int[], 3)
    @test [s.index for s in out] == collect(1:12)
    @test [s.value for s in out] == seeds .* 2
    @test all(s -> s.success, out)
    # fail_fast: the first failure is rethrown after the queue drains.
    specf = SCamp.MonteCarloSpec(seeds = seeds, threads = 1, fail_fast = true)
    @test_throws Exception SCamp._run_monte_carlo_mixed(
        x -> (x == 105 ? error("boom") : x), seeds, specf, Int[], 2)
    # Nothing to consume the queue is a caller bug, not a silent no-op.
    @test_throws ArgumentError SCamp._run_monte_carlo_mixed(identity, seeds, spec, Int[], 0)
end

@testset "MonteCarloResult carries the route it ran" begin
    r = SCamp.run_monte_carlo(identity, 1:3; threads = 1)
    @test r.route === :none
    @test r.local_slots == 0
    if Base.Threads.nthreads() >= 2
        r2 = SCamp.run_monte_carlo(identity, 1:4; threads = 2)
        @test r2.route === :threads
        @test r2.local_slots == 0
    end
    # The three-argument constructor keeps working for callers that predate the fields.
    e = SCamp.MonteCarloResult(SCamp.MonteCarloSampleResult[], 0.0, 0)
    @test e.route === :none && e.local_slots == 0
end

@testset "adopt_process_workers! registers without spawning" begin
    pool = SpaceAGORA.ProcessPool(Base.active_project())
    @test SpaceAGORA.adopt_process_workers!(pool, [7, 9]) == [7, 9]
    @test SpaceAGORA.adopt_process_workers!(pool, [9, 11]) == [7, 9, 11]   # idempotent on 9
end
