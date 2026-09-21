using Test
using SpaceAGORA

# The R7 predictive planner as the campaign runner sees it:
# `run_monte_carlo(threads=:auto)` under SPACEAGORA_CAMPAIGN_PLANNER=predictive,
# with trivial sample functions.
#
# Every campaign here is routed with a tuning whose pool is one worker, so the
# process route is never offered and no `Distributed` worker is ever spawned:
# these tests are about the campaign-level plan, the result's shape and the
# feedback that must NOT be recorded, none of which need a real pool. The
# mixed-dispatch path (pool workers beside coordinator slots, and the guard
# that watches them) is exercised by the reduced-scale harness runs, which have
# one.

const SCamp = SpaceAGORA.SimulationCampaigns
const PPr = SpaceAGORA.ParallelProfiles

_feat(; samples) = SCamp.campaign_route_features(
    samples = samples, n_sats = 1, density_family = "exponential", mission_time_s = 3600.0)

# One worker means no pool: `predictive` sees no :process plan, and the bandit
# path cannot route to one either, so the two are compared on equal ground.
_tuning() = PPr.OuterRouteTuning(process_max_workers = 1, process_workers_resident = 0,
                                 memory_aware = false, mixed_dispatch = false)

_predictive(f, seeds; kwargs...) = withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "predictive") do
    SCamp.run_monte_carlo(f, seeds; threads = :auto, kwargs...)
end

@testset "a predictive campaign runs every sample once, in seed order" begin
    seeds = collect(101:140)
    state = PPr.OuterRouteState()
    r = _predictive(x -> x * 2, seeds;
                    route_features = _feat(samples = length(seeds)),
                    route_state = state, route_tuning = _tuning())
    @test length(r.samples) == length(seeds)
    @test [s.index for s in r.samples] == collect(1:length(seeds))
    @test [s.seed for s in r.samples] == seeds
    @test [s.value for s in r.samples] == seeds .* 2
    @test isempty(r.failed)
    @test r.elapsed_s >= 0.0
end

@testset "the result carries the plan the campaign ran" begin
    T = Base.Threads.nthreads()
    state = PPr.OuterRouteState()
    r = _predictive(identity, 1:32;
                    route_features = _feat(samples = 32), route_state = state,
                    route_tuning = _tuning())
    @test r.local_slots == 0        # no pool, so never mixed dispatch
    if T > 1
        @test r.route === :threads
        @test r.threads == min(32, T)
    else
        @test r.route === :none
        @test r.threads == 1
    end
    # A single sample is the serial plan on any machine.
    one = _predictive(identity, 1:1;
                      route_features = _feat(samples = 1), route_state = state,
                      route_tuning = _tuning())
    @test one.route === :none && one.threads == 1 && one.local_slots == 0
    @test length(one.samples) == 1
    # An empty campaign is still empty.
    @test isempty(_predictive(identity, Int[];
                              route_features = _feat(samples = 0), route_state = state,
                              route_tuning = _tuning()).samples)
end

@testset "the predictive path records no route feedback and the bandit path does" begin
    features = _feat(samples = 24)
    signature = PPr.outer_route_signature(features)
    predictive_state = PPr.OuterRouteState()
    _predictive(identity, 1:24; route_features = features, route_state = predictive_state,
                route_tuning = _tuning())
    # R7 has no learner: nothing is credited to any arm, so a second campaign
    # of the same shape gets the same plan for the same reason.
    @test isempty(PPr.outer_route_stats_snapshot(predictive_state, signature))
    bandit_state = PPr.OuterRouteState()
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => nothing) do
        SCamp.run_monte_carlo(identity, 1:24; threads = :auto, route_features = features,
                              route_state = bandit_state, route_tuning = _tuning())
    end
    @test !isempty(PPr.outer_route_stats_snapshot(bandit_state, signature))
end

@testset "the predictive path neither loads nor saves persisted route state" begin
    dir = mktempdir()
    path = joinpath(dir, "outer_route_state.toml")
    features = _feat(samples = 16)
    state = PPr.OuterRouteState()
    withenv("SPACEAGORA_OUTER_ROUTE_STATE_PATH" => path,
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1") do
        SCamp.reset_campaign_route_state_persistence!()
        _predictive(identity, 1:16; route_features = features, route_state = state,
                    route_tuning = _tuning())
        @test !isfile(path)
        SCamp.reset_campaign_route_state_persistence!()
    end
end

@testset "nested under an outer split the predictive path runs serially" begin
    state = PPr.OuterRouteState()
    r = withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
        _predictive(identity, 1:12; route_features = _feat(samples = 12),
                    route_state = state, route_tuning = _tuning())
    end
    @test r.route === :none
    @test r.threads == 1
    @test length(r.samples) == 12
    @test isempty(PPr.outer_route_stats_snapshot(state, PPr.outer_route_signature(_feat(samples = 12))))
end

@testset "sample failures are reported, not swallowed, and fail_fast still throws" begin
    f = x -> (x == 7 ? error("boom") : x)
    r = _predictive(f, 1:16; route_features = _feat(samples = 16),
                    route_state = PPr.OuterRouteState(), route_tuning = _tuning())
    @test length(r.failed) == 1
    @test r.failed[1].seed == 7
    @test length(r.successful) == 15
    @test_throws Exception _predictive(f, 1:16; fail_fast = true,
                                       route_features = _feat(samples = 16),
                                       route_state = PPr.OuterRouteState(),
                                       route_tuning = _tuning())
end

@testset "an unrecognized planner name fails the campaign rather than choosing one" begin
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "r7") do
        @test_throws ArgumentError SCamp.run_monte_carlo(identity, 1:4; threads = :auto,
                                                         route_features = _feat(samples = 4),
                                                         route_state = PPr.OuterRouteState(),
                                                         route_tuning = _tuning())
    end
end

@testset "the mixed dispatcher tags each sample with the class that ran it" begin
    # The guard's whole input. Without a pool only the local class appears,
    # which is enough to prove the tagging lines up with the sample indices.
    seeds = collect(1:12)
    spec = SCamp.MonteCarloSpec(seeds = seeds, threads = 1)
    sink = fill(:unset, length(seeds))
    taken = zeros(Float64, length(seeds))
    out = SCamp._run_monte_carlo_mixed(x -> x * 3, seeds, spec, Int[], 3;
                                       class_sink = sink, take_sink = taken)
    @test [s.value for s in out] == seeds .* 3
    @test all(c -> c === :local, sink)
    # Every sample was taken before it finished, and the bracket is its own.
    @test all(t -> t > 0.0, taken)
    for s in out
        @test s.finished_ns > taken[s.index]
        @test (s.finished_ns - taken[s.index]) / 1e9 >= s.elapsed_s * 0.9
    end
    @test_throws ArgumentError SCamp._run_monte_carlo_mixed(
        identity, seeds, spec, Int[], 2; take_sink = zeros(Float64, 3))
    # A sink of the wrong length is a caller bug, not a silently partial record.
    @test_throws ArgumentError SCamp._run_monte_carlo_mixed(
        identity, seeds, spec, Int[], 2; class_sink = fill(:unset, 3))
    # No sink is the bandit path, and it behaves exactly as before.
    @test [s.index for s in SCamp._run_monte_carlo_mixed(identity, seeds, spec, Int[], 2)] ==
        collect(1:12)
end

@testset "the bandit path is still the default and still races widths" begin
    # Neutrality check at the campaign level: with the switch unset the plan
    # comes from `_campaign_route_plan`, which the predictive path never calls.
    features = _feat(samples = 64)
    state = PPr.OuterRouteState()
    tuning = _tuning()
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => nothing) do
        @test SCamp.campaign_planner_mode() === :bandit
        plan = SCamp._campaign_route_plan(features, 64; state = state, tuning = tuning)
        @test plan.record
        @test hasproperty(plan, :split_race)
        r = SCamp.run_monte_carlo(identity, 1:64; threads = :auto, route_features = features,
                                  route_state = state, route_tuning = tuning)
        @test length(r.samples) == 64
    end
end
