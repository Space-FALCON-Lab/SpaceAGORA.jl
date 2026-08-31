using Test
using SpaceAGORA
using TOML

const PPr = SpaceAGORA.ParallelProfiles
const SCamp = SpaceAGORA.SimulationCampaigns

function _feat(; samples = 64)
    return SCamp.campaign_route_features(
        samples = samples, n_sats = 1,
        density_family = "exponential", mission_time_s = 3600.0,
    )
end

function _record!(state, feat, route, mean_s; n = 64, reps = 5)
    for _ in 1:reps
        PPr.record_outer_route_feedback!(
            state, feat; route = route, successes = n, failures = 0,
            elapsed_success_s = n * mean_s, elapsed_success_sq_sum_s = n * mean_s^2,
        )
    end
end

@testset "Outer-route state survives a process boundary" begin
    # save/load_outer_route_state existed and were exported from the start and
    # nothing in src/ called them, so the bandit began cold in every process and
    # re-paid for exploration it had already done. This asserts the wiring, not
    # the serialiser.
    mktempdir() do dir
        path = joinpath(dir, "route.toml")
        feat = _feat()
        sig = PPr.outer_route_signature(feat)

        withenv("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1",
                "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => path) do
            learned = PPr.OuterRouteState()
            _record!(learned, feat, :process, 0.10)
            _record!(learned, feat, :threads, 0.90)
            _record!(learned, feat, :none, 1.50)
            SCamp.save_campaign_route_state(learned)
            @test isfile(path)

            fresh = PPr.OuterRouteState()
            SCamp.reset_campaign_route_state_persistence!()
            SCamp.ensure_campaign_route_state_loaded!(fresh)
            snap = PPr.outer_route_stats_snapshot(fresh, sig)
            @test Set(keys(snap)) == Set([:none, :threads, :process])
            @test snap[:process].mean_s ≈ 0.10 atol = 1e-9
            @test snap[:threads].mean_s ≈ 0.90 atol = 1e-9

            # The point of persisting: a fresh state exploits immediately rather
            # than spending trials rediscovering what the last run measured.
            @test PPr.select_outer_route!(
                fresh, feat; machine_class = :large, threads_available = true) === :process
        end
    end
end

@testset "Persistence is opt-in and never throws" begin
    mktempdir() do dir
        path = joinpath(dir, "route.toml")
        feat = _feat()
        # Disabled: nothing is written, and loading is a no-op.
        withenv("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
                "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
                "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => path) do
            st = PPr.OuterRouteState()
            _record!(st, feat, :process, 0.1)
            SCamp.save_campaign_route_state(st)
            @test !isfile(path)
        end
        # A malformed state file must degrade to cold start, not break a run.
        bad = joinpath(dir, "bad.toml")
        write(bad, "not valid toml {{{")
        withenv("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1",
                "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => bad) do
            SCamp.reset_campaign_route_state_persistence!()
            fresh = PPr.OuterRouteState()
            @test SCamp.ensure_campaign_route_state_loaded!(fresh) === nothing
        end
    end
end

@testset "A proven default is not re-explored" begin
    # _under_sampled_candidate guarantees every candidate a trial before the
    # selector will exploit. Cold that is right; with restored history it is a
    # regression -- a signature carrying hundreds of observations of :process at
    # 0.1 s would still divert to an unsampled :none, which cold start would
    # never have chosen because default_outer_route answers :process outright.
    feat = _feat()
    partial = PPr.OuterRouteState()
    _record!(partial, feat, :process, 0.10)
    _record!(partial, feat, :threads, 0.90)
    # :none deliberately absent.
    @test PPr.select_outer_route!(
        partial, feat; machine_class = :large, threads_available = true) === :process

    # But an unproven default still gets explored: one sample is not evidence.
    thin = PPr.OuterRouteState()
    PPr.record_outer_route_feedback!(
        thin, feat; route = :process, successes = 1, failures = 0,
        elapsed_success_s = 0.1, elapsed_success_sq_sum_s = 0.01)
    chosen_thin = PPr.select_outer_route!(
        thin, feat; machine_class = :large, threads_available = true)
    @test chosen_thin in (:none, :threads, :process)

    # And a default that is genuinely beaten must not be treated as proven.
    beaten = PPr.OuterRouteState()
    _record!(beaten, feat, :process, 0.90)
    _record!(beaten, feat, :threads, 0.10)
    @test PPr._route_is_proven(
        Symbol[:none, :threads, :process],
        PPr.outer_route_stats_snapshot(beaten, PPr.outer_route_signature(feat)),
        :process, 2) == false
end

@testset "A default tested against nothing is not proven" begin
    # The guard must not fire before any alternative has been tried. "No sampled
    # route beats the default" is vacuously true when the default is the only
    # sampled route, and stopping exploration there would mean the router never
    # learns anything -- it would lock in whatever `default_outer_route` guessed
    # on the first campaign. The integration suite caught exactly this.
    feat = _feat()
    only_default = PPr.OuterRouteState()
    _record!(only_default, feat, :process, 0.10)
    snap = PPr.outer_route_stats_snapshot(only_default, PPr.outer_route_signature(feat))
    @test PPr._route_is_proven(Symbol[:none, :threads, :process], snap, :process, 2) == false

    # One sampled alternative that loses is enough to make it proven.
    _record!(only_default, feat, :threads, 0.90)
    snap2 = PPr.outer_route_stats_snapshot(only_default, PPr.outer_route_signature(feat))
    @test PPr._route_is_proven(Symbol[:none, :threads, :process], snap2, :process, 2) == true
end

@testset "Outer/inner split: candidate ladder" begin
    T = PPr.OuterRouteTuning(process_max_workers = 12)
    # Geometric, and always including the widest split -- the widest is what the
    # previous arithmetic rule produced, so the selector must be able to
    # reproduce it exactly.
    @test PPr.outer_split_candidates(:threads; budget = 12, n_units = 64, tuning = T) == [1, 2, 4, 8, 12]
    # Capped by the work available, not just the budget: four samples cannot
    # occupy twelve workers.
    @test PPr.outer_split_candidates(:threads; budget = 12, n_units = 4, tuning = T) == [1, 2, 4]
    # A serial route has exactly one width.
    @test PPr.outer_split_candidates(:none; budget = 12, n_units = 64, tuning = T) == [1]
    # Process width is bounded by process_max_workers rather than the thread pool.
    @test PPr.outer_split_candidates(:process; budget = 2, n_units = 64, tuning = T) == [1, 2, 4, 8, 12]
end

@testset "Outer/inner split: cold behaviour matches the old arithmetic" begin
    # With no history the selector must return the widest split, because that is
    # exactly what the fixed `min(n_samples, nthreads)` rule produced. An
    # uncalibrated machine therefore behaves as it did before the adaptation
    # existed, which is what makes this change unable to regress a cold run.
    feat = _feat()
    cold = PPr.OuterRouteState()
    for route in (:threads, :process)
        w = PPr.select_outer_split!(cold, feat; route = route, budget = 12,
                                    n_units = 64, tuning = PPr.OuterRouteTuning())
        @test w == PPr.outer_split_candidates(route; budget = 12, n_units = 64,
                                              tuning = PPr.OuterRouteTuning())[end]
    end
    # A route with a single feasible width is returned without consulting history.
    @test PPr.select_outer_split!(cold, feat; route = :none, budget = 12, n_units = 64) == 1
end

@testset "Outer/inner split: history moves the choice" begin
    feat = _feat()
    st = PPr.OuterRouteState()
    T = PPr.OuterRouteTuning()
    # Teach it that a narrow split is much faster than a wide one -- the case
    # where each sample can use several threads and starving it to widen the
    # queue is the wrong trade.
    for w in PPr.outer_split_candidates(:threads; budget = 12, n_units = 64, tuning = T)
        mean_s = w == 4 ? 0.10 : 0.90
        for _ in 1:5
            PPr.record_outer_split_feedback!(
                st, feat; route = :threads, workers = w, successes = 64, failures = 0,
                elapsed_success_s = 64 * mean_s, elapsed_success_sq_sum_s = 64 * mean_s^2)
        end
    end
    @test PPr.select_outer_split!(st, feat; route = :threads, budget = 12,
                                  n_units = 64, tuning = T) == 4

    # Split arms must not disturb route selection: the two share a bucket, and
    # each selector may only score the arms it enumerated.
    _record!(st, feat, :process, 0.10)
    _record!(st, feat, :threads, 0.90)
    _record!(st, feat, :none, 1.50)
    @test PPr.select_outer_route!(st, feat; machine_class = :large,
                                  threads_available = true) === :process
end

@testset "Outer/inner split: arms survive persistence" begin
    # save_outer_route_state enumerated a fixed three route symbols, which would
    # have silently dropped every split arm -- so a restored state would carry
    # route history but no split history, and the split bandit would re-explore
    # from cold on every run while the route bandit did not.
    mktempdir() do dir
        path = joinpath(dir, "split.toml")
        feat = _feat()
        st = PPr.OuterRouteState()
        for _ in 1:4
            PPr.record_outer_split_feedback!(
                st, feat; route = :threads, workers = 4, successes = 64, failures = 0,
                elapsed_success_s = 64 * 0.1, elapsed_success_sq_sum_s = 64 * 0.01)
        end
        PPr.save_outer_route_state(st, path)

        back = PPr.OuterRouteState()
        PPr.load_outer_route_state!(back, path)
        snap = PPr.outer_route_stats_snapshot(
            back, "split|" * PPr.outer_route_signature(feat))
        @test haskey(snap, Symbol("split_threads_w4"))
        @test snap[Symbol("split_threads_w4")].mean_s ≈ 0.1 atol = 1e-9
    end
end

@testset "Split arms cannot make a route look proven" begin
    # The route and split selectors share a per-signature bucket, and
    # _route_is_proven must only weigh the arms its own selector enumerated.
    # Scanning the whole bucket meant one campaign's split arm -- which carries
    # that same campaign's timing under a second name -- counted as an
    # independent alternative, so the route was declared proven after a single
    # observation and exploration stopped early.
    feat = _feat()
    st = PPr.OuterRouteState()
    _record!(st, feat, :threads, 0.50)
    PPr.record_outer_split_feedback!(
        st, feat; route = :threads, workers = 4, successes = 64, failures = 0,
        elapsed_success_s = 64 * 0.50, elapsed_success_sq_sum_s = 64 * 0.25)
    # Split arms are not in the route signature's bucket at all.
    snap = PPr.outer_route_stats_snapshot(st, PPr.outer_route_signature(feat))
    @test !haskey(snap, Symbol("split_threads_w4"))
    # Route arm only: 5 reps x 64 samples. Before the namespace split this
    # summed to 384, because the split arm's 64 were counted again here.
    @test sum(i.samples for i in values(snap)) == 320
    split_snap = PPr.outer_route_stats_snapshot(
        st, "split|" * PPr.outer_route_signature(feat))
    @test haskey(split_snap, Symbol("split_threads_w4"))
    # :threads has no route-level alternative measured, so it is NOT proven.
    @test PPr._route_is_proven(Symbol[:none, :threads, :process], snap, :threads, 2) == false
end

@testset "Observation counters survive the round trip" begin
    # Both counters gate behaviour and neither was serialised, so a restored
    # arm read as never-observed (campaigns) and had its whole history replaced
    # by the second live timing (observations, via the cold-eviction branch).
    mktempdir() do dir
        path = joinpath(dir, "counters.toml")
        feat = _feat()
        sig = PPr.outer_route_signature(feat)

        learned = PPr.OuterRouteState()
        _record!(learned, feat, :process, 0.10)   # reps = 5
        _record!(learned, feat, :threads, 0.90)
        PPr.save_outer_route_state(learned, path)

        fresh = PPr.OuterRouteState()
        PPr.load_outer_route_state!(fresh, path)
        snap = PPr.outer_route_stats_snapshot(fresh, sig)
        # Four rounds' worth, not five: the cold-eviction branch replaces the
        # arm wholesale on its second observation, so the first is discarded.
        # Both counters still count all five, which is the point -- they record
        # how often the arm was observed, not how much of that survived into
        # the average.
        @test snap[:process].samples == 4 * 64

        restored = fresh.history[sig][:process]
        @test restored.campaigns == 5
        @test restored.observations == 5

        # The load-then-observe sequence that used to discard everything: with
        # observations restored at 5 the eviction branch (observations == 2)
        # cannot fire, so the restored weight survives.
        before_samples = restored.samples
        PPr.record_outer_route_feedback!(
            fresh, feat; route = :process, successes = 64, failures = 0,
            elapsed_success_s = 64 * 0.10, elapsed_success_sq_sum_s = 64 * 0.10^2,
        )
        PPr.record_outer_route_feedback!(
            fresh, feat; route = :process, successes = 64, failures = 0,
            elapsed_success_s = 64 * 0.10, elapsed_success_sq_sum_s = 64 * 0.10^2,
        )
        @test fresh.history[sig][:process].samples == before_samples + 2 * 64
        @test fresh.history[sig][:process].campaigns == 7
    end
end

@testset "Schema-2 files still load" begin
    # The bump to schema 3 records what a file contains; it must not reject one
    # written before the counters existed.
    mktempdir() do dir
        path = joinpath(dir, "legacy.toml")
        feat = _feat()
        sig = PPr.outer_route_signature(feat)
        legacy = Dict{String, Any}(
            "schema_version" => 2,
            "history" => [Dict{String, Any}(
                "signature" => sig,
                "route" => "process",
                "stats" => Dict{String, Any}(
                    "samples" => 64,
                    "successes" => 64,
                    "failures" => 0,
                    "elapsed_sum_s" => 6.4,
                    "elapsed_sq_sum_s" => 0.64,
                ),
            )],
        )
        open(path, "w") do io
            TOML.print(io, legacy)
        end

        st = PPr.OuterRouteState()
        result = PPr.load_outer_route_state!(st, path)
        @test result.rows == 1
        stats = st.history[sig][:process]
        @test stats.samples == 64
        # One, not zero: a row with samples > 0 ran at least once, and nothing
        # in a schema-2 file says how many times.
        @test stats.campaigns == 1
        @test stats.observations == 1
    end
end
