using Test
using SpaceAGORA

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
    @test PPr._route_is_proven(snap, :process, 2) == false

    # One sampled alternative that loses is enough to make it proven.
    _record!(only_default, feat, :threads, 0.90)
    snap2 = PPr.outer_route_stats_snapshot(only_default, PPr.outer_route_signature(feat))
    @test PPr._route_is_proven(snap2, :process, 2) == true
end
