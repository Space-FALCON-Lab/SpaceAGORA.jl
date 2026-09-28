using Test
using Distributed

const _CPR_REPO_ROOT = isdefined(Main, :REPO_ROOT) ? Main.REPO_ROOT : normpath(joinpath(@__DIR__, "..", ".."))

# `using SpaceAGORA` (a real Pkg-loaded root package), not a raw `include`:
# the :process route round-trips SpaceAGORA-defined types (MonteCarloSampleResult)
# through Distributed, and deserializing a value built on a worker (which is
# always Pkg-loaded via `using SpaceAGORA`, see `_bootstrap_process_worker!`)
# requires the *coordinator* side to resolve the same `SpaceAGORA` PkgId as a
# root module too -- a raw-included module was never registered that way and
# deserialization fails with `KeyError: ... PkgId(..., "SpaceAGORA") not found`.
using SpaceAGORA

const SC = SpaceAGORA.SimulationCampaigns
const PP = SpaceAGORA.ParallelProfiles

@testset "Campaign Process-Route Coverage Debt Probes" begin
    @testset "_run_campaign_with_route_env process dispatch" begin
        # `spec.seeds` is a `Vector` here (not a bare range) to match how
        # `_run_campaign_adaptive` always calls this: it `collect`s seeds
        # before building the spec, and `_run_monte_carlo_process` requires
        # a `Vector` (it indexes `seeds` while draining the shared job channel).
        spec = SC.MonteCarloSpec(seeds=collect(1:6), threads=2, fail_fast=false)
        plan = (route=:process, threads=2, inner_thread_budget=1, record=true)
        result = SC._run_campaign_with_route_env(seed -> seed^2, spec, plan)
        @test result isa SC.MonteCarloResult
        @test length(result.samples) == 6
        @test all(s -> s.success, result.samples)
        @test [s.value for s in result.samples] == [1, 4, 9, 16, 25, 36]
        @test result.threads == 2

        # worker_count collapses to 1 when seeds are scarcer than requested
        # threads, so this bypasses the process pool entirely.
        # threads=1, not e.g. 4: run_monte_carlo(f, spec)'s own thread-count
        # validation checks spec.threads against Threads.nthreads() even on
        # this fallback path, so an oversized spec.threads would throw here
        # regardless of the worker_count collapse this test is exercising.
        spec_one = SC.MonteCarloSpec(seeds=collect(1:1), threads=1, fail_fast=false)
        plan_one = (route=:process, threads=1, inner_thread_budget=1, record=true)
        result_one = SC._run_campaign_with_route_env(x -> x, spec_one, plan_one)
        @test result_one.threads == 1
        @test result_one.samples[1].value == 1

        # fail_fast surfaces the first sample failure on the process route too.
        spec_fail = SC.MonteCarloSpec(seeds=collect(1:4), threads=2, fail_fast=true)
        plan_fail = (route=:process, threads=2, inner_thread_budget=1, record=true)
        @test_throws ErrorException SC._run_campaign_with_route_env(
            seed -> (seed == 2 ? error("probe-injected sample failure") : seed), spec_fail, plan_fail
        )
    end

    @testset "adaptive campaign selects :threads for analytic-density Monte Carlo" begin
        # This testset used to be "selects :process on a forced-large machine",
        # and it asserted the behaviour that made R4/R5 lose every Monte Carlo
        # case in the benchmark catalog by +33% to +201% against the best static
        # route. The Monte Carlo route rule preferred :process on any medium or
        # large machine once the sample count or the SIMULATED mission time
        # cleared a threshold, neither of which is a proxy for per-sample
        # compute. Measured at 64 samples on 12 threads, process was 1.41x to
        # 2.90x SLOWER than threads across five Monte Carlo cases spanning
        # 0.038 s to 3.179 s of compute per sample -- no crossover anywhere in
        # that range.
        #
        # An analytic-density campaign now routes to threads.
        tuning = PP.OuterRouteTuning(process_max_workers=2)
        state = PP.OuterRouteState()
        features = SC.campaign_route_features(samples=16, density_family="exponential", mission_time_s=1000.0)
        result = withenv("SPACEAGORA_PERF_HARDWARE_CLASS" => "large") do
            SC._run_campaign_adaptive(seed -> seed, 1:16; fail_fast=false, features=features, state=state, tuning=tuning)
        end
        @test result isa SC.MonteCarloResult
        @test length(result.successful) == 16

        sig = PP.outer_route_signature(features)
        snap = PP.outer_route_stats_snapshot(state, sig)
        # Which non-process route is taken depends on the thread count this
        # suite happens to run under: with worker threads it is :threads, and on
        # a single-threaded run threads_available is false and the rule falls
        # through to :none. Both are correct and neither is :process, which is
        # the property under test. Asking default_outer_route for the expected
        # answer keeps this green at any thread count without weakening it.
        expected = PP.default_outer_route(
            features; tuning=tuning, machine_class=:large,
            threads_available=(Threads.nthreads() > 1), parallel_enabled=true
        )
        @test expected in (:threads, :none)
        @test haskey(snap, expected)
        @test snap[expected].samples == 16
        @test snap[expected].success_rate == 1.0
        @test !haskey(snap, :process)
        # The process route must remain a CANDIDATE even though it is no longer
        # the default, or the bandit could never rediscover it where it wins.
        @test :process in PP.outer_route_candidates(
            features; tuning=tuning, machine_class=:large,
            threads_available=true, parallel_enabled=true
        )
    end

    @testset "adaptive campaign still selects :process for native GRAM density" begin
        # Native GRAM is not thread-safe, so there the process route is a
        # correctness requirement rather than a performance preference, and
        # _is_native_gram_point_density decides it before the Monte Carlo rule
        # above is consulted. This keeps the process path covered end-to-end now
        # that analytic-density campaigns no longer exercise it.
        tuning = PP.OuterRouteTuning(process_max_workers=2)
        state = PP.OuterRouteState()
        features = SC.campaign_route_features(samples=16, density_family="gram_point", mission_time_s=1000.0)
        result = withenv("SPACEAGORA_PERF_HARDWARE_CLASS" => "large") do
            SC._run_campaign_adaptive(seed -> seed, 1:16; fail_fast=false, features=features, state=state, tuning=tuning)
        end
        @test result isa SC.MonteCarloResult
        @test length(result.successful) == 16
        @test result.threads <= 2

        sig = PP.outer_route_signature(features)
        snap = PP.outer_route_stats_snapshot(state, sig)
        @test haskey(snap, :process)
        @test snap[:process].samples == 16
        @test snap[:process].success_rate == 1.0
    end

    @testset "_record_campaign_route_feedback! short-circuits on empty result" begin
        state = PP.OuterRouteState()
        features = SC.campaign_route_features(samples=0)
        empty_result = SC.MonteCarloResult(SC.MonteCarloSampleResult[], 0.0, 0)
        @test SC._record_campaign_route_feedback!(state, features, :process, empty_result; tuning=PP.OuterRouteTuning()) === nothing
        @test isempty(state.history)
    end
end
