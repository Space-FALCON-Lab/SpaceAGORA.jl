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

    @testset "adaptive campaign selects :process on a forced-large machine" begin
        tuning = PP.OuterRouteTuning(process_max_workers=2)
        state = PP.OuterRouteState()
        features = SC.campaign_route_features(samples=16, density_family="exponential", mission_time_s=1000.0)
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
