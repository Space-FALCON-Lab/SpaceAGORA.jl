using Test
using SpaceAGORA

# The bookkeeping the R7 predictive campaign runner keeps around a dispatch:
# the guard's running per-class observation, the mean work that feeds the
# sample-time correction, which pool workers count as warm, folding a finished
# campaign into the corrections file, and the trace a campaign prints when
# SPACEAGORA_CAMPAIGN_DISPATCH_TRACE is on.
#
# None of it needs a process pool: the guard and the fold are fed by hand, and
# the traced campaign runs with a one-worker tuning, so the process route is
# never offered (see predictive_campaign_tests.jl).

const BkCamp = SpaceAGORA.SimulationCampaigns
const BkPP = SpaceAGORA.ParallelProfiles

_bk_sample(index; elapsed_s, success = true, finished_ns = NaN) =
    BkCamp.MonteCarloSampleResult(index = index, seed = index, success = success,
                                  elapsed_s = elapsed_s, finished_ns = finished_ns)

@testset "the guard averages work and occupancy per class and skips each worker's first sample" begin
    g = BkCamp._PredictiveGuardState(2, 1)
    @test length(g.worker_seen) == 2 && !any(g.worker_seen)
    @test length(g.local_seen) == 1 && !any(g.local_seen)
    @test !g.decided[]
    @test g.completed == 0 && g.verdict === nothing
    # A negative width is an empty class, not an error.
    @test isempty(BkCamp._PredictiveGuardState(-1, 0).worker_seen)

    # Samples were taken at 1 s on the coordinator clock (ns) and landed later.
    taken = fill(1.0e9, 4)
    # Worker 1, first sample: counted in work and occupancy, not in steady.
    BkCamp._predictive_guard_observe!(g, _bk_sample(1; elapsed_s = 0.2, finished_ns = 1.5e9),
                                      :worker, 1, taken)
    @test g.worker_seen == [true, false]
    @test g.worker_n == 1 && g.worker_steady_n == 0
    @test g.worker_occupancy_s ≈ 0.5
    # Worker 1 again: now it is a steady sample.
    BkCamp._predictive_guard_observe!(g, _bk_sample(2; elapsed_s = 0.4, finished_ns = 1.3e9),
                                      :worker, 1, taken)
    @test g.worker_steady_n == 1
    @test g.worker_steady_occupancy_s ≈ 0.3
    # A local slot that failed: its work counts, and so does the failure.
    BkCamp._predictive_guard_observe!(g, _bk_sample(3; elapsed_s = 0.1, success = false,
                                                    finished_ns = 1.25e9), :local, 1, taken)
    @test g.local_seen == [true]
    @test g.failures == 1
    @test g.completed == 3

    @test BkCamp._predictive_guard_mean(g.worker_work_s, g.worker_n) ≈ 0.3
    @test BkCamp._predictive_guard_mean(g.local_occupancy_s, g.local_n) ≈ 0.25
    # Nothing landed in a class: its mean is unknown, not zero.
    @test isnan(BkCamp._predictive_guard_mean(0.0, 0))
end

@testset "mean pool work is over the successful samples only" begin
    samples = [_bk_sample(1; elapsed_s = 1.0), _bk_sample(2; elapsed_s = 3.0),
               _bk_sample(3; elapsed_s = 100.0, success = false)]
    @test BkCamp._predictive_mean_work(samples) ≈ 2.0
    @test isnan(BkCamp._predictive_mean_work([_bk_sample(1; elapsed_s = 1.0, success = false)]))
    @test isnan(BkCamp._predictive_mean_work(BkCamp.MonteCarloSampleResult[]))
end

@testset "a pool smaller than the plan's worker count is cold" begin
    pool = SpaceAGORA.campaign_process_pool()
    # This process has spawned no pool workers, so any plan that needs them is
    # cold, and marking the (empty) pool warm cannot change that.
    if isempty(pool.workers)
        @test BkCamp._predictive_pool_cold(2)
        @test BkCamp._predictive_mark_pool_warm!() === nothing
        @test BkCamp._predictive_pool_cold(2)
        # A plan with no pool has nothing to warm up.
        @test !BkCamp._predictive_pool_cold(0)
    end
end

@testset "folding a campaign updates the corrections and saves them only in :on mode" begin
    dir = mktempdir()
    path = joinpath(dir, "campaign_corrections_fold.toml")
    rules = BkCamp.CampaignCorrectionRules(step = 0.05, prior_share = 0.25, stale_campaigns = 20)
    consts = BkCamp.PredictiveCampaignConstants(round_tail = 0.5)
    fold!(c; trace) = BkCamp._predictive_fold_and_save!(c, rules, consts;
        signature = "sig", shape_key = "sig|n=16", final_plan = "process@w4+l2",
        worker_sample_s = 0.8, tail_observed = 0.4, trace = trace)

    withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS" => "on",
            "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => path) do
        c = BkCamp.CampaignCorrections("fp", "tok")
        @test fold!(c; trace = false) === nothing
        @test c.campaigns == 1
        @test c.last_plan["sig|n=16"] == "process@w4+l2"
        @test haskey(c.sample_time_s, "sig")
        @test c.round_tail !== nothing
        @test isfile(path)
        back = BkCamp.load_campaign_corrections(path; fingerprint = "fp", code_token = "tok")
        @test back.campaigns == 1
        @test back.last_plan == c.last_plan

        # The trace line reports the folded state and the mode.
        out = mktemp() do tmp, io
            redirect_stdout(io) do
                fold!(c; trace = true)
            end
            close(io)
            read(tmp, String)
        end
        @test occursin("[predictive] corrections campaigns=2", out)
        @test occursin("mode=on", out)
        @test BkCamp.load_campaign_corrections(path; fingerprint = "fp",
                                               code_token = "tok").campaigns == 2
    end

    read_path = joinpath(dir, "campaign_corrections_read.toml")
    withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS" => "read",
            "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => read_path) do
        c = BkCamp.CampaignCorrections("fp", "tok")
        fold!(c; trace = false)
        # Read mode still learns in memory but never writes the file.
        @test c.campaigns == 1
        @test !isfile(read_path)
    end
end

@testset "a traced predictive campaign prints its shape, its plan and its dispatch" begin
    seeds = collect(1:8)
    features = BkCamp.campaign_route_features(samples = length(seeds), n_sats = 1,
                                              density_family = "exponential", mission_time_s = 3600.0)
    tuning = BkPP.OuterRouteTuning(process_max_workers = 1, process_workers_resident = 0,
                                   memory_aware = false, mixed_dispatch = false)
    result = nothing
    out = mktemp() do tmp, io
        redirect_stdout(io) do
            withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "predictive",
                    "SPACEAGORA_CAMPAIGN_DISPATCH_TRACE" => "1",
                    "SPACEAGORA_PREDICTIVE_INNER_CURVE" => "0") do
                result = BkCamp.run_monte_carlo(x -> x + 1, seeds; threads = :auto,
                                                route_features = features,
                                                route_state = BkPP.OuterRouteState(),
                                                route_tuning = tuning)
            end
        end
        close(io)
        read(tmp, String)
    end
    @test [s.value for s in result.samples] == seeds .+ 1
    @test occursin("[predictive] shape n=8 threads=$(Base.Threads.nthreads()) pool=0", out)
    @test occursin("[predictive] terms heap_scale=", out)
    @test occursin("inner_curve=none", out)
    @test occursin("[predictive]   candidate ", out)
    @test occursin("[predictive] chosen ", out)
    # No pool means no local slots to guard: the campaign runs unguarded.
    @test occursin("n=8 failures=0 unguarded", out)
end
