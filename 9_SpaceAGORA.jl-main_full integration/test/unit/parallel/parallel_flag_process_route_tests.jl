using Test
using SpaceAGORA

# run_monte_carlo(...; parallel=true) sees only the sample function. Without
# route_features it cannot know whether that function (a closure defined in a
# script's Main, as here) exists in a worker process, so the flag must never
# send it there on its own; a caller that describes the workload gets the
# process route back as a candidate.

const PRF_SC = SpaceAGORA.SimulationCampaigns
const PRF_PC = SpaceAGORA.SimulationModel.ParallelCost
const PRF_TMP = mktempdir()

@testset "without route_features the flag withholds the process route" begin
    withheld = PRF_SC._monte_carlo_flag_tuning(nothing, nothing)
    @test withheld.process_max_workers == 1
    base = PRF_SC._campaign_route_tuning()
    for name in fieldnames(typeof(base))
        name === :process_max_workers && continue
        @test getfield(withheld, name) == getfield(base, name)
    end
    # A caller's tuning is kept apart from the withheld route.
    custom = SpaceAGORA.ParallelProfiles.OuterRouteTuning(mc_process_min_samples = 2, process_max_workers = 8)
    kept = PRF_SC._monte_carlo_flag_tuning(nothing, custom)
    @test kept.mc_process_min_samples == 2 && kept.process_max_workers == 1
    # Described workloads are planned with the caller's tuning, unchanged.
    features = PRF_SC.campaign_route_features(samples = 64, mission_time_s = 7200.0)
    @test PRF_SC._monte_carlo_flag_tuning(features, custom) === custom
    @test PRF_SC._monte_carlo_flag_tuning(features, nothing) === nothing
    # The planner's pool width under the withheld tuning is below a pool.
    @test SpaceAGORA.ParallelProfiles.effective_process_workers(
        PRF_SC.campaign_route_features(samples = 64), withheld) < 2
end

withenv("SPACEAGORA_COST_CONSTANTS_PATH" => joinpath(PRF_TMP, "constants.toml"),
        "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => joinpath(PRF_TMP, "policy_state.toml"),
        "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(PRF_TMP, "rhs_calibration.toml"),
        "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => joinpath(PRF_TMP, "outer_route_state.toml"),
        "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => joinpath(PRF_TMP, "campaign_corrections.toml"),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
    # Settle the constants path without measuring anything.
    @test_logs (:warn,) (:info,) match_mode = :any PRF_PC.ensure_machine_constants!(
        path = joinpath(PRF_TMP, "constants.toml"), calibrate = () -> error("no calibration in this test"))

    @testset "a script-defined sample function runs in this process" begin
        pool_before = length(SpaceAGORA.ParallelProcess.campaign_process_pool().workers)
        offset = 1000   # captured script state a worker process would not have
        result = run_monte_carlo(1:32; parallel = true) do seed
            seed + offset
        end
        @test result.route !== :process
        @test length(result.successful) == 32
        @test isempty(result.failed)
        @test sort([s.value for s in result.successful]) == collect(1001:1032)
        @test length(SpaceAGORA.ParallelProcess.campaign_process_pool().workers) == pool_before
    end
end
