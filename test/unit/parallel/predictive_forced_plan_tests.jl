using Test
using SpaceAGORA

# SPACEAGORA_PREDICTIVE_FORCE_PLAN, the measurement hook that runs a named
# campaign plan in place of the chosen one. Unset it is inert; set it must name
# a plan the host can run; and a forced campaign must leave the corrections file
# alone, or it would steer later unforced campaigns of the same shape.

const FPCamp = SpaceAGORA.SimulationCampaigns
const FPPP = SpaceAGORA.ParallelProfiles
const FP_VAR = "SPACEAGORA_PREDICTIVE_FORCE_PLAN"

# n = 64 samples, T = 8 threads, a pool of 8 workers, at most 7 local slots.
fp_cfg() = FPCamp.PredictivePlannerConfig(
    margin = 0.15, guard_factor = 1.5, local_slots_max = 7, remote_overhead = 0.0,
    route_switch = false, heap_model = :locals, local_thrash = 3.0)
fp_planning(; split = nothing, unsplit = split) = FPCamp.predictive_plan(
    n_samples = 64, threads = 8, process_workers = 8, threads_candidate = true,
    local_slots_cap = 7, constants = nothing, config = fp_cfg(),
    inner_curve = split, inner_curve_unsplit = unsplit)
fp_force(value, planning = fp_planning(); split = nothing, unsplit = nothing,
         allowed_routes = [:none, :threads, :process]) =
    withenv(FP_VAR => value) do
        FPCamp._predictive_forced_plan(planning, 64, fp_cfg(), nothing, FPCamp.PredictiveCostTerms(),
                                       split, unsplit; threads = 8, pool_workers = 8, local_cap = 7,
                                       allowed_routes = allowed_routes)
    end
fp_message(value) = try
    fp_force(value); ""
catch err
    err isa ArgumentError ? err.msg : rethrow()
end

@testset "unset or blank forces nothing" begin
    @test fp_force(nothing) === nothing
    @test fp_force("") === nothing
    @test fp_force("   ") === nothing
end

@testset "a malformed key throws" begin
    @test_throws ArgumentError fp_force("bogus")
    @test_throws ArgumentError fp_force("Process@w8+l1")
end

@testset "an infeasible key throws, naming the variable and the bound" begin
    @test occursin("affordable pool of 8", fp_message("process@w9+l0"))
    @test occursin(FP_VAR, fp_message("process@w9+l0"))
    @test occursin("workers must be >= 1", fp_message("process@w0+l0"))
    @test occursin("at least 2 workers", fp_message("process@w1+l0"))
    @test occursin("at least 2 workers", fp_message("process@w1+l1"))
    @test occursin("planner's bound of 7", fp_message("process@w8+l8"))
    @test occursin("b=16", fp_message("threads@w1+l0+b16"))
    @test occursin("whose bound is 0", fp_message("none@w1+l3"))
    @test occursin("runs 1 worker", fp_message("none@w5+l0"))
    @test occursin("workers * b = 16", fp_message("threads@w8+l0+b2"))
end

@testset "forced plans preserve workload route restrictions" begin
    features = FPCamp.campaign_route_features(samples = 64, density_family = "gram")
    tuning = FPPP.OuterRouteTuning(process_max_workers = 8, memory_aware = false)
    routes = FPPP.outer_route_candidates(features; tuning, threads_available = true)
    @test routes == [:none, :process]
    planning = FPCamp.predictive_plan(n_samples = 64, threads = 8, process_workers = 8,
        threads_candidate = false, local_slots_cap = 7, constants = nothing, config = fp_cfg())
    @test_throws ArgumentError fp_force("threads@w2+l0", planning; allowed_routes = routes)
    # Also enforce restrictions before returning a previously enumerated plan.
    @test_throws ArgumentError fp_force("threads@w8+l0"; allowed_routes = routes)
    @test fp_force("none@w1+l0", planning; allowed_routes = routes).route == :none
end

@testset "serial dispatch applies its inner ceiling and restores the environment" begin
    seeds = collect(1:8)
    budget = min(2, Threads.nthreads())
    plan = FPCamp._predictive_plan(:none, 1, 0, length(seeds), false, nothing, 0.0;
                                  inner_budget = budget)
    no_ceiling = FPCamp._predictive_plan(:none, 1, 0, length(seeds), false, nothing, 0.0)
    tuning = FPPP.OuterRouteTuning(memory_aware = false)
    measure(p) = FPCamp._predictive_dispatch(
        _ -> FPCamp.ParallelPolicy.effective_inner_thread_budget(), seeds, p, tuning; fail_fast = true)
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
        @test all(s -> s.value == budget, measure(plan).samples)
        @test !haskey(ENV, "SPACEAGORA_INNER_THREAD_BUDGET")
        @test all(s -> s.value == Threads.nthreads(), measure(no_ceiling).samples)
    end
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "1") do
        @test all(s -> s.value == 1, measure(plan).samples)
        @test all(s -> s.value == 1, measure(no_ceiling).samples)
        @test ENV["SPACEAGORA_INNER_THREAD_BUDGET"] == "1"
    end
end

@testset "a valid key returns the forced plan" begin
    planning = fp_planning()
    enumerated = only(filter(p -> FPCamp.predictive_plan_key(p) == "process@w8+l2", planning.plans))
    @test fp_force("process@w8+l2", planning) === enumerated
    @test fp_force("  process@w8+l2 ", planning) === enumerated
    # Not enumerated without a curve, but feasible: priced as an off-list plan.
    off = fp_force("threads@w4+l0+b2", planning)
    @test FPCamp.predictive_plan_key(off) == "threads@w4+l0+b2"
    @test !off.static_equivalent
    @test !any(p -> p === off, planning.plans)
end

@testset "the key is compared after parsing, not as a string" begin
    planning = fp_planning()
    threads = only(filter(p -> FPCamp.predictive_plan_key(p) == "threads@w8+l0", planning.plans))
    @test fp_force("threads@w8+l0+b1", planning) === threads
end

@testset "an off-list serial plan is priced from the unsplit curve" begin
    split = FPCamp.InnerSpeedupCurve([1.0, 1.1, 1.1, 1.2, 1.2, 1.2, 1.2, 1.2])
    unsplit = FPCamp.InnerSpeedupCurve([1.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0, 8.0])
    # Planned without the unsplit curve, so none@w1+l0+b8 is not enumerated.
    planning = fp_planning(split = split, unsplit = nothing)
    forced = fp_force("none@w1+l0+b8", planning; split = split, unsplit = unsplit)
    expected = FPCamp._predictive_plan(:none, 1, 0, 64, false, nothing, 0.0;
                                       terms = FPCamp.PredictiveCostTerms(),
                                       inner_budget = 8, inner_time = 1.0 / 8.0)
    @test forced.makespan == expected.makespan
end

@testset "a forced campaign neither folds nor saves the corrections" begin
    path = joinpath(mktempdir(), "campaign_corrections_forced.toml")
    seeds = collect(1:8)
    # A pool of two is affordable, so the campaign loads the corrections; the
    # forced serial plan then runs without ever spawning a worker.
    features = FPCamp.campaign_route_features(samples = length(seeds), n_sats = 1,
                                              density_family = "exponential", mission_time_s = 1.0e6)
    tuning = FPPP.OuterRouteTuning(process_max_workers = 2, process_workers_resident = 0,
                                   memory_aware = false, mixed_dispatch = false,
                                   mc_route_by_core_budget = true)
    result = nothing
    FPCamp.reset_campaign_corrections!()
    try
        out = mktemp() do tmp, io
            redirect_stdout(io) do
                withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "predictive",
                        "SPACEAGORA_CAMPAIGN_DISPATCH_TRACE" => "1",
                        "SPACEAGORA_PREDICTIVE_INNER_CURVE" => "0",
                        "SPACEAGORA_CAMPAIGN_CORRECTIONS" => "on",
                        "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => path,
                        FP_VAR => "none@w1+l0") do
                    result = FPCamp.run_monte_carlo(x -> x + 1, seeds; threads = :auto,
                                                    route_features = features,
                                                    route_state = FPPP.OuterRouteState(),
                                                    route_tuning = tuning)
                end
            end
            close(io)
            read(tmp, String)
        end
        @test [s.value for s in result.samples] == seeds .+ 1
        # The precondition: the shape had a pool, so corrections were in play.
        @test occursin("pool=2", out)
        @test occursin("corrections=on", out)
        @test occursin("[predictive] forced none@w1+l0", out)
        @test !occursin("[predictive] corrections campaigns=", out)
        @test !isfile(path)
        withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS" => "on",
                "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => path) do
            @test FPCamp.campaign_corrections().campaigns == 0
        end
    finally
        FPCamp.reset_campaign_corrections!()
    end
end
