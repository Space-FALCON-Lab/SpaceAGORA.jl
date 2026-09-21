using Test
using SpaceAGORA

# The R7 predictive campaign planner, as pure logic: no simulations, no
# dispatch, no machine of any particular size. Every plan-space case is driven
# by passing the shape explicitly, so the P3-P5 points below are checked on
# this machine exactly as they would be on the TRX50 they were measured on.
#
# The TRX50 numbers quoted in the comments are medians in seconds from the
# cold-11 run (`paper_benchmarks_trx50_cold11`). They are used here as the
# expected DIRECTION of a decision -- which plan must win, which must not be
# chosen -- and never as a target the model is fitted to.

const SCamp = SpaceAGORA.SimulationCampaigns
const PCost = SpaceAGORA.SimulationModel.ParallelCost

_curve() = PCost.RateCurve([0.0, 10.0], [1.0, 1.0])

# A MachineConstants carrying nothing but the two USL terms the planner reads.
# The other fields are required by the struct and are not consulted here.
function _constants(; alpha = 0.05, beta_alloc = 0.004)
    return PCost.MachineConstants(
        simd_lane = _curve(), coeff_touch = _curve(), parallel_speedup = _curve(),
        ns_per_scalar_item = 1.0, ns_per_queue_node = 1.0,
        dispatch_pool_ns_base = 1.0, dispatch_pool_ns_per_worker = 1.0,
        dispatch_batch_ns_base = 1.0, dispatch_batch_ns_per_worker = 1.0,
        ns_per_atomic = 1.0, reference_fma_ns = 1.0, reference_mem_ns = 1.0,
        usl_alpha_base = alpha, usl_beta_alloc = beta_alloc,
    )
end

_cfg(; margin = 0.15, guard_factor = 1.5, local_slots_max = 64) =
    SCamp.PredictivePlannerConfig(margin = margin, guard_factor = guard_factor,
                                  local_slots_max = local_slots_max)

# Shorthand for "the plan the planner chose for this shape", with the machine
# supplied rather than inherited.
function _chosen(; n, threads, pool, local_cap = max(0, threads - 1),
                 constants = nothing, config = _cfg(local_slots_max = max(0, threads - 1)),
                 threads_candidate = threads > 1)
    return SCamp.predictive_plan(
        n_samples = n, threads = threads, process_workers = pool,
        threads_candidate = threads_candidate, local_slots_cap = local_cap,
        constants = constants, config = config).chosen
end

@testset "campaign_planner_mode is the switch, and bandit is the default" begin
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => nothing) do
        @test SCamp.campaign_planner_mode() === :bandit
    end
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "") do
        @test SCamp.campaign_planner_mode() === :bandit
    end
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "bandit") do
        @test SCamp.campaign_planner_mode() === :bandit
    end
    for value in ("predictive", "PREDICTIVE", " Predictive ")
        withenv("SPACEAGORA_CAMPAIGN_PLANNER" => value) do
            @test SCamp.campaign_planner_mode() === :predictive
        end
    end
    # An unrecognized value must not fall back to a planner the caller did not
    # ask for: a typo in a profile's env would otherwise run the other one.
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "r7") do
        err = try
            SCamp.campaign_planner_mode(); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("bandit", sprint(showerror, err))
        @test occursin("predictive", sprint(showerror, err))
    end
end

@testset "PredictivePlannerConfig reads env, keywords override, bad values throw" begin
    withenv("SPACEAGORA_PREDICTIVE_MARGIN" => nothing,
            "SPACEAGORA_PREDICTIVE_GUARD_FACTOR" => nothing,
            "SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX" => nothing) do
        c = SCamp.PredictivePlannerConfig()
        @test c.margin == 0.15
        @test c.guard_factor == 1.5
        @test c.local_slots_max == max(0, Base.Threads.nthreads() - 1)
    end
    withenv("SPACEAGORA_PREDICTIVE_MARGIN" => "0.4",
            "SPACEAGORA_PREDICTIVE_GUARD_FACTOR" => "2.5",
            "SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX" => "3") do
        c = SCamp.PredictivePlannerConfig()
        @test c.margin == 0.4
        @test c.guard_factor == 2.5
        @test c.local_slots_max == 3
        @test SCamp.PredictivePlannerConfig(margin = 0.05).margin == 0.05
    end
    withenv("SPACEAGORA_PREDICTIVE_MARGIN" => "soon") do
        @test_throws ArgumentError SCamp.PredictivePlannerConfig()
    end
    @test_throws ArgumentError SCamp.PredictivePlannerConfig(margin = -0.1)
    @test_throws ArgumentError SCamp.PredictivePlannerConfig(guard_factor = 0.5)
    @test_throws ArgumentError SCamp.PredictivePlannerConfig(local_slots_max = -1)
end

@testset "makespan is a greedy list schedule in units of one sample" begin
    # Equal consumers: the round count.
    @test SCamp.predictive_makespan(0, [1.0]) == 0.0
    @test SCamp.predictive_makespan(8, [1.0]) == 8.0
    @test SCamp.predictive_makespan(8, fill(1.0, 4)) == 2.0
    @test SCamp.predictive_makespan(9, fill(1.0, 4)) == 3.0
    @test SCamp.predictive_makespan(8, fill(2.0, 4)) == 4.0
    # No consumer can finish nothing.
    @test SCamp.predictive_makespan(4, Float64[]) == Inf
    # Heterogeneous: one fast consumer and one three times slower. Eight
    # samples go 1,2,3 to the fast one before the slow one finishes its first,
    # which is what the one-queue dispatchers do.
    @test SCamp.predictive_makespan(1, [1.0, 3.0]) == 1.0
    @test SCamp.predictive_makespan(2, [1.0, 3.0]) == 3.0
    @test SCamp.predictive_makespan(4, [1.0, 3.0]) == 3.0
    @test SCamp.predictive_makespan(8, [1.0, 3.0]) == 6.0
    # A slower consumer can never raise the makespan above what the fast ones
    # alone would need.
    @test SCamp.predictive_makespan(12, [1.0, 1.0, 5.0]) <= SCamp.predictive_makespan(12, [1.0, 1.0])
end

@testset "heap slowdown is 1 without constants and grows with width with them" begin
    @test SCamp.predictive_heap_slowdown(nothing, 1) == 1.0
    @test SCamp.predictive_heap_slowdown(nothing, 32) == 1.0
    mc = _constants()
    @test SCamp.predictive_heap_slowdown(mc, 1) == 1.0
    s2 = SCamp.predictive_heap_slowdown(mc, 2)
    s8 = SCamp.predictive_heap_slowdown(mc, 8)
    s32 = SCamp.predictive_heap_slowdown(mc, 32)
    @test 1.0 <= s2 <= s8 <= s32
    @test s32 > 1.0
    # Zeroed USL terms are an uncalibrated machine in all but name.
    @test SCamp.predictive_heap_slowdown(_constants(alpha = 0.0, beta_alloc = 0.0), 16) == 1.0
end

@testset "the plan space is the three routes and nothing narrower" begin
    # One sample: serial, with the whole budget.
    p = _chosen(n = 1, threads = 8, pool = 8)
    @test p.route === :none && p.consumers == 1 && p.static_equivalent
    # A single-threaded coordinator with no pool has only the serial plan.
    @test _chosen(n = 16, threads = 1, pool = 0).route === :none
    # Threads route is min(n, T), never a narrower rung.
    @test _chosen(n = 4, threads = 8, pool = 0).workers == 4
    @test _chosen(n = 64, threads = 8, pool = 0).workers == 8
    # Every enumerated process plan is the same pool with a different L.
    planning = SCamp.predictive_plan(
        n_samples = 64, threads = 8, process_workers = 8, threads_candidate = true,
        local_slots_cap = 7, constants = nothing, config = _cfg(local_slots_max = 7))
    process_plans = [p for p in planning.plans if p.route === :process]
    @test all(p -> p.workers == 8, process_plans)
    @test sort([p.local_slots for p in process_plans]) == collect(0:7)
    @test count(p -> p.static_equivalent, process_plans) == 1
    @test all(p -> p.consumers == p.workers + p.local_slots, process_plans)
    # Every concurrent plan declares an inner budget of one thread.
    @test all(p -> p.consumers <= 1 || p.inner_thread_budget == 1, planning.plans)
    # The local-slot cap is the tightest of the three bounds.
    capped = SCamp.predictive_plan(
        n_samples = 64, threads = 8, process_workers = 8, threads_candidate = true,
        local_slots_cap = 2, constants = nothing, config = _cfg(local_slots_max = 7))
    @test maximum(p.local_slots for p in capped.plans if p.route === :process) == 2
    knob = SCamp.predictive_plan(
        n_samples = 64, threads = 8, process_workers = 8, threads_candidate = true,
        local_slots_cap = 7, constants = nothing, config = _cfg(local_slots_max = 1))
    @test maximum(p.local_slots for p in knob.plans if p.route === :process) == 1
    # n - W_p bounds it too: there is no sample left for a local slot to run.
    tight = SCamp.predictive_plan(
        n_samples = 8, threads = 8, process_workers = 8, threads_candidate = true,
        local_slots_cap = 7, constants = nothing, config = _cfg(local_slots_max = 7))
    @test all(p -> p.local_slots == 0, tight.plans)
end

@testset "the static-equivalent plan is the default and the margin is the toll" begin
    # 64 samples, 8 workers, 7 local slots: 8 rounds against 5. A 37% predicted
    # gain clears any sane margin.
    wide = SCamp.predictive_plan(
        n_samples = 64, threads = 8, process_workers = 8, threads_candidate = true,
        local_slots_cap = 7, constants = nothing, config = _cfg(margin = 0.15))
    @test wide.reason === :predicted_gain
    @test wide.chosen.route === :process && wide.chosen.local_slots > 0
    @test wide.gain > 0.15
    # The same shape with a margin larger than the gain returns the static plan.
    held = SCamp.predictive_plan(
        n_samples = 64, threads = 8, process_workers = 8, threads_candidate = true,
        local_slots_cap = 7, constants = nothing, config = _cfg(margin = 0.9))
    @test held.reason === :margin_not_met
    @test held.chosen.static_equivalent
    @test held.chosen.local_slots == 0
    # The deviation taken is the SMALLEST one that buys the predicted gain:
    # among plans that tie on makespan the ranking prefers fewer local slots,
    # so the planner ends as close to the static plan as the gain allows. 17
    # samples over 8 workers is three rounds; one local slot already makes it
    # two, and seven would not make it one.
    @test _chosen(n = 17, threads = 8, pool = 8, local_cap = 7,
                  config = _cfg(local_slots_max = 7)).local_slots == 1
    # 64 samples over 8 workers: five rounds needs 13 consumers, so five slots.
    @test wide.chosen.local_slots == 5
end

@testset "among static-equivalent plans the pool is the tie-break, not an override" begin
    # 32 samples, 8 workers, 8 threads, no local slots allowed: both static
    # routes need four rounds, and the pool takes the tie.
    tie = _chosen(n = 32, threads = 8, pool = 8, local_cap = 0,
                  config = _cfg(local_slots_max = 0))
    @test tie.route === :process && tie.local_slots == 0
    # A pool too small to keep up loses outright: 2 workers against 8 threads
    # is 8 rounds against 2, far past any tie.
    @test _chosen(n = 16, threads = 8, pool = 2, local_cap = 0,
                  config = _cfg(local_slots_max = 0)).route === :threads
end

# ── The TRX50 P3-P5 points, as directions ────────────────────────────────────

@testset "P3 independent_1sat_1hr n=256 at (W_p=2, T=2): mixed 2+1" begin
    # Measured: R6 mixed 3.576, outer_threads 4.735, outer_process 5.072.
    # Three consumers against two is a third off the makespan, which clears the
    # margin; the planner must deviate to mixed.
    p = _chosen(n = 256, threads = 2, pool = 2, local_cap = 1,
                config = _cfg(local_slots_max = 1))
    @test p.route === :process
    @test p.workers == 2 && p.local_slots == 1
    # And with constants, where a single local slot is uncontended by
    # definition (s_heap(1) = 1), the same answer.
    q = _chosen(n = 256, threads = 2, pool = 2, local_cap = 1,
                constants = _constants(), config = _cfg(local_slots_max = 1))
    @test q.route === :process && q.local_slots == 1
end

@testset "P3 at (32, 32): mixed or capped mixed, never the threads route" begin
    # Measured: R6 mixed 32+31 1.069, outer_process 1.124, outer_threads 2.743.
    for constants in (nothing, _constants())
        p = _chosen(n = 256, threads = 32, pool = 32, local_cap = 31,
                    constants = constants, config = _cfg(local_slots_max = 31))
        @test p.route === :process
        @test p.workers == 32
    end
    # Uncalibrated: five rounds against eight, bought with the fewest slots
    # that reach five (52 consumers, so 20 of them).
    uncal = _chosen(n = 256, threads = 32, pool = 32, local_cap = 31,
                    config = _cfg(local_slots_max = 31))
    @test uncal.local_slots == 20
    @test uncal.makespan == 5.0
end

@testset "P4 montecarlo_heavy_aerobraking n=32 at (32, 32): process at L=0" begin
    # Measured: outer_process 1.678, outer_threads 1.856, R6 1.924. n - W_p = 0
    # leaves no sample for a local slot, so the plan space is the two static
    # routes and the pool takes the one-round tie.
    p = _chosen(n = 32, threads = 32, pool = 32, local_cap = 31,
                config = _cfg(local_slots_max = 31))
    @test p.route === :process
    @test p.local_slots == 0 && p.static_equivalent
end

@testset "P4 at (8, 8): mixed 8+7" begin
    # Measured: R6 mixed 2.757, outer_process 3.781. 32 samples over 15
    # consumers is three rounds against four.
    p = _chosen(n = 32, threads = 8, pool = 8, local_cap = 7,
                config = _cfg(local_slots_max = 7))
    @test p.route === :process
    @test p.local_slots > 0
    @test !p.static_equivalent
    @test p.makespan == 3.0          # against four for either static route
end

@testset "P5 mcgrid_8sat_16mc n=16 at (W_p=4, T=8): the static-equivalent plan" begin
    # Measured: outer_threads 1.384, R6 mixed 4+7 1.500. Mixed at 4+7 and the
    # threads route at W=8 both take two rounds, so the deviation buys nothing
    # the margin will pay for and the static plan stands.
    planning = SCamp.predictive_plan(
        n_samples = 16, threads = 8, process_workers = 4, threads_candidate = true,
        local_slots_cap = 7, constants = nothing, config = _cfg(local_slots_max = 7))
    @test planning.chosen.static_equivalent
    @test planning.chosen.route === :threads && planning.chosen.workers == 8
    @test planning.chosen.inner_thread_budget == 1
    # Static-equivalent and best outright, so the margin is never consulted.
    @test planning.reason === :static_equivalent_best
    @test planning.chosen.makespan == 2.0
    @test minimum(p.makespan for p in planning.plans) == 2.0
end

@testset "P5 with no pool: the threads route at min(n, T), budget 1" begin
    # mcgrid_16sat_8mc n=8 at (1, 32): measured outer_threads 1.896, R6 2.104.
    p = _chosen(n = 8, threads = 32, pool = 1, local_cap = 0,
                config = _cfg(local_slots_max = 31))
    @test p.route === :threads && p.workers == 8 && p.inner_thread_budget == 1
    # mcgrid_8sat_16mc n=16 at (1, 32): measured outer_threads 1.016; R6 threw
    # in 10 of 11 campaigns.
    q = _chosen(n = 16, threads = 32, pool = 1, local_cap = 0,
                config = _cfg(local_slots_max = 31))
    @test q.route === :threads && q.workers == 16 && q.inner_thread_budget == 1
end

@testset "P5 at (32, 1): the pool, with no local slots to offer" begin
    # Measured: outer_process 0.728 against 1.467. A one-thread coordinator has
    # no slot to spare and no threads route to take.
    p = _chosen(n = 16, threads = 1, pool = 32, local_cap = 0,
                threads_candidate = false, config = _cfg(local_slots_max = 0))
    @test p.route === :process
    @test p.local_slots == 0 && p.static_equivalent
end

@testset "the guard only ever reduces the local slots" begin
    cfg = _cfg(guard_factor = 1.5)
    plan = _chosen(n = 64, threads = 8, pool = 8, local_cap = 7,
                   config = _cfg(local_slots_max = 7))
    @test plan.local_slots > 1
    # The observation matches the prediction: nothing changes.
    v = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.10, local_mean_s = 0.10 * plan.heap_slowdown, failures = 0)
    @test !v.replan && v.local_slots == plan.local_slots
    @test v.reason === :observation_matches
    @test isapprox(v.ratio, 1.0; atol = 1e-9)
    # Local slots twice as slow as predicted: halve them.
    v2 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.10, local_mean_s = 0.20 * plan.heap_slowdown, failures = 0)
    @test v2.replan && v2.local_slots == plan.local_slots ÷ 2
    @test v2.reason === :local_slots_slower
    # Far past the factor: drop them entirely, which is the static plan.
    v3 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.10, local_mean_s = 0.50 * plan.heap_slowdown, failures = 0)
    @test v3.replan && v3.local_slots == 0
    @test v3.reason === :local_slots_far_slower
    # A failure is the strongest signal there is.
    v4 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.10, local_mean_s = 0.10, failures = 1)
    @test v4.replan && v4.local_slots == 0 && v4.reason === :sample_failure
    # Faster than predicted is never a reason to widen.
    v5 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.10, local_mean_s = 0.01, failures = 0)
    @test !v5.replan && v5.local_slots == plan.local_slots
    # An unobservable class leaves the plan alone rather than guessing.
    v6 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = NaN, local_mean_s = 0.1, failures = 0)
    @test !v6.replan && v6.reason === :not_observed
    # A static plan has nothing the guard can reduce.
    static = _chosen(n = 32, threads = 32, pool = 32, local_cap = 31,
                     config = _cfg(local_slots_max = 31))
    v7 = SCamp.predictive_guard_verdict(static, cfg;
        worker_mean_s = 0.1, local_mean_s = 10.0, failures = 3)
    @test !v7.replan && v7.reason === :nothing_to_reduce
end

@testset "a re-plan keeps the route and re-prices the remainder" begin
    plan = _chosen(n = 64, threads = 8, pool = 8, local_cap = 7,
                   config = _cfg(local_slots_max = 7))
    @test plan.local_slots >= 3
    replanned = SCamp.predictive_replan(plan, 3, 56, nothing)
    @test replanned.route === plan.route
    @test replanned.workers == plan.workers
    @test replanned.local_slots == 3
    @test replanned.consumers == plan.workers + 3
    @test !replanned.static_equivalent
    # To zero is to the static plan.
    @test SCamp.predictive_replan(plan, 0, 56, nothing).static_equivalent
    # Never upward, whatever it is asked for.
    @test SCamp.predictive_replan(plan, 99, 56, nothing).local_slots == plan.local_slots
end

@testset "machine constants load once and are a legal absence" begin
    SCamp.reset_predictive_machine_constants!()
    withenv("SPACEAGORA_COST_CONSTANTS_PATH" => joinpath(mktempdir(), "no_such_file.toml")) do
        @test SCamp.predictive_machine_constants() === nothing
    end
    SCamp.reset_predictive_machine_constants!()
    # Whatever this machine has, asking twice gives the same object and the
    # planner reports whether it got one.
    first_load = SCamp.predictive_machine_constants()
    @test SCamp.predictive_machine_constants() === first_load
    planning = SCamp.predictive_plan(n_samples = 4, threads = 2, process_workers = 0,
                                     threads_candidate = true, local_slots_cap = 0,
                                     constants = nothing, config = _cfg())
    @test planning.constants_loaded == false
    planning2 = SCamp.predictive_plan(n_samples = 4, threads = 2, process_workers = 0,
                                      threads_candidate = true, local_slots_cap = 0,
                                      constants = _constants(), config = _cfg())
    @test planning2.constants_loaded == true
end
