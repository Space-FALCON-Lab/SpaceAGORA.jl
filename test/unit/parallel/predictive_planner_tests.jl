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

_cfg(; margin = 0.15, guard_factor = 1.5, local_slots_max = 64, remote_overhead = 0.0,
     route_switch = true, heap_model = :none) =
    SCamp.PredictivePlannerConfig(margin = margin, guard_factor = guard_factor,
                                  local_slots_max = local_slots_max,
                                  remote_overhead = remote_overhead,
                                  route_switch = route_switch,
                                  heap_model = heap_model)

# The TRX50's own calibration, from job 20260921-202331-989854 (fingerprint
# f48f5e44b83aeac4). SOURCED: these are the numbers the machine measured, and
# the two points below are what the planner did with them.
_trx50_constants() = _constants(alpha = 0.156, beta_alloc = 0.00605)

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
        @test c.remote_overhead == 0.0
        @test c.route_switch == false          # built, measured, shipped off
        @test c.heap_model === :none           # measured and refuted, see below
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
    @test_throws ArgumentError SCamp.PredictivePlannerConfig(remote_overhead = -0.1)
    withenv("SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD" => "0.6") do
        @test SCamp.PredictivePlannerConfig().remote_overhead == 0.6
    end
    withenv("SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH" => "on") do
        @test SCamp.PredictivePlannerConfig().route_switch
    end
    withenv("SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH" => "maybe") do
        @test_throws ArgumentError SCamp.PredictivePlannerConfig()
    end
    withenv("SPACEAGORA_PREDICTIVE_HEAP_MODEL" => "usl") do
        @test SCamp.PredictivePlannerConfig().heap_model === :usl
    end
    withenv("SPACEAGORA_PREDICTIVE_HEAP_MODEL" => "NONE") do
        @test SCamp.PredictivePlannerConfig().heap_model === :none
    end
    withenv("SPACEAGORA_PREDICTIVE_HEAP_MODEL" => "amdahl") do
        @test_throws ArgumentError SCamp.PredictivePlannerConfig()
    end
    @test_throws ArgumentError SCamp.PredictivePlannerConfig(heap_model = :usl2)
end

@testset "a declared remote overhead prices the pool against the threads route" begin
    # Measured on this repo's workstation, independent_1sat_1hr, 64 samples of
    # ~38 ms, 8 pool workers against 8 coordinator threads: outer_process
    # 0.722 s against outer_threads 0.427 s, i.e. a worker class roughly 1.6x a
    # thread rather than the 1.0x the default assumes. With the overhead left
    # at its default the planner prefers the pool here; declared, it does not.
    shape = (n_samples = 64, threads = 8, process_workers = 8, threads_candidate = true,
             local_slots_cap = 4, constants = nothing)
    default = SCamp.predictive_plan(; shape..., config = _cfg(local_slots_max = 7))
    @test default.chosen.route === :process
    declared = SCamp.predictive_plan(; shape...,
                                     config = _cfg(local_slots_max = 7, remote_overhead = 0.6))
    @test declared.chosen.route === :threads
    @test declared.chosen.workers == 8
    # The worker class is priced as declared, and a re-plan keeps that price.
    pool_plan = first(p for p in declared.plans if p.route === :process)
    @test pool_plan.worker_slowdown == 1.6
    @test SCamp.predictive_replan(pool_plan, 0, 32, nothing).worker_slowdown == 1.6
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

# ── The two TRX50 points the USL heap model got wrong ────────────────────────
#
# Cold 11-repeat run, job 20260921-202331-989854, tree 197c53c39, machine
# constants LOADED (usl_alpha_base = 0.156, usl_beta_alloc = 0.00605). R7
# passed 34 of 36 points and failed these two, both P5 mcgrid_8sat_16mc:
#
#   2 workers x 16 threads   chose mixed, ended w2+l10   1.524 s
#                            against pinned threads      0.957 s
#   4 workers x 8 threads    chose mixed, ended w4+l2    1.846 s
#                            against pinned threads      1.382 s
#
# The cause is the contention term, not the ranking: s_heap(16) = 4.8 and
# s_heap(8) = 3.0 at those constants, which prices the pinned-threads plan at
# several times one round and hands the ranking to anything with pool workers,
# however small the pool. The measured threads route is not several times a
# round. The alloc-kernel USL fit does not transfer to whole samples in
# magnitude -- the mapping the contract flagged as ASSUMED -- so the heap model
# is off by default and these points are checked with constants PRESENT.

@testset "P5 mcgrid_8sat_16mc at 2 workers x 16 threads: the pinned threads plan" begin
    planning = SCamp.predictive_plan(
        n_samples = 16, threads = 16, process_workers = 2, threads_candidate = true,
        local_slots_cap = 15, constants = _trx50_constants(),
        config = _cfg(local_slots_max = 15))
    # The constants are loaded and reported; the default heap model declines to
    # charge contention with them.
    @test planning.constants_loaded
    @test planning.chosen.route === :threads
    @test planning.chosen.workers == 16
    @test planning.chosen.local_slots == 0
    @test planning.chosen.static_equivalent
    @test planning.chosen.heap_slowdown == 1.0
    @test planning.reason === :static_equivalent_best
end

@testset "P5 mcgrid_8sat_16mc at 4 workers x 8 threads: the pinned threads plan" begin
    planning = SCamp.predictive_plan(
        n_samples = 16, threads = 8, process_workers = 4, threads_candidate = true,
        local_slots_cap = 7, constants = _trx50_constants(),
        config = _cfg(local_slots_max = 7))
    @test planning.constants_loaded
    @test planning.chosen.route === :threads
    @test planning.chosen.workers == 8
    @test planning.chosen.local_slots == 0
    @test planning.chosen.static_equivalent
    @test planning.reason === :static_equivalent_best
end

@testset "the usl heap model still produces its old ranking when asked for" begin
    # Same shape, same constants, the two models. `:none` ranks by round count
    # and takes the pinned threads plan; `:usl` prices sixteen tasks on one
    # heap at 4.8x a sample and hands the point to the pool. Both are still
    # reachable; only the default changed.
    shape = (n_samples = 16, threads = 16, process_workers = 8,
             threads_candidate = true, local_slots_cap = 15,
             constants = _trx50_constants())
    off = SCamp.predictive_plan(; shape..., config = _cfg(local_slots_max = 15))
    on  = SCamp.predictive_plan(; shape...,
                                config = _cfg(local_slots_max = 15, heap_model = :usl))
    @test off.chosen.route === :threads && off.chosen.workers == 16
    @test on.chosen.route === :process
    @test on.constants_loaded && off.constants_loaded
    # The term itself is unchanged and still steep at these widths; what
    # changed is whether the planner is handed it.
    mc = _trx50_constants()
    @test SCamp.predictive_heap_slowdown(mc, 16) > 4.0
    @test SCamp.predictive_heap_slowdown(mc, 8) > 2.0
    threads_off = first(p for p in off.plans if p.route === :threads)
    threads_on = first(p for p in on.plans if p.route === :threads)
    @test threads_off.heap_slowdown == 1.0
    @test threads_on.heap_slowdown == SCamp.predictive_heap_slowdown(mc, 16)
    @test threads_on.makespan > 4.0 * threads_off.makespan
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

# ── The guard's second direction: which side of the machine the campaign is on ──
#
# Round one measures each class twice: the work inside the sample (`mean_s`)
# and what the consumer occupied to deliver it (`occupancy_s`, work plus the
# round trip and serialization around it). The first ratio says whether the
# local slots cost what the contention model charged them; the second says
# whether the pool is paying for itself at all. Both verdicts end on a
# static-equivalent plan.
#
# The occupancy numbers below are from this repo's workstation,
# independent_1sat_1hr at 64 samples of ~38 ms over 8 pool workers and 3 local
# slots, and from the TRX50 cold-11 ratios quoted for W <= 8. They are used as
# the expected DIRECTION of the decision, never as a fitted target.

# A mixed plan with a known shape, priced the way the planner would price it
# on an uncalibrated machine.
_guard_plan(; local_slots = 3, workers = 8, n = 64, remote_overhead = 0.0) =
    SCamp._predictive_plan(:process, workers, local_slots, n, local_slots == 0,
                           nothing, remote_overhead)

@testset "workers occupying far more than local slots move the remainder to threads" begin
    cfg = _cfg(guard_factor = 1.5)
    plan = _guard_plan()
    @test plan.route === :process && plan.local_slots == 3 && plan.workers == 8
    v = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.038, local_mean_s = 0.046,
        worker_occupancy_s = 0.333, local_occupancy_s = 0.050,
        remaining = 53, threads = 8, threads_candidate = true, constants = nothing)
    @test v.replan
    @test v.route === :threads
    @test v.workers == 8                 # the width budget it was handed
    @test v.local_slots == 0
    @test v.reason === :workers_occupying_more_than_threads
    @test v.occupancy_ratio > cfg.guard_factor
    @test v.threads_s < v.continue_s
    # The work-only ratio alone would have said "nothing is wrong": 0.046/0.038
    # is inside the factor. That is the point of the second measurement.
    @test v.ratio < cfg.guard_factor
    # And the re-plan the verdict asks for is the pinned threads plan.
    replanned = SCamp.predictive_replan(plan, v.local_slots, 53, nothing;
                                        route = v.route, workers = v.workers)
    @test replanned.route === :threads
    @test replanned.workers == 8 && replanned.local_slots == 0
    @test replanned.consumers == 8
    @test replanned.static_equivalent
    @test replanned.inner_thread_budget == 1
end

@testset "a slow-looking pool that would still finish sooner is left alone" begin
    # Ratio past the factor, but only two threads to move to: the pool's eight
    # workers finish the remainder first even at twice a local slot's cost.
    cfg = _cfg(guard_factor = 1.5)
    plan = _guard_plan()
    v = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.09, local_mean_s = 0.045,
        worker_occupancy_s = 0.10, local_occupancy_s = 0.05,
        remaining = 53, threads = 2, threads_candidate = true, constants = nothing)
    @test !v.replan
    @test v.route === plan.route && v.local_slots == plan.local_slots
    @test v.reason === :threads_no_better
    @test v.threads_s > v.continue_s
end

@testset "the TRX50 W<=8 occupancy band is the no-change band" begin
    # 1.07 to 1.11 there; inside the factor, so the guard does nothing and the
    # campaign stays on the plan it was given.
    cfg = _cfg(guard_factor = 1.5)
    plan = _guard_plan()
    for ratio in (1.0, 1.07, 1.11, 1.49)
        v = SCamp.predictive_guard_verdict(plan, cfg;
            worker_mean_s = 0.10, local_mean_s = 0.10,
            worker_occupancy_s = 0.10 * ratio, local_occupancy_s = 0.10,
            remaining = 53, threads = 8, threads_candidate = true, constants = nothing)
        @test !v.replan
        @test v.route === :process && v.local_slots == plan.local_slots
        @test v.reason === :observation_matches
    end
end

@testset "the two directions coexist, and the route change wins a tie" begin
    cfg = _cfg(guard_factor = 1.5)
    plan = _guard_plan()
    # Healthy occupancies, slow local work: direction one, as before.
    v1 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.05, local_mean_s = 0.20,
        worker_occupancy_s = 0.05, local_occupancy_s = 0.05,
        remaining = 53, threads = 8, threads_candidate = true, constants = nothing)
    @test v1.replan && v1.route === :process && v1.local_slots == 0
    @test v1.reason === :local_slots_far_slower
    # Both qualify: being on the wrong side of the machine outranks holding too
    # many local slots on the right one.
    v2 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.05, local_mean_s = 0.20,
        worker_occupancy_s = 0.60, local_occupancy_s = 0.10,
        remaining = 53, threads = 8, threads_candidate = true, constants = nothing)
    @test v2.replan && v2.route === :threads
    @test v2.reason === :workers_occupying_more_than_threads
    # A failure still outranks both, and still lands on :process@L=0.
    v3 = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.05, local_mean_s = 0.05,
        worker_occupancy_s = 0.60, local_occupancy_s = 0.10,
        remaining = 53, threads = 8, threads_candidate = true, constants = nothing,
        failures = 1)
    @test v3.replan && v3.route === :process && v3.local_slots == 0
    @test v3.reason === :sample_failure
end

@testset "the guard cannot move to a threads route that is not a candidate" begin
    # Native GRAM point density withdraws the threads route entirely; a
    # single-threaded coordinator has none to move to either.
    cfg = _cfg(guard_factor = 1.5)
    plan = _guard_plan()
    obs = (worker_mean_s = 0.05, local_mean_s = 0.05,
           worker_occupancy_s = 0.60, local_occupancy_s = 0.10, remaining = 53)
    v_none = SCamp.predictive_guard_verdict(plan, cfg; obs...,
        threads = 8, threads_candidate = false, constants = nothing)
    @test !v_none.replan && v_none.route === :process
    v_one = SCamp.predictive_guard_verdict(plan, cfg; obs...,
        threads = 1, threads_candidate = true, constants = nothing)
    @test !v_one.replan && v_one.route === :process
    # Nothing left to run is nothing to re-plan.
    v_empty = SCamp.predictive_guard_verdict(plan, cfg;
        worker_mean_s = 0.05, local_mean_s = 0.05,
        worker_occupancy_s = 0.60, local_occupancy_s = 0.10,
        remaining = 0, threads = 8, threads_candidate = true, constants = nothing)
    @test !v_empty.replan
end

@testset "widening from local slots to threads tasks is charged, never credited" begin
    # The threads plan is priced from the observed local class. Going from 3
    # slots to 8 tasks on one heap costs the USL ratio when the machine has
    # constants, and nothing when it does not -- but never less than nothing.
    cfg = _cfg(guard_factor = 1.5)
    plan = _guard_plan()
    obs = (worker_mean_s = 0.05, local_mean_s = 0.05,
           worker_occupancy_s = 0.60, local_occupancy_s = 0.10,
           remaining = 53, threads = 8, threads_candidate = true)
    usl = _cfg(guard_factor = 1.5, heap_model = :usl)
    flat = SCamp.predictive_guard_verdict(plan, usl; obs..., constants = nothing)
    charged = SCamp.predictive_guard_verdict(plan, usl; obs..., constants = _constants())
    @test charged.threads_s >= flat.threads_s
    # A machine contended enough can make the move not worth it at all.
    steep = SCamp.predictive_guard_verdict(plan, usl; obs...,
        constants = _constants(alpha = 0.5, beta_alloc = 0.5))
    @test steep.threads_s > charged.threads_s
    # Under the default heap model the same constants are not charged at all,
    # so the guard prices the move exactly as an uncalibrated machine would.
    @test SCamp.predictive_guard_verdict(plan, cfg; obs...,
        constants = _constants()).threads_s == flat.threads_s
end

@testset "the reachable threads width is the local slots already running" begin
    # The caller hands the verdict the width it can actually reach WITHOUT a
    # barrier, which is the local slots already consuming the queue: closing
    # the pool class leaves them, and starting new consumers mid-dispatch is
    # the widening the guard is not allowed to do. Priced at a width it could
    # not run, the comparison would be against a plan that does not exist.
    cfg = _cfg(guard_factor = 1.5)
    plan = _guard_plan(local_slots = 3)
    obs = (worker_mean_s = 0.038, local_mean_s = 0.046,
           worker_occupancy_s = 0.333, local_occupancy_s = 0.050,
           remaining = 53, threads_candidate = true, constants = nothing)
    reachable = SCamp.predictive_guard_verdict(plan, cfg; obs..., threads = plan.local_slots)
    wide = SCamp.predictive_guard_verdict(plan, cfg; obs..., threads = 8)
    # Three slots at 50 ms cannot beat eight workers at 333 ms plus those same
    # three slots, so at the reachable width this verdict does not fire --
    # which is the honest answer, and the optimistic one would not have been.
    @test !reachable.replan
    @test reachable.reason === :threads_no_better
    @test reachable.threads_s > reachable.continue_s
    # A verdict that does not move reports the plan it is leaving alone, not a
    # width it declined to use.
    @test reachable.workers == plan.workers
    @test reachable.route === plan.route
    # The same observation at a width the caller cannot reach would have fired,
    # and the difference is entirely the pricing: eight consumers against three.
    @test wide.reason === :workers_occupying_more_than_threads
    @test wide.workers == 8
    @test reachable.threads_s > wide.threads_s
    @test reachable.continue_s == wide.continue_s
end

@testset "the route switch is off by default, and the evidence is still reported" begin
    # The capability exists and is measured; the default is what measured
    # better. With the switch off the guard computes the same comparison and
    # says so, but leaves the plan alone.
    off = _cfg(guard_factor = 1.5, route_switch = false)
    plan = _guard_plan()
    obs = (worker_mean_s = 0.038, local_mean_s = 0.046,
           worker_occupancy_s = 0.333, local_occupancy_s = 0.050,
           remaining = 53, threads = 8, threads_candidate = true, constants = nothing)
    v = SCamp.predictive_guard_verdict(plan, off; obs...)
    @test !v.replan
    @test v.route === :process && v.local_slots == plan.local_slots
    @test v.reason === :route_switch_disabled
    # The measurement that would have driven it is still there to read.
    @test v.occupancy_ratio > off.guard_factor
    @test v.threads_s < v.continue_s
    # And with the switch on, the same observation moves the remainder.
    on = _cfg(guard_factor = 1.5, route_switch = true)
    @test SCamp.predictive_guard_verdict(plan, on; obs...).route === :threads
    # The local-slot direction is unaffected by the switch either way.
    slow_locals = (worker_mean_s = 0.05, local_mean_s = 0.20,
                   worker_occupancy_s = 0.05, local_occupancy_s = 0.05,
                   remaining = 53, threads = 8, threads_candidate = true,
                   constants = nothing)
    for cfg in (off, on)
        v2 = SCamp.predictive_guard_verdict(plan, cfg; slow_locals...)
        @test v2.replan && v2.route === :process && v2.local_slots == 0
    end
end

@testset "predictive_replan refuses anything that is not a move toward static" begin
    plan = _guard_plan()
    @test_throws ArgumentError SCamp.predictive_replan(plan, 0, 53, nothing; route = :none,
                                                        workers = 1)
    @test_throws ArgumentError SCamp.predictive_replan(plan, 2, 53, nothing; route = :threads,
                                                        workers = 8)
    # The reduce-L path is unchanged and still cannot widen.
    @test SCamp.predictive_replan(plan, 99, 53, nothing).local_slots == plan.local_slots
    # Pricing after a re-plan. The reduce-L path keeps the pool pricing the
    # round was planned with; the threads path does not carry it and must not,
    # because a plan with no pool workers has no pool-worker class -- its cost
    # is entirely `heap_slowdown`, and `worker_slowdown` is the route's own 1.0.
    priced = SCamp._predictive_plan(:process, 8, 3, 64, false, nothing, 0.6)
    @test priced.worker_slowdown == 1.6
    @test SCamp.predictive_replan(priced, 1, 53, nothing).worker_slowdown == 1.6
    moved = SCamp.predictive_replan(priced, 0, 53, nothing; route = :threads, workers = 4)
    @test moved.route === :threads
    @test moved.worker_slowdown == 1.0
    @test moved.heap_slowdown == SCamp.predictive_heap_slowdown(nothing, 4)
    @test moved.makespan == SCamp.predictive_makespan(53, fill(moved.heap_slowdown, 4))
end
