using Test
using SpaceAGORA

# The predictive planner's campaign terms (the final round's tail and the
# cold pool's start-up), the online corrections to its cost-model parameters, and the leash.
# Pure logic: no simulations and no pool.
#
# The shapes and numbers are the TRX50's, from the archived targeted run
# trx50_targeted_cold_20260923_140152 (P3 at 32, P4 at 8 and 16; 11 repeats
# each) and the machine's calibration (job 20260921-202331-989854). They fix
# the DIRECTION a decision must take; the model is not fitted to them.

const SCamp = SpaceAGORA.SimulationCampaigns
const PCost = SpaceAGORA.SimulationModel.ParallelCost

_curve() = PCost.RateCurve([0.0, 10.0], [1.0, 1.0])

function _constants(; alpha = 0.156, beta_alloc = 0.00605)
    return PCost.MachineConstants(
        simd_lane = _curve(), coeff_touch = _curve(), parallel_speedup = _curve(),
        ns_per_scalar_item = 1.0, ns_per_queue_node = 1.0,
        dispatch_pool_ns_base = 1.0, dispatch_pool_ns_per_worker = 1.0,
        dispatch_batch_ns_base = 1.0, dispatch_batch_ns_per_worker = 1.0,
        ns_per_atomic = 1.0, reference_fma_ns = 1.0, reference_mem_ns = 1.0,
        usl_alpha_base = alpha, usl_beta_alloc = beta_alloc,
        schema_version = PCost.CALIBRATION_SCHEMA_VERSION,
    )
end

_cfg(; threads, heap_model = :locals) = SCamp.PredictivePlannerConfig(
    margin = 0.15, guard_factor = 1.5, local_slots_max = threads - 1, remote_overhead = 0.0,
    route_switch = false, heap_model = heap_model, local_thrash = 3.0)

# The round tail extracted from the archived run: the median over its 55
# pure-process rows (scripts/extract_campaign_cost_terms.py).
const TRX50_TAIL = 0.2261

_plan(; n, threads, pool = threads, heap_model = :locals, terms = SCamp.PredictiveCostTerms(),
      constants = _constants()) =
    SCamp.predictive_plan(n_samples = n, threads = threads, process_workers = pool,
                          threads_candidate = true, local_slots_cap = threads - 1,
                          constants = constants, config = _cfg(threads = threads, heap_model = heap_model),
                          terms = terms)

_find(planning, L) = only(filter(p -> p.route === :process && p.local_slots == L, planning.plans))

_rules() = SCamp.CampaignCorrectionRules(step = 0.05, prior_share = 0.25, stale_campaigns = 20)

@testset "default terms price every plan exactly as the model without them" begin
    for (n, T) in ((32, 8), (32, 16), (256, 32), (64, 8), (7, 4))
        with = _plan(n = n, threads = T)
        without = SCamp.predictive_plan(n_samples = n, threads = T, process_workers = T,
                                        threads_candidate = true, local_slots_cap = T - 1,
                                        constants = _constants(), config = _cfg(threads = T))
        @test [p.makespan for p in with.plans] == [p.makespan for p in without.plans]
        @test SCamp.predictive_plan_key(with.chosen) == SCamp.predictive_plan_key(without.chosen)
    end
    s = [1.0, 1.0, 1.54, 1.54]
    @test SCamp.predictive_plan_makespan(13, 2, s, SCamp.PredictiveCostTerms()) ==
          SCamp.predictive_makespan(13, s)
    @test_throws ArgumentError SCamp.PredictiveCostTerms(heap_scale = 0.0)
    @test_throws ArgumentError SCamp.PredictiveCostTerms(round_tail = -0.1)
    @test_throws ArgumentError SCamp.PredictiveCostTerms(startup = NaN)
end

@testset "each class pays the tail of its final round" begin
    @test SCamp.predictive_round_excess(1) == 0.0
    @test SCamp.predictive_round_excess(2) == 0.5
    @test SCamp.predictive_round_excess(16) ≈ sum(1 / i for i in 2:16)
    t = SCamp.PredictiveCostTerms(round_tail = 0.5)
    # Pure pool, 32 samples on 16 workers: the final round is sixteen wide.
    @test SCamp.predictive_plan_makespan(32, 16, fill(1.0, 16), t) ≈
          2.0 + 0.5 * SCamp.predictive_round_excess(16)
    # Sixteen workers and fifteen local slots: the pool's final round is one
    # sample and has no tail; the local slots' only round is fifteen wide.
    @test SCamp.predictive_plan_makespan(32, 16, fill(1.0, 31), t) ≈
          max(2.0, 1.0 + 0.5 * SCamp.predictive_round_excess(15))
    # The threads route (no pool) pays the same tail as the pool: it is the
    # samples' spread, not the pool's.
    @test SCamp.predictive_plan_makespan(32, 0, fill(1.0, 16), t) ≈
          SCamp.predictive_plan_makespan(32, 16, fill(1.0, 16), t)
    # A single consumer has no tail.
    @test SCamp.predictive_plan_makespan(5, 0, [1.0], t) == 5.0
    # A cold pool's workers start late, and the local slots take the first
    # samples meanwhile.
    cold = SCamp.PredictiveCostTerms(startup = 0.5)
    @test SCamp.predictive_plan_makespan(4, 2, fill(1.0, 4), cold) == 1.5
    @test SCamp.predictive_plan_makespan(2, 2, fill(1.0, 4), cold) == 1.0
end

@testset "the round tail breaks the P4-at-16 tie in favor of the mixed plans" begin
    # Without the term, every plan up to six local slots is priced at exactly
    # two rounds and the static plan wins the tie -- what R7 did on the TRX50.
    flat = _plan(n = 32, threads = 16)
    @test all(L -> _find(flat, L).makespan == 2.0, 0:6)
    @test flat.chosen.local_slots == 0
    @test flat.reason === :static_equivalent_best

    terms = SCamp.PredictiveCostTerms(round_tail = TRX50_TAIL)
    tailed = _plan(n = 32, threads = 16, terms = terms)
    spans = [_find(tailed, L).makespan for L in 0:6]
    @test issorted(spans; rev = true) && allunique(spans)
    @test !first(tailed.plans).static_equivalent
    # Under the calibrated heap term the best mixed plan (six slots) is
    # predicted 4% ahead, which is inside the margin: ranked ahead, not taken.
    @test tailed.reason === :margin_not_met
    @test tailed.chosen.local_slots == 0
    # The threads plan pays the same tail and still loses the tie to the pool.
    @test tailed.chosen.route === :process

    # With the heap term at what the archived traces measured for fifteen
    # local slots beside sixteen workers (local work 1.53x a worker's, where
    # the calibrated curve says 4.45x), the ranking is the measured one and
    # the gain clears the margin: 16 workers + 15 slots, R6's plan, 1.859 s
    # against the pure pool's 2.302 s.
    measured = SCamp.predictive_heap_slowdown(_constants(), 15)
    scaled = SCamp.PredictiveCostTerms(heap_scale = 1.53 / measured, round_tail = TRX50_TAIL)
    mixed = _plan(n = 32, threads = 16, terms = scaled)
    @test mixed.chosen.route === :process
    @test mixed.chosen.local_slots == 15
    @test mixed.reason === :predicted_gain

    # The same shape with no heap term at all ranks the same way.
    none = _plan(n = 32, threads = 16, heap_model = :none, terms = terms)
    @test none.chosen.local_slots == 15
    @test none.gain >= 0.15
end

@testset "the round tail puts the pure pool last at P4 at 8, as measured" begin
    # Measured (s): w8+l7 2.620 < w8+l4 3.002 < w8+l0 3.742.
    tailed = _plan(n = 32, threads = 8, terms = SCamp.PredictiveCostTerms(round_tail = TRX50_TAIL))
    @test _find(tailed, 0).makespan > _find(tailed, 7).makespan
    @test _find(tailed, 0).makespan > _find(tailed, 4).makespan
    # Without it the pure pool ranked ahead of seven slots.
    flat = _plan(n = 32, threads = 8)
    @test _find(flat, 0).makespan < _find(flat, 7).makespan
end

@testset "predictive_round_tail_observation is the extraction formula" begin
    # One archived P4-at-16 row: wall 2.28 s, mean sample 0.8956 s.
    @test SCamp.predictive_round_tail_observation(2.28, 0.8956, 32, 16) ≈
          (2.28 / 0.8956 - 2) / SCamp.predictive_round_excess(16)
    @test SCamp.predictive_round_tail_observation(1.5, 0.5, 20, 16) ≈
          (3.0 - 2) / SCamp.predictive_round_excess(4)
    # A final round of one sample has no tail to measure.
    @test isnan(SCamp.predictive_round_tail_observation(2.25, 1.0, 17, 16))
    @test isnan(SCamp.predictive_round_tail_observation(1.0, 0.0, 32, 16))
end

@testset "a heap correction moves toward the observation by a bounded step" begin
    rules = _rules()
    p = nothing
    previous = SCamp.correction_value(p, 1.0, 0, rules)
    @test previous == 1.0
    for k in 1:12
        p = SCamp.correction_observe(p, 1.0, 0.8, k, rules)
        value = SCamp.correction_value(p, 1.0, k, rules)
        @test value <= previous
        @test previous - value <= rules.step + 1e-12
        @test value >= 0.25 * 1.0 + 0.75 * 0.8 - 1e-12
        previous = value
    end
    @test p.evidence ≈ 0.8
    @test SCamp.correction_value(p, 1.0, 12, rules) ≈ 0.85
    @test p.observations == 12

    # One campaign cannot lock a verdict in: a single wild observation moves
    # the value by one step.
    one = SCamp.correction_observe(nothing, 1.0, 0.1, 1, rules)
    @test SCamp.correction_value(one, 1.0, 1, rules) ≈ 1.0 - 0.75 * 0.05

    # Through the campaign fold, with the guard's observation.
    c = SCamp.CampaignCorrections("fp", "token")
    consts = SCamp.PredictiveCampaignConstants(round_tail = TRX50_TAIL)
    for _ in 1:3
        SCamp.predictive_fold_campaign!(c, rules, consts; signature = "sig", shape_key = "shape",
                                        final_plan = "process@w8+l4", heap_scale_observed = 0.8)
    end
    priced = SCamp.predictive_cost_terms(c, consts, rules; signature = "sig", pool_cold = false)
    @test priced.terms.heap_scale ≈ 0.25 + 0.75 * 0.85
    @test priced.terms.round_tail == TRX50_TAIL
    @test c.last_plan["shape"] == "process@w8+l4"
end

@testset "a stale correction decays to the prior" begin
    rules = _rules()
    c = SCamp.CampaignCorrections("fp", "token")
    consts = SCamp.PredictiveCampaignConstants(round_tail = 0.5, pool_startup_s = 0.8)
    SCamp.predictive_fold_campaign!(c, rules, consts; signature = "sig", shape_key = "s",
                                    final_plan = "process@w8+l0", heap_scale_observed = 0.5,
                                    tail_observed = 0.3, worker_sample_s = 0.4)
    corrected = SCamp.predictive_cost_terms(c, consts, rules; signature = "sig", pool_cold = true)
    @test corrected.terms.heap_scale < 1.0
    @test corrected.terms.round_tail < 0.5
    @test corrected.sample_time_s == 0.4
    @test corrected.terms.startup ≈ 0.8 / 0.4
    for _ in 1:rules.stale_campaigns
        SCamp.predictive_fold_campaign!(c, rules, consts; signature = "other", shape_key = "s",
                                        final_plan = "process@w8+l0")
    end
    @test c.heap_scale !== nothing              # exactly at the limit: still held
    SCamp.predictive_fold_campaign!(c, rules, consts; signature = "other", shape_key = "s",
                                    final_plan = "process@w8+l0")
    @test c.heap_scale === nothing
    @test c.round_tail === nothing
    @test !haskey(c.sample_time_s, "sig")
    back = SCamp.predictive_cost_terms(c, consts, rules; signature = "sig", pool_cold = true)
    @test back.terms.heap_scale == 1.0
    @test back.terms.round_tail == 0.5
    @test back.terms.startup == 0.0            # sample time unknown again
    # Without a prior in the constants file the tail term is off and is never
    # learned.
    bare = SCamp.PredictiveCampaignConstants()
    d = SCamp.CampaignCorrections("fp", "token")
    SCamp.predictive_fold_campaign!(d, rules, bare; signature = "sig", shape_key = "s",
                                    final_plan = "process@w8+l0", tail_observed = 0.9)
    @test d.round_tail === nothing
    @test SCamp.predictive_cost_terms(d, bare, rules; signature = "sig", pool_cold = false).terms.round_tail == 0.0
end

@testset "the leash allows one slot per campaign and refuses a two-slot jump" begin
    # P4 at 8 with the round tail: the model's plan is four local slots, and
    # it prices three slots above two (4.158 against 4.113 sample times).
    planning = _plan(n = 32, threads = 8, terms = SCamp.PredictiveCostTerms(round_tail = TRX50_TAIL))
    @test planning.chosen.local_slots == 4
    @test planning.reason === :predicted_gain
    @test _find(planning, 3).makespan > _find(planning, 2).makespan

    plan, why = SCamp.predictive_leash(planning, nothing)
    @test why === :no_previous && plan === planning.chosen
    plan, why = SCamp.predictive_leash(planning, "process@w8+l3")
    @test why === :within_leash && plan === planning.chosen
    plan, why = SCamp.predictive_leash(planning, "process@w8+l5")
    @test why === :within_leash && plan === planning.chosen
    # Two slots away and the step between is priced worse than where the shape
    # already is: the jump is refused and the shape stays where it was.
    plan, why = SCamp.predictive_leash(planning, "process@w8+l2")
    @test why === :leash_held
    @test plan.local_slots == 2
    # Three away with a step the model favors: one step toward the plan.
    plan, why = SCamp.predictive_leash(planning, "process@w8+l1")
    @test why === :leashed && plan.local_slots == 2
    # From the static plan every walk starts at one slot, whichever route the
    # static plan was on.
    plan, why = SCamp.predictive_leash(planning, "process@w8+l0")
    @test why === :leashed && plan.local_slots == 1
    plan, why = SCamp.predictive_leash(planning, "threads@w8+l0")
    @test why === :leashed && plan.local_slots == 1
    # Downward too.
    plan, why = SCamp.predictive_leash(planning, "process@w8+l7")
    @test why === :leashed && plan.local_slots == 6

    # A static choice is never held back.
    static = _plan(n = 32, threads = 16)
    @test static.chosen.static_equivalent
    plan, why = SCamp.predictive_leash(static, "process@w16+l9")
    @test why === :static && plan === static.chosen
end

@testset "P3 at 32 still runs the pure pool after the corrections P4 induces" begin
    rules = _rules()
    consts = SCamp.PredictiveCampaignConstants(round_tail = TRX50_TAIL, pool_startup_s = 0.8146)
    c = SCamp.CampaignCorrections("fp", "token")
    # Eleven P4-at-8 campaigns: guard observed/predicted 0.71-0.85 on the
    # archived run (0.8 here), and a round tail below the prior.
    for _ in 1:11
        SCamp.predictive_fold_campaign!(c, rules, consts; signature = "p4", shape_key = "p4@8",
                                        final_plan = "process@w8+l4", heap_scale_observed = 0.8,
                                        tail_observed = 0.17, worker_sample_s = 0.85)
    end
    priced = SCamp.predictive_cost_terms(c, consts, rules; signature = "p3", pool_cold = false)
    @test priced.terms.heap_scale < 1.0
    @test priced.terms.round_tail < TRX50_TAIL
    p3 = _plan(n = 256, threads = 32, terms = priced.terms)
    @test p3.chosen.route === :process
    @test p3.chosen.local_slots == 0
    @test SCamp.predictive_plan_key(first(SCamp.predictive_leash(p3, "process@w32+l0"))) == "process@w32+l0"
    # On a cold pool the start-up is charged only where this signature's sample
    # time is known: P4's is, P3's is not.
    cold_p3 = SCamp.predictive_cost_terms(c, consts, rules; signature = "p3", pool_cold = true)
    @test cold_p3.terms.startup == 0.0
    cold_p4 = SCamp.predictive_cost_terms(c, consts, rules; signature = "p4", pool_cold = true)
    @test cold_p4.terms.startup ≈ 0.8146 / 0.85
end

@testset "the corrections file round-trips and is keyed by machine and code" begin
    dir = mktempdir()
    path = joinpath(dir, "campaign_corrections_test.toml")
    rules = _rules()
    consts = SCamp.PredictiveCampaignConstants(round_tail = 0.5)
    c = SCamp.CampaignCorrections("fp", "token-a")
    SCamp.predictive_fold_campaign!(c, rules, consts; signature = "a|b=c", shape_key = "a|b=c|n=32",
                                    final_plan = "process@w8+l5", heap_scale_observed = 0.9,
                                    tail_observed = 0.4, worker_sample_s = 0.85)
    SCamp.save_campaign_corrections(c, path)
    back = SCamp.load_campaign_corrections(path; fingerprint = "fp", code_token = "token-a")
    @test back.campaigns == 1
    @test back.heap_scale == c.heap_scale
    @test back.round_tail == c.round_tail
    @test back.sample_time_s == c.sample_time_s
    @test back.last_plan == c.last_plan
    # Another code token, another machine, another schema: cold.
    @test SCamp.load_campaign_corrections(path; fingerprint = "fp", code_token = "token-b").campaigns == 0
    @test SCamp.load_campaign_corrections(path; fingerprint = "other", code_token = "token-a").campaigns == 0
    write(path, replace(read(path, String), "schema_version = 1" => "schema_version = 99"))
    @test SCamp.load_campaign_corrections(path; fingerprint = "fp", code_token = "token-a").campaigns == 0
    @test SCamp.load_campaign_corrections(joinpath(dir, "absent.toml");
                                          fingerprint = "fp", code_token = "token-a").campaigns == 0

    withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS" => nothing) do
        @test SCamp.campaign_corrections_mode() === :on
    end
    withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS" => "0") do
        @test SCamp.campaign_corrections_mode() === :off
        @test SCamp.campaign_corrections() === nothing
    end
    withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS" => "read") do
        @test SCamp.campaign_corrections_mode() === :read
    end
    withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS" => "sometimes") do
        @test_throws ArgumentError SCamp.campaign_corrections_mode()
    end
    withenv("SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => path) do
        @test SCamp.campaign_corrections_path() == path
    end
    @test_throws ArgumentError SCamp.CampaignCorrectionRules(step = 0.0)
    @test_throws ArgumentError SCamp.CampaignCorrectionRules(prior_share = 0.0)
    @test_throws ArgumentError SCamp.CampaignCorrectionRules(stale_campaigns = 0)
end

@testset "the [campaign] table is read, and survives a re-calibration" begin
    dir = mktempdir()
    path = joinpath(dir, "cost_constants_test.toml")
    PCost.save_machine_constants(_constants(), path)
    @test SCamp.load_predictive_campaign_constants(path).round_tail === nothing
    open(path, "a") do io
        println(io, "\n[campaign]\nround_tail = 0.2261\npool_startup_s = 0.8146\nsource = \"archived run\"")
    end
    loaded = SCamp.load_predictive_campaign_constants(path)
    @test loaded.round_tail == 0.2261
    @test loaded.pool_startup_s == 0.8146
    @test loaded.source == "archived run"
    PCost.save_machine_constants(_constants(alpha = 0.2), path)
    again = SCamp.load_predictive_campaign_constants(path)
    @test again.round_tail == 0.2261
    @test again.pool_startup_s == 0.8146
    @test PCost.load_machine_constants(path).usl_alpha_base == 0.2
    @test SCamp.load_predictive_campaign_constants(joinpath(dir, "absent.toml")).pool_startup_s === nothing
end
