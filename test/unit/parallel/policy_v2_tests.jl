using Test
using SpaceAGORA

# SPACEAGORA_PARALLEL_POLICY_V2 / profile R6: the revised inner-policy and
# outer-route behaviours, each asserted against the shipped behaviour it
# replaces so the switch is demonstrably the only difference.

const PP = SpaceAGORA.SimulationModel.ParallelPolicy
const PPr = SpaceAGORA.ParallelProfiles
const SCamp = SpaceAGORA.SimulationCampaigns

function _v2_feat(; samples = 64)
    return SCamp.campaign_route_features(
        samples = samples, n_sats = 1,
        density_family = "exponential", mission_time_s = 3600.0,
    )
end

function _v2_record!(state, feat, route, mean_s; n = 64, reps = 5)
    for _ in 1:reps
        PPr.record_outer_route_feedback!(
            state, feat; route = route, successes = n, failures = 0,
            elapsed_success_s = n * mean_s, elapsed_success_sq_sum_s = n * mean_s^2,
        )
    end
end

@testset "R6 is R5 plus the V2 switch" begin
    @test PPr.parse_parallel_profile("R6") === PPr.R6
    @test PPr.parse_parallel_profile("policy_v2") === PPr.R6
    @test PPr.parallel_profile_name(PPr.R6) == "R6"
    r5 = PPr.profile_config(PPr.R5)
    r6 = PPr.profile_config(PPr.R6)
    @test r6.policy_v2 && !r5.policy_v2
    # Everything else identical, by construction.
    for name in fieldnames(PPr.ParallelProfileConfig)
        name in (:profile, :label, :policy_v2) && continue
        @test getfield(r5, name) == getfield(r6, name)
    end
    pairs5 = Dict(PPr.profile_env_pairs(PPr.R5; preserve_existing = false))
    pairs6 = Dict(PPr.profile_env_pairs(PPr.R6; preserve_existing = false))
    @test pairs5["SPACEAGORA_PARALLEL_POLICY_V2"] == "0"
    @test pairs6["SPACEAGORA_PARALLEL_POLICY_V2"] == "1"
    @test pairs6["SPACEAGORA_PARALLEL_PROFILE"] == "R6"
    # The switch is the only env-level difference.
    diffs = [k for k in keys(pairs6) if pairs5[k] != pairs6[k]]
    @test Set(diffs) == Set(["SPACEAGORA_PARALLEL_POLICY_V2", "SPACEAGORA_PARALLEL_PROFILE"])
end

@testset "The switch is snapshotted, default off" begin
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => nothing) do
        @test !PP.policy_v2_enabled()
        @test !PP.snapshot_policy_decision_env().policy_v2
        @test !PPr.OuterRouteTuning().mc_route_by_core_budget
    end
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1") do
        @test PP.policy_v2_enabled()
        @test PP.snapshot_policy_decision_env().policy_v2
        @test PPr.OuterRouteTuning().mc_route_by_core_budget
    end
end

@testset "Hint confidence width scales with the arm's spread under V2" begin
    # Four samples: 0.8, 1.2, 0.8, 1.2 ms. Mean 1.0 ms, std 0.2 ms.
    spread = PP.AdaptiveChoiceStats(
        samples = 4, successes = 4, failures = 0,
        elapsed_sum_ns = 4.0e6, elapsed_sq_sum_ns = 2 * 0.64e12 + 2 * 1.44e12,
    )
    flat = PP.AdaptiveChoiceStats(
        samples = 4, successes = 4, failures = 0,
        elapsed_sum_ns = 4.0e6, elapsed_sq_sum_ns = 4 * 1.0e12,
    )
    total = Int64(8)
    mean0, w0 = PP._hint_mean_and_width(spread, total, 1.5)
    mean1, w1 = PP._hint_mean_and_width(spread, total, 1.5; scaled = true)
    @test mean0 ≈ 1.0e6 && mean1 ≈ 1.0e6
    # Shipped: a few nanoseconds against a millisecond mean -- the bound is the
    # mean, so the chooser was greedy.
    @test w0 < 10.0
    @test w1 ≈ w0 * 2.0e5
    _, wflat = PP._hint_mean_and_width(flat, total, 1.5; scaled = true)
    @test wflat == 0.0
    # The chooser threads the flag through; both must still return a candidate.
    sig = "profile=test|machine=test|src=v2_probe|items=5_8|thr=1|budget=5_8|outer=0|heavy_only=0|heavy=1"
    cands = PP._hint_candidate_allotments(8, 8)
    a = PP._hint_choose_allotment(sig, cands)
    b = PP._hint_choose_allotment(sig, cands; scaled_width = true)
    @test a.allotment in cands && b.allotment in cands
end

@testset "Hint gate fails closed on an unmeasured source under V2" begin
    env_off = withenv("SPACEAGORA_PARALLEL_POLICY_V2" => nothing) do
        PP.snapshot_policy_decision_env()
    end
    env_on = withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1") do
        PP.snapshot_policy_decision_env()
    end
    @test env_off.hint_work_ratio > 0.0
    PP.with_policy_context() do
        # No observation yet for this source in this context.
        @test PP._hint_layer_pays(:v2_gate_probe, env_off)
        @test !PP._hint_layer_pays(:v2_gate_probe, env_on)
        # One large observation opens the gate either way.
        PP.record_policy_observation!(
            :v2_gate_probe; mode = :auto, num_items = 8, use_threads = false,
            elapsed_ns = 10^9, env = env_on)
        @test PP._hint_layer_pays(:v2_gate_probe, env_on)
        @test PP._hint_layer_pays(:v2_gate_probe, env_off)
    end
end

@testset "Observations from worker tasks reach the scoped context when told where it is" begin
    env = PP.snapshot_policy_decision_env()
    PP.with_policy_context() do
        scoped = PP._active_policy_context()
        # Task-local storage does not cross task boundaries: a spawned task sees
        # the global context, not this one. This is the mechanism behind the
        # misrouted per-satellite observations.
        @test fetch(Threads.@spawn PP._active_policy_context()) !== scoped

        fetch(Threads.@spawn PP.record_policy_observation!(
            :v2_ctx_explicit; mode = :auto, num_items = 4, use_threads = false,
            elapsed_ns = 5_000, env = env, ctx = scoped))
        @test haskey(scoped.adaptive_state, :v2_ctx_explicit)
        @test scoped.adaptive_state[:v2_ctx_explicit].elapsed_ema_ns ≈ 5_000.0

        fetch(Threads.@spawn PP.record_policy_observation!(
            :v2_ctx_implicit; mode = :auto, num_items = 4, use_threads = false,
            elapsed_ns = 5_000, env = env))
        @test !haskey(scoped.adaptive_state, :v2_ctx_implicit)

        # The buffer-side accessor: nothing without buffers, the context with.
        @test PP.policy_context_hint(nothing) === nothing
        @test PP.policy_context_hint((shared_buffers = (x = 1,),)) === nothing
        buffers = (policy_context = Ref{Union{Nothing, PP.AbstractPolicyContext}}(scoped),)
        @test PP.policy_context_hint((shared_buffers = buffers,)) === scoped
        buffers[:policy_context][] = nothing
        @test PP.policy_context_hint((shared_buffers = buffers,)) === nothing
    end
end

@testset "A beaten default exploits the proven winner instead of exploring" begin
    # Monte Carlo default is :threads. Threads measured slow, process measured
    # fast, :none never tried. Shipped: the default is beaten so it is unproven,
    # exploration re-enables, and the under-sampled :none is chosen -- serial,
    # on a campaign where a 9x-faster route has five campaigns of evidence.
    feat = _v2_feat()
    beaten = PPr.OuterRouteState()
    _v2_record!(beaten, feat, :threads, 0.90)
    _v2_record!(beaten, feat, :process, 0.10)
    shipped = PPr.OuterRouteTuning(explore_until_any_proven = false)
    revised = PPr.OuterRouteTuning(explore_until_any_proven = true)
    @test PPr.default_outer_route(feat; tuning = shipped, machine_class = :large,
                                  threads_available = true) === :threads
    @test PPr.select_outer_route!(beaten, feat; tuning = shipped, machine_class = :large,
                                  threads_available = true) === :none
    @test PPr.select_outer_route!(beaten, feat; tuning = revised, machine_class = :large,
                                  threads_available = true) === :process

    # Cold, nothing is proven, so the revised guard still explores.
    cold = PPr.OuterRouteState()
    _v2_record!(cold, feat, :threads, 0.10; reps = 1)
    snap = PPr.outer_route_stats_snapshot(cold, PPr.outer_route_signature(feat))
    @test !PPr._any_candidate_proven(Symbol[:none, :threads, :process], snap, 2)
    chosen = PPr.select_outer_route!(cold, feat; tuning = revised, machine_class = :large,
                                     threads_available = true)
    @test chosen in (:none, :threads, :process)
end

@testset "Monte Carlo default follows the wider outer axis under V2" begin
    feat = _v2_feat()
    # (threads, workers) -> route. Every row is a measured point of the
    # 2026-09-02 B12/B13 grids; see _priority_outer_route_montecarlo.
    rows = [(12, 1, :threads), (6, 2, :threads), (4, 3, :threads),
            (3, 4, :process), (2, 6, :process), (1, 12, :process),
            (12, 12, :process), (8, 12, :process)]
    for (t, w, expect) in rows
        revised = PPr.OuterRouteTuning(mc_route_by_core_budget = true,
                                       process_max_workers = w, outer_thread_budget = t)
        shipped = PPr.OuterRouteTuning(mc_route_by_core_budget = false,
                                       process_max_workers = w, outer_thread_budget = t)
        @test PPr.default_outer_route(feat; tuning = revised, machine_class = :medium,
                                      threads_available = t > 1) === expect
        @test PPr.default_outer_route(feat; tuning = shipped, machine_class = :medium,
                                      threads_available = t > 1) === (t > 1 ? :threads : :process)
    end
    wide = PPr.OuterRouteTuning(mc_route_by_core_budget = true,
                                process_max_workers = 12, outer_thread_budget = 4)
    # A small machine never affords the process route.
    @test PPr.default_outer_route(feat; tuning = wide, machine_class = :small,
                                  threads_available = true) === :threads
    # Nor does a campaign below the size threshold.
    tiny = SCamp.campaign_route_features(samples = 8, n_sats = 1,
                                         density_family = "exponential", mission_time_s = 600.0)
    @test PPr.default_outer_route(tiny; tuning = wide, machine_class = :medium,
                                  threads_available = true) === :threads
end

@testset "The process worker cap honours SPACEAGORA_PERF_PROCS under V2" begin
    cores = PPr.usable_core_budget()
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1", "SPACEAGORA_PERF_PROCS" => "2") do
        @test PPr.OuterRouteTuning().process_max_workers == min(2, cores)
    end
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1", "SPACEAGORA_PERF_PROCS" => "junk") do
        @test PPr.OuterRouteTuning().process_max_workers == cores
    end
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => nothing, "SPACEAGORA_PERF_PROCS" => "2") do
        @test PPr.OuterRouteTuning().process_max_workers == cores
    end
end

@testset "Under V2 the other parallel arm is measured before any arm is exploited" begin
    feat = _v2_feat()
    st = PPr.OuterRouteState()
    revised = PPr.OuterRouteTuning(explore_until_any_proven = true, mc_route_by_core_budget = true,
                                   process_max_workers = 12, outer_thread_budget = 12)
    @test PPr.default_outer_route(feat; tuning = revised, machine_class = :large,
                                  threads_available = true) === :process
    # Five campaigns of the default, nothing else tried: threads is a
    # must-measure arm, so it is explored before process can be called proven.
    _v2_record!(st, feat, :process, 0.10)
    @test PPr.select_outer_route!(st, feat; tuning = revised, machine_class = :large,
                                  threads_available = true) === :threads
    _v2_record!(st, feat, :threads, 0.90; reps = 2)
    # Threads measured and slower: process is proven. Serial is never required.
    @test PPr.select_outer_route!(st, feat; tuning = revised, machine_class = :large,
                                  threads_available = true) === :process
    snap = PPr.outer_route_stats_snapshot(st, PPr.outer_route_signature(feat))
    @test !haskey(snap, :none)
    # The shipped guard would have called process proven without ever trying
    # threads, had serial been sampled first.
    cands = Symbol[:none, :threads, :process]
    only_default = PPr.OuterRouteState()
    _v2_record!(only_default, feat, :process, 0.10)
    _v2_record!(only_default, feat, :none, 1.50)
    snap2 = PPr.outer_route_stats_snapshot(only_default, PPr.outer_route_signature(feat))
    @test PPr._route_is_proven(cands, snap2, :process, 2)
    @test !PPr._route_is_proven(cands, snap2, :process, 2; must_measure = Symbol[:threads, :process])
end

@testset "V2 takes the static width on the per-call path" begin
    if Threads.nthreads() >= 4
        env_v2 = withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1",
                         "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
                         "SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
            PP.snapshot_policy_decision_env()
        end
        env_r5 = withenv("SPACEAGORA_PARALLEL_POLICY_V2" => nothing,
                         "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
                         "SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
            PP.snapshot_policy_decision_env()
        end
        PP.with_policy_context() do
            d = PP.thread_policy_decision(8; mode = :auto, threshold = 1,
                                          source = :v2_static_probe, env = env_v2)
            @test d.adaptive_enabled
            @test d.use_threads
            @test d.allotment == env_v2.inner_thread_budget
            # The hint path was not entered: no signature was built.
            @test PP.policy_telemetry_snapshot().last_signature == ""
        end
        PP.with_policy_context() do
            PP.thread_policy_decision(8; mode = :auto, threshold = 1,
                                      source = :r5_probe, env = env_r5)
            @test !isempty(PP.policy_telemetry_snapshot().last_signature)
        end
    else
        @test_skip "needs >= 4 Julia threads"
    end
end

@testset "V2 honours a cached heuristic verdict on a long solve" begin
    SE = SpaceAGORA.SimulationEngine
    sig = "v6|machine=v2test|budget=8|sats=257p|effs=1|harm=1|eff=v2test|dens=NoAtmosphereModel|outer=0"
    withenv("SPACEAGORA_RHS_CALIBRATE" => "auto",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => tempname() * ".toml") do
        setsolve = ns -> lock(SE._rhs_calib_lock) do
            SE._rhs_calib_cache[sig]["solve_ns"] = ns
        end
        SE._rhs_calib_store_heuristic!(sig, 1.0)
        setsolve(5.0e9)
        @test SE._rhs_calib_solve_exceeds_threshold(sig)
        @test SE._rhs_calib_cached_verdict(sig, false) === nothing
        @test SE._rhs_calib_cached_verdict(sig, true) === :heuristic
        setsolve(0.5e9)
        @test SE._rhs_calib_cached_verdict(sig, false) === :heuristic
        # A cached PLAN is still re-swept on a long solve, switch or not.
        SE._rhs_calib_store!(sig, SE._make_calib_satellite_batch_plan(), 1.0)
        setsolve(5.0e9)
        @test SE._rhs_calib_cached_verdict(sig, true) === nothing
        @test SE._rhs_calib_cached_verdict(sig, false) === nothing
        setsolve(0.5e9)
        @test SE._rhs_calib_cached_verdict(sig, true) !== nothing
        lock(SE._rhs_calib_lock) do
            delete!(SE._rhs_calib_cache, sig)
        end
    end
end

@testset "Density callback width ladder and no-regret chooser" begin
    SE = SpaceAGORA.SimulationEngine
    @test SE._callback_width_candidates(12) == [1, 2, 4, 8, 12]
    @test SE._callback_width_candidates(8) == [1, 2, 4, 8]
    @test SE._callback_width_candidates(1) == [1]
    @test SE._callback_width_candidates(3) == [1, 2, 3]
    # Static width 8 measured at 50; width 4 at 40 clears a 10% margin.
    @test SE._choose_callback_width(Dict(1 => 100.0, 2 => 60.0, 4 => 40.0, 8 => 50.0), 8, 0.10) == 4
    # Inside the margin: static retained (0 = no override).
    @test SE._choose_callback_width(Dict(1 => 100.0, 2 => 60.0, 4 => 46.0, 8 => 50.0), 8, 0.10) == 0
    # Static fastest: retained.
    @test SE._choose_callback_width(Dict(1 => 100.0, 4 => 60.0, 8 => 50.0), 8, 0.10) == 0
    # Static unmeasurable: the best measured candidate is taken.
    @test SE._choose_callback_width(Dict(1 => 100.0, 4 => 60.0, 8 => Inf), 8, 0.10) == 4
    # Ties resolve to the narrower width.
    @test SE._choose_callback_width(Dict(2 => 40.0, 4 => 40.0, 8 => 50.0), 8, 0.10) == 2
end

@testset "Split race batch sizing" begin
    # Three rounds per worker: 64 samples over [4, 8, 12] is one dispatch
    # round per width and not a measurement, so it does not race.
    @test SCamp._split_race_batch(64, [4, 8, 12]) == 0
    @test SCamp._split_race_batch(256, [4, 8, 12]) == 36   # max(36, cld(256, 12)); 1 + 108 + 12 <= 256
    @test SCamp._split_race_batch(40, [2, 4]) == 12        # max(12, cld(40, 8)); 1 + 24 + 4 <= 40
    @test SCamp._split_race_batch(16, [4, 8, 12]) == 0     # too small to race
    @test SCamp._split_race_batch(64, [12]) == 0           # nothing to race
    @test SCamp._split_race_batch(1, [2, 4]) == 0
    # A process-route warm-up batch touches every worker of the widest split.
    @test SCamp._split_race_batch(256, [4, 8, 12]; warm = 12) == 36   # 12 + 108 + 12 <= 256
    @test SCamp._split_race_batch(130, [4, 8, 12]; warm = 12) == 0    # 12 + 108 + 12 > 130
    @test SCamp._split_race_warm_count((route = :process,), [4, 8, 12], 64) == 12
    @test SCamp._split_race_warm_count((route = :process,), [4, 8, 12], 6) == 6
    @test SCamp._split_race_warm_count((route = :threads,), [4, 8, 12], 64) == 1
end

@testset "In-campaign split race on the threads route" begin
    if Threads.nthreads() >= 4
        feat = _v2_feat()
        st = PPr.OuterRouteState()
        tuning = PPr.OuterRouteTuning(split_race = true, explore_until_any_proven = true,
                                      mc_route_by_core_budget = false, process_max_workers = 1)
        seeds = collect(101:140)
        result = withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
                         "SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
            SCamp.run_monte_carlo(seeds; threads = :auto, route_features = feat,
                                  route_state = st, route_tuning = tuning) do seed
                seed * 2
            end
        end
        cands = PPr.outer_split_candidates(:threads; budget = Threads.nthreads(), n_units = 40, tuning = tuning)
        @test length(cands) >= 2
        @test length(result.samples) == 40
        @test [s.index for s in result.samples] == collect(1:40)
        @test [s.seed for s in result.samples] == seeds
        @test all(s -> s.success && s.value == 2 * s.seed, result.samples)
        @test result.threads in cands
        # Every raced width was recorded with min-samples weight, so the next
        # campaign can exploit without re-trying each width.
        sig = PPr._split_signature_chain(feat)[1]
        snap = PPr._outer_route_stats_snapshot_internal(st, sig)
        for w in cands
            arm = PPr._split_arm(:threads, w)
            @test haskey(snap, arm)
            @test snap[arm].campaigns == tuning.adaptive_min_samples
        end
        @test PPr.outer_split_history_present(st, feat, :threads)
        # Second campaign: history present, no race; the selector exploits.
        result2 = withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
                          "SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
            SCamp.run_monte_carlo(seeds; threads = :auto, route_features = feat,
                                  route_state = st, route_tuning = tuning) do seed
                seed * 2
            end
        end
        @test length(result2.samples) == 40
        @test result2.threads in cands
    else
        @test_skip "needs >= 4 Julia threads"
    end
end

@testset "Under V2 the route selector does not manufacture evidence" begin
    feat = _v2_feat()
    v2 = PPr.OuterRouteTuning(explore_routes = false, explore_until_any_proven = true,
                              mc_route_by_core_budget = true, process_max_workers = 12,
                              outer_thread_budget = 12)
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1") do
        @test !PPr.OuterRouteTuning().explore_routes
    end
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => nothing) do
        @test PPr.OuterRouteTuning().explore_routes
    end
    # Only the default measured: no trial of threads, the default stands.
    st = PPr.OuterRouteState()
    _v2_record!(st, feat, :process, 0.10; reps = 1)
    @test PPr.select_outer_route!(st, feat; tuning = v2, machine_class = :large,
                                  threads_available = true) === :process
    _v2_record!(st, feat, :process, 0.10; reps = 3)
    @test PPr.select_outer_route!(st, feat; tuning = v2, machine_class = :large,
                                  threads_available = true) === :process
    # History it already holds still decides: threads measured faster wins.
    _v2_record!(st, feat, :threads, 0.02)
    @test PPr.select_outer_route!(st, feat; tuning = v2, machine_class = :large,
                                  threads_available = true) === :threads
end
