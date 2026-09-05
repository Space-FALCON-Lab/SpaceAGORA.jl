@inline function thread_policy_decision(
    num_items::Int;
    mode::Symbol,
    threshold::Int,
    heavy_work::Bool=true,
    heavy_only::Bool=false,
    outer_active::Bool=false,
    allow_with_outer::Bool=true,
    source::Symbol=:other,
    env::Union{Nothing, PolicyDecisionEnvConfig}=nothing
)
    budget = env === nothing ? effective_inner_thread_budget() : env.inner_thread_budget
    min_auto_budget = env === nothing ? auto_thread_min_budget(source) : _snapshot_auto_min_budget(env, source)
    auto_budget_allowed = mode != :auto || budget >= min_auto_budget
    adaptive_enabled = (mode == :auto) && (env === nothing ? adaptive_policy_enabled() : env.adaptive_enabled)
    # Short-circuit when the decision is already forced.
    #
    # `use_threads` below opens with `budget <= 1 || num_items <= 1 -> false`,
    # unconditionally and before anything the adaptive branch computes is
    # consulted, and `allotted` then collapses to 1. So in those two cases the
    # adaptive branch cannot change the answer, yet it still built a workload
    # signature string, performed a persistent-hint lookup, and took two lock
    # acquisitions -- once per decision.
    #
    # That is the dominant differential cost of adaptive routing on
    # process-routed campaigns, where every worker is launched with one thread
    # and so has budget 1 on every call. Measured on independent_1sat_1hr at 256
    # samples: R5 ran 58% slower per sample than the same workload on a pinned
    # process route, and disabling the policy engine alone recovered 96% of it.
    #
    # Only work is skipped, never a decision: use_threads and allotment are
    # identical either way, verified across 1728 input combinations.
    # adaptive_enabled itself is left as the caller configured it, so telemetry
    # still reports whether adaptive routing was *on*, distinct from whether it
    # had anything to decide.
    # Every way `use_threads` below resolves to false without consulting the
    # adaptive branch, not just the two budget/item cases it started with.
    #
    # `use_threads` already short-circuits on an under-threshold item count, on
    # an outer split that disallows inner threading, and on the heavy-only guard
    # rejecting light work -- but all three were tested *after* the adaptive
    # branch had built a signature string, taken a persistent-hint lock, done a
    # hint lookup and taken the telemetry lock. The result was discarded every
    # time.
    #
    # The heavy-only case is the expensive one in practice. A single spacecraft
    # with two cheap effectors puts num_items at 2 -- above the forced floor --
    # so the full adaptive machinery ran once per RHS call to decide whether to
    # thread two effectors that the heavy-work guard was always going to reject.
    # Measured on the montecarlo_heavy_aerobraking shape (1 satellite,
    # inverse-square + aero), turning the effector policy off was worth 48% of
    # wall time and disabling adaptation entirely 45%, on a workload where
    # neither can ever enable threading.
    #
    # `mode == :off` needs no entry: adaptive_enabled already requires :auto.
    # `mode == :on` is deliberately absent -- it forces use_threads true, but the
    # allotment still comes from the adaptive branch, so that decision is not
    # forced. It cannot reach here anyway for the same :auto reason.
    #
    # As before, only work is skipped and never a decision: use_threads and
    # allotment are identical either way, because each disjunct below
    # independently pins use_threads false and allotted collapses to 1.
    decision_forced =
        budget <= 1 || num_items <= 1 ||
        !auto_budget_allowed ||
        num_items < max(1, threshold) ||
        (outer_active && !allow_with_outer) ||
        (heavy_only && !heavy_work)
    # V2: callback and effector widths come from the static rule,
    # min(items, budget), not from the per-call hint store or AIMD.
    #
    # Nothing in the 2026-09-02 paper run shows the per-call layer winning. The
    # AIMD branch never reads elapsed time (its window score is a fill ratio,
    # observation_tracking.jl), so under R4 it is the static width plus lock and
    # signature overhead; the R5 hint chooser was a two-sample greedy argmin
    # (see _hint_mean_and_width). Where R5 wins large (stack256_e3..e5, -26% to
    # -46%) the plan telemetry credits the sweep-pinned satellite_batch route,
    # and on the one case dominated by a threaded callback
    # (atmo256_gram_surrogate) the static inner_only route is best at both
    # thread counts, with R5 +25% and +8% behind it. The hint and AIMD paths
    # stay reachable with the switch off; R6 simply does not take them.
    static_width_only = env !== nothing && env.policy_v2
    adaptive_active = adaptive_enabled && !decision_forced && !static_width_only
    # The hint layer is consulted only when it pays for itself on this machine;
    # see _hint_layer_pays and hint_work_ratio.
    hint_layer_active = adaptive_active && _hint_layer_pays(source, env)
    measured_reward = hint_layer_active &&
        (env === nothing ? persistent_hints_enabled() : env.persistent_hints) &&
        (env === nothing ? adaptive_measured_reward_enabled() : env.adaptive_measured_reward)
    signature = ""
    hint_allotment = Int64(1)
    hint_confidence = 0.0
    hint_regret_ns = 0.0
    hints_loaded = false
    hints_entries = 0
    desire = 1
    allotment = 1
    if adaptive_active
        signature = _hint_workload_signature(
            source,
            num_items,
            threshold,
            budget,
            outer_active,
            heavy_only,
            heavy_work
        )
        if hint_layer_active
            hint = _hint_choose_allotment(
                signature,
                _hint_candidate_allotments(num_items, budget);
                scaled_width=(env !== nothing && env.policy_v2)
            )
            hint_allotment = hint.allotment
            hint_confidence = hint.confidence
            hint_regret_ns = hint.regret_ns
            lock(_persistent_hint_lock) do
                hints_loaded = _persistent_hint_state[].loaded
                hints_entries = _hint_entry_count(_persistent_hint_state[])
            end
        end
        # Measured-reward mode drives the width from the hint layer; otherwise
        # it is the static answer, min(items, budget), the same one the
        # non-adaptive branch below gives. The AIMD controller that used to sit
        # here was removed: its window score was a fill ratio that never read
        # elapsed time, and since the width seed at the cap it had been a fixed
        # point at full width.
        lock(_active_policy_context().lock) do
            st = _adaptive_state_for(source)
            st.desire = measured_reward ? max(1, min(Int(hint_allotment), budget)) : max(1, budget)
            desire = st.desire
        end
        allotment = max(1, min(desire, budget))
    else
        if mode == :on
            desire = max(2, budget)
        elseif mode == :auto
            desire = max(1, budget)
        end
        allotment = if mode == :off
            1
        else
            max(1, min(max(1, desire), budget))
        end
    end

    use_threads =
        if budget <= 1 || num_items <= 1
            false
        elseif !auto_budget_allowed
            false
        elseif mode == :off
            false
        elseif mode == :on
            true
        elseif outer_active && !allow_with_outer
            false
        elseif heavy_only && !heavy_work
            false
        elseif adaptive_active
            desire > 1 && num_items >= max(1, threshold)
        else
            num_items >= max(1, threshold)
        end

    allotted = use_threads ? allotment : 1
    _record_policy_decision!(
        source,
        mode,
        threshold,
        num_items,
        budget,
        adaptive_enabled,
        desire,
        allotted,
        outer_active,
        allow_with_outer,
        heavy_only,
        heavy_work,
        use_threads,
        signature,
        hint_allotment,
        hint_confidence,
        hint_regret_ns,
        hints_loaded,
        hints_entries,
        env
    )
    return (
        use_threads=use_threads,
        allotment=allotted,
        budget=budget,
        mode=mode,
        threshold=max(1, threshold),
        num_items=max(0, num_items),
        adaptive_enabled=adaptive_enabled,
        desire=max(1, desire)
    )
end

@inline function use_threads_policy(
    num_items::Int;
    mode::Symbol,
    threshold::Int,
    heavy_work::Bool=true,
    heavy_only::Bool=false,
    outer_active::Bool=false,
    allow_with_outer::Bool=true,
    source::Symbol=:other
)::Bool
    return thread_policy_decision(
        num_items;
        mode=mode,
        threshold=threshold,
        heavy_work=heavy_work,
        heavy_only=heavy_only,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=source
    ).use_threads
end
