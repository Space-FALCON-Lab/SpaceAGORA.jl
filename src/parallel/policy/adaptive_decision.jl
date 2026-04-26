@inline function thread_policy_decision(
    num_items::Int;
    mode::Symbol,
    threshold::Int,
    heavy_work::Bool=true,
    heavy_only::Bool=false,
    outer_active::Bool=false,
    allow_with_outer::Bool=true,
    source::Symbol=:other
)
    budget = effective_inner_thread_budget()
    adaptive_enabled = (mode == :auto) && adaptive_policy_enabled()
    measured_reward = adaptive_enabled && persistent_hints_enabled() && adaptive_measured_reward_enabled()
    bootstrap_threads = adaptive_enabled && adaptive_bootstrap_threads()
    control_tail_guard = adaptive_enabled && adaptive_control_tail_guard()
    signature = ""
    hint_allotment = Int64(1)
    hint_confidence = 0.0
    hint_regret_ns = 0.0
    hints_loaded = false
    hints_entries = 0
    desire = 1
    allotment = 1
    if adaptive_enabled
        signature = _hint_workload_signature(
            source,
            num_items,
            threshold,
            budget,
            outer_active,
            heavy_only,
            heavy_work
        )
        hint = _hint_choose_allotment(signature, _hint_candidate_allotments(num_items, budget))
        hint_allotment = hint.allotment
        hint_confidence = hint.confidence
        hint_regret_ns = hint.regret_ns
        lock(_persistent_hint_lock) do
            hints_loaded = _persistent_hint_state[].loaded
            hints_entries = _hint_entry_count(_persistent_hint_state[])
        end
        ρ = adaptive_rho()
        desire_cap = _adaptive_desire_cap(budget, ρ)
        lock(_policy_telemetry_lock) do
            st = _adaptive_state_for(source)
            st.desire = min(max(1, st.desire), desire_cap)
            if measured_reward
                # Measured-reward mode: drive desire directly from elapsed-time hint choice.
                st.desire = min(desire_cap, max(1, Int(hint_allotment)))
            else
                # Tail guard: avoid a cold-start serial window on obviously parallel workloads.
                if bootstrap_threads && st.desire == 1 && budget > 1 && num_items >= max(1, threshold)
                    st.desire = min(desire_cap, 2)
                end
                if control_tail_guard && source == :control_callback && budget > 1 && num_items >= max(1, threshold)
                    stable_desire = min(desire_cap, min(budget, max(2, num_items)))
                    st.desire = max(st.desire, stable_desire)
                end
                if hint_allotment > 1
                    # Blend persisted hint with live AIMD state; this reuses past wins without hard pinning.
                    blended = max(st.desire, Int(hint_allotment))
                    st.desire = min(desire_cap, max(1, blended))
                end
            end
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
        elseif mode == :off
            false
        elseif mode == :on
            true
        elseif outer_active && !allow_with_outer
            false
        elseif heavy_only && !heavy_work
            false
        elseif adaptive_enabled
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
        hints_entries
    )
    lock(_policy_telemetry_lock) do
        ctx = _active_policy_context()
        ctx.decision_signature[source] = signature
        ctx.decision_allotment[source] = Int64(allotted)
    end
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
