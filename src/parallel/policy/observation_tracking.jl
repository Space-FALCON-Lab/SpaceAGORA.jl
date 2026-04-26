function record_policy_observation!(
    source::Symbol;
    mode::Symbol,
    num_items::Int,
    use_threads::Bool,
    elapsed_ns::Integer
)
    budget = effective_inner_thread_budget()
    elapsed_ns_i64 = try
        Int64(elapsed_ns)
    catch
        typemax(Int64)
    end
    elapsed_ns_clamped = max(Int64(0), elapsed_ns_i64)
    adaptive_enabled = (mode == :auto) && adaptive_policy_enabled()
    measured_reward = adaptive_enabled && persistent_hints_enabled() && adaptive_measured_reward_enabled()
    hint_signature = ""
    hint_allotment = Int64(1)

    lock(_policy_telemetry_lock) do
        ctx = _active_policy_context()
        t = ctx.telemetry
        t.observations_total += 1
        t.last_elapsed_ns = elapsed_ns_clamped
        t.elapsed_ns_total += elapsed_ns_clamped
        if use_threads
            t.threaded_elapsed_ns_total += elapsed_ns_clamped
        else
            t.serial_elapsed_ns_total += elapsed_ns_clamped
        end

        hint_signature = get(ctx.decision_signature, source, "")
        hint_allotment = get(ctx.decision_allotment, source, use_threads ? Int64(max(1, min(budget, num_items))) : Int64(1))
        t.last_signature = hint_signature

        if !adaptive_enabled
            return nothing
        end

        st = _adaptive_state_for(source)
        if measured_reward
            st.desire = max(Int64(1), hint_allotment)
            st.last_classification = :measured_reward
            st.last_utilization = use_threads ? 1.0 : 0.0
            st.window_calls = 0
            st.window_allotment_sum = 0
            st.window_useful_sum = 0.0
            st.window_deprived_calls = 0

            t.adaptation_updates_total += 1
            t.last_classification = st.last_classification
            t.last_utilization = st.last_utilization
            t.last_desire = st.desire
            t.quantum_length = 1
            t.trim_quanta_budget = 0
            t.quantums_total += 1
            t.quantums_accounted_proxy += 1
            return nothing
        end

        ρ = adaptive_rho()
        δ = adaptive_delta()
        L = adaptive_window_size()
        trim_quanta = adaptive_trim_quanta_budget()

        st.desire = min(max(1, st.desire), _adaptive_desire_cap(budget, ρ))
        allotment = use_threads ? max(1, min(st.desire, budget)) : 1
        useful = min(max(0, num_items), allotment)

        st.window_calls += 1
        st.window_allotment_sum += allotment
        st.window_useful_sum += useful
        st.window_deprived_calls += (allotment < st.desire) ? 1 : 0

        if st.window_calls < L
            return nothing
        end

        utilization = st.window_useful_sum / max(1, st.window_allotment_sum)
        efficient = utilization >= δ
        if !efficient
            st.last_classification = :inefficient
            st.desire = max(1, floor(Int, st.desire / ρ))
        elseif st.window_deprived_calls == 0
            st.last_classification = :efficient_satisfied
            st.desire = min(_adaptive_desire_cap(budget, ρ), max(1, ceil(Int, st.desire * ρ)))
        else
            st.last_classification = :efficient_deprived
        end

        st.last_utilization = utilization
        st.window_calls = 0
        st.window_allotment_sum = 0
        st.window_useful_sum = 0.0
        st.window_deprived_calls = 0

        t.adaptation_updates_total += 1
        t.last_classification = st.last_classification
        t.last_utilization = st.last_utilization
        t.last_desire = st.desire
        t.quantum_length = L
        t.trim_quanta_budget = trim_quanta
        t.quantums_total += 1
        if st.last_classification == :inefficient
            t.quantums_inefficient += 1
            t.quantums_deductible_proxy += 1
        elseif st.last_classification == :efficient_satisfied
            t.quantums_efficient_satisfied += 1
            t.quantums_deductible_proxy += 1
        elseif st.last_classification == :efficient_deprived
            t.quantums_efficient_deprived += 1
            t.quantums_accounted_proxy += 1
        else
            t.quantums_deductible_proxy += 1
        end
    end

    if adaptive_enabled && persistent_hints_enabled() && !isempty(hint_signature)
        _hint_record_observation!(
            hint_signature,
            max(Int64(1), hint_allotment),
            elapsed_ns_clamped;
            success=(elapsed_ns_clamped > 0)
        )
        hints_loaded = false
        hints_entries = 0
        hints_path = ""
        lock(_persistent_hint_lock) do
            hints_loaded = _persistent_hint_state[].loaded
            hints_entries = _hint_entry_count(_persistent_hint_state[])
            hints_path = _persistent_hint_state[].path
        end
        lock(_policy_telemetry_lock) do
            t = _active_policy_context().telemetry
            t.persistent_hints_updates += 1
            t.persistent_hints_loaded = hints_loaded
            t.persistent_hints_entries = hints_entries
            t.persistent_hints_path = hints_path
        end
    end
    return nothing
end
