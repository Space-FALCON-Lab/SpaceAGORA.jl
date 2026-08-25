@inline function _adaptive_state_for(source::Symbol)::AdaptiveControllerState
    ctx = _active_policy_context()
    return get!(ctx.adaptive_state, source) do
        AdaptiveControllerState()
    end
end

function _record_policy_decision!(
    source::Symbol,
    mode::Symbol,
    threshold::Int,
    num_items::Int,
    budget::Int,
    adaptive_enabled::Bool,
    desire::Int,
    allotment::Int,
    outer_active::Bool,
    allow_with_outer::Bool,
    heavy_only::Bool,
    heavy_work::Bool,
    use_threads::Bool,
    signature::String,
    hint_allotment::Int64,
    hint_confidence::Float64,
    hint_regret_ns::Float64,
    hints_loaded::Bool,
    hints_entries::Int64
)
    lock(_policy_telemetry_lock) do
        t = _active_policy_context().telemetry
        t.decisions_total += 1
        t.threads_enabled_total += use_threads ? 1 : 0
        t.policy_threading_proposed_total += use_threads ? 1 : 0
        t.adaptive_decisions_total += adaptive_enabled ? 1 : 0

        bucket = _telemetry_bucket(source)
        if bucket == :density
            t.density_decisions += 1
            t.density_threads_enabled += use_threads ? 1 : 0
        elseif bucket == :control
            t.control_decisions += 1
            t.control_threads_enabled += use_threads ? 1 : 0
        elseif bucket == :multibody
            t.multibody_decisions += 1
            t.multibody_threads_enabled += use_threads ? 1 : 0
        else
            t.other_decisions += 1
            t.other_threads_enabled += use_threads ? 1 : 0
        end

        t.last_source = source
        t.last_mode = mode
        t.last_threshold = max(1, threshold)
        t.last_num_items = max(0, num_items)
        t.last_budget = max(1, budget)
        t.last_adaptive_enabled = adaptive_enabled
        t.last_desire = max(1, desire)
        t.last_allotment = max(1, allotment)
        t.last_outer_active = outer_active
        t.last_allow_with_outer = allow_with_outer
        t.last_heavy_only = heavy_only
        t.last_heavy_work = heavy_work
        t.last_use_threads = use_threads
        t.persistent_hints_enabled = persistent_hints_enabled()
        t.persistent_hints_loaded = hints_loaded
        t.persistent_hints_entries = hints_entries
        t.persistent_hints_path = _persistent_hint_state[].path
        t.last_signature = signature
        t.last_hint_allotment = max(1, hint_allotment)
        t.last_hint_confidence = max(0.0, hint_confidence)
        t.last_hint_regret_ns = max(0.0, hint_regret_ns)
    end
    return nothing
end

function reset_policy_telemetry!()
    lock(_policy_telemetry_lock) do
        ctx = _active_policy_context()
        ctx.telemetry = PolicyTelemetry()
        empty!(ctx.adaptive_state)
        empty!(ctx.decision_signature)
        empty!(ctx.decision_allotment)
    end
    return nothing
end

function policy_telemetry_snapshot()
    lock(_policy_telemetry_lock) do
        t = _active_policy_context().telemetry
        quantums_effective = max(0, t.quantums_total - min(t.quantums_total, t.trim_quanta_budget))
        trimmed_accounted = min(t.quantums_accounted_proxy, quantums_effective)
        accounted_fraction_proxy = t.quantums_total == 0 ? 0.0 : t.quantums_accounted_proxy / t.quantums_total
        trimmed_accounted_fraction_proxy = quantums_effective == 0 ? 0.0 : trimmed_accounted / quantums_effective
        return (
            decisions_total=t.decisions_total,
            threads_enabled_total=t.threads_enabled_total,
            policy_threading_proposed_total=t.policy_threading_proposed_total,
            policy_threading_dispatched_total=t.policy_threading_dispatched_total,
            policy_discarded_by_route_total=t.policy_discarded_by_route_total,
            adaptive_decisions_total=t.adaptive_decisions_total,
            density_decisions=t.density_decisions,
            density_threads_enabled=t.density_threads_enabled,
            control_decisions=t.control_decisions,
            control_threads_enabled=t.control_threads_enabled,
            multibody_decisions=t.multibody_decisions,
            multibody_threads_enabled=t.multibody_threads_enabled,
            other_decisions=t.other_decisions,
            other_threads_enabled=t.other_threads_enabled,
            last_source=String(t.last_source),
            last_mode=String(t.last_mode),
            last_threshold=t.last_threshold,
            last_num_items=t.last_num_items,
            last_budget=t.last_budget,
            last_adaptive_enabled=t.last_adaptive_enabled,
            last_desire=t.last_desire,
            last_allotment=t.last_allotment,
            last_outer_active=t.last_outer_active,
            last_allow_with_outer=t.last_allow_with_outer,
            last_heavy_only=t.last_heavy_only,
            last_heavy_work=t.last_heavy_work,
            last_use_threads=t.last_use_threads,
            observations_total=t.observations_total,
            adaptation_updates_total=t.adaptation_updates_total,
            last_elapsed_ns=t.last_elapsed_ns,
            elapsed_ns_total=t.elapsed_ns_total,
            threaded_elapsed_ns_total=t.threaded_elapsed_ns_total,
            serial_elapsed_ns_total=t.serial_elapsed_ns_total,
            last_classification=String(t.last_classification),
            last_utilization=t.last_utilization,
            quantum_length=t.quantum_length,
            trim_quanta_budget=t.trim_quanta_budget,
            quantums_total=t.quantums_total,
            quantums_inefficient=t.quantums_inefficient,
            quantums_efficient_satisfied=t.quantums_efficient_satisfied,
            quantums_efficient_deprived=t.quantums_efficient_deprived,
            quantums_accounted_proxy=t.quantums_accounted_proxy,
            quantums_deductible_proxy=t.quantums_deductible_proxy,
            persistent_hints_enabled=t.persistent_hints_enabled,
            persistent_hints_loaded=t.persistent_hints_loaded,
            persistent_hints_updates=t.persistent_hints_updates,
            persistent_hints_entries=t.persistent_hints_entries,
            persistent_hints_path=t.persistent_hints_path,
            last_signature=t.last_signature,
            last_hint_allotment=t.last_hint_allotment,
            last_hint_confidence=t.last_hint_confidence,
            last_hint_regret_ns=t.last_hint_regret_ns,
            rhs_plan_source=String(t.rhs_plan_source),
            rhs_plan_mode=String(t.rhs_plan_mode),
            rhs_plan_allotment=t.rhs_plan_allotment,
            accounted_fraction_proxy=accounted_fraction_proxy,
            trimmed_accounted_fraction_proxy=trimmed_accounted_fraction_proxy
        )
    end
end

"""
    record_rhs_plan_selection!(source::Symbol, mode::Symbol, allotment::Integer)

Record which RHS execution plan pre-solve calibration installed, so the choice
survives the solve and can be read back through
[`policy_telemetry_snapshot`](@ref). `source` is `:cache` when the plan came from
the calibration cache and `:sweep` when it was measured by a fresh route sweep.

Accounting only — no policy path reads these fields back.
"""
function record_rhs_plan_selection!(source::Symbol, mode::Symbol, allotment::Integer)
    lock(_policy_telemetry_lock) do
        t = _active_policy_context().telemetry
        t.rhs_plan_source = source
        t.rhs_plan_mode = mode
        t.rhs_plan_allotment = Int64(max(0, allotment))
    end
    return nothing
end
