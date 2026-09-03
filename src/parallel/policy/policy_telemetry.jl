# Does consulting the persistent-hint layer pay for itself for this source?
#
# Compares the measured cost of one consultation against the measured work the
# decision guards. Both come from this machine; see hint_work_ratio in
# env_config.jl. A source with no observation yet is allowed through, so cold
# start behaves as it always did and the estimate converges within a few calls.
#
# Takes the context lock briefly and releases it before the caller touches the
# hint store, preserving the lock order (hint store before context) that the
# adaptive branch already relies on.
# `env` is the run-scoped snapshot when the caller has one.
#
# hint_work_ratio() is strip(get(ENV, ...)) plus a float parse -- 92.4 ns and an
# allocation of the 106.6 ns this predicate costs -- and it runs on every
# adaptive decision and every observation. Reading it live here reintroduced
# exactly the defect that "Forced-region decision: -30% time, -40% allocation"
# removed from _record_policy_decision! earlier on this branch. The snapshot
# already resolves every other per-decision knob once per run; this one now
# travels with them.
@inline function _hint_layer_pays(
    source::Symbol,
    env::Union{Nothing, PolicyDecisionEnvConfig}=nothing;
    ctx::Union{Nothing, PolicyContext}=nothing
)::Bool
    ratio = env === nothing ? hint_work_ratio() : env.hint_work_ratio
    ratio <= 0.0 && return true
    c = ctx === nothing ? _active_policy_context() : ctx
    work_ns = lock(c.lock) do
        st = get(c.adaptive_state, source, nothing)
        st === nothing ? 0.0 : st.elapsed_ema_ns
    end
    if work_ns <= 0.0
        # No work estimate yet. The shipped answer is "consult anyway", which
        # was meant to cover cold start and instead covered forever: a source
        # whose observations never reach this context (see policy_context_hint)
        # reads zero on every call and pays the lookup on every call. Under V2
        # an unmeasured region fails CLOSED -- one skipped consultation until
        # the first observation lands, against one paid consultation per call
        # until it never does.
        return !(env !== nothing && env.policy_v2)
    end
    return work_ns >= ratio * hint_overhead_ns()
end

@inline function _adaptive_state_for(ctx::PolicyContext, source::Symbol)::AdaptiveControllerState
    return get!(ctx.adaptive_state, source) do
        AdaptiveControllerState()
    end
end

@inline function _adaptive_state_for(source::Symbol)::AdaptiveControllerState
    return _adaptive_state_for(_active_policy_context(), source)
end

# `env` carries the run-scoped PolicyDecisionEnvConfig when the caller has one.
#
# This function used to call persistent_hints_enabled() directly, which is
# parse_bool_env -> lowercase(strip(get(ENV, ...))) -- a process-global ENV
# lookup that ALLOCATES a String, on every policy decision, including every
# forced-region decision the short-circuit above exists to make cheap. The
# snapshot that already resolves this exact knob once per run was threaded
# through thread_policy_decision and then not used by the one caller on the hot
# path that still read ENV.
#
# `ctx_signature`/`ctx_allotment` fold in the context stamp that
# thread_policy_decision used to perform in a SECOND lock(_policy_telemetry_lock)
# block immediately after this call returned. Two uncontended acquisitions per
# decision where one will do; the writes are to the same context object under
# the same lock, so merging them changes nothing observable.
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
    hints_entries::Int64,
    env::Union{Nothing, PolicyDecisionEnvConfig}=nothing
)
    ctx = _active_policy_context()
    lock(ctx.lock) do
        t = ctx.telemetry
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
        t.persistent_hints_enabled = (env === nothing || !policy_telemetry_uses_snapshot()) ?
            persistent_hints_enabled() : env.persistent_hints
        t.persistent_hints_loaded = hints_loaded
        t.persistent_hints_entries = hints_entries
        t.persistent_hints_path = _persistent_hint_state[].path
        t.last_signature = signature
        t.last_hint_allotment = max(1, hint_allotment)
        t.last_hint_confidence = max(0.0, hint_confidence)
        t.last_hint_regret_ns = max(0.0, hint_regret_ns)

        # Folded in from thread_policy_decision's second lock block.
        if policy_telemetry_uses_snapshot()
            ctx.decision_signature[source] = signature
            ctx.decision_allotment[source] = Int64(allotment)
        end
    end
    return nothing
end

function reset_policy_telemetry!()
    ctx = _active_policy_context()
    lock(ctx.lock) do
        ctx.telemetry = PolicyTelemetry()
        empty!(ctx.adaptive_state)
        empty!(ctx.decision_signature)
        empty!(ctx.decision_allotment)
    end
    return nothing
end

function policy_telemetry_snapshot()
    ctx = _active_policy_context()
    lock(ctx.lock) do
        t = ctx.telemetry
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
            rhs_plan_scheduler=String(t.rhs_plan_scheduler),
            accounted_fraction_proxy=accounted_fraction_proxy,
            trimmed_accounted_fraction_proxy=trimmed_accounted_fraction_proxy
        )
    end
end

"""
    record_rhs_plan_selection!(source::Symbol, mode::Symbol, allotment::Integer,
                              scheduler::Symbol=:none)

Record which RHS execution plan pre-solve calibration installed, so the choice
survives the solve and can be read back through
[`policy_telemetry_snapshot`](@ref). `source` is `:cache` when the plan came from
the calibration cache and `:sweep` when it was measured by a fresh route sweep.

Accounting only — no policy path reads these fields back.
"""
function record_rhs_plan_selection!(
    source::Symbol,
    mode::Symbol,
    allotment::Integer,
    scheduler::Symbol=:none
)
    ctx = _active_policy_context()
    lock(ctx.lock) do
        t = ctx.telemetry
        t.rhs_plan_source = source
        t.rhs_plan_mode = mode
        t.rhs_plan_allotment = Int64(max(0, allotment))
        t.rhs_plan_scheduler = scheduler
    end
    return nothing
end
