function record_route_discard!()
    ctx = _active_policy_context()
    lock(ctx.lock) do
        ctx.telemetry.policy_discarded_by_route_total += 1
    end
    return nothing
end

function record_policy_observation!(
    source::Symbol;
    mode::Symbol,
    num_items::Int,
    use_threads::Bool,
    elapsed_ns::Integer,
    env::Union{Nothing, PolicyDecisionEnvConfig}=nothing,
    # The context the matching decision ran in, when the caller can supply it
    # (`policy_context_hint(p)`). `nothing` resolves to the task's active
    # context, which is the shipped behaviour and is wrong on a worker task --
    # see policy_context_hint for why.
    ctx::Union{Nothing, PolicyContext}=nothing
)
    budget = env === nothing ? effective_inner_thread_budget() : env.inner_thread_budget
    ctx = ctx === nothing ? _active_policy_context() : ctx
    elapsed_ns_i64 = try
        Int64(elapsed_ns)
    catch
        typemax(Int64)
    end
    elapsed_ns_clamped = max(Int64(0), elapsed_ns_i64)
    adaptive_enabled = (mode == :auto) &&
        (env === nothing ? adaptive_policy_enabled() : env.adaptive_enabled)

    # Skip the adaptive branch when nothing can consume what it would produce.
    #
    # This is the observation-side counterpart to the short-circuit in
    # thread_policy_decision (adaptive_decision.jl), which was added on its own
    # and left this path unguarded. Everything below the `!adaptive_active`
    # return exists to move `AdaptiveControllerState.desire`, and `desire` is
    # read in exactly one place: the adaptive branch of the *next*
    # thread_policy_decision for the same source. Two cases where that read
    # cannot happen:
    #
    #   decision_forced -- budget <= 1 or num_items <= 1 makes the next decision
    #     take the forced path, which never reads desire. Every Distributed
    #     worker in a process-routed campaign runs --threads=1, so this is the
    #     common case there, not an edge case.
    #
    #   empty signature -- thread_policy_decision writes ctx.decision_signature
    #     only from its adaptive branch, so an empty entry means no adaptive
    #     decision has ever been made for this source in this context. That is
    #     the whole of a vacuum-gravity constellation solve: the RHS flat/batch
    #     routes call this once per RHS call with plan.policy_applied set, but
    #     the plan came from pre-solve calibration and no callback ever consults
    #     the policy (measured on gravity_1024sat_l50_vacuum_1hr:
    #     policy_decisions_total = 0, adaptive_decisions_total = 0, against
    #     thousands of observations per second).
    #
    # The telemetry writes above the guard are unconditional as before, so
    # observations_total and the elapsed accumulators are unchanged. Only the
    # controller update -- a dict lookup plus ~8 field writes, and for R4 four
    # more ENV reads -- is skipped, and only when its result is unreadable.
    decision_forced = budget <= 1 || num_items <= 1
    # V2 takes the static width (adaptive_decision.jl), so there is no
    # controller state to feed; the telemetry writes below still happen.
    adaptive_possible = adaptive_enabled && !decision_forced &&
        !(env !== nothing && env.policy_v2)
    # Consulted only when it pays for itself on this machine; see
    # _hint_layer_pays and hint_work_ratio.
    hints_enabled = adaptive_possible &&
        _hint_layer_pays(source, env; ctx=ctx) &&
        (env === nothing ? persistent_hints_enabled() : env.persistent_hints)
    measured_reward = hints_enabled &&
        (env === nothing ? adaptive_measured_reward_enabled() : env.adaptive_measured_reward)

    hint_signature = ""
    hint_allotment = Int64(1)
    adaptive_active = false

    lock(ctx.lock) do
        t = ctx.telemetry
        # Keep the per-source work estimate current even when the adaptive
        # branch below is skipped: the hint layer's gate reads it, and a source
        # whose estimate froze at its first observation could never be
        # reclassified. One multiply-add on a path that already holds the lock.
        st_work = get!(ctx.adaptive_state, source) do
            AdaptiveControllerState()
        end
        st_work.elapsed_ema_ns = st_work.elapsed_ema_ns <= 0.0 ?
            Float64(elapsed_ns_clamped) :
            0.8 * st_work.elapsed_ema_ns + 0.2 * Float64(elapsed_ns_clamped)
        t.observations_total += 1
        t.last_elapsed_ns = elapsed_ns_clamped
        t.elapsed_ns_total += elapsed_ns_clamped
        if use_threads
            t.threaded_elapsed_ns_total += elapsed_ns_clamped
            t.policy_threading_dispatched_total += 1
        else
            t.serial_elapsed_ns_total += elapsed_ns_clamped
        end

        hint_signature = get(ctx.decision_signature, source, "")
        hint_allotment = get(ctx.decision_allotment, source, use_threads ? Int64(max(1, min(budget, num_items))) : Int64(1))
        t.last_signature = hint_signature

        adaptive_active = adaptive_possible && !isempty(hint_signature)
        if !adaptive_active
            return nothing
        end

        if measured_reward
            # Measured-reward mode: the hint layer drives the width.
            st = _adaptive_state_for(ctx, source)
            st.desire = max(Int64(1), hint_allotment)
            t.adaptation_updates_total += 1
            t.last_desire = st.desire
        end
        # No AIMD: without measured reward the width is static (see
        # thread_policy_decision) and there is nothing to update.
        return nothing
    end

    # `adaptive_active` already implies a non-empty signature, which in turn
    # implies the decision that produced it ran its adaptive branch -- so this
    # keeps recording exactly the hint samples it recorded before. The forced
    # region drops out for free: thread_policy_decision stores "" there, so the
    # old `!isempty(hint_signature)` test already excluded it.
    if adaptive_active && hints_enabled
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
        ctx2 = ctx
        lock(ctx2.lock) do
            t = ctx2.telemetry
            t.persistent_hints_updates += 1
            t.persistent_hints_loaded = hints_loaded
            t.persistent_hints_entries = hints_entries
            t.persistent_hints_path = hints_path
        end
    end
    return nothing
end
