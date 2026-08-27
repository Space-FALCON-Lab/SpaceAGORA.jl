Base.@kwdef mutable struct AdaptiveControllerState
    desire::Int64 = 1
    window_calls::Int64 = 0
    window_allotment_sum::Int64 = 0
    window_useful_sum::Float64 = 0.0
    window_deprived_calls::Int64 = 0
    last_classification::Symbol = :none
    last_utilization::Float64 = 1.0
end

Base.@kwdef mutable struct PolicyTelemetry
    decisions_total::Int64 = 0
    threads_enabled_total::Int64 = 0
    # Separate proposed (policy output) from dispatched (actually threaded) so that
    # route overrides — e.g. satellite_batch forcing effector threads off — are visible.
    policy_threading_proposed_total::Int64 = 0
    policy_threading_dispatched_total::Int64 = 0
    policy_discarded_by_route_total::Int64 = 0
    adaptive_decisions_total::Int64 = 0
    density_decisions::Int64 = 0
    density_threads_enabled::Int64 = 0
    control_decisions::Int64 = 0
    control_threads_enabled::Int64 = 0
    multibody_decisions::Int64 = 0
    multibody_threads_enabled::Int64 = 0
    other_decisions::Int64 = 0
    other_threads_enabled::Int64 = 0
    last_source::Symbol = :none
    last_mode::Symbol = :none
    last_threshold::Int64 = 0
    last_num_items::Int64 = 0
    last_budget::Int64 = 0
    last_outer_active::Bool = false
    last_allow_with_outer::Bool = true
    last_heavy_only::Bool = false
    last_heavy_work::Bool = true
    last_use_threads::Bool = false
    last_adaptive_enabled::Bool = false
    last_desire::Int64 = 1
    last_allotment::Int64 = 1
    observations_total::Int64 = 0
    adaptation_updates_total::Int64 = 0
    last_elapsed_ns::Int64 = 0
    elapsed_ns_total::Int64 = 0
    threaded_elapsed_ns_total::Int64 = 0
    serial_elapsed_ns_total::Int64 = 0
    last_classification::Symbol = :none
    last_utilization::Float64 = 1.0
    quantum_length::Int64 = 0
    trim_quanta_budget::Int64 = 0
    quantums_total::Int64 = 0
    quantums_inefficient::Int64 = 0
    quantums_efficient_satisfied::Int64 = 0
    quantums_efficient_deprived::Int64 = 0
    quantums_accounted_proxy::Int64 = 0
    quantums_deductible_proxy::Int64 = 0
    persistent_hints_enabled::Bool = false
    persistent_hints_loaded::Bool = false
    persistent_hints_updates::Int64 = 0
    persistent_hints_entries::Int64 = 0
    persistent_hints_path::String = ""
    last_signature::String = ""
    last_hint_allotment::Int64 = 1
    last_hint_confidence::Float64 = 0.0
    last_hint_regret_ns::Float64 = 0.0
    # Which RHS execution plan pre-solve calibration installed via
    # `rhs_plan_override`, and how it got there. Recorded by
    # `_record_rhs_plan_selection!` (simulation/engine/rhs_calibration.jl) purely
    # so the choice is observable after the fact: the plan lives in per-solve
    # `shared_buffers` and is otherwise unreachable once the solve returns, which
    # makes "the router picked the wrong route for this workload" impossible to
    # attribute from outside. Read-only accounting; nothing dispatches on these.
    #   rhs_plan_source  :none (never calibrated) | :cache | :sweep
    #   rhs_plan_mode    the plan's own mode, e.g. :satellite_batch or :flat
    rhs_plan_source::Symbol = :none
    rhs_plan_mode::Symbol = :none
    rhs_plan_allotment::Int64 = 0
    #   rhs_plan_scheduler  the plan's inner scheduler (:static/:dynamic), or
    #                       :auto when it defers to the env var. Recorded since
    #                       the scheduler became part of the calibrated plan
    #                       rather than a process-global setting.
    rhs_plan_scheduler::Symbol = :none
end

Base.@kwdef mutable struct AdaptiveChoiceStats
    samples::Int64 = 0
    successes::Int64 = 0
    failures::Int64 = 0
    elapsed_sum_ns::Float64 = 0.0
    elapsed_sq_sum_ns::Float64 = 0.0
end

Base.@kwdef mutable struct _HintLayerStatsAccumulator
    signatures::Set{String} = Set{String}()
    choice_count::Int64 = 0
    samples_total::Int64 = 0
    successes_total::Int64 = 0
    failures_total::Int64 = 0
    elapsed_sum_ns::Float64 = 0.0
    elapsed_sq_sum_ns::Float64 = 0.0
    confidence_sum::Float64 = 0.0
    regret_sum_ns::Float64 = 0.0
    signature_metric_count::Int64 = 0
end

Base.@kwdef mutable struct PolicyContext
    telemetry::PolicyTelemetry = PolicyTelemetry()
    adaptive_state::Dict{Symbol, AdaptiveControllerState} = Dict{Symbol, AdaptiveControllerState}()
    decision_signature::Dict{Symbol, String} = Dict{Symbol, String}()
    decision_allotment::Dict{Symbol, Int64} = Dict{Symbol, Int64}()
    # Per-context, not process-global, and that distinction is the whole point.
    #
    # Every field above belongs to ONE PolicyContext, and with_policy_context
    # gives each run_simulation its own (execution.jl). They were nonetheless all
    # protected by a single process-global ReentrantLock, taken on every policy
    # decision and every policy observation -- so under a thread-routed Monte
    # Carlo campaign all N concurrently running samples serialised on one lock,
    # once per RHS call each, to guard state no two of them shared.
    #
    # That was the entire residual regret of R5 against the best static route on
    # thread-routed Monte Carlo once the outer-route default was fixed. Measured
    # on independent_1sat_1hr, 64 samples, 12 threads: the ablation arm that
    # stops entering the decision/observation path at all
    # (full_smart_innermodes_off) ran -3.8% against the best static where
    # full_smart ran +22.4%, i.e. it recovered 117% of the gap, and the arm that
    # disables the adaptive policy outright recovered 142%.
    #
    # A lock per context makes the contention structural rather than incidental:
    # two samples can only contend if they genuinely share a context, which is
    # exactly when they share the state it guards.
    lock::ReentrantLock = ReentrantLock()
end

# Retained for the global fallback context and for callers outside a scoped
# context. Prefer `lock(_active_policy_context().lock)`, which is what every hot
# path now does; this exists so that the handful of places holding no context
# still have something to take.
const _policy_telemetry_lock = ReentrantLock()
const _policy_context_tls_key = :spaceagora_parallel_policy_context
const _global_policy_context = Ref{PolicyContext}(PolicyContext())
const _persistent_foreach_lock = ReentrantLock()

Base.@kwdef mutable struct _PersistentHintState
    loaded::Bool = false
    dirty::Bool = false
    path::String = ""
    history::Dict{String, Dict{Int64, AdaptiveChoiceStats}} = Dict{String, Dict{Int64, AdaptiveChoiceStats}}()
end

const _persistent_hint_lock = ReentrantLock()
const _persistent_hint_state = Ref{_PersistentHintState}(_PersistentHintState())
const _persistent_hint_atexit_registered = Ref(false)

Base.@kwdef mutable struct _PersistentForeachPool
    workers::Int
    request_channels::Vector{Channel{Any}}
    done_channel::Channel{Any}
    run_lock::ReentrantLock = ReentrantLock()
end

const _persistent_foreach_pools = Dict{Tuple{UInt, Symbol}, _PersistentForeachPool}()
const _persistent_foreach_worker_pools = Dict{Tuple{UInt, Symbol}, _PersistentForeachPool}()

# Spin-barrier pool: workers poll an atomic generation counter instead of sleeping
# on a Channel. Dispatch latency drops from ~1-5 µs (condvar wake) to ~10-50 ns
# (spin-poll check), allowing fine-grained parallelism at 32-128+ threads.
# Trade-off: worker threads burn CPU while idle. Use only for very-high-frequency
# kernels where the compute time per invocation is sub-microsecond.
mutable struct _SpinBarrierPool
    workers::Int
    # Per-worker generation: coordinator increments to signal "new work ready".
    worker_gen::Vector{Threads.Atomic{Int}}
    # Number of workers that have finished the current dispatch round.
    done_count::Threads.Atomic{Int}
    # True once the pool should stop.
    stop::Threads.Atomic{Bool}
    # Shared request payload written by coordinator before bumping generations.
    request::Base.RefValue{Any}
    # Per-worker error slot: each worker writes only its own index before
    # incrementing done_count; the coordinator reads after the barrier.
    errors::Vector{Union{Nothing, Base.CapturedException}}
    run_lock::ReentrantLock
end

function _SpinBarrierPool(workers::Int)
    return _SpinBarrierPool(
        workers,
        [Threads.Atomic{Int}(0) for _ in 1:workers],
        Threads.Atomic{Int}(0),
        Threads.Atomic{Bool}(false),
        Ref{Any}(nothing),
        Union{Nothing, Base.CapturedException}[nothing for _ in 1:workers],
        ReentrantLock(),
    )
end

const _spin_barrier_pools = Dict{Tuple{UInt, Symbol}, _SpinBarrierPool}()
const _spin_barrier_lock  = ReentrantLock()
