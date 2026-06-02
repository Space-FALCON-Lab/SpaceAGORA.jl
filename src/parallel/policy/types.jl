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
end

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
    run_lock::ReentrantLock
end

function _SpinBarrierPool(workers::Int)
    return _SpinBarrierPool(
        workers,
        [Threads.Atomic{Int}(0) for _ in 1:workers],
        Threads.Atomic{Int}(0),
        Threads.Atomic{Bool}(false),
        Ref{Any}(nothing),
        ReentrantLock(),
    )
end

const _spin_barrier_pools = Dict{Tuple{UInt, Symbol}, _SpinBarrierPool}()
const _spin_barrier_lock  = ReentrantLock()
