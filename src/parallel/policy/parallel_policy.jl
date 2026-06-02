module ParallelPolicy

using Base.Threads
using TOML

export parse_bool_env, parse_parallel_mode_env, parse_thread_threshold_env
export outer_parallel_active, effective_inner_thread_budget, use_threads_policy
export thread_policy_decision, threaded_foreach, threaded_reduce, threaded_foreach_persistent, threaded_foreach_worker_persistent, with_policy_context
export threaded_foreach_worker_spin, harmonics_batch_spin_barrier_enabled
export auto_thread_min_budget
export thread_worker_count
export reset_policy_telemetry!, policy_telemetry_snapshot, record_policy_observation!, record_route_discard!
export reset_persistent_hint_state!, persistent_hints_state_reset_requested
export hint_layer_stats_snapshot

include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "context.jl"))
include(joinpath(@__DIR__, "env_config.jl"))
include(joinpath(@__DIR__, "persistent_hints.jl"))
include(joinpath(@__DIR__, "policy_telemetry.jl"))
include(joinpath(@__DIR__, "adaptive_decision.jl"))
include(joinpath(@__DIR__, "thread_execution.jl"))
include(joinpath(@__DIR__, "observation_tracking.jl"))

end # module ParallelPolicy
