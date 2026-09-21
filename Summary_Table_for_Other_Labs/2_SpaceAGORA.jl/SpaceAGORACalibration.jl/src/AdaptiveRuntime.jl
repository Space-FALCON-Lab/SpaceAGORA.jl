module AdaptiveRuntime

using Dates
using TOML
import SpaceAGORA

using ..Spec: CalibrationSpec
using ..Backend: BackendEvaluation

export RuntimeDecision, RuntimeController
export full_auto_enabled, init_runtime_controller
export choose_runtime_decision!, record_runtime_feedback!, save_runtime_controller!

Base.@kwdef struct RuntimeDecision
    workers::Int = 1
    batch_size::Int = 1
    outer_backend::Symbol = :none
    inner_thread_budget::Int = 1
end

Base.@kwdef mutable struct RuntimeStats
    count::Int = 0
    success::Int = 0
    fail::Int = 0
    fidelity_fail::Int = 0
    runtime_s::Float64 = 0.0
end

Base.@kwdef mutable struct RuntimeController
    enabled::Bool = false
    profile_name::String = "R0"
    machine_key::String = "machine"
    cache_path::String = ""
    warmup_batches::Int = 6
    rebalance_window::Int = 4
    fail_rate_limit::Float64 = 0.35
    fidelity_fail_rate_limit::Float64 = 0.35
    exploration_weight::Float64 = 0.35
    policy_space::Vector{RuntimeDecision} = RuntimeDecision[]
    stage_batches::Dict{String, Int} = Dict{String, Int}()
    stage_selected::Dict{String, RuntimeDecision} = Dict{String, RuntimeDecision}()
    stats::Dict{String, RuntimeStats} = Dict{String, RuntimeStats}()
    batches_recorded::Int = 0
end

@inline function full_auto_enabled(profile::SpaceAGORA.ParallelProfile)::Bool
    token = try
        lowercase(strip(SpaceAGORA.parallel_profile_name(profile)))
    catch
        ""
    end
    return token in ("r5", "r5_full_auto", "r5fullauto", "r4_full_auto", "r4fullauto")
end

@inline function _parse_positive_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return max(1, default)
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be a positive integer, got '$raw'"))
    end
    return max(1, parsed)
end

@inline function _parse_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return default
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a float, got '$raw'"))
    end
    return parsed
end

@inline function _machine_key()::String
    label = strip(get(ENV, "SPACEAGORA_CALIBRATION_MACHINE_LABEL", ""))
    isempty(label) && (label = strip(get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", "")))
    isempty(label) && (label = strip(get(ENV, "HOSTNAME", "")))
    isempty(label) && (label = strip(get(ENV, "COMPUTERNAME", "")))
    isempty(label) && (label = "machine")
    token = replace(lowercase(label), r"[^a-z0-9._-]+" => "_")
    return isempty(token) ? "machine" : token
end

@inline function _stats_key(stage::String, policy_sig::String)::String
    return string(stage, "|", policy_sig)
end

@inline function _policy_signature(policy::RuntimeDecision)::String
    return string(
        "w", policy.workers,
        "_q", policy.batch_size,
        "_b", String(policy.outer_backend),
        "_t", policy.inner_thread_budget
    )
end

@inline function _serialize_policy(policy::RuntimeDecision)::Dict{String, Any}
    return Dict(
        "workers" => policy.workers,
        "batch_size" => policy.batch_size,
        "outer_backend" => String(policy.outer_backend),
        "inner_thread_budget" => policy.inner_thread_budget
    )
end

@inline function _deserialize_policy(row)::RuntimeDecision
    return RuntimeDecision(
        workers=Int(get(row, "workers", 1)),
        batch_size=Int(get(row, "batch_size", 1)),
        outer_backend=Symbol(String(get(row, "outer_backend", "none"))),
        inner_thread_budget=Int(get(row, "inner_thread_budget", 1))
    )
end

@inline function _thread_capacity()::Int
    return max(1, Threads.nthreads())
end

function _build_policy_space(base_parallel::Int, base_batch_size::Int)::Vector{RuntimeDecision}
    max_workers = max(1, base_parallel)
    workers_opts = sort!(unique([
        1,
        min(max_workers, 2),
        min(max_workers, 4),
        max_workers
    ]))

    batch_opts_base = sort!(unique([
        1,
        min(max_workers, 2),
        min(max_workers, 4),
        max(1, base_batch_size)
    ]))

    thread_cap = _thread_capacity()
    raw = RuntimeDecision[]
    for workers in workers_opts
        backends = workers <= 1 ? (:none,) : (:process, :threads)
        batch_opts = sort!(unique(min(workers, q) for q in batch_opts_base))
        inner_cap = max(1, fld(thread_cap, workers))
        inner_opts = sort!(unique([1, min(inner_cap, 2), min(inner_cap, 4), inner_cap]))
        for backend in backends
            for batch_size in batch_opts
                for inner_budget in inner_opts
                    push!(raw, RuntimeDecision(
                        workers=workers,
                        batch_size=batch_size,
                        outer_backend=backend,
                        inner_thread_budget=inner_budget
                    ))
                end
            end
        end
    end

    out = RuntimeDecision[]
    seen = Set{String}()
    for policy in raw
        sig = _policy_signature(policy)
        sig in seen && continue
        push!(seen, sig)
        push!(out, policy)
    end
    return out
end

@inline function _cache_path(spec::CalibrationSpec, profile_name::String, machine_key::String)::String
    root = joinpath(spec.output_root, "runtime_policy_cache")
    mkpath(root)
    return joinpath(root, string(machine_key, "_", lowercase(profile_name), ".toml"))
end

function _load_cache!(controller::RuntimeController)::Nothing
    isfile(controller.cache_path) || return nothing
    parsed = try
        TOML.parsefile(controller.cache_path)
    catch
        return nothing
    end

    cached_profile = String(get(parsed, "profile_name", ""))
    cached_machine = String(get(parsed, "machine_key", ""))
    if cached_profile != controller.profile_name || cached_machine != controller.machine_key
        return nothing
    end

    rows = get(parsed, "stats", Any[])
    for row in rows
        row isa AbstractDict || continue
        stage = String(get(row, "stage", ""))
        sig = String(get(row, "policy_signature", ""))
        isempty(stage) && continue
        isempty(sig) && continue

        key = _stats_key(stage, sig)
        controller.stats[key] = RuntimeStats(
            count=max(0, Int(get(row, "count", 0))),
            success=max(0, Int(get(row, "success", 0))),
            fail=max(0, Int(get(row, "fail", 0))),
            fidelity_fail=max(0, Int(get(row, "fidelity_fail", 0))),
            runtime_s=max(0.0, Float64(get(row, "runtime_s", 0.0)))
        )
    end
    return nothing
end

function save_runtime_controller!(controller::RuntimeController)::Nothing
    controller.enabled || return nothing
    isempty(controller.cache_path) && return nothing

    rows = Dict{String, Any}[]
    for (key, stats) in sort(collect(controller.stats); by=first)
        sep = findfirst('|', key)
        sep === nothing && continue
        stage = key[begin:prevind(key, sep)]
        sig = key[nextind(key, sep):end]
        push!(rows, Dict(
            "stage" => stage,
            "policy_signature" => sig,
            "count" => stats.count,
            "success" => stats.success,
            "fail" => stats.fail,
            "fidelity_fail" => stats.fidelity_fail,
            "runtime_s" => stats.runtime_s
        ))
    end

    payload = Dict{String, Any}(
        "schema_version" => 1,
        "profile_name" => controller.profile_name,
        "machine_key" => controller.machine_key,
        "updated_utc" => string(now(UTC)),
        "warmup_batches" => controller.warmup_batches,
        "rebalance_window" => controller.rebalance_window,
        "fail_rate_limit" => controller.fail_rate_limit,
        "fidelity_fail_rate_limit" => controller.fidelity_fail_rate_limit,
        "exploration_weight" => controller.exploration_weight,
        "policy_space" => [_serialize_policy(policy) for policy in controller.policy_space],
        "stats" => rows
    )

    mkpath(dirname(controller.cache_path))
    open(controller.cache_path, "w") do io
        TOML.print(io, payload)
    end
    return nothing
end

function init_runtime_controller(
    spec::CalibrationSpec,
    profile::SpaceAGORA.ParallelProfile;
    base_parallel::Int,
    base_batch_size::Int=1,
    full_auto_requested::Bool=false
)::Union{Nothing, RuntimeController}
    (full_auto_requested || full_auto_enabled(profile)) || return nothing

    profile_name = full_auto_requested ? "R5" : SpaceAGORA.parallel_profile_name(profile)
    machine_key = _machine_key()
    policy_space = _build_policy_space(base_parallel, base_batch_size)
    isempty(policy_space) && return nothing

    controller = RuntimeController(
        enabled=true,
        profile_name=profile_name,
        machine_key=machine_key,
        cache_path=_cache_path(spec, profile_name, machine_key),
        warmup_batches=_parse_positive_int_env("SPACEAGORA_CALIBRATION_AUTO_WARMUP_BATCHES", min(length(policy_space), 6)),
        rebalance_window=_parse_positive_int_env("SPACEAGORA_CALIBRATION_AUTO_REBALANCE_WINDOW", 4),
        fail_rate_limit=clamp(_parse_float_env("SPACEAGORA_CALIBRATION_AUTO_FAIL_RATE_LIMIT", 0.35), 0.0, 1.0),
        fidelity_fail_rate_limit=clamp(_parse_float_env("SPACEAGORA_CALIBRATION_AUTO_FIDELITY_FAIL_RATE_LIMIT", 0.35), 0.0, 1.0),
        exploration_weight=max(0.0, _parse_float_env("SPACEAGORA_CALIBRATION_AUTO_EXPLORATION_WEIGHT", 0.35)),
        policy_space=policy_space
    )
    _load_cache!(controller)
    return controller
end

@inline function _policy_stats(controller::RuntimeController, stage::String, policy::RuntimeDecision)::RuntimeStats
    return get(controller.stats, _stats_key(stage, _policy_signature(policy)), RuntimeStats())
end

function _stage_total_observations(controller::RuntimeController, stage::String)::Int
    prefix = string(stage, "|")
    total = 0
    for (key, stats) in controller.stats
        startswith(key, prefix) || continue
        total += stats.count
    end
    return total
end

function _policy_score(
    controller::RuntimeController,
    stage::String,
    policy::RuntimeDecision,
    total_obs::Int
)::Float64
    stats = _policy_stats(controller, stage, policy)
    if stats.count == 0
        # Prefer unexplored options during and after warmup.
        return 1.0e6
    end

    fail_rate = stats.fail / max(stats.count, 1)
    fidelity_fail_rate = stats.fidelity_fail / max(stats.success, 1)
    runtime_s = max(stats.runtime_s, 1.0e-6 * stats.count)
    throughput = stats.success / runtime_s

    constrained_ok = fail_rate <= controller.fail_rate_limit &&
                     fidelity_fail_rate <= controller.fidelity_fail_rate_limit
    constrained_ok || (throughput *= 0.25)

    explore = controller.exploration_weight * sqrt(log(total_obs + 1.0) / (stats.count + 1.0))
    return throughput + explore
end

function _select_policy(controller::RuntimeController, stage::String)::RuntimeDecision
    total_obs = _stage_total_observations(controller, stage)
    best_policy = controller.policy_space[1]
    best_score = -Inf
    for policy in controller.policy_space
        score = _policy_score(controller, stage, policy, total_obs)
        if score > best_score
            best_score = score
            best_policy = policy
        end
    end
    return best_policy
end

function choose_runtime_decision!(
    controller::RuntimeController,
    stage::String,
    remaining::Int
)::RuntimeDecision
    controller.enabled || return RuntimeDecision()
    remaining = max(1, remaining)

    stage_batch = get(controller.stage_batches, stage, 0) + 1
    controller.stage_batches[stage] = stage_batch

    chosen = if stage_batch <= controller.warmup_batches
        idx = ((stage_batch - 1) % length(controller.policy_space)) + 1
        controller.policy_space[idx]
    else
        if !haskey(controller.stage_selected, stage) ||
           ((stage_batch - controller.warmup_batches - 1) % max(1, controller.rebalance_window) == 0)
            controller.stage_selected[stage] = _select_policy(controller, stage)
        end
        controller.stage_selected[stage]
    end

    workers = max(1, min(chosen.workers, remaining))
    batch_size = max(1, min(chosen.batch_size, remaining))
    outer_backend = workers <= 1 ? :none : chosen.outer_backend
    return RuntimeDecision(
        workers=workers,
        batch_size=batch_size,
        outer_backend=outer_backend,
        inner_thread_budget=max(1, chosen.inner_thread_budget)
    )
end

function record_runtime_feedback!(
    controller::RuntimeController,
    stage::String,
    decision::RuntimeDecision,
    evaluations::Vector{BackendEvaluation}
)::Nothing
    controller.enabled || return nothing
    isempty(evaluations) && return nothing

    key = _stats_key(stage, _policy_signature(decision))
    stats = get!(controller.stats, key, RuntimeStats())
    for ev in evaluations
        ok = ev.success && isfinite(ev.score)
        stats.count += 1
        stats.success += ok ? 1 : 0
        stats.fail += ok ? 0 : 1
        fidelity_ok = get(ev.metrics, "all_pass", 1.0) >= 0.5
        if ok && !fidelity_ok
            stats.fidelity_fail += 1
        end
        stats.runtime_s += max(ev.runtime_s, 1.0e-6)
    end

    controller.batches_recorded += 1
    if controller.batches_recorded % max(1, controller.rebalance_window) == 0
        save_runtime_controller!(controller)
    end
    return nothing
end

end # module AdaptiveRuntime
