module ParallelProfiles

using Dates
using TOML

export ParallelProfile, ParallelProfileConfig
export parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
export OuterRouteFeatures, OuterRouteTuning, OuterRouteState
export reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
export default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
export load_outer_route_state!, save_outer_route_state

@enum ParallelProfile begin
    R0
    R1_a
    R1_b
    R2
    R3
    R4
    R5
end

# Backward-compatible alias for historical profile naming.
const R4_full_auto = R5

Base.@kwdef struct ParallelProfileConfig
    profile::ParallelProfile
    label::String
    outer_backend::Symbol
    inner_adaptive::Bool
    outer_route_adaptive::Bool
    density_mode::String
    control_mode::String
    thermal_mode::String
    multibody_mode::String
    effector_mode::String
end

@inline function parallel_profile_name(profile::ParallelProfile)::String
    if profile == R0
        return "R0"
    elseif profile == R1_a
        return "R1_a"
    elseif profile == R1_b
        return "R1_b"
    elseif profile == R2
        return "R2"
    elseif profile == R3
        return "R3"
    elseif profile == R4
        return "R4"
    end
    return "R5"
end

@inline function _normalize_profile_token(raw::AbstractString)::String
    token = lowercase(strip(String(raw)))
    token = replace(token, "-" => "_")
    return replace(token, " " => "")
end

function parse_parallel_profile(raw::AbstractString)::ParallelProfile
    token = _normalize_profile_token(raw)
    if token in ("r0", "serial", "r0_true_serial", "true_serial")
        return R0
    elseif token in ("r1_a", "r1a", "r1_a_outer_only", "outer_only", "threads", "outer_only_threads")
        return R1_a
    elseif token in ("r1_b", "r1b", "r1_b_outer_only_process", "outer_only_process", "process_outer")
        return R1_b
    elseif token in ("r2", "r2_inner_only", "inner_only")
        return R2
    elseif token in ("r3", "r3_outer_inner_static", "outer_inner_static", "process_static")
        return R3
    elseif token in ("r4", "r4_outer_inner_adaptive", "outer_inner_adaptive", "auto_adaptive")
        return R4
    elseif token in (
        "r5",
        "r5_full_auto",
        "r5fullauto",
        "r5_calibration_full_auto",
        "full_smart",
        # Legacy aliases:
        "r4_full_auto",
        "r4fullauto",
        "r4_calibration_full_auto"
    )
        return R5
    end
    throw(ArgumentError(
        "Unsupported parallel profile '$raw'. Use one of: R0, R1_a, R1_b, R2, R3, R4, R5."
    ))
end

parse_parallel_profile(raw::Symbol)::ParallelProfile = parse_parallel_profile(String(raw))
parse_parallel_profile(profile::ParallelProfile)::ParallelProfile = profile

function profile_config(profile_in)::ParallelProfileConfig
    profile = parse_parallel_profile(profile_in)
    if profile == R0
        return ParallelProfileConfig(
            profile=profile,
            label="r0_true_serial",
            outer_backend=:none,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="off",
            control_mode="off",
            thermal_mode="off",
            multibody_mode="off",
            effector_mode="off"
        )
    elseif profile == R1_a
        return ParallelProfileConfig(
            profile=profile,
            label="r1_a_outer_only",
            outer_backend=:threads,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="off",
            control_mode="off",
            thermal_mode="off",
            multibody_mode="off",
            effector_mode="off"
        )
    elseif profile == R1_b
        return ParallelProfileConfig(
            profile=profile,
            label="r1_b_outer_only_process",
            outer_backend=:process,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="off",
            control_mode="off",
            thermal_mode="off",
            multibody_mode="off",
            effector_mode="off"
        )
    elseif profile == R2
        return ParallelProfileConfig(
            profile=profile,
            label="r2_inner_only",
            outer_backend=:none,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="auto",
            control_mode="auto",
            thermal_mode="auto",
            multibody_mode="auto",
            effector_mode="auto"
        )
    elseif profile == R3
        return ParallelProfileConfig(
            profile=profile,
            label="r3_outer_inner_static",
            outer_backend=:auto,
            inner_adaptive=false,
            outer_route_adaptive=false,
            density_mode="auto",
            control_mode="auto",
            thermal_mode="auto",
            multibody_mode="auto",
            effector_mode="auto"
        )
    elseif profile == R4
        return ParallelProfileConfig(
            profile=profile,
            label="r4_outer_inner_adaptive",
            outer_backend=:auto,
            inner_adaptive=true,
            outer_route_adaptive=true,
            density_mode="auto",
            control_mode="auto",
            thermal_mode="auto",
            multibody_mode="auto",
            effector_mode="auto"
        )
    end
    return ParallelProfileConfig(
        profile=profile,
        label="r5",
        outer_backend=:auto,
        inner_adaptive=true,
        outer_route_adaptive=true,
        density_mode="auto",
        control_mode="auto",
        thermal_mode="auto",
        multibody_mode="auto",
        effector_mode="auto"
    )
end

@inline function _coerce_env_bool(v::Bool)::String
    return v ? "1" : "0"
end

@inline function _profile_outer_backend_token(backend::Symbol)::String
    if backend == :none
        return "none"
    elseif backend == :threads
        return "threads"
    elseif backend == :process
        return "process"
    elseif backend == :auto
        return "auto"
    end
    throw(ArgumentError("Unsupported outer backend '$backend'."))
end

@inline function _env_or_default(name::String, fallback::String; preserve_existing::Bool)::String
    if !preserve_existing
        return fallback
    end
    existing = strip(get(ENV, name, ""))
    return isempty(existing) ? fallback : existing
end

function profile_env_pairs(
    profile_in;
    preserve_existing::Bool=true,
    outer_parallel_active::Bool=false
)::Vector{Pair{String, String}}
    cfg = profile_config(profile_in)
    return Pair{String, String}[
        "SPACEAGORA_PARALLEL_PROFILE" => _env_or_default(
            "SPACEAGORA_PARALLEL_PROFILE",
            parallel_profile_name(cfg.profile);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => _env_or_default(
            "SPACEAGORA_PERF_PARALLEL_BACKEND",
            _profile_outer_backend_token(cfg.outer_backend);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE" => _env_or_default(
            "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE",
            _coerce_env_bool(cfg.outer_route_adaptive);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _env_or_default(
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE",
            _coerce_env_bool(outer_parallel_active);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE",
            _coerce_env_bool(cfg.inner_adaptive);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => _env_or_default(
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL",
            cfg.density_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => _env_or_default(
            "SPACEAGORA_CONTROL_CALLBACK_PARALLEL",
            cfg.control_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => _env_or_default(
            "SPACEAGORA_THERMAL_CALLBACK_PARALLEL",
            cfg.thermal_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_MULTIBODY_PARALLEL" => _env_or_default(
            "SPACEAGORA_MULTIBODY_PARALLEL",
            cfg.multibody_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_EFFECTOR_PARALLEL" => _env_or_default(
            "SPACEAGORA_EFFECTOR_PARALLEL",
            cfg.effector_mode;
            preserve_existing=preserve_existing
        )
    ]
end

function with_parallel_profile(
    f::Function,
    profile_in;
    preserve_existing::Bool=true,
    outer_parallel_active::Bool=false
)
    env_pairs = profile_env_pairs(
        profile_in;
        preserve_existing=preserve_existing,
        outer_parallel_active=outer_parallel_active
    )
    return withenv(env_pairs...) do
        f()
    end
end

function with_parallel_profile(
    profile_in,
    f::Function;
    preserve_existing::Bool=true,
    outer_parallel_active::Bool=false
)
    return with_parallel_profile(
        f,
        profile_in;
        preserve_existing=preserve_existing,
        outer_parallel_active=outer_parallel_active
    )
end

Base.@kwdef struct OuterRouteFeatures
    category::String = "deterministic"
    n_sats::Int = 1
    n_links::Int = 1
    mission_time_s::Float64 = 0.0
    has_nbody::Bool = false
    has_srp::Bool = false
    harmonics_degree::Int = 0
    has_control::Bool = false
    orientation_on::Bool = false
    montecarlo_samples::Int = 0
end

Base.@kwdef struct OuterRouteTuning
    inner_sat_threshold::Int = 8
    inner_link_threshold::Int = 12
    outer_light_sat_threshold::Int = 2
    outer_light_link_threshold::Int = 4
    outer_light_mission_threshold_s::Float64 = 14_400.0
    spice_constellation_process_enabled::Bool = true
    spice_constellation_min_sats::Int = 4
    adaptive_enabled::Bool = true
    adaptive_min_samples::Int = 2
    failure_penalty_s::Float64 = 120.0
    mc_process_min_samples::Int = 16
    mc_process_min_mission_s::Float64 = 3600.0
    trace::Bool = false
end

Base.@kwdef mutable struct OuterRouteStats
    samples::Int = 0
    successes::Int = 0
    failures::Int = 0
    elapsed_sum_s::Float64 = 0.0
end

Base.@kwdef mutable struct OuterRouteState
    lock::ReentrantLock = ReentrantLock()
    history::Dict{String, Dict{Symbol, OuterRouteStats}} = Dict{String, Dict{Symbol, OuterRouteStats}}()
end

function reset_outer_route_state!(state::OuterRouteState)
    lock(state.lock) do
        empty!(state.history)
    end
    return nothing
end

@inline function _route_stats_payload(stats::OuterRouteStats)::Dict{String, Any}
    return Dict{String, Any}(
        "samples" => max(0, Int(stats.samples)),
        "successes" => max(0, Int(stats.successes)),
        "failures" => max(0, Int(stats.failures)),
        "elapsed_sum_s" => max(0.0, Float64(stats.elapsed_sum_s))
    )
end

@inline function _route_payload_stats(payload)::Union{Nothing, OuterRouteStats}
    payload isa AbstractDict || return nothing
    samples = try
        max(0, Int(get(payload, "samples", 0)))
    catch
        0
    end
    successes = try
        max(0, Int(get(payload, "successes", 0)))
    catch
        0
    end
    failures = try
        max(0, Int(get(payload, "failures", 0)))
    catch
        0
    end
    elapsed_sum_s = try
        max(0.0, Float64(get(payload, "elapsed_sum_s", 0.0)))
    catch
        0.0
    end
    samples <= 0 && return nothing
    successes = min(samples, successes)
    failures = min(samples - successes, failures)
    elapsed_sum_s = max(0.0, elapsed_sum_s)
    return OuterRouteStats(
        samples=samples,
        successes=successes,
        failures=failures,
        elapsed_sum_s=elapsed_sum_s
    )
end

function save_outer_route_state(
    state::OuterRouteState,
    path::AbstractString;
    metadata::AbstractDict=Dict{String, Any}()
)::NamedTuple{(:path, :signatures, :rows), Tuple{String, Int, Int}}
    path_s = String(path)
    rows = Dict{String, Any}[]
    signatures = 0
    lock(state.lock) do
        for signature in sort!(collect(keys(state.history)))
            bucket = state.history[signature]
            isempty(bucket) && continue
            signature_rows = 0
            for route in (:none, :threads, :process)
                stats = get(bucket, route, nothing)
                stats isa OuterRouteStats || continue
                stats.samples > 0 || continue
                push!(rows, Dict{String, Any}(
                    "signature" => signature,
                    "route" => String(route),
                    "stats" => _route_stats_payload(stats)
                ))
                signature_rows += 1
            end
            signature_rows > 0 && (signatures += 1)
        end
    end

    metadata_out = Dict{String, Any}()
    for (k, v) in metadata
        key = String(k)
        metadata_out[key] = if v isa Number || v isa Bool
            v
        else
            String(v)
        end
    end

    payload = Dict{String, Any}(
        "schema_version" => 1,
        "updated_utc" => string(now(UTC)),
        "metadata" => metadata_out,
        "history" => rows
    )

    mkpath(dirname(path_s))
    open(path_s, "w") do io
        TOML.print(io, payload)
    end
    return (path=path_s, signatures=signatures, rows=length(rows))
end

function load_outer_route_state!(
    state::OuterRouteState,
    path::AbstractString;
    replace::Bool=true
)::NamedTuple{(:path, :signatures, :rows), Tuple{String, Int, Int}}
    path_s = String(path)
    isfile(path_s) || return (path=path_s, signatures=0, rows=0)
    parsed = TOML.parsefile(path_s)
    rows_in = get(parsed, "history", Any[])
    rows_in isa AbstractVector || return (path=path_s, signatures=0, rows=0)

    loaded_rows = 0
    loaded_signatures = Set{String}()
    lock(state.lock) do
        replace && empty!(state.history)
        for row in rows_in
            row isa AbstractDict || continue
            signature = strip(String(get(row, "signature", "")))
            isempty(signature) && continue
            route = Symbol(String(get(row, "route", "")))
            route in (:none, :threads, :process) || continue
            stats = _route_payload_stats(get(row, "stats", nothing))
            stats === nothing && continue

            bucket = get!(state.history, signature) do
                Dict{Symbol, OuterRouteStats}()
            end
            existing = get!(bucket, route) do
                OuterRouteStats()
            end
            existing.samples += stats.samples
            existing.successes += stats.successes
            existing.failures += stats.failures
            existing.elapsed_sum_s += stats.elapsed_sum_s
            loaded_rows += 1
            push!(loaded_signatures, signature)
        end
    end
    return (path=path_s, signatures=length(loaded_signatures), rows=loaded_rows)
end

@inline function _threads_or_none(threads_available::Bool)::Symbol
    return threads_available ? :threads : :none
end

@inline function _route_sat_bucket(n_sats::Int)::String
    if n_sats <= 1
        return "1"
    elseif n_sats <= 2
        return "2"
    elseif n_sats <= 4
        return "3_4"
    end
    return "5p"
end

@inline function _route_link_bucket(n_links::Int)::String
    if n_links <= 1
        return "1"
    elseif n_links <= 4
        return "2_4"
    elseif n_links <= 8
        return "5_8"
    end
    return "9p"
end

@inline function _route_mission_bucket(mission_time_s::Float64)::String
    if mission_time_s <= 1800.0
        return "short"
    elseif mission_time_s <= 7200.0
        return "medium"
    end
    return "long"
end

@inline function _route_harmonics_bucket(L::Int)::String
    if L <= 0
        return "0"
    elseif L <= 10
        return "1_10"
    elseif L <= 20
        return "11_20"
    end
    return "21p"
end

@inline function outer_route_signature(f::OuterRouteFeatures)::String
    return join((
        "cat=$(f.category)",
        "sat=$(_route_sat_bucket(f.n_sats))",
        "links=$(_route_link_bucket(f.n_links))",
        "mission=$(_route_mission_bucket(f.mission_time_s))",
        "nbody=$(f.has_nbody ? "1" : "0")",
        "srp=$(f.has_srp ? "1" : "0")",
        "harm=$(_route_harmonics_bucket(f.harmonics_degree))",
        "ctrl=$(f.has_control ? "1" : "0")",
        "orient=$(f.orientation_on ? "1" : "0")"
    ), "|")
end

function outer_route_stats_snapshot(
    state::OuterRouteState,
    signature::String
)::Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}}
    lock(state.lock) do
        entry = get(state.history, signature, nothing)
        snap = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}}()
        if entry === nothing
            return snap
        end
        for (route, stats) in entry
            if stats.samples <= 0
                continue
            end
            success_rate = stats.successes / max(1, stats.samples)
            snap[route] = (samples=stats.samples, mean_s=stats.elapsed_sum_s / stats.samples, success_rate=success_rate)
        end
        return snap
    end
end

@inline function _feature_is_lightweight(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool
    return f.n_sats <= t.outer_light_sat_threshold &&
        f.n_links <= t.outer_light_link_threshold &&
        f.mission_time_s <= t.outer_light_mission_threshold_s &&
        !f.has_nbody &&
        !f.has_control &&
        !f.orientation_on &&
        f.harmonics_degree == 0
end

@inline function _feature_heavy_for_process(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool
    if _feature_is_lightweight(f, t)
        return false
    end
    return f.has_nbody || f.harmonics_degree >= 20 || f.mission_time_s > t.outer_light_mission_threshold_s
end

@inline function _priority_outer_route_montecarlo(
    f::OuterRouteFeatures,
    t::OuterRouteTuning;
    machine_class::Symbol,
    threads_available::Bool,
    parallel_enabled::Bool
)::Symbol
    if !parallel_enabled
        return :none
    end
    if f.montecarlo_samples <= 1
        return :none
    end
    if machine_class in (:large, :medium) &&
       (f.montecarlo_samples >= t.mc_process_min_samples ||
        f.mission_time_s >= t.mc_process_min_mission_s)
        return :process
    end
    return _threads_or_none(threads_available)
end

function default_outer_route(
    f::OuterRouteFeatures;
    tuning::OuterRouteTuning=OuterRouteTuning(),
    machine_class::Symbol=:small,
    threads_available::Bool=true,
    parallel_enabled::Bool=true
)::Symbol
    if !parallel_enabled
        return :none
    end

    if lowercase(strip(f.category)) == "montecarlo"
        return _priority_outer_route_montecarlo(
            f,
            tuning;
            machine_class=machine_class,
            threads_available=threads_available,
            parallel_enabled=parallel_enabled
        )
    end

    if _feature_is_lightweight(f, tuning)
        return :none
    end

    if f.n_sats >= tuning.inner_sat_threshold || f.n_links >= tuning.inner_link_threshold
        return :none
    end

    if tuning.spice_constellation_process_enabled &&
       f.n_sats >= tuning.spice_constellation_min_sats &&
       (f.has_nbody || f.has_srp)
        return machine_class in (:large, :medium) ? :process : _threads_or_none(threads_available)
    end

    if lowercase(strip(f.category)) == "satellite_scaling" && f.n_sats >= 4
        return _threads_or_none(threads_available)
    end

    if f.has_nbody || f.harmonics_degree >= 20
        return machine_class in (:large, :medium) ? :process : _threads_or_none(threads_available)
    end

    if machine_class in (:large, :medium)
        return :process
    end
    return _threads_or_none(threads_available)
end

function outer_route_candidates(
    f::OuterRouteFeatures;
    tuning::OuterRouteTuning=OuterRouteTuning(),
    machine_class::Symbol=:small,
    threads_available::Bool=true,
    parallel_enabled::Bool=true
)::Vector{Symbol}
    candidates = Symbol[:none]
    if threads_available
        push!(candidates, :threads)
    end
    allow_process = if lowercase(strip(f.category)) == "montecarlo"
        _priority_outer_route_montecarlo(
            f,
            tuning;
            machine_class=machine_class,
            threads_available=threads_available,
            parallel_enabled=parallel_enabled
        ) == :process
    else
        _feature_heavy_for_process(f, tuning)
    end
    if allow_process
        push!(candidates, :process)
    end
    return unique(candidates)
end

@inline function _route_ranked_candidates(candidates::Vector{Symbol}, default_route::Symbol)::Vector{Symbol}
    ranked = Symbol[]
    if default_route in candidates
        push!(ranked, default_route)
    end
    for route in (:threads, :none, :process)
        if route in candidates && !(route in ranked)
            push!(ranked, route)
        end
    end
    for route in candidates
        if !(route in ranked)
            push!(ranked, route)
        end
    end
    return ranked
end

@inline function _under_sampled_candidate(
    candidates::Vector{Symbol},
    snapshot::Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}},
    default_route::Symbol,
    min_samples::Int
)::Union{Nothing, Symbol}
    ranked = _route_ranked_candidates(candidates, default_route)
    for route in ranked
        info = get(snapshot, route, (samples=0, mean_s=Inf, success_rate=0.0))
        if info.samples < max(1, min_samples)
            return route
        end
    end
    return nothing
end

@inline function _best_candidate(
    candidates::Vector{Symbol},
    snapshot::Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}}
)::Union{Nothing, Symbol}
    best_route = nothing
    best_mean = Inf
    best_success_rate = -Inf
    for route in candidates
        info = get(snapshot, route, (samples=0, mean_s=Inf, success_rate=0.0))
        if info.samples <= 0 || !isfinite(info.mean_s)
            continue
        end
        if info.mean_s < best_mean - 1e-12
            best_route = route
            best_mean = info.mean_s
            best_success_rate = info.success_rate
        elseif isapprox(info.mean_s, best_mean; atol=1e-12, rtol=0.0) && info.success_rate > best_success_rate
            best_route = route
            best_success_rate = info.success_rate
        end
    end
    return best_route
end

function select_outer_route!(
    state::OuterRouteState,
    f::OuterRouteFeatures;
    tuning::OuterRouteTuning=OuterRouteTuning(),
    machine_class::Symbol=:small,
    threads_available::Bool=true,
    parallel_enabled::Bool=true
)::Symbol
    default_route = default_outer_route(
        f;
        tuning=tuning,
        machine_class=machine_class,
        threads_available=threads_available,
        parallel_enabled=parallel_enabled
    )
    if !tuning.adaptive_enabled
        return default_route
    end

    candidates = outer_route_candidates(
        f;
        tuning=tuning,
        machine_class=machine_class,
        threads_available=threads_available,
        parallel_enabled=parallel_enabled
    )
    isempty(candidates) && return :none
    if !(default_route in candidates)
        default_route = first(candidates)
    end

    signature = outer_route_signature(f)
    snapshot = outer_route_stats_snapshot(state, signature)
    chosen = default_route
    reason = "default"
    if !isempty(snapshot)
        explore = _under_sampled_candidate(candidates, snapshot, default_route, tuning.adaptive_min_samples)
        if !(explore === nothing)
            chosen = explore
            reason = "explore"
        else
            best = _best_candidate(candidates, snapshot)
            if !(best === nothing)
                chosen = best
                reason = chosen == default_route ? "default_best" : "exploit_best"
            end
        end
    end
    if tuning.trace
        println(
            "[outer-route] signature=$(signature) default=$(default_route) chosen=$(chosen) " *
            "reason=$(reason) candidates=$(join(string.(candidates), ','))"
        )
    end
    return chosen
end

function record_outer_route_feedback!(
    state::OuterRouteState,
    f::OuterRouteFeatures;
    route::Symbol,
    successes::Int,
    failures::Int,
    elapsed_success_s::Float64=0.0,
    tuning::OuterRouteTuning=OuterRouteTuning()
)::Nothing
    route in (:none, :threads, :process) || return nothing
    success_count = max(0, successes)
    failure_count = max(0, failures)
    samples = success_count + failure_count
    samples <= 0 && return nothing
    elapsed_sum_s = max(0.0, elapsed_success_s) + failure_count * max(0.0, tuning.failure_penalty_s)
    signature = outer_route_signature(f)
    lock(state.lock) do
        bucket = get!(state.history, signature) do
            Dict{Symbol, OuterRouteStats}()
        end
        stats = get!(bucket, route) do
            OuterRouteStats()
        end
        stats.samples += samples
        stats.successes += success_count
        stats.failures += failure_count
        stats.elapsed_sum_s += elapsed_sum_s
    end
    return nothing
end

end # module ParallelProfiles
