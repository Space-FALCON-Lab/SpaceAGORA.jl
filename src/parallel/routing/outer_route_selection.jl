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

@inline function _route_max_link_bucket(n_links::Int)::String
    if n_links <= 1
        return "1"
    elseif n_links <= 2
        return "2"
    elseif n_links <= 4
        return "3_4"
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

@inline function _route_density_bucket(family::AbstractString)::String
    token = lowercase(strip(String(family)))
    if token in ("none", "vacuum", "noatmosphere")
        return "none"
    elseif token in ("gram_point", "gram")
        return "gram_pt"
    elseif token in ("gram_surrogate", "gram_offline_surrogate")
        return "gram_srg"
    elseif token in ("exponential", "exp")
        return "exp"
    elseif token in ("polynomialfit", "polyfit", "poly")
        return "poly"
    elseif token in ("nrlmsise00", "nrl")
        return "nrl"
    elseif isempty(token)
        return "unknown"
    end
    return replace(token, "|" => "_")
end

@inline function _route_solver_bucket(mode::AbstractString)::String
    token = lowercase(strip(String(mode)))
    if token in ("auto", "auto_stiff", "autostiff")
        return "auto"
    elseif token in ("rodas5p", "rodas")
        return "rodas"
    elseif token in ("tsit5",)
        return "tsit5"
    elseif startswith(token, "split_imex")
        return "split"
    elseif token in ("multirate",)
        return "mrate"
    elseif isempty(token)
        return "auto"
    end
    return replace(token, "|" => "_")
end

@inline function _route_interval_bucket(interval_s::Float64)::String
    if !isfinite(interval_s) || interval_s <= 0.0
        return "na"
    elseif interval_s <= 0.2
        return "ultra"
    elseif interval_s <= 1.0
        return "fast"
    elseif interval_s <= 10.0
        return "med"
    end
    return "slow"
end

@inline function _route_count_bucket(v::Int)::String
    if v <= 0
        return "0"
    elseif v == 1
        return "1"
    elseif v <= 3
        return "2_3"
    elseif v <= 6
        return "4_6"
    end
    return "7p"
end

@inline function _route_effector_cost_bucket(raw::AbstractString)::String
    token = lowercase(strip(String(raw)))
    if token in ("light", "medium", "heavy")
        return token
    elseif isempty(token)
        return "unknown"
    end
    return replace(token, "|" => "_")
end

@inline function _compat_outer_route_signature(f::OuterRouteFeatures)::String
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

@inline function outer_route_signature(f::OuterRouteFeatures)::String
    return join((
        "cat=$(f.category)",
        "sat=$(_route_sat_bucket(f.n_sats))",
        "links=$(_route_link_bucket(f.n_links))",
        "maxlinks=$(_route_max_link_bucket(f.max_links_per_sat))",
        "mission=$(_route_mission_bucket(f.mission_time_s))",
        "nbody=$(f.has_nbody ? "1" : "0")",
        "srp=$(f.has_srp ? "1" : "0")",
        "harm=$(_route_harmonics_bucket(f.harmonics_degree))",
        "ctrl=$(f.has_control ? "1" : "0")",
        "orient=$(f.orientation_on ? "1" : "0")",
        "dens=$(_route_density_bucket(f.density_family))",
        "solver=$(_route_solver_bucket(f.solver_mode))",
        "dt=$(_route_interval_bucket(f.dt_max_orbit_s))",
        "ctrl_rate=$(_route_interval_bucket(f.control_rate_s))",
        "guid_rate=$(_route_interval_bucket(f.guidance_rate_s))",
        "nav_rate=$(_route_interval_bucket(f.navigation_rate_s))",
        "gram_srg=$(f.gram_surrogate_enabled ? "1" : "0")",
        "gram_grid=$(f.gram_static_grid_enabled ? "1" : "0")",
        "ctrl_eff=$(_route_count_bucket(f.control_effector_count))",
        "thermal=$(f.thermal_enabled ? "1" : "0")",
        "eff_cnt=$(_route_count_bucket(f.dynamic_effector_count))",
        "eff_cost=$(_route_effector_cost_bucket(f.effector_cost_class))"
    ), "|")
end

@inline function _outer_route_signature_hierarchy(f::OuterRouteFeatures)::Vector{String}
    full = outer_route_signature(f)
    mid = join((
        "cat=$(f.category)",
        "sat=$(_route_sat_bucket(f.n_sats))",
        "links=$(_route_link_bucket(f.n_links))",
        "maxlinks=$(_route_max_link_bucket(f.max_links_per_sat))",
        "mission=$(_route_mission_bucket(f.mission_time_s))",
        "nbody=$(f.has_nbody ? "1" : "0")",
        "srp=$(f.has_srp ? "1" : "0")",
        "harm=$(_route_harmonics_bucket(f.harmonics_degree))",
        "ctrl=$(f.has_control ? "1" : "0")",
        "orient=$(f.orientation_on ? "1" : "0")",
        "dens=$(_route_density_bucket(f.density_family))",
        "solver=$(_route_solver_bucket(f.solver_mode))",
        "thermal=$(f.thermal_enabled ? "1" : "0")",
        "eff_cost=$(_route_effector_cost_bucket(f.effector_cost_class))"
    ), "|")
    legacy = _compat_outer_route_signature(f)
    return unique(String[full, mid, legacy])
end

@inline function _route_elapsed_stats(
    stats::OuterRouteStats
)::NamedTuple{(:mean_s, :std_s), Tuple{Float64, Float64}}
    samples = max(1, stats.samples)
    mean_s = stats.elapsed_sum_s / samples
    mean_sq_s = stats.elapsed_sq_sum_s / samples
    variance_s = max(0.0, mean_sq_s - mean_s^2)
    return (mean_s=mean_s, std_s=sqrt(variance_s))
end

function _outer_route_stats_snapshot_internal(
    state::OuterRouteState,
    signature::String
)::Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate, :std_s), Tuple{Int, Float64, Float64, Float64}}}
    lock(state.lock) do
        entry = get(state.history, signature, nothing)
        snap = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate, :std_s), Tuple{Int, Float64, Float64, Float64}}}()
        if entry === nothing
            return snap
        end
        for (route, stats) in entry
            if stats.samples <= 0
                continue
            end
            success_rate = stats.successes / max(1, stats.samples)
            elapsed = _route_elapsed_stats(stats)
            snap[route] = (
                samples=stats.samples,
                mean_s=elapsed.mean_s,
                success_rate=success_rate,
                std_s=elapsed.std_s
            )
        end
        return snap
    end
end

function outer_route_stats_snapshot(
    state::OuterRouteState,
    signature::String
)::Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}}
    base = _outer_route_stats_snapshot_internal(state, signature)
    snap = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}}()
    for (route, info) in base
        snap[route] = (samples=info.samples, mean_s=info.mean_s, success_rate=info.success_rate)
    end
    return snap
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

@inline function _is_native_gram_point_density(f::OuterRouteFeatures)::Bool
    return _route_density_bucket(f.density_family) == "gram_pt" &&
        !f.gram_surrogate_enabled &&
        !f.gram_static_grid_enabled
end

@inline function _feature_heavy_for_process(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool
    if _is_native_gram_point_density(f)
        # Native GRAM point calls are lock-limited; prefer outer process isolation.
        return true
    end
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

    if _is_native_gram_point_density(f)
        return :process
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
    if !parallel_enabled
        return Symbol[:none]
    end
    if _is_native_gram_point_density(f)
        return Symbol[:none, :process]
    end

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
    snapshot,
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
    snapshot
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

@inline function _candidate_confidence_width(
    std_s::Float64,
    samples::Int,
    total_samples::Int,
    exploration_c::Float64
)::Float64
    samples <= 0 && return 0.0
    exploration = max(0.0, exploration_c)
    scaled = std_s * sqrt(log(max(2.0, Float64(total_samples) + 1.0)) / max(1.0, Float64(samples)))
    width = exploration * scaled
    return isfinite(width) ? max(0.0, width) : 0.0
end

@inline function _best_candidate_confidence(
    candidates::Vector{Symbol},
    snapshot,
    default_route::Symbol,
    exploration_c::Float64
)::NamedTuple{(:route, :confidence_s, :regret_s), Tuple{Union{Nothing, Symbol}, Float64, Float64}}
    total_samples = 0
    for route in candidates
        info = get(snapshot, route, nothing)
        info === nothing && continue
        total_samples += max(0, info.samples)
    end
    total_samples = max(1, total_samples)

    best_route = nothing
    best_score = Inf
    best_mean = Inf
    best_width = 0.0
    best_success_rate = -Inf
    best_observed_mean = Inf

    ranked = _route_ranked_candidates(candidates, default_route)
    for route in ranked
        info = get(snapshot, route, nothing)
        info === nothing && continue
        if info.samples <= 0 || !isfinite(info.mean_s)
            continue
        end
        best_observed_mean = min(best_observed_mean, info.mean_s)
        width = _candidate_confidence_width(info.std_s, info.samples, total_samples, exploration_c)
        score = info.mean_s - width
        if score < best_score - 1e-12
            best_route = route
            best_score = score
            best_mean = info.mean_s
            best_width = width
            best_success_rate = info.success_rate
        elseif isapprox(score, best_score; atol=1e-12, rtol=0.0) && info.success_rate > best_success_rate + 1e-12
            best_route = route
            best_mean = info.mean_s
            best_width = width
            best_success_rate = info.success_rate
        end
    end

    if best_route === nothing
        return (route=nothing, confidence_s=0.0, regret_s=0.0)
    end
    regret_s = isfinite(best_observed_mean) ? max(0.0, best_mean - best_observed_mean) : 0.0
    return (route=best_route, confidence_s=best_width, regret_s=regret_s)
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

    signature_chain = _outer_route_signature_hierarchy(f)
    signature = first(signature_chain)
    snapshot = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate, :std_s), Tuple{Int, Float64, Float64, Float64}}}()
    signature_used = signature
    for candidate_sig in signature_chain
        snap = _outer_route_stats_snapshot_internal(state, candidate_sig)
        if !isempty(snap)
            snapshot = snap
            signature_used = candidate_sig
            break
        end
    end
    chosen = default_route
    reason = "default"
    confidence_s = 0.0
    regret_s = 0.0
    if !isempty(snapshot)
        explore = _under_sampled_candidate(candidates, snapshot, default_route, tuning.adaptive_min_samples)
        if !(explore === nothing)
            chosen = explore
            reason = "explore_hier"
        else
            best = _best_candidate_confidence(
                candidates,
                snapshot,
                default_route,
                tuning.adaptive_exploration_c
            )
            if !(best.route === nothing)
                chosen = best.route
                confidence_s = best.confidence_s
                regret_s = best.regret_s
                reason = chosen == default_route ? "default_ucb_hier" : "exploit_ucb_hier"
            end
        end
    end
    if tuning.trace
        println(
            "[outer-route] signature=$(signature) default=$(default_route) chosen=$(chosen) " *
            "reason=$(reason) signature_used=$(signature_used) candidates=$(join(string.(candidates), ',')) " *
            "confidence_s=$(round(confidence_s; digits=6)) regret_s=$(round(regret_s; digits=6))"
        )
    end
    return chosen
end

