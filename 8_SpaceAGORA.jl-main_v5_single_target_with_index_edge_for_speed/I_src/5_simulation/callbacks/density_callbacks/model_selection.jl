@inline function _density_model_for_sat(
    density_models::AbstractVector{<:AbstractDensityModel},
    fallback_model::AbstractDensityModel,
    sat_idx::Int
)
    if sat_idx <= length(density_models)
        return density_models[sat_idx]
    end
    return fallback_model
end

@inline function _density_model_for_sat(p, sat_idx::Int)
    return _density_model_for_sat(
        p.shared_buffers.density_models,
        p.args.environment_model.density_model,
        sat_idx,
    )
end

@inline function _density_batch_model_for_callback(
    density_models::AbstractVector{<:AbstractDensityModel},
    fallback_model::AbstractDensityModel,
    num_sats::Int
)
    if isempty(density_models)
        return fallback_model
    end
    if length(density_models) < num_sats
        return nothing
    end
    first_model = density_models[1]
    @inbounds for i in 2:num_sats
        density_models[i] === first_model || return nothing
    end
    return first_model
end

@inline function _density_batch_model_for_callback(p, num_sats::Int)
    return _density_batch_model_for_callback(
        p.shared_buffers.density_models,
        p.args.environment_model.density_model,
        num_sats,
    )
end

@inline function _gram_isolated_pool_batch_model_for_callback(
    density_models::AbstractVector{<:AbstractDensityModel},
    fallback_model::AbstractDensityModel,
    num_sats::Int
)
    _gram_isolated_pool_enabled(num_sats) || return nothing
    isempty(density_models) || return nothing
    model = fallback_model
    return model isa EnvironmentModels.GRAMAtmosphereModel ? model : nothing
end

@inline function _gram_isolated_pool_batch_model_for_callback(
    env::CallbackEnvConfig,
    density_models::AbstractVector{<:AbstractDensityModel},
    fallback_model::AbstractDensityModel,
    num_sats::Int
)
    _gram_isolated_pool_enabled(env, num_sats) || return nothing
    isempty(density_models) || return nothing
    model = fallback_model
    return model isa EnvironmentModels.GRAMAtmosphereModel ? model : nothing
end

@inline function _gram_isolated_pool_batch_model_for_callback(p, num_sats::Int)
    return _gram_isolated_pool_batch_model_for_callback(
        _callback_env_config(p),
        p.shared_buffers.density_models,
        p.args.environment_model.density_model,
        num_sats,
    )
end

@inline function _gram_batch_elapsed_time(el_time::Float64, idx::Int)::Float64
    return el_time
end

@inline function _gram_batch_elapsed_time(el_time::AbstractVector{<:Real}, idx::Int)::Float64
    return Float64(el_time[idx])
end

@inline function _gram_isolated_pool_density_state(
    model::EnvironmentModels.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p,
    model_lock::ReentrantLock
)::Tuple{Float64, Float64, SVector{3, Float64}}
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0
    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return EnvironmentModels.density_polyfit(h, p)
    end
    return EnvironmentModels._gram_core_density_state(
        model.core,
        h,
        lat,
        lon,
        el_time,
        wind,
        model_lock,
        p.args.environment_model.planet.T_ref
    )
end

@inline function _gram_density_cache_for_sat!(
    caches::Vector{Union{Nothing, GramTrackCache}},
    sat_idx::Int
)::GramTrackCache
    if sat_idx <= length(caches)
        cache = @inbounds caches[sat_idx]
        if cache === nothing
            cache = GramTrackCache()
            @inbounds caches[sat_idx] = cache
        end
        return cache
    end
    return GramTrackCache()
end

function _ensure_gram_isolated_pool!(
    p,
    template_model::EnvironmentModels.GRAMAtmosphereModel,
    workers::Int
)::Tuple{Vector{EnvironmentModels.GRAMAtmosphereModel}, Vector{ReentrantLock}}
    shared = p.shared_buffers
    models = shared.gram_isolated_pool_models
    locks = shared.gram_isolated_pool_locks
    if workers <= 0
        return models, locks
    end
    rebuild = length(models) != workers || length(locks) != workers
    if rebuild
        empty!(models)
        empty!(locks)
        sizehint!(models, workers)
        sizehint!(locks, workers)
        @inbounds for _ in 1:workers
            push!(models, deepcopy(template_model))
            push!(locks, ReentrantLock())
        end
    end
    return models, locks
end

@inline function _gram_isolated_pool_batch_eval!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    density_model,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p;
    allotment_hint::Int=max(1, Threads.nthreads())
)::Bool
    return false
end

function _gram_isolated_pool_batch_eval!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    density_model::EnvironmentModels.GRAMAtmosphereModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p;
    allotment_hint::Int=max(1, Threads.nthreads())
)::Bool
    n = length(hs)
    env = _callback_env_config(p)
    _gram_isolated_pool_enabled(env, n) || return false
    length(rhos) == n || return false
    length(Ts) == n || return false
    length(winds) == n || return false
    length(lats) == n || return false
    length(lons) == n || return false
    if el_time isa AbstractVector{<:Real}
        length(el_time) == n || return false
    end

    max_allotment = min(max(1, allotment_hint), env.gram_isolated_pool_max_workers)
    workers = ParallelPolicy.thread_worker_count(n, max_allotment)
    workers > 1 || return false

    models, locks = _ensure_gram_isolated_pool!(p, density_model, workers)
    ParallelPolicy.threaded_foreach_worker_persistent(:rhs_gram_batch, n, workers) do worker_id, idx
        model_i = models[worker_id]::EnvironmentModels.GRAMAtmosphereModel
        lock_i = locks[worker_id]
        h = Float64(hs[idx])
        lat = Float64(lats[idx])
        lon = Float64(lons[idx])
        t_i = _gram_batch_elapsed_time(el_time, idx)
        rho_i, T_i, wind_i = _gram_isolated_pool_density_state(model_i, h, lat, lon, t_i, wind, p, lock_i)
        @inbounds begin
            rhos[idx] = rho_i
            Ts[idx] = T_i
            winds[idx] = wind_i
        end
    end
    return true
end
