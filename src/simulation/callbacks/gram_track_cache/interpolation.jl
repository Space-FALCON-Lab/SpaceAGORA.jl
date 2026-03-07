@inline function _gram_track_cache_enabled(cfg::GramTrackCacheConfig, density_model)::Bool
    cfg.mode == :off && return false
    return _is_gram_density_model(density_model)
end

@inline function _gram_track_cache_profile(cfg::GramTrackCacheConfig, p, alt::Float64)
    EI_m = p.args.environment_model.EI * 1e3
    if alt <= EI_m + cfg.transition_band_m
        return cfg.entry_horizon_s, cfg.entry_alt_tol_m, cfg.entry_ang_tol_rad, cfg.entry_points
    end
    return cfg.orbit_horizon_s, cfg.orbit_alt_tol_m, cfg.orbit_ang_tol_rad, cfg.orbit_points
end

@inline _lerp(a::Float64, b::Float64, x::Float64) = a + x * (b - a)
@inline _angdiff_rad(a::Float64, b::Float64) = abs(atan(sin(a - b), cos(a - b)))
@inline function _lerp_angle_rad(a::Float64, b::Float64, x::Float64)
    δ = atan(sin(b - a), cos(b - a))
    return atan(sin(a + x * δ), cos(a + x * δ))
end

@inline function _gram_track_cache_segment(
    cache::GramTrackCache,
    t::Float64
)::Union{Nothing, Tuple{Int, Float64}}
    n = length(cache.times)
    if !cache.valid || n < 2
        return nothing
    end
    ignore_time_window = _gram_track_cache_ignore_time_window()
    if !ignore_time_window && (t < cache.t0 || t > cache.t1)
        return nothing
    end
    tq = ignore_time_window ? clamp(t, cache.t0, cache.t1) : t
    idx = clamp(cache.index_hint, 1, n - 1)
    if tq < cache.times[idx] || tq > cache.times[idx + 1]
        idx = clamp(searchsortedlast(cache.times, tq), 1, n - 1)
    end
    t_lo = cache.times[idx]
    t_hi = cache.times[idx + 1]
    if tq < t_lo || tq > t_hi
        return nothing
    end
    x = t_hi == t_lo ? 0.0 : (tq - t_lo) / (t_hi - t_lo)
    cache.index_hint = idx
    return idx, x
end

@inline function _gram_track_cache_ready(
    cache::GramTrackCache,
    t::Float64,
    alt::Float64,
    lat::Float64,
    lon::Float64,
    alt_tol_m::Float64,
    ang_tol_rad::Float64
)::Union{Nothing, Tuple{Int, Float64}}
    seg = _gram_track_cache_segment(cache, t)
    seg === nothing && return nothing
    idx, x = seg
    alt_interp = _lerp(cache.alts[idx], cache.alts[idx + 1], x)
    lat_interp = _lerp_angle_rad(cache.lats[idx], cache.lats[idx + 1], x)
    lon_interp = _lerp_angle_rad(cache.lons[idx], cache.lons[idx + 1], x)
    if abs(alt - alt_interp) <= alt_tol_m &&
       _angdiff_rad(lat, lat_interp) <= ang_tol_rad &&
       _angdiff_rad(lon, lon_interp) <= ang_tol_rad
        return idx, x
    end
    return nothing
end

@inline function _gram_track_cache_eval(cache::GramTrackCache, idx::Int, x::Float64)
    ρ = _lerp(cache.rhos[idx], cache.rhos[idx + 1], x)
    T = _lerp(cache.Ts[idx], cache.Ts[idx + 1], x)
    wind = cache.winds[idx] + x * (cache.winds[idx + 1] - cache.winds[idx])
    return ρ, T, wind
end
