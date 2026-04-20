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

# Fast signed angular difference in (-π, π].  For most consecutive orbital
# cache points |b-a| ≪ π, so we first try the cheap branch before falling
# back to the full atan2 path.
@inline function _angdiff_rad(a::Float64, b::Float64)::Float64
    d = b - a
    # Cheap path: |d| already in (-π, π]
    if -π < d ≤ π
        return abs(d)
    end
    return abs(atan(sin(d), cos(d)))
end

# Linear angle interpolation.  Uses the same short-arc wrapping but avoids
# a second atan2 call by computing the result as a + x*δ and re-normalising
# only when the result drifts outside (-π, π].
@inline function _lerp_angle_rad(a::Float64, b::Float64, x::Float64)::Float64
    d = b - a
    # Wrap the difference to the short arc without trig when possible.
    if d > π
        d -= 2π
    elseif d < -π
        d += 2π
    end
    r = a + x * d
    # Re-normalise cheaply (no trig needed — just a conditional subtract/add).
    if r > π
        r -= 2π
    elseif r ≤ -π
        r += 2π
    end
    return r
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

    # Walk the hint forward before falling back to binary search.
    # Integrators advance time monotonically; the common case is that tq sits
    # in the current hint segment or the one immediately after it.
    idx = clamp(cache.index_hint, 1, n - 1)
    times = cache.times
    @inbounds if tq < times[idx]
        # Time went backward (e.g. rejected step); fall through to binary search.
        idx = clamp(searchsortedlast(times, tq), 1, n - 1)
    elseif tq > times[idx + 1]
        # Advance hint one step first (O(1) common case).
        next = idx + 1
        if next < n && tq <= times[next + 1]
            idx = next
        else
            idx = clamp(searchsortedlast(times, tq), 1, n - 1)
        end
    end
    # No else needed: tq ∈ [times[idx], times[idx+1]] — already in segment.

    @inbounds t_lo = times[idx]
    @inbounds t_hi = times[idx + 1]
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
    # Avoid a subtraction+add of full SVectors; use muladd element-wise via
    # the two-argument form: w = w0 + x*(w1 - w0)  →  w = (1-x)*w0 + x*w1
    @inbounds w0 = cache.winds[idx]
    @inbounds w1 = cache.winds[idx + 1]
    mx = 1.0 - x
    wind = SVector{3, Float64}(
        muladd(x, w1[1], mx * w0[1]),
        muladd(x, w1[2], mx * w0[2]),
        muladd(x, w1[3], mx * w0[3])
    )
    return ρ, T, wind
end
