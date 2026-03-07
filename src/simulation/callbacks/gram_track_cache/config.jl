@inline GramTrackCache() = GramTrackCache(
    false,
    0.0,
    0.0,
    1,
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    SVector{3, Float64}[]
)

@inline function _parse_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    return parsed
end

@inline function _parse_float_env_optional(name::String)::Union{Nothing, Float64}
    haskey(ENV, name) || return nothing
    raw = strip(ENV[name])
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    return parsed
end

@inline function _parse_int_env_optional(name::String)::Union{Nothing, Int}
    haskey(ENV, name) || return nothing
    raw = strip(ENV[name])
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer value, got '$raw'"))
    end
    return parsed
end

struct GramTrackCacheConfig
    mode::Symbol
    entry_horizon_s::Float64
    entry_alt_tol_m::Float64
    entry_ang_tol_rad::Float64
    entry_points::Int
    orbit_horizon_s::Float64
    orbit_alt_tol_m::Float64
    orbit_ang_tol_rad::Float64
    orbit_points::Int
    transition_band_m::Float64
end

@inline function _gram_track_cache_mode()::Symbol
    # Benchmark note (2026-02-27):
    # Track-cache refresh can dominate runtime due frequent large GRAM trajectory builds.
    # Example entry case: point-to-point ~0.45 s vs cache+Allen-Eggers ~42.8 s.
    # Keep point-to-point as primary default; enable cache explicitly via env.
    raw = haskey(ENV, "SPACEAGORA_GRAM_TRACK_CACHE") ?
        ENV["SPACEAGORA_GRAM_TRACK_CACHE"] :
        get(ENV, "SPACEAGORA_GRAM_SEGMENT_CACHE", "off")
    mode = lowercase(strip(raw))
    if mode in ("off", "none", "0", "false", "no")
        return :off
    elseif mode in ("on", "1", "true", "yes")
        return :on
    elseif mode == "auto"
        return :auto
    end
    throw(ArgumentError("Unsupported GRAM track-cache mode '$mode'. Use one of: off, auto, on."))
end

function _gram_track_cache_config()::GramTrackCacheConfig
    mode = _gram_track_cache_mode()

    # Backward-compatible global knobs still apply unless regime-specific values are provided.
    compat_horizon = _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_HORIZON_S")
    compat_alt_tol = _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ALT_TOL_M")
    compat_ang_tol_deg = _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ANG_TOL_DEG")
    compat_points = _parse_int_env_optional("SPACEAGORA_GRAM_TRACK_CACHE_NPOS")
    if compat_points === nothing
        compat_points = _parse_int_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_POINTS")
    end

    entry_horizon_s = max(
        1e-3,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_HORIZON_S"),
            something(compat_horizon, 1.0)
        )
    )
    orbit_horizon_s = max(
        1e-3,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_HORIZON_S"),
            something(compat_horizon, 8.0)
        )
    )
    entry_points = max(
        2,
        something(
            _parse_int_env_optional("SPACEAGORA_GRAM_TRACK_CACHE_ENTRY_NPOS"),
            something(
                _parse_int_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_POINTS"),
                something(compat_points, 16)
            )
        )
    )
    orbit_points = max(
        2,
        something(
            _parse_int_env_optional("SPACEAGORA_GRAM_TRACK_CACHE_ORBIT_NPOS"),
            something(
                _parse_int_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_POINTS"),
                something(compat_points, 48)
            )
        )
    )
    entry_alt_tol_m = max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_ALT_TOL_M"),
            something(compat_alt_tol, 500.0)
        )
    )
    orbit_alt_tol_m = max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_ALT_TOL_M"),
            something(compat_alt_tol, 3000.0)
        )
    )
    entry_ang_tol_rad = deg2rad(max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_ANG_TOL_DEG"),
            something(compat_ang_tol_deg, 0.6)
        )
    ))
    orbit_ang_tol_rad = deg2rad(max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_ANG_TOL_DEG"),
            something(compat_ang_tol_deg, 4.0)
        )
    ))
    transition_band_m = max(0.0, _parse_float_env("SPACEAGORA_GRAM_SEGMENT_CACHE_TRANSITION_BAND_M", 20e3))

    return GramTrackCacheConfig(
        mode,
        entry_horizon_s,
        entry_alt_tol_m,
        entry_ang_tol_rad,
        entry_points,
        orbit_horizon_s,
        orbit_alt_tol_m,
        orbit_ang_tol_rad,
        orbit_points,
        transition_band_m
    )
end

# Backward-compatibility shims for older cached method bodies/tests that still
# reference the previous "segment cache" helper names.
@inline _gram_segment_cache_mode() = _gram_track_cache_mode()
@inline _gram_segment_cache_config() = _gram_track_cache_config()
@inline _gram_segment_cache_enabled(cfg::GramTrackCacheConfig, density_model) = _gram_track_cache_enabled(cfg, density_model)
@inline _gram_segment_cache_profile(cfg::GramTrackCacheConfig, p, alt::Float64) = _gram_track_cache_profile(cfg, p, alt)
@inline _gram_segment_cache_ready(cache::GramTrackCache, t::Float64, alt::Float64, lat::Float64, lon::Float64, alt_tol_m::Float64, ang_tol_rad::Float64) =
    _gram_track_cache_ready(cache, t, alt, lat, lon, alt_tol_m, ang_tol_rad)
@inline _gram_segment_cache_eval(cache::GramTrackCache, idx::Int, x::Float64) = _gram_track_cache_eval(cache, idx, x)
@inline _gram_segment_cache_refresh!(cache::GramTrackCache, density_model, p, pos::SVector{3, Float64}, vel::SVector{3, Float64}, alt::Float64, lat::Float64, lon::Float64, t::Float64, horizon_s::Float64, n_points::Int) =
    _gram_track_cache_refresh!(cache, density_model, p, pos, vel, alt, lat, lon, t, horizon_s, n_points)
