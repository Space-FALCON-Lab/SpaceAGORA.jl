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

@inline _gram_track_cache_periapsis_split_enabled() = _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT", true)

@inline function _gram_track_cache_max_npos()::Int
    raw = strip(get(ENV, "SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS", "512"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS must be an integer value, got '$raw'"))
    end
    return max(2, parsed)
end

@inline _angle_delta_rad(a::Float64, b::Float64) = atan(sin(b - a), cos(b - a))

@inline function _gram_expected_track_length_m(
    planet,
    alt0::Float64,
    lat0::Float64,
    lon0::Float64,
    alt1::Float64,
    lat1::Float64,
    lon1::Float64
)::Tuple{Float64, Float64}
    radius_m = max(1.0, planet.Rp_m + 0.5 * (alt0 + alt1))
    cos_angle = clamp(
        sin(lat0) * sin(lat1) + cos(lat0) * cos(lat1) * cos(_angle_delta_rad(lon0, lon1)),
        -1.0,
        1.0
    )
    central_angle = acos(cos_angle)
    horizontal_m = radius_m * central_angle
    vertical_m = alt1 - alt0
    length_m = hypot(horizontal_m, vertical_m)
    return max(length_m, 1.0), radius_m
end

@inline function _gram_track_cache_target_spacing_m(
    alt_tol_m::Float64,
    ang_tol_rad::Float64,
    radius_m::Float64
)::Float64
    # Keep interpolation spacing tighter than cache acceptance tolerances.
    ang_tol_m = max(1.0, radius_m * max(ang_tol_rad, 1e-9))
    tol_scale_m = min(max(1.0, alt_tol_m), ang_tol_m)
    return max(1.0, 0.5 * tol_scale_m)
end

@inline function _solve_kepler_elliptic(M::Float64, e::Float64)::Float64
    M2π = mod(M, 2pi)
    E = e < 0.8 ? M2π : pi
    @inbounds for _ in 1:20
        f = E - e * sin(E) - M2π
        fp = 1.0 - e * cos(E)
        abs(fp) < 1e-14 && break
        dE = f / fp
        E -= dE
        abs(dE) <= 1e-12 && break
    end
    return E
end

@inline function _true_to_eccentric_anomaly(ν::Float64, e::Float64)::Float64
    E = atan(sqrt(max(0.0, 1.0 - e^2)) * sin(ν), e + cos(ν))
    return mod(E, 2pi)
end

@inline function _eccentric_to_true_anomaly(E::Float64, e::Float64)::Float64
    ν = atan(sqrt(max(0.0, 1.0 - e^2)) * sin(E), cos(E) - e)
    return mod(ν, 2pi)
end

@inline function _gram_j2_rates(a::Float64, e::Float64, i::Float64, n::Float64, planet)::Tuple{Float64, Float64}
    if !isfinite(planet.J2) || planet.J2 == 0.0
        return 0.0, 0.0
    end
    p = a * (1.0 - e^2)
    if !isfinite(p) || p <= 0.0
        return 0.0, 0.0
    end
    scale = planet.J2 * (planet.Rp_e / p)^2
    Ωdot = -1.5 * n * scale * cos(i)
    ωdot = 0.75 * n * scale * (5.0 * cos(i)^2 - 1.0)
    return Ωdot, ωdot
end

function _gram_kepler_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64;
    include_j2::Bool = true
)::Union{Nothing, Tuple{Float64, Float64, Float64}}
    try
        oe = rvtoorbitalelement(pos, vel, planet)
        a, e, i, Ω, ω, ν = Float64(oe[1]), Float64(oe[2]), Float64(oe[3]), Float64(oe[4]), Float64(oe[5]), Float64(oe[6])
        if !isfinite(dt) || dt <= 1e-6 || !isfinite(a) || !isfinite(e) || !isfinite(ν) || a <= 0.0 || e < 0.0 || e >= 1.0
            return nothing
        end
        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            return nothing
        end
        E0 = _true_to_eccentric_anomaly(ν, e)
        M0 = E0 - e * sin(E0)
        E1 = _solve_kepler_elliptic(M0 + n * dt, e)
        ν1 = _eccentric_to_true_anomaly(E1, e)

        Ω1 = Ω
        ω1 = ω
        if include_j2
            Ωdot, ωdot = _gram_j2_rates(a, e, i, n, planet)
            Ω1 += Ωdot * dt
            ω1 += ωdot * dt
        end

        oe_target = SVector{7, Float64}(a, e, i, Ω1, ω1, ν1, 0.0)
        r_target_vec, v_target_vec = orbitalelemtorv(oe_target, planet)
        r_target_i = SVector{3, Float64}(Float64.(r_target_vec))
        v_target_i = SVector{3, Float64}(Float64.(v_target_vec))
        r_target_p, _ = r_intor_p!(r_target_i, v_target_i, planet)
        alt_target, lat_target, lon_target = rtolatlong(r_target_p, planet)
        return alt_target, lat_target, lon_target
    catch
        return nothing
    end
end

function _gram_periapsis_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet;
    include_j2::Bool = true
)::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}}
    try
        oe = rvtoorbitalelement(pos, vel, planet)
        a, e, ν = Float64(oe[1]), Float64(oe[2]), Float64(oe[6])
        if !isfinite(a) || !isfinite(e) || !isfinite(ν) || a <= 0.0 || e < 0.0 || e >= 1.0
            return nothing
        end

        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            return nothing
        end
        E = _true_to_eccentric_anomaly(ν, e)
        M = mod(E - e * sin(E), 2pi)
        dt_peri = (2pi - M) / n
        if !isfinite(dt_peri) || dt_peri <= 1e-6
            return nothing
        end
        target = _gram_kepler_target(pos, vel, planet, dt_peri; include_j2=include_j2)
        target === nothing && return nothing
        alt_peri, lat_peri, lon_peri = target
        return dt_peri, alt_peri, lat_peri, lon_peri
    catch
        return nothing
    end
end

@inline _gram_descending_periapsis_target(pos::SVector{3, Float64}, vel::SVector{3, Float64}, planet; include_j2::Bool = true) =
    _gram_periapsis_target(pos, vel, planet; include_j2=include_j2)

function _gram_orbit_period_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet;
    include_j2::Bool = true
)::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}}
    try
        oe = rvtoorbitalelement(pos, vel, planet)
        a, e = Float64(oe[1]), Float64(oe[2])
        if !isfinite(a) || !isfinite(e) || a <= 0.0 || e < 0.0 || e >= 1.0
            return nothing
        end
        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            return nothing
        end
        dt_orbit = 2pi / n
        if !isfinite(dt_orbit) || dt_orbit <= 1e-6
            return nothing
        end
        target = _gram_kepler_target(pos, vel, planet, dt_orbit; include_j2=include_j2)
        target === nothing && return nothing
        alt_end, lat_end, lon_end = target
        return dt_orbit, alt_end, lat_end, lon_end
    catch
        return nothing
    end
end

@inline function _gram_linear_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64
)::Tuple{Float64, Float64, Float64}
    pos_target = pos + vel * dt
    rp_target, _ = r_intor_p!(pos_target, vel, planet)
    target_alt, target_lat, target_lon = rtolatlong(rp_target, planet)
    return target_alt, target_lat, target_lon
end

@inline function _gram_entry_reference_area_m2(p, sat_index::Int)::Float64
    try
        spacecraft = p.args.dynamics_model.spacecraft[sat_index]
        area = 0.0
        @inbounds for link in spacecraft.links
            a = Float64(link.ref_area)
            if isfinite(a) && a > 0.0
                area += a
            end
        end
        if area > 0.0
            return area
        end
        a0 = Float64(spacecraft.root.ref_area)
        return (isfinite(a0) && a0 > 0.0) ? a0 : 1.0
    catch
        return 1.0
    end
end

@inline function _gram_entry_mass_kg(p, sat_index::Int, current_mass_kg::Float64)::Float64
    if isfinite(current_mass_kg) && current_mass_kg > 0.0
        return current_mass_kg
    end
    try
        spacecraft = p.args.dynamics_model.spacecraft[sat_index]
        m = Float64(spacecraft.dry_mass + spacecraft.prop_mass)
        return (isfinite(m) && m > 0.0) ? m : 100.0
    catch
        return 100.0
    end
end

@inline function _gram_entry_reference_density(planet, h::Float64)::Float64
    H = max(1.0, planet.H)
    ρ = planet.ρ_ref * exp((planet.h_ref - h) / H)
    return (isfinite(ρ) && ρ > 0.0) ? ρ : 0.0
end

function _gram_entry_target_allen_eggers(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64,
    mass_kg::Float64,
    ref_area_m2::Float64;
    cd::Float64 = _gram_entry_target_cd()
)::Union{Nothing, Tuple{Float64, Float64, Float64}}
    if !isfinite(dt) || dt <= 1e-6 || !isfinite(mass_kg) || mass_kg <= 0.0 || !isfinite(ref_area_m2) || ref_area_m2 <= 0.0
        return nothing
    end
    if !isfinite(planet.ρ_ref) || planet.ρ_ref <= 0.0
        return nothing
    end

    rp, vp = r_intor_p!(pos, vel, planet)
    alt0, lat0, lon0 = rtolatlong(rp, planet)
    uD, uN, uE = latlongtoNED(SVector{3, Float64}(alt0, lat0, lon0))

    vN0 = dot(vp, uN)
    vE0 = dot(vp, uE)
    vU0 = -dot(vp, uD)
    vh0 = hypot(vN0, vE0)
    v = norm(vp)
    if !isfinite(v) || v <= 1e-3
        return nothing
    end

    γ = atan(vU0, max(vh0, 1e-9))
    χ = atan(vE0, vN0)
    h = alt0
    lat = lat0
    lon = lon0
    β = mass_kg / max(1e-6, cd * ref_area_m2)

    n_steps = Int(ceil(dt / _gram_entry_target_dt()))
    n_steps = clamp(n_steps, 2, _gram_entry_target_max_steps())
    Δτ = dt / n_steps

    @inbounds for _ in 1:n_steps
        r = max(1.0, planet.Rp_e + h)
        cosγ = cos(γ)
        sinγ = sin(γ)
        cosχ = cos(χ)
        sinχ = sin(χ)
        coslat = max(1e-6, abs(cos(lat)))
        tanlat = clamp(tan(lat), -100.0, 100.0)

        ρ = _gram_entry_reference_density(planet, h)
        g = planet.μ / (r * r)
        drag_acc = 0.5 * ρ * v * v / β

        v_dot = -drag_acc - g * sinγ
        γ_dot = (v / r - g / max(v, 1e-3)) * cosγ
        χ_dot = v * cosγ * sinχ * tanlat / r
        h_dot = v * sinγ
        lat_dot = v * cosγ * cosχ / r
        lon_dot = v * cosγ * sinχ / (r * coslat) - planet.ω[3]

        v = max(1e-3, v + v_dot * Δτ)
        γ = clamp(γ + γ_dot * Δτ, -deg2rad(89.0), deg2rad(89.0))
        χ = atan(sin(χ + χ_dot * Δτ), cos(χ + χ_dot * Δτ))
        h += h_dot * Δτ
        lat = clamp(lat + lat_dot * Δτ, -deg2rad(89.9), deg2rad(89.9))
        lon = atan(sin(lon + lon_dot * Δτ), cos(lon + lon_dot * Δτ))
    end

    return h, lat, lon
end

function _gram_track_cache_fill_from_trajectory!(
    cache::GramTrackCache,
    trajectory
)
    n = length(trajectory)
    n >= 2 || throw(ArgumentError("GRAM trajectory must contain at least 2 points, got $n."))
    if length(cache.times) != n
        resize!(cache.times, n)
        resize!(cache.alts, n)
        resize!(cache.lats, n)
        resize!(cache.lons, n)
        resize!(cache.rhos, n)
        resize!(cache.Ts, n)
        resize!(cache.winds, n)
    end

    @inbounds for k in 1:n
        pt = trajectory[k]
        pos = pt.position
        dyn = pt.dynamics
        wnd = pt.winds

        cache.times[k] = Float64(pos.elapsedTime)
        cache.alts[k] = Float64(pos.height) * 1e3
        cache.lats[k] = deg2rad(Float64(pos.latitude))
        cache.lons[k] = atan(sin(deg2rad(Float64(pos.longitude))), cos(deg2rad(Float64(pos.longitude))))
        cache.rhos[k] = Float64(dyn.density)
        cache.Ts[k] = Float64(dyn.temperature)
        cache.winds[k] = SVector{3, Float64}(
            Float64(wnd.perturbedEWWind),
            Float64(wnd.perturbedNSWind),
            Float64(wnd.perturbedVerticalWind)
        )
    end

    return nothing
end

function _gram_track_cache_refresh!(
    cache::GramTrackCache,
    density_model,
    p,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    alt::Float64,
    lat::Float64,
    lon::Float64,
    t::Float64,
    horizon_s::Float64,
    n_points::Int,
    alt_tol_m::Float64 = 500.0,
    ang_tol_rad::Float64 = deg2rad(1.0),
    transition_band_m::Float64 = 0.0,
    segment_end_t::Float64 = NaN,
    target_include_j2::Bool = true,
    sat_index::Int = 1,
    current_mass_kg::Float64 = NaN
)
    stats_enabled = _gram_runtime_stats_enabled()
    refresh_t0_ns = stats_enabled ? time_ns() : 0
    try
        planet = p.args.environment_model.planet
        base_horizon_s = max(1e-3, horizon_s)
        include_j2 = target_include_j2
        EI_m = p.args.environment_model.EI * 1e3
        in_atm_band = alt <= EI_m + max(0.0, transition_band_m)
        mission_is_orbit = p.args.mission_configuration.mission_type == MissionOrbits

        target_for_dt = function (dt::Float64)
            target = _gram_kepler_target(pos, vel, planet, dt; include_j2=include_j2)
            return target === nothing ? _gram_linear_target(pos, vel, planet, dt) : target
        end

        dt_segment = base_horizon_s
        target_alt, target_lat, target_lon = target_for_dt(dt_segment)

        if _gram_track_cache_periapsis_split_enabled() && in_atm_band
            # Drag passage: build the cache segment to periapsis when possible.
            peri_target = _gram_periapsis_target(pos, vel, planet; include_j2=include_j2)
            if peri_target !== nothing
                dt_segment, target_alt, target_lat, target_lon = peri_target
            else
                # Entry-like/open trajectories: use Allen-Eggers-type endpoint targeting.
                used_entry_model = false
                if _gram_entry_target_mode() == :allen_eggers
                    dt_entry = if isfinite(segment_end_t)
                        max(1e-3, segment_end_t - t)
                    else
                        base_horizon_s
                    end
                    if isfinite(dt_entry) && dt_entry > 1e-6
                        mass_kg = _gram_entry_mass_kg(p, sat_index, current_mass_kg)
                        area_m2 = _gram_entry_reference_area_m2(p, sat_index)
                        entry_target = _gram_entry_target_allen_eggers(
                            pos,
                            vel,
                            planet,
                            dt_entry,
                            mass_kg,
                            area_m2;
                            cd=_gram_entry_target_cd()
                        )
                        if entry_target !== nothing
                            dt_segment = dt_entry
                            target_alt, target_lat, target_lon = entry_target
                            used_entry_model = true
                        end
                    end
                end

                # If entry model is disabled/unavailable, use solver endpoint propagation.
                if !used_entry_model && isfinite(segment_end_t)
                    dt_to_end = segment_end_t - t
                    if dt_to_end > 1e-6
                        dt_segment = max(1e-3, dt_to_end)
                        target_alt, target_lat, target_lon = target_for_dt(dt_segment)
                    end
                end
            end
        elseif mission_is_orbit
            # Orbit mission: build one segment over a full orbital period.
            orbit_target = _gram_orbit_period_target(pos, vel, planet; include_j2=include_j2)
            if orbit_target !== nothing
                dt_segment, target_alt, target_lat, target_lon = orbit_target
            end
        else
            # Time mission: propagate to the remaining solver time span.
            used_tspan_endpoint = false
            if isfinite(segment_end_t)
                dt_to_end = segment_end_t - t
                if dt_to_end > 1e-6
                    dt_segment = max(1e-3, dt_to_end)
                    target_alt, target_lat, target_lon = target_for_dt(dt_segment)
                    used_tspan_endpoint = true
                end
            end
            if !used_tspan_endpoint
                orbit_target = _gram_orbit_period_target(pos, vel, planet; include_j2=include_j2)
                if orbit_target !== nothing
                    dt_segment, target_alt, target_lat, target_lon = orbit_target
                end
            end
        end

        base_n = max(2, n_points)
        dt_per_sample = max(1e-3, base_horizon_s / max(1, base_n - 1))
        n_from_time = Int(ceil(dt_segment / dt_per_sample)) + 1
        expected_length_m, radius_m = _gram_expected_track_length_m(
            planet,
            alt,
            lat,
            lon,
            target_alt,
            target_lat,
            target_lon
        )
        target_spacing_m = _gram_track_cache_target_spacing_m(alt_tol_m, ang_tol_rad, radius_m)
        n_from_length = Int(ceil(expected_length_m / target_spacing_m)) + 1
        n = max(base_n, n_from_time, n_from_length)
        n = min(n, _gram_track_cache_max_npos())

        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.refresh_calls += 1
                s.refresh_points_total += n
                s.refresh_points_max = max(s.refresh_points_max, n)
            end)
        end

        if length(cache.times) != n
            resize!(cache.times, n)
            resize!(cache.alts, n)
            resize!(cache.lats, n)
            resize!(cache.lons, n)
            resize!(cache.rhos, n)
            resize!(cache.Ts, n)
            resize!(cache.winds, n)
        end

        gram_traj_built = false
        if _gram_track_trajectory_supported(density_model)
            gram_driver = getproperty(density_model, :gram)
            gram_atmosphere = getproperty(density_model, :gram_atmosphere)
            denom = max(1, n - 1)
            Δalt_km = (target_alt - alt) * 1e-3 / denom
            Δlat_deg = rad2deg(_angle_delta_rad(lat, target_lat)) / denom
            Δlon_deg = rad2deg(_angle_delta_rad(lon, target_lon)) / denom
            Δt = dt_segment / denom

            gen_track = () -> Base.invokelatest(
                getproperty(gram_driver, :generate_trajectory),
                gram_atmosphere;
                initial_height=alt * 1e-3,
                initial_latitude=rad2deg(lat),
                initial_longitude=rad2deg(lon),
                initial_elapsed_time=t,
                delta_height=Δalt_km,
                delta_latitude=Δlat_deg,
                delta_longitude=Δlon_deg,
                delta_elapsed_time=Δt,
                n_points=n,
                update_initial_perturbations=true
            )

            trajectory = if EnvironmentModels._gram_use_global_lock()
                lock(GRAM_LOCK) do
                    gen_track()
                end
            else
                gen_track()
            end
            _gram_track_cache_fill_from_trajectory!(cache, trajectory)
            gram_traj_built = true
        end

        if !gram_traj_built
            @inbounds for k in 1:n
                x = (k - 1) / (n - 1)
                t_k = t + x * dt_segment
                pos_k = pos + vel * (t_k - t)
                rp_k, _ = r_intor_p!(pos_k, vel, planet)
                alt_k, lat_k, lon_k = rtolatlong(rp_k, planet)
                cache.times[k] = t_k
                cache.alts[k] = alt_k
                cache.lats[k] = lat_k
                cache.lons[k] = lon_k
            end
            if !_gram_isolated_pool_batch_eval!(
                cache.rhos,
                cache.Ts,
                cache.winds,
                density_model,
                cache.alts,
                cache.lats,
                cache.lons,
                cache.times,
                true,
                p
            )
                getDensityBatch!(
                    cache.rhos,
                    cache.Ts,
                    cache.winds,
                    density_model,
                    cache.alts,
                    cache.lats,
                    cache.lons,
                    cache.times,
                    true,
                    p
                )
            end
        end
        cache.valid = true
        cache.t0 = cache.times[1]
        cache.t1 = cache.times[end]
        cache.index_hint = 1
        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.refresh_elapsed_s += (time_ns() - refresh_t0_ns) * 1e-9
            end)
        end
        return cache.rhos[1], cache.Ts[1], cache.winds[1]
    catch err
        cache.valid = false
        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.refresh_failures += 1
                s.refresh_elapsed_s += (time_ns() - refresh_t0_ns) * 1e-9
            end)
        end
        if !_gram_track_cache_warning_emitted[]
            _gram_track_cache_warning_emitted[] = true
            @warn "GRAM track cache refresh failed; falling back to direct GRAM sampling." exception=err
        end
        return getDensity(density_model, alt, lat, lon, t, true, p)
    end
end

