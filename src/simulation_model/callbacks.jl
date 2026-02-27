module SimulationCallbacks
include("../utils/Reference_system.jl") # Get the reference system types for the callback

using DifferentialEquations
using LinearAlgebra
using SPICE
using Dates
using ..SimulationModel: SPICE_LOCK, GRAM_LOCK
using ..EnvironmentModels
using ..EnvironmentModels: getDensity, getHeatRate, NoAtmosphereModel
using ..DynamicEffectors: BaseThrusterModel, AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
using ..ConfigTypes: SaveData
using ..ControlEffectors: calcControlEffect!
using ..GuidanceEffectors: calcGuidanceEffect!
using ..SimConfig: SimulationConfiguration, MissionOrbits
export get_callbacks

@inline callback_verbose(integrator) = integrator.p.args.simulation_settings.verbose
@inline callback_use_invokelatest() = get(ENV, "SPACEAGORA_DEV_HOT_RELOAD", "0") == "1"
const _gram_track_cache_warning_emitted = Ref(false)

@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _density_callback_parallel_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_DENSITY_CALLBACK_PARALLEL", "auto")))
    if mode in ("off", "none", "serial", "0", "false", "no")
        return :off
    elseif mode in ("on", "thread", "threads", "1", "true", "yes")
        return :on
    elseif mode == "auto"
        return :auto
    end
    throw(ArgumentError("Unsupported SPACEAGORA_DENSITY_CALLBACK_PARALLEL='$mode'. Use one of: off, auto, on."))
end

@inline function _density_callback_thread_threshold()::Int
    raw = strip(get(ENV, "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD", "8"))
    threshold = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD must be an integer, got '$raw'"))
    end
    return max(1, threshold)
end

@inline function _control_callback_parallel_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_CONTROL_CALLBACK_PARALLEL", "auto")))
    if mode in ("off", "none", "serial", "0", "false", "no")
        return :off
    elseif mode in ("on", "thread", "threads", "1", "true", "yes")
        return :on
    elseif mode == "auto"
        return :auto
    end
    throw(ArgumentError("Unsupported SPACEAGORA_CONTROL_CALLBACK_PARALLEL='$mode'. Use one of: off, auto, on."))
end

@inline function _control_callback_thread_threshold()::Int
    raw = strip(get(ENV, "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD", "8"))
    threshold = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD must be an integer, got '$raw'"))
    end
    return max(1, threshold)
end

# Extend this for custom user density models as needed:
# SimulationModel.SimulationCallbacks.density_model_threadsafe(::MyDensityModel) = true
@inline density_model_threadsafe(::AbstractDensityModel)::Bool = false
@inline density_model_threadsafe(::NoAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.ExponentialAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.PolynomialFitAtmosphereModel)::Bool = true
# GRAM C-wrapper calls are serialized inside getDensity via SimulationModel.GRAM_LOCK.
@inline density_model_threadsafe(::EnvironmentModels.GRAMAtmosphereModel)::Bool = true

@inline function _density_callback_use_threads(args::SimulationConfiguration, num_sats::Int)::Bool
    if Threads.nthreads() <= 1 || num_sats <= 1
        return false
    end
    mode = _density_callback_parallel_mode()
    mode == :off && return false

    model = args.environment_model.density_model
    model_threadsafe = density_model_threadsafe(model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_ASSUME_THREADSAFE", false)
        return false
    end
    if mode == :on
        return true
    end
    return num_sats >= _density_callback_thread_threshold()
end

# Extend this for custom user control models as needed:
# SimulationModel.SimulationCallbacks.control_model_threadsafe(::MyControlModel) = true
@inline control_model_threadsafe(::Any)::Bool = false
@inline control_model_threadsafe(::BaseThrusterModel)::Bool = true

@inline function _control_callback_use_threads(control_model, num_sats::Int, use_invokelatest::Bool)::Bool
    if use_invokelatest || Threads.nthreads() <= 1 || num_sats <= 1
        return false
    end
    mode = _control_callback_parallel_mode()
    mode == :off && return false

    model_threadsafe = control_model_threadsafe(control_model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE", false)
        return false
    end
    if mode == :on
        return true
    end
    return num_sats >= _control_callback_thread_threshold()
end

@inline function _density_model_for_sat(p, sat_idx::Int)
    density_models = p.shared_buffers.density_models
    if sat_idx <= length(density_models)
        return density_models[sat_idx]
    end
    return p.args.environment_model.density_model
end

mutable struct GramTrackCache
    valid::Bool
    t0::Float64
    t1::Float64
    index_hint::Int
    times::Vector{Float64}
    alts::Vector{Float64}
    lats::Vector{Float64}
    lons::Vector{Float64}
    rhos::Vector{Float64}
    Ts::Vector{Float64}
    winds::Vector{SVector{3, Float64}}
end

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
    raw = haskey(ENV, "SPACEAGORA_GRAM_TRACK_CACHE") ?
        ENV["SPACEAGORA_GRAM_TRACK_CACHE"] :
        get(ENV, "SPACEAGORA_GRAM_SEGMENT_CACHE", "auto")
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
    legacy_horizon = _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_HORIZON_S")
    legacy_alt_tol = _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ALT_TOL_M")
    legacy_ang_tol_deg = _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ANG_TOL_DEG")
    legacy_points = _parse_int_env_optional("SPACEAGORA_GRAM_TRACK_CACHE_NPOS")
    if legacy_points === nothing
        legacy_points = _parse_int_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_POINTS")
    end

    entry_horizon_s = max(
        1e-3,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_HORIZON_S"),
            something(legacy_horizon, 1.0)
        )
    )
    orbit_horizon_s = max(
        1e-3,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_HORIZON_S"),
            something(legacy_horizon, 8.0)
        )
    )
    entry_points = max(
        2,
        something(
            _parse_int_env_optional("SPACEAGORA_GRAM_TRACK_CACHE_ENTRY_NPOS"),
            something(
                _parse_int_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_POINTS"),
                something(legacy_points, 16)
            )
        )
    )
    orbit_points = max(
        2,
        something(
            _parse_int_env_optional("SPACEAGORA_GRAM_TRACK_CACHE_ORBIT_NPOS"),
            something(
                _parse_int_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_POINTS"),
                something(legacy_points, 48)
            )
        )
    )
    entry_alt_tol_m = max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_ALT_TOL_M"),
            something(legacy_alt_tol, 500.0)
        )
    )
    orbit_alt_tol_m = max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_ALT_TOL_M"),
            something(legacy_alt_tol, 3000.0)
        )
    )
    entry_ang_tol_rad = deg2rad(max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ENTRY_ANG_TOL_DEG"),
            something(legacy_ang_tol_deg, 0.6)
        )
    ))
    orbit_ang_tol_rad = deg2rad(max(
        0.0,
        something(
            _parse_float_env_optional("SPACEAGORA_GRAM_SEGMENT_CACHE_ORBIT_ANG_TOL_DEG"),
            something(legacy_ang_tol_deg, 4.0)
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
    return density_model isa EnvironmentModels.GRAMAtmosphereModel
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
    if !cache.valid || n < 2 || t < cache.t0 || t > cache.t1
        return nothing
    end
    idx = clamp(cache.index_hint, 1, n - 1)
    if t < cache.times[idx] || t > cache.times[idx + 1]
        idx = clamp(searchsortedlast(cache.times, t), 1, n - 1)
    end
    t_lo = cache.times[idx]
    t_hi = cache.times[idx + 1]
    if t < t_lo || t > t_hi
        return nothing
    end
    x = t_hi == t_lo ? 0.0 : (t - t_lo) / (t_hi - t_lo)
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

function _gram_descending_periapsis_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet
)::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}}
    try
        oe = rvtoorbitalelement(pos, vel, planet)
        a, e, i, Ω, ω, ν = Float64(oe[1]), Float64(oe[2]), Float64(oe[3]), Float64(oe[4]), Float64(oe[5]), Float64(oe[6])
        if !isfinite(a) || !isfinite(e) || !isfinite(ν) || a <= 0.0 || e < 0.0 || e >= 1.0
            return nothing
        end
        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            return nothing
        end
        E = atan(sqrt(max(0.0, 1.0 - e^2)) * sin(ν), e + cos(ν))
        E = mod(E, 2pi)
        M = mod(E - e * sin(E), 2pi)
        dt_peri = (2pi - M) / n
        if !isfinite(dt_peri) || dt_peri <= 1e-6
            return nothing
        end

        oe_peri = SVector{7, Float64}(a, e, i, Ω, ω, 0.0, 0.0)
        r_peri_vec, v_peri_vec = orbitalelemtorv(oe_peri, planet)
        r_peri_i = SVector{3, Float64}(Float64.(r_peri_vec))
        v_peri_i = SVector{3, Float64}(Float64.(v_peri_vec))
        r_peri_p, _ = r_intor_p!(r_peri_i, v_peri_i, planet)
        alt_peri, lat_peri, lon_peri = rtolatlong(r_peri_p, planet)
        return dt_peri, alt_peri, lat_peri, lon_peri
    catch
        return nothing
    end
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
    transition_band_m::Float64 = 0.0
)
    try
        planet = p.args.environment_model.planet
        base_horizon_s = max(1e-3, horizon_s)
        dt_segment = base_horizon_s
        pos_target = pos + vel * dt_segment
        rp_target, _ = r_intor_p!(pos_target, vel, planet)
        target_alt, target_lat, target_lon = rtolatlong(rp_target, planet)

        if _gram_track_cache_periapsis_split_enabled() &&
           alt <= p.args.environment_model.EI * 1e3 + max(0.0, transition_band_m) &&
           dot(pos, vel) < 0.0
            peri_target = _gram_descending_periapsis_target(pos, vel, planet)
            if peri_target !== nothing
                dt_segment, target_alt, target_lat, target_lon = peri_target
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
        if density_model isa EnvironmentModels.GRAMAtmosphereModel &&
           isdefined(density_model.gram, :generate_trajectory)
            denom = max(1, n - 1)
            Δalt_km = (target_alt - alt) * 1e-3 / denom
            Δlat_deg = rad2deg(_angle_delta_rad(lat, target_lat)) / denom
            Δlon_deg = rad2deg(_angle_delta_rad(lon, target_lon)) / denom
            Δt = dt_segment / denom

            gen_track = () -> density_model.gram.generate_trajectory(
                density_model.gram_atmosphere;
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
                ρ_k, T_k, wind_k = getDensity(density_model, alt_k, lat_k, lon_k, t_k, true, p)
                cache.times[k] = t_k
                cache.alts[k] = alt_k
                cache.lats[k] = lat_k
                cache.lons[k] = lon_k
                cache.rhos[k] = ρ_k
                cache.Ts[k] = T_k
                cache.winds[k] = wind_k
            end
        end
        cache.valid = true
        cache.t0 = cache.times[1]
        cache.t1 = cache.times[end]
        cache.index_hint = 1
        return cache.rhos[1], cache.Ts[1], cache.winds[1]
    catch err
        cache.valid = false
        if !_gram_track_cache_warning_emitted[]
            _gram_track_cache_warning_emitted[] = true
            @warn "GRAM track cache refresh failed; falling back to direct GRAM sampling." exception=err
        end
        return getDensity(density_model, alt, lat, lon, t, true, p)
    end
end

@inline function _uses_atmospheric_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
            return true
        end
    end
    return false
end

@inline function _requires_density_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _uses_atmospheric_dynamic_effector(effectors) || !(args.environment_model.density_model isa NoAtmosphereModel)
end

@inline function _requires_orbit_end_callback(args::SimulationConfiguration)::Bool
    return args.mission_configuration.mission_type == MissionOrbits
end

@inline function _requires_drag_state_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    if !_requires_density_callback(effectors, args)
        return false
    end
    tol = args.integration_tolerances
    return tol.dt_max_atmosphere != tol.dt_max_orbit ||
           tol.reltol_atmosphere != tol.reltol_orbit ||
           tol.abstol_atmosphere != tol.abstol_orbit
end

@inline function _requires_thermal_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    # Thermal-rate evaluation requires atmospheric state (ρ, T, wind), populated by density callback.
    return _requires_density_callback(effectors, args)
end

@inline function _resolved_component_tolerance(component_tol::Float64, baseline_tol::Float64)::Float64
    return component_tol == 0.0 ? baseline_tol : component_tol
end

function _callback_tolerances_for_phase(template_reltol, template_abstol, args::SimulationConfiguration, in_atmosphere::Bool)
    tol = args.integration_tolerances
    baseline_reltol = in_atmosphere ? tol.reltol_atmosphere : tol.reltol_orbit
    baseline_abstol = in_atmosphere ? tol.abstol_atmosphere : tol.abstol_orbit
    if template_reltol isa Number && template_abstol isa Number
        return baseline_reltol, baseline_abstol
    end

    reltol_mass = _resolved_component_tolerance(tol.reltol_mass, baseline_reltol)
    abstol_mass = _resolved_component_tolerance(tol.abstol_mass, baseline_abstol)
    reltol_heat = _resolved_component_tolerance(tol.reltol_heat_load, baseline_reltol)
    abstol_heat = _resolved_component_tolerance(tol.abstol_heat_load, baseline_abstol)
    reltol_ω = _resolved_component_tolerance(tol.reltol_angular_rate, baseline_reltol)
    abstol_ω = _resolved_component_tolerance(tol.abstol_angular_rate, baseline_abstol)

    reltol_new = copy(template_reltol)
    abstol_new = copy(template_abstol)
    reltol_new .= baseline_reltol
    abstol_new .= baseline_abstol

    @inbounds for i in eachindex(reltol_new.sc)
        reltol_new.sc[i].mass = reltol_mass
        abstol_new.sc[i].mass = abstol_mass
        reltol_new.sc[i].heat_loads .= reltol_heat
        abstol_new.sc[i].heat_loads .= abstol_heat
        if hasproperty(reltol_new.sc[i], :q)
            reltol_new.sc[i].q .= tol.reltol_quaternion
            abstol_new.sc[i].q .= tol.abstol_quaternion
        end
        if hasproperty(reltol_new.sc[i], :ω)
            reltol_new.sc[i].ω .= reltol_ω
            abstol_new.sc[i].ω .= abstol_ω
        end
    end
    return reltol_new, abstol_new
end

function get_callbacks(num_sats::Int, effectors::Tuple, args::SimulationConfiguration)::CallbackSet
    callbacks = Any[
        get_impact_callback(num_sats),
        update_planet_frame_callback()
    ]

    if _requires_density_callback(effectors, args)
        push!(callbacks, get_density_callback(num_sats, args))
    end

    if _requires_thermal_callback(effectors, args)
        push!(callbacks, get_thermal_callback(num_sats, args))
    end

    if _requires_orbit_end_callback(args)
        push!(callbacks, get_orbit_end_callback(num_sats))
    end

    if _requires_drag_state_callback(effectors, args)
        push!(callbacks, get_drag_state_callback(num_sats))
    end

    append!(callbacks, get_control_callbacks(num_sats, args))
    append!(callbacks, get_guidance_callbacks(num_sats, args))

    return CallbackSet(callbacks...)
end
function update_planet_frame_callback()
    condition(u, t, integrator) = true
    et_start = Ref{Float64}(0.0) # Store the start epoch in a Ref so it can be updated in the affect! function
    function affect!(integrator)
        p = integrator.p
        planet = p.args.environment_model.planet
        et = et_start[] + integrator.t
        lock(SPICE_LOCK) do
            planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(planet.name)", et)) * planet.J2000_to_pci'
        end
    end

    function init_affect!(cb, u, t, integrator)
        p = integrator.p
        initial_time = p.args.initial_time
        start_epoch = from_utc(DateTime(
            initial_time.year,
            initial_time.month,
            initial_time.day,
            initial_time.hour,
            initial_time.minute,
            initial_time.second
        ))
        lock(SPICE_LOCK) do
            et_start[] = utc2et(to_utc(start_epoch))
        end
        affect!(integrator) # Call the affect! function at the start of the simulation to initialize the planet frame
    end

    return DiscreteCallback(condition, affect!, initialize=init_affect!)
end
# Factory function to build the callback
function get_density_callback(num_sats::Int, args::SimulationConfiguration)
    use_threads = _density_callback_use_threads(args, num_sats)
    cache_cfg = _gram_track_cache_config()
    function update_density_sat!(i::Int, p, u, t)
        pos = SVector{3, Float64}(u.sc[i].pos)
        vel = SVector{3, Float64}(u.sc[i].vel)
        rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet)
        alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet)
        density_model = _density_model_for_sat(p, i)
        ρ, T, wind_vec = if _gram_track_cache_enabled(cache_cfg, density_model)
            cache_horizon_s, cache_alt_tol_m, cache_ang_tol_rad, cache_points = _gram_track_cache_profile(cache_cfg, p, alt)
            caches = p.shared_buffers.gram_density_cache
            cache = if i <= length(caches) && caches[i] isa GramTrackCache
                caches[i]::GramTrackCache
            else
                c = GramTrackCache()
                if i <= length(caches)
                    caches[i] = c
                end
                c
            end
            seg = _gram_track_cache_ready(cache, t, alt, lat, lon, cache_alt_tol_m, cache_ang_tol_rad)
            if seg !== nothing
                idx, x = seg
                _gram_track_cache_eval(cache, idx, x)
            else
                _gram_track_cache_refresh!(
                    cache,
                    density_model,
                    p,
                    pos,
                    vel,
                    alt,
                    lat,
                    lon,
                    t,
                    cache_horizon_s,
                    cache_points,
                    cache_alt_tol_m,
                    cache_ang_tol_rad,
                    cache_cfg.transition_band_m
                )
            end
        else
            getDensity(density_model, alt, lat, lon, t, true, p)
        end
        p.shared_buffers.densities[i] = ρ
        p.shared_buffers.temperatures[i] = T
        p.shared_buffers.winds[i] = wind_vec
        return nothing
    end
    
    # Logic for when the callback should trigger
    # Returning true means it runs at every integration step
    condition(u, t, integrator) = true

    # Logic for what happens when it triggers
    function affect!(integrator)
        p = integrator.p
        u = integrator.u

        if use_threads
            Threads.@threads for i in 1:num_sats
                @inbounds update_density_sat!(i, p, u, integrator.t)
            end
        else
            @inbounds for i in 1:num_sats
                update_density_sat!(i, p, u, integrator.t)
            end
        end
    end

    return DiscreteCallback(condition, affect!, initialize=(cb, u, t, integrator) -> affect!(integrator))
end

function get_thermal_callback(num_sats::Int, args::SimulationConfiguration)
    use_threads = _density_callback_use_threads(args, num_sats)

    function update_thermal_sat!(i::Int, p, u)
        links = p.args.dynamics_model.spacecraft[i].links
        n_links = length(links)
        heat_rates = p.shared_buffers.heat_rates[i]
        if length(heat_rates) != n_links
            resize!(heat_rates, n_links)
        end
        fill!(heat_rates, 0.0)

        ρ = p.shared_buffers.densities[i]
        T = p.shared_buffers.temperatures[i]
        wind = p.shared_buffers.winds[i]
        if !isfinite(ρ) || !isfinite(T) || ρ <= 0.0 || T <= 0.0
            return nothing
        end

        planet = p.args.environment_model.planet
        thermal_model = p.args.environment_model.thermal_model
        pos_ii = SVector{3, Float64}(u.sc[i].pos)
        vel_ii = SVector{3, Float64}(u.sc[i].vel)
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, planet)
        alt_lat_lon = rtolatlong(pos_pp, planet)
        uD, uN, uE = latlongtoNED(alt_lat_lon)
        wE, wN, wU = wind
        wind_pp = wN * uN + wE * uE - wU * uD
        vel_pp_rw = vel_pp + wind_pp
        v = norm(vel_pp_rw)
        sound_velocity = sqrt(planet.γ * planet.R * T)
        if !isfinite(v) || !isfinite(sound_velocity) || v <= 0.0 || sound_velocity <= 0.0
            return nothing
        end
        mach = v / sound_velocity
        S = sqrt(planet.γ * 0.5) * mach
        @inbounds for j in eachindex(links)
            α = links[j].α
            if !isfinite(α)
                continue
            end
            qdot = getHeatRate(thermal_model, S, T, ρ, v, α)
            heat_rates[j] = (isfinite(qdot) && qdot > 0.0) ? qdot : 0.0
        end
        return nothing
    end

    condition(u, t, integrator) = true

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        if use_threads
            Threads.@threads for i in 1:num_sats
                @inbounds update_thermal_sat!(i, p, u)
            end
        else
            @inbounds for i in 1:num_sats
                update_thermal_sat!(i, p, u)
            end
        end
    end

    return DiscreteCallback(condition, affect!, initialize=(cb, u, t, integrator) -> affect!(integrator))
end

function get_impact_callback(num_sats::Int)
    function condition(u, t, integrator)
        p = integrator.p
        # Check for impacts based on the current state
        @inbounds for i in 1:num_sats
            if p.is_active[i] && abs(norm(u.sc[i].pos) - p.args.environment_model.planet.Rp_e) < 50e3 # Example condition for impact (e.g., altitude below 80 km)
                return true
            end
        end
        return false
    end

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        
        # Check for impacts and log them (placeholder logic)
        @inbounds for i in 1:num_sats
            if p.is_active[i] && abs(norm(u.sc[i].pos) - p.args.environment_model.planet.Rp_e) < 50e3 # Example condition for impact (e.g., altitude below 80 km)
                if callback_verbose(integrator)
                    println("Impact detected for satellite $i at time $(integrator.t) seconds!")
                end
                p.is_active[i] = false # Mark the satellite as inactive after impact
                if all(p.is_active .== false) # If all satellites are inactive, we can stop the simulation
                    if callback_verbose(integrator)
                        println("All satellites have impacted. Stopping simulation.")
                    end
                    terminate!(integrator)
                end
                # Log impact details to a file or data structure as needed
            end
        end
    end

    return DiscreteCallback(condition, affect!)
end

# Only call if simulating a single satellite
function get_orbit_end_callback(num_sats::Int)
    function condition!(out, u, t, integrator)
        # Use radial-velocity root events for orbit bookkeeping.
        # At apoapsis, dot(r,v) crosses + -> -, so -dot(r,v) crosses - -> + (upcrossing),
        # which matches the single affect! handler below.
        @inbounds for i in 1:num_sats
            pos = SVector{3, Float64}(u.sc[i].pos)
            vel = SVector{3, Float64}(u.sc[i].vel)
            out[i] = -dot(pos, vel)
        end
    end

    function affect!(integrator, idx::Int64)
        p = integrator.p

        p.orbit_counter[idx] += 1 # Increment the orbit counter in the shared buffers
        completed_orbits = p.orbit_counter[idx] - 1
        target_orbits = p.args.mission_configuration.number_of_orbits
        # Check for orbit completion and log it (placeholder logic)
        if callback_verbose(integrator)
            println("Orbit $(completed_orbits) completed by Satellite $idx at time $(integrator.t) seconds!")
        end

        # MissionOrbits should end when every active satellite reaches the requested
        # number of completed orbits, analogous to drag-passage event termination.
        if completed_orbits >= target_orbits
            all_active_reached_target = true
            @inbounds for sat_idx in eachindex(p.orbit_counter)
                if p.is_active[sat_idx] && (p.orbit_counter[sat_idx] - 1) < target_orbits
                    all_active_reached_target = false
                    break
                end
            end
            if all_active_reached_target
                if callback_verbose(integrator)
                    println("Target orbit count reached for all active satellites. Stopping simulation.")
                end
                if applicable(terminate!, integrator)
                    terminate!(integrator)
                end
            end
        end
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
end

# Only call if simulating single satellite
# Callback to adjust the max timestep depending on whether the satellite is in the atmosphere
function get_drag_state_callback(num_sats::Int)
    # out = zeros(num_sats) # Output array to store the condition for each satellite
    condition!(out, u, t, integrator) = begin
        @inbounds for i in 1:length(u.sc)
            alt = norm(u.sc[i].pos) - integrator.p.args.environment_model.planet.Rp_e
            # println("Satellite $i altitude: $(alt) meters at time $(integrator.t) seconds")
            out[i] = alt - integrator.p.args.environment_model.EI*1e3 # Positive when above the atmosphere, negative when in the atmosphere
        end
        # return out
    end
    function affect_upcrossing!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u
        if callback_verbose(integrator)
            println("Switching to space integration at time $(integrator.t) seconds!")
        end
        # sleep(3.0)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_orbit # Increase the maximum timestep when exiting the atmosphere
        reltol_new, abstol_new = _callback_tolerances_for_phase(
            integrator.opts.reltol,
            integrator.opts.abstol,
            p.args,
            false
        )
        integrator.opts.reltol = reltol_new # Adjust tolerances when exiting the atmosphere
        integrator.opts.abstol = abstol_new
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u
        if callback_verbose(integrator)
            println("Switching to atmosphere integration at time $(integrator.t) seconds!")
        end
        # sleep(3.0)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_atmosphere # Decrease the maximum timestep when entering the atmosphere
        reltol_new, abstol_new = _callback_tolerances_for_phase(
            integrator.opts.reltol,
            integrator.opts.abstol,
            p.args,
            true
        )
        integrator.opts.reltol = reltol_new # Adjust tolerances when entering the atmosphere
        integrator.opts.abstol = abstol_new
    end

    return VectorContinuousCallback(condition!, affect_upcrossing!, affect_downcrossing!, num_sats)
end

function get_data_saving_callback(num_sats::Int)
    # Implement a callback to save the state of the simulation at regular intervals defined by num_steps_to_save in the mission configuration
    # This will involve storing the state in the shared buffers and writing to a file when the buffer is full, then clearing the buffer for the next set of data
    # The condition can be based on the number of steps taken or the simulation time, depending on how you want to structure it
    saved_values = SaveValues(Float64, SaveData)
    function save_func(u, t, integrator)
        p = integrator.p

        # Save the current state to the SaveData struct

    end
    return SavingCallback(save_func, saved_values) # Placeholder, implement the actual logic for the data saving callback
end

function get_guidance_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    # Implement a callback to calculate guidance commands at each time step based on the current state and the guidance model defined in the simulation configuration
    guidance_models = args.guidance_model.guidance_effectors
    guidance_rates = args.guidance_model.guidance_rates
    use_invokelatest = callback_use_invokelatest()
    callbacks = Vector{DiscreteCallback}(undef, length(guidance_models))
    for i in eachindex(guidance_models)
        guidance_model = guidance_models[i]
        guidance_rate = guidance_rates[i]
        # Implement a callback for this guidance model that triggers at the specified guidance rate and calculates the guidance commands based on the current state and the guidance model
        # The calculated guidance commands should be stored in the shared buffers for use in the dynamics calculations
        guidance_func = (integrator) -> begin
            if use_invokelatest
                # Dev mode: keep Revise/hot-reload workflows free of world-age errors.
                Base.invokelatest(calcGuidanceEffect!, guidance_model, integrator.u, integrator.p, integrator.t, i)
            else
                # Production mode: direct dispatch avoids invokelatest overhead.
                calcGuidanceEffect!(guidance_model, integrator.u, integrator.p, integrator.t, i)
            end
        end
        callbacks[i] = PeriodicCallback(guidance_func, guidance_rate)
    end
    return callbacks
end

function get_control_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    # Perform the control effects' calculations at specific rates given by the control_rates field in the ControlModel
    control_models = args.control_model.control_effectors
    control_rates = args.control_model.control_rates
    use_invokelatest = callback_use_invokelatest()
    callbacks = Vector{DiscreteCallback}(undef, length(control_models))
    for i in eachindex(control_models)
        control_model = control_models[i]
        control_rate = control_rates[i]
        use_threads = _control_callback_use_threads(control_model, num_sats, use_invokelatest)
        if control_model isa BaseThrusterModel
            n_slots = length(control_model.thrust)
            if n_slots != num_sats
                throw(ArgumentError(
                    "BaseThrusterModel vector length ($n_slots) must match number of spacecraft ($num_sats). " *
                    "Use one shared model with per-spacecraft vectors."
                ))
            end
        end
        # Each control effector callback runs at its own rate and updates
        # all spacecraft states. The spacecraft index is passed explicitly
        # to avoid conflating effector-index with spacecraft-index.
        apply_control! = if use_invokelatest
            # Dev mode: keep Revise/hot-reload workflows free of world-age errors.
            (integrator, sat_idx) -> Base.invokelatest(calcControlEffect!, control_model, integrator.u, integrator.p, integrator.t, sat_idx)
        else
            # Production mode: direct dispatch avoids invokelatest overhead.
            (integrator, sat_idx) -> calcControlEffect!(control_model, integrator.u, integrator.p, integrator.t, sat_idx)
        end
        control_func = (integrator) -> begin
            if use_threads
                Threads.@threads for sat_idx in 1:num_sats
                    @inbounds apply_control!(integrator, sat_idx)
                end
            else
                @inbounds for sat_idx in 1:num_sats
                    apply_control!(integrator, sat_idx)
                end
            end
        end
        callbacks[i] = PeriodicCallback(control_func, control_rate)
    end
    return callbacks
end

function get_periapsis_save_callback(num_sats::Int)
    # Implement a callback to save the state of the simulation at periapsis for each orbit, which can be useful for analyzing the changes in the orbit after each pass through the atmosphere
    function condition!(out, u, t, integrator)
        @inbounds for i in 1:num_sats
            OE = rvtoorbitalelement(SVector{3, Float64}(u.sc[i].pos), SVector{3, Float64}(u.sc[i].vel), integrator.p.args.environment_model.planet)
            out[i] = OE[6] # Return the true anomaly (ν) which is zero at periapsis
        end
    end

    function affect!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u

        if callback_verbose(integrator)
            println("Periapsis reached for Satellite $idx at time $(integrator.t) seconds!")
        end
        # Save the state at periapsis to a file or data structure as needed
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
end
end # module
