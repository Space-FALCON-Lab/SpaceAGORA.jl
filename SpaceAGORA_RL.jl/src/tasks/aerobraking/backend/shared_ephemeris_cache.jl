struct SpaceAGORARLSharedEphemerisCache
    nbody::Any
    srp_sun::Any
    planet_frame::Any
    et_start::Float64
    et_end::Float64
    dt_s::Float64
    sample_count::Int
end

const _SPACEAGORA_RL_SHARED_EPHEMERIS_LOCK = ReentrantLock()
const _SPACEAGORA_RL_SHARED_EPHEMERIS_CACHES = Dict{Any,SpaceAGORARLSharedEphemerisCache}()

function _spaceagora_rl_bool_env(name::AbstractString, default::Bool)
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    raw in ("1", "true", "yes", "on") && return true
    raw in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be a boolean value, got '$(get(ENV, name, ""))'."))
end

_spaceagora_rl_shared_ephemeris_enabled() =
    _spaceagora_rl_bool_env("SPACEAGORA_RL_SHARED_EPHEMERIS_CACHE", true)

function _spaceagora_rl_shared_ephemeris_dt_s()
    raw = strip(get(ENV, "SPACEAGORA_RL_SHARED_EPHEMERIS_DT_S", "30.0"))
    dt_s = tryparse(Float64, raw)
    dt_s !== nothing && isfinite(dt_s) && dt_s > 0.0 ||
        throw(ArgumentError("SPACEAGORA_RL_SHARED_EPHEMERIS_DT_S must be a positive number, got '$raw'."))
    return dt_s::Float64
end

function _spaceagora_rl_shared_ephemeris_max_samples()
    raw = strip(get(ENV, "SPACEAGORA_RL_SHARED_EPHEMERIS_MAX_SAMPLES", "1000000"))
    max_samples = tryparse(Int, raw)
    max_samples !== nothing && max_samples >= 2 ||
        throw(ArgumentError("SPACEAGORA_RL_SHARED_EPHEMERIS_MAX_SAMPLES must be an integer >= 2, got '$raw'."))
    return max_samples::Int
end

function _spaceagora_rl_ephemeris_epoch_window(config::AerobrakingScenarioConfig)
    rand_cfg = config.randomization_config
    if rand_cfg.nominal
        return (start=config.nominal_epoch, latest_start=config.nominal_epoch)
    elseif rand_cfg.initial_date_start !== nothing && rand_cfg.initial_date_days > 0
        start = DateTime(rand_cfg.initial_date_start::Date)
        return (start=start, latest_start=start + Day(rand_cfg.initial_date_days))
    end
    return (start=config.nominal_epoch, latest_start=config.nominal_epoch + Day(28))
end

function _spaceagora_rl_ephemeris_coverage(config::AerobrakingScenarioConfig,
                                           max_passes_per_campaign::Integer;
                                           dt_s::Real=_spaceagora_rl_shared_ephemeris_dt_s())
    window = _spaceagora_rl_ephemeris_epoch_window(config)
    pass_cap = min(max(1, Int(max_passes_per_campaign)),
                   max(1, config.termination_config.max_passes))
    max_apoapsis_radius_m =
        config.initial_apoapsis_radius_m + abs(config.randomization_config.apoapsis_jitter_m)
    max_periapsis_altitude_m = 180.0e3
    max_period_s = orbital_period_s(
        config,
        max_apoapsis_radius_m,
        max_periapsis_altitude_m,
    )
    mission_duration_s = max(2max_period_s, 1.25 * (pass_cap + 1) * max_period_s)
    calendar_span_s = Dates.value(window.latest_start - window.start) / 1.0e3
    total_span_s = calendar_span_s + mission_duration_s
    sample_count = max(2, Int(ceil(total_span_s / Float64(dt_s))) + 1)
    return (
        start=window.start,
        latest_start=window.latest_start,
        mission_duration_s=mission_duration_s,
        total_span_s=total_span_s,
        pass_cap=pass_cap,
        dt_s=Float64(dt_s),
        sample_count=sample_count,
    )
end

function _spaceagora_rl_build_planet_frame_cache_latest(SM,
                                                         planet,
                                                         ephemerides_model,
                                                         ets::Vector{Float64})
    planet_frame_lpi = getproperty(SM, :planet_frame_lpi)
    dcm_to_quaternion = getproperty(SM, :dcm_to_quaternion)
    first_quaternion = dcm_to_quaternion(planet_frame_lpi(planet, ets[1], ephemerides_model))
    quaternions = Vector{typeof(first_quaternion)}(undef, length(ets))
    quaternions[1] = first_quaternion
    @inbounds for idx in 2:length(ets)
        quaternions[idx] = dcm_to_quaternion(
            planet_frame_lpi(planet, ets[idx], ephemerides_model),
        )
    end
    return getproperty(SM, :PlanetFrameEphemerisCache)(ets, quaternions)
end

function _spaceagora_rl_build_shared_ephemeris_cache(config::AerobrakingScenarioConfig,
                                                      coverage)
    max_samples = _spaceagora_rl_shared_ephemeris_max_samples()
    coverage.sample_count <= max_samples ||
        throw(ArgumentError(
            "RL shared ephemeris cache requires $(coverage.sample_count) samples, exceeding " *
            "SPACEAGORA_RL_SHARED_EPHEMERIS_MAX_SAMPLES=$max_samples."
        ))

    spaceagora = _load_spaceagora!(; load_gramsuite=_spaceagora_live_needs_gramsuite(config))
    SM = getproperty(spaceagora, :SimulationModel)
    engine = getproperty(spaceagora, :SimulationEngine)
    runtime_services = getproperty(spaceagora, :RuntimeServices)
    planet = deepcopy(_spaceagora_mars(spaceagora, _spaceagora_spice_path()))
    ephemerides_model = Base.invokelatest(getproperty(SM, :SpiceEphemeridesModel))
    et_start = Base.invokelatest(
        getproperty(SM, :ephemerides_time_seconds),
        _initial_time_from_datetime(coverage.start; spaceagora=spaceagora),
        ephemerides_model,
    )
    latest_start_et = Base.invokelatest(
        getproperty(SM, :ephemerides_time_seconds),
        _initial_time_from_datetime(coverage.latest_start; spaceagora=spaceagora),
        ephemerides_model,
    )
    et_end = latest_start_et + coverage.mission_duration_s
    span_s = et_end - et_start
    sample_count = max(2, Int(ceil(span_s / coverage.dt_s)) + 1)
    sample_count <= max_samples ||
        throw(ArgumentError(
            "RL shared ephemeris cache requires $sample_count ET samples after time conversion, " *
            "exceeding SPACEAGORA_RL_SHARED_EPHEMERIS_MAX_SAMPLES=$max_samples."
        ))

    @printf(
        "precomputing RL shared ephemeris cache start=%s latest_campaign_start=%s span_days=%.2f dt_s=%.1f samples=%d\n",
        string(coverage.start),
        string(coverage.latest_start),
        span_s / 86400.0,
        coverage.dt_s,
        sample_count,
    )
    flush(stdout)
    started = time()

    nbody_cache, planet_frame_cache = lock(getproperty(runtime_services, :SPICE_LOCK)) do
        nbody = Base.invokelatest(
            getproperty(engine, :_build_nbody_ephemeris_cache),
            "mars_barycenter",
            ["sun"],
            et_start,
            span_s,
            coverage.dt_s,
        )
        frame = Base.invokelatest(
            _spaceagora_rl_build_planet_frame_cache_latest,
            SM,
            planet,
            ephemerides_model,
            nbody.ets,
        )
        return nbody, frame
    end

    sun_positions = copy(vec(nbody_cache.positions_j2000_m[:, 1]))
    srp_cache = Base.invokelatest(
        getproperty(SM, :SRPSunEphemerisCache),
        nbody_cache.ets,
        sun_positions,
    )
    cache = SpaceAGORARLSharedEphemerisCache(
        nbody_cache,
        srp_cache,
        planet_frame_cache,
        et_start,
        et_end,
        coverage.dt_s,
        length(nbody_cache.ets),
    )
    @printf(
        "RL shared ephemeris cache ready elapsed=%.1fs samples=%d\n",
        time() - started,
        cache.sample_count,
    )
    flush(stdout)
    return cache
end

function _spaceagora_rl_shared_ephemeris_cache(config::AerobrakingScenarioConfig,
                                                max_passes_per_campaign::Integer;
                                                build_if_missing::Bool=true)
    _spaceagora_rl_shared_ephemeris_enabled() || return nothing
    coverage = _spaceagora_rl_ephemeris_coverage(config, max_passes_per_campaign)
    key = (
        coverage.start,
        coverage.latest_start,
        coverage.mission_duration_s,
        coverage.dt_s,
        normpath(_spaceagora_spice_path()),
    )
    return lock(_SPACEAGORA_RL_SHARED_EPHEMERIS_LOCK) do
        if build_if_missing
            return get!(_SPACEAGORA_RL_SHARED_EPHEMERIS_CACHES, key) do
                _spaceagora_rl_build_shared_ephemeris_cache(config, coverage)
            end
        end
        return get(_SPACEAGORA_RL_SHARED_EPHEMERIS_CACHES, key, nothing)
    end
end

function _prewarm_spaceagora_rl_shared_ephemeris_cache!(
    config::AerobrakingScenarioConfig,
    max_passes_per_campaign::Integer,
)
    _is_spaceagora_live_backend(config.backend_mode) || return nothing
    return _spaceagora_rl_shared_ephemeris_cache(config, max_passes_per_campaign)
end

function _install_spaceagora_rl_shared_ephemeris_cache!(p,
                                                         cache::SpaceAGORARLSharedEphemerisCache)
    mission_end_et = p.shared_buffers.et_start[] +
                     Float64(p.args.mission_configuration.mission_time)
    tolerance_s = max(1.0e-6, 16eps(max(abs(cache.et_end), 1.0)))
    if p.shared_buffers.et_start[] < cache.et_start - tolerance_s ||
       mission_end_et > cache.et_end + tolerance_s
        throw(ArgumentError(
            "RL shared ephemeris cache does not cover campaign ET interval " *
            "[$(p.shared_buffers.et_start[]), $mission_end_et]; cache interval is " *
            "[$(cache.et_start), $(cache.et_end)]."
        ))
    end
    p.shared_buffers.nbody_ephemeris_cache[] = cache.nbody
    p.shared_buffers.srp_sun_ephemeris_cache[] = cache.srp_sun
    p.shared_buffers.planet_frame_ephemeris_cache[] = cache.planet_frame
    growth_state = p.shared_buffers.ephemeris_cache_growth_state
    growth_state.enabled = false
    growth_state.next_growth_t_s = Inf
    return nothing
end

function _spaceagora_rl_shared_ephemeris_callback(spaceagora,
                                                   config::AerobrakingScenarioConfig,
                                                   max_passes_per_campaign::Integer)
    cache = _spaceagora_rl_shared_ephemeris_cache(
        config,
        max_passes_per_campaign;
        build_if_missing=false,
    )
    cache === nothing && return nothing
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    discrete_callback = getproperty(callbacks, :DiscreteCallback)
    condition(u, t, integrator) = false
    affect!(integrator) = nothing
    initialize = (cb, u, t, integrator) ->
        _install_spaceagora_rl_shared_ephemeris_cache!(integrator.p, cache)
    return Base.invokelatest(discrete_callback, condition, affect!; initialize=initialize)
end

function _spaceagora_rl_core_ephemeris_env_pairs()
    _spaceagora_rl_shared_ephemeris_enabled() || return Pair{String,String}[]
    return Pair{String,String}[
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "0",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE" => "0",
        "SPACEAGORA_PLANET_FRAME_CACHE" => "0",
    ]
end
