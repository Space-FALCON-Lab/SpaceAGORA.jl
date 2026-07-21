struct MarsOdysseyPhaseConstants
    phase::String
    final_apoapsis_radius_m::Float64
    r_norm_m::Float64
    initial_apoapsis_radius_m::Float64
    nominal_periapsis_altitude_m::Float64
    nominal_inclination_deg::Float64
    nominal_raan_deg::Float64
    nominal_argument_of_periapsis_deg::Float64
    nominal_epoch::DateTime
    mars_radius_m::Float64
end

function mars_odyssey_phase_constants(phase::AbstractString)
    key = lowercase(String(phase))
    mars_radius = 3.3962e6
    if key == "walkout"
        return MarsOdysseyPhaseConstants("Walkout", 3900e3, 5000e3, 4906e3, 86e3,
                                         93.6, 115.0, 66.0, DateTime(2002, 1, 3), mars_radius)
    elseif key == "endgame"
        return MarsOdysseyPhaseConstants("Endgame", 4906e3, 10100e3, 6000e3, 86e3,
                                         93.6, 115.0, 66.0, DateTime(2001, 12, 26), mars_radius)
    elseif key == "campaign"
        return MarsOdysseyPhaseConstants("Campaign", 3900e3, 30000e3, 28533e3, 87.5e3,
                                         93.6, 114.0, 109.0, DateTime(2001, 11, 6), mars_radius)
    elseif key == "main"
        return MarsOdysseyPhaseConstants("Main", 3900e3, 10100e3, 10038e3, 92e3,
                                         93.6, 115.0, 89.0, DateTime(2001, 12, 18), mars_radius)
    else
        throw(ArgumentError("unknown aerobraking phase: $phase"))
    end
end

struct AerobrakingScenarioConfig
    phase::String
    final_apoapsis_radius_m::Float64
    r_norm_m::Float64
    initial_apoapsis_radius_m::Float64
    nominal_periapsis_altitude_m::Float64
    nominal_inclination_rad::Float64
    nominal_raan_rad::Float64
    nominal_argument_of_periapsis_rad::Float64
    nominal_epoch::DateTime
    mars_radius_m::Float64
    j2::Float64
    mu_m3_s2::Float64
    rho_ref_kg_m3::Float64
    h_ref_m::Float64
    scale_height_m::Float64
    gas_constant_j_kg_k::Float64
    gamma::Float64
    temperature_k::Float64
    g_ref_m_s2::Float64
    base_apoapsis_decay_m::Float64
    heat_decay_gain_m::Float64
    heat_nominal_altitude_m::Float64
    heat_nominal_velocity_mps::Float64
    heat_nominal_w_cm2::Float64
    backend_mode::Symbol
    reward_config::RewardConfig
    termination_config::TerminationConfig
    randomization_config::AerobrakingRandomizationConfig
    normalization_bounds::NormalizationBounds
    training::Bool
end

function default_aerobraking_config(; phase::AbstractString="Main",
                                    nominal::Bool=true,
                                    max_passes::Int=80,
                                    backend_mode::Symbol=:paper_surrogate,
                                    training::Bool=true,
                                    reward_config::RewardConfig=RewardConfig(),
                                    termination_config::Union{Nothing,TerminationConfig}=nothing,
                                    randomization_config::Union{Nothing,AerobrakingRandomizationConfig}=nothing,
                                    nominal_inclination_deg::Union{Nothing,Real}=nothing,
                                    nominal_raan_deg::Union{Nothing,Real}=nothing,
                                    nominal_argument_of_periapsis_deg::Union{Nothing,Real}=nothing)
    constants = mars_odyssey_phase_constants(phase)
    term = termination_config === nothing ? TerminationConfig(max_passes=max_passes) : termination_config
    rand_cfg = randomization_config === nothing ? AerobrakingRandomizationConfig(nominal=nominal) : randomization_config
    return AerobrakingScenarioConfig(
        constants.phase,
        constants.final_apoapsis_radius_m,
        constants.r_norm_m,
        constants.initial_apoapsis_radius_m,
        constants.nominal_periapsis_altitude_m,
        deg2rad(Float64(something(nominal_inclination_deg, constants.nominal_inclination_deg))),
        deg2rad(Float64(something(nominal_raan_deg, constants.nominal_raan_deg))),
        deg2rad(Float64(something(nominal_argument_of_periapsis_deg, constants.nominal_argument_of_periapsis_deg))),
        constants.nominal_epoch,
        constants.mars_radius_m,
        1.96045e-3,
        4.2828e13,
        8.748923102971180e-7,
        90e3,
        6.308278108290950e3,
        188.92,
        1.33,
        150.0,
        3.71,
        20e3,
        130e3,
        92e3,
        3450.0,
        0.15,
        backend_mode,
        reward_config,
        term,
        rand_cfg,
        paper_normalization_bounds(constants.phase),
        training,
    )
end

Base.@kwdef struct AerobrakingDecisionState
    pass_index::Int = 0
    apoapsis_radius_m::Float64
    periapsis_altitude_m::Float64
    inclination_rad::Float64
    raan_rad::Float64
    argument_of_periapsis_rad::Float64
    epoch::DateTime
    mission_elapsed_s::Float64 = 0.0
    total_delta_v_mps::Float64 = 0.0
    maneuver_count::Int = 0
    previous_density_kg_m3::Float64 = 0.0
    previous_heat_rate_w_cm2::Float64 = 0.0
    last_drag_passage_time_s::Float64 = 400.0
    drag_coefficient_scale::Float64 = 1.0
    lift_coefficient_scale::Float64 = 1.0
    marsgram_seed::Int = 1001
    marsgram_prediction_seed::Int = 2001
end

struct AerobrakingStepResult
    state::AerobrakingDecisionState
    action::AerobrakingAction
    raw_observation::PaperObservation
    normalized_observation::Vector{Float32}
    reward::Float64
    flags::TerminationFlags
    metrics::AerobrakingPassMetrics
    simulation_config
end

struct SpaceAGORAAerobrakingBackend <: AbstractRLBackend
    config::AerobrakingScenarioConfig
    adapter::SpaceAGORACoreAdapter
end

SpaceAGORAAerobrakingBackend(config::AerobrakingScenarioConfig=default_aerobraking_config()) =
    SpaceAGORAAerobrakingBackend(config, SpaceAGORACoreAdapter(config.backend_mode))

function observe_state(config::AerobrakingScenarioConfig, state::AerobrakingDecisionState)
    return PaperObservation(
        state.last_drag_passage_time_s,
        state.apoapsis_radius_m,
        state.periapsis_altitude_m,
        state.argument_of_periapsis_rad,
        state.raan_rad,
        state.inclination_rad,
        Float64(python_ordinal(state.epoch)),
        state.previous_density_kg_m3,
        state.previous_heat_rate_w_cm2,
    )
end

function reset_scenario(config::AerobrakingScenarioConfig, rng::AbstractRNG)
    rand_cfg = config.randomization_config
    if rand_cfg.nominal
        ra = config.initial_apoapsis_radius_m + uniform_jitter(rng, rand_cfg.apoapsis_jitter_m)
        hp = config.nominal_periapsis_altitude_m + uniform_jitter(rng, rand_cfg.periapsis_jitter_m)
        inc = config.nominal_inclination_rad + deg2rad(uniform_jitter(rng, rand_cfg.angle_jitter_deg))
        raan = config.nominal_raan_rad + deg2rad(uniform_jitter(rng, rand_cfg.angle_jitter_deg))
        omega = config.nominal_argument_of_periapsis_rad + deg2rad(uniform_jitter(rng, rand_cfg.angle_jitter_deg))
        epoch = config.nominal_epoch
    else
        ra = config.initial_apoapsis_radius_m + uniform_jitter(rng, rand_cfg.apoapsis_jitter_m)
        hp = config.nominal_periapsis_altitude_m + uniform_jitter(rng, rand_cfg.periapsis_jitter_m)
        inc = deg2rad(rand_cfg.nonnominal_inclination_low_deg +
                      rand(rng) * (rand_cfg.nonnominal_inclination_high_deg -
                                   rand_cfg.nonnominal_inclination_low_deg))
        raan = deg2rad(rand_cfg.nonnominal_raan_low_deg +
                       rand(rng) * (rand_cfg.nonnominal_raan_high_deg -
                                    rand_cfg.nonnominal_raan_low_deg))
        omega = deg2rad(rand_cfg.nonnominal_aop_low_deg +
                        rand(rng) * (rand_cfg.nonnominal_aop_high_deg -
                                     rand_cfg.nonnominal_aop_low_deg))
        epoch = config.nominal_epoch + Day(rand(rng, 0:27)) + Hour(rand(rng, 0:23))
    end
    drag_scale, lift_scale = paper_aerodynamic_coefficient_scales(config, rng)
    marsgram_seed, marsgram_prediction_seed = paper_marsgram_campaign_seeds(config, rng)
    return AerobrakingDecisionState(
        apoapsis_radius_m = ra,
        periapsis_altitude_m = hp,
        inclination_rad = inc,
        raan_rad = raan,
        argument_of_periapsis_rad = omega,
        epoch = epoch,
        drag_coefficient_scale = drag_scale,
        lift_coefficient_scale = lift_scale,
        marsgram_seed = marsgram_seed,
        marsgram_prediction_seed = marsgram_prediction_seed,
    )
end

reset_scenario(backend::SpaceAGORAAerobrakingBackend, rng::AbstractRNG) =
    reset_scenario(backend.config, rng)

function paper_aerodynamic_coefficient_scales(config::AerobrakingScenarioConfig, rng::AbstractRNG)
    rand_cfg = config.randomization_config
    if rand_cfg.aerodynamic_coefficient_dispersion && rand_cfg.aerodynamic_coefficient_span > 0
        span = rand_cfg.aerodynamic_coefficient_span
        return 1.0 + uniform_jitter(rng, span), 1.0 + uniform_jitter(rng, span)
    end
    return 1.0, 1.0
end

function paper_marsgram_campaign_seeds(config::AerobrakingScenarioConfig, rng::AbstractRNG)
    base = config.randomization_config.marsgram_seed_base
    actual = base + rand(rng, 0:typemax(Int32) - base - 1)
    predicted = base + rand(rng, 0:typemax(Int32) - base - 1)
    predicted == actual && (predicted += 1)
    return actual, predicted
end

function orbital_period_s(config, apoapsis_radius_m::Real, periapsis_altitude_m::Real)
    rp = config.mars_radius_m + Float64(periapsis_altitude_m)
    a = (Float64(apoapsis_radius_m) + rp) / 2
    return 2pi * sqrt(a^3 / config.mu_m3_s2)
end

function j2_angle_update(config, state, apoapsis_radius_m::Real, periapsis_altitude_m::Real, elapsed_s::Real)
    rp = config.mars_radius_m + Float64(periapsis_altitude_m)
    ra = Float64(apoapsis_radius_m)
    a = (ra + rp) / 2
    e = (ra - rp) / (ra + rp)
    denom = max((1 - e^2) * a^(7 / 2), eps(Float64))
    raan_rate = -(1.5 * sqrt(config.mu_m3_s2) * config.j2 * config.mars_radius_m^2 / denom) *
                cos(state.inclination_rad)
    arg_rate = abs(cos(state.inclination_rad)) < 1e-8 ? 0.0 :
               raan_rate * (((5 / 2) * sin(state.inclination_rad)^2 - 2) / cos(state.inclination_rad))
    return state.raan_rad + raan_rate * Float64(elapsed_s),
           state.argument_of_periapsis_rad + arg_rate * Float64(elapsed_s)
end

function periapsis_after_action_m(config, state, action::AerobrakingAction)
    return clamp(apply_apoapsis_maneuver(config, state, action), 50e3, 180e3)
end

function periapsis_velocity_mps(config, state, periapsis_altitude_m::Real)
    rp = config.mars_radius_m + Float64(periapsis_altitude_m)
    a_before_drag = (state.apoapsis_radius_m + rp) / 2
    return sqrt(config.mu_m3_s2 * (2 / rp - 1 / a_before_drag))
end

function apply_density_process_noise(config, density::Real, heat_rate::Real, rng::AbstractRNG)
    density_out = Float64(density)
    heat_out = Float64(heat_rate)
    if config.randomization_config.process_noise && config.randomization_config.process_noise_scale > 0
        density_multiplier = exp(config.randomization_config.process_noise_scale * randn(rng) -
                                 0.5 * config.randomization_config.process_noise_scale^2)
        density_out *= density_multiplier
        heat_out *= sqrt(density_multiplier)
    end
    return max(0.0, density_out), max(0.0, heat_out)
end

function density_heat_for_paper_surrogate(config, state, action::AerobrakingAction, rng::AbstractRNG)
    periapsis_after_maneuver = periapsis_after_action_m(config, state, action)
    periapsis_velocity = periapsis_velocity_mps(config, state, periapsis_after_maneuver)
    density = exponential_mars_density(config, periapsis_after_maneuver)
    heat_rate = paper_heat_rate_from_density_w_cm2(config, density, periapsis_velocity)
    density, heat_rate = apply_density_process_noise(config, density, heat_rate, rng)
    return periapsis_after_maneuver, density, heat_rate
end

function pass_transition_from_atmosphere(config::AerobrakingScenarioConfig,
                                         state::AerobrakingDecisionState,
                                         action::AerobrakingAction,
                                         rng::AbstractRNG,
                                         periapsis_after_maneuver::Real,
                                         density::Real,
                                         heat_rate::Real)
    hp = Float64(periapsis_after_maneuver)
    drag_time = drag_passage_duration_s(hp)

    heat_scale = clamp((max(heat_rate, 1e-6) / config.heat_nominal_w_cm2)^1.15, 0.05, 4.0)
    apoapsis_drop = (config.base_apoapsis_decay_m + config.heat_decay_gain_m * heat_scale) *
                    state.drag_coefficient_scale
    if config.randomization_config.process_noise && config.backend_mode != :spaceagora_marsgram
        apoapsis_drop *= 1 + config.randomization_config.process_noise_scale * randn(rng)
    end
    next_apoapsis = state.apoapsis_radius_m - max(0.0, apoapsis_drop)

    rp = config.mars_radius_m + hp
    period = orbital_period_s(config, max(next_apoapsis, rp + 1), hp)
    elapsed_s = max(period, drag_time)
    raan, omega = j2_angle_update(config, state, max(next_apoapsis, rp + 1), hp, elapsed_s)
    next_epoch = state.epoch + Millisecond(round(Int, elapsed_s * 1000))

    total_delta_v = state.total_delta_v_mps + action.magnitude_mps
    maneuver_count = state.maneuver_count + (action.magnitude_mps > 1e-9 ? 1 : 0)

    next_state = AerobrakingDecisionState(
        pass_index = state.pass_index + 1,
        apoapsis_radius_m = next_apoapsis,
        periapsis_altitude_m = hp,
        inclination_rad = state.inclination_rad,
        raan_rad = raan,
        argument_of_periapsis_rad = omega,
        epoch = next_epoch,
        mission_elapsed_s = state.mission_elapsed_s + elapsed_s,
        total_delta_v_mps = total_delta_v,
        maneuver_count = maneuver_count,
        previous_density_kg_m3 = Float64(density),
        previous_heat_rate_w_cm2 = Float64(heat_rate),
        last_drag_passage_time_s = drag_time,
        drag_coefficient_scale = state.drag_coefficient_scale,
        lift_coefficient_scale = state.lift_coefficient_scale,
        marsgram_seed = state.marsgram_seed,
        marsgram_prediction_seed = state.marsgram_prediction_seed,
    )
    return next_state
end

function paper_surrogate_pass(config::AerobrakingScenarioConfig, state::AerobrakingDecisionState,
                              action::AerobrakingAction, rng::AbstractRNG)
    periapsis_after_maneuver, density, heat_rate =
        density_heat_for_paper_surrogate(config, state, action, rng)
    return pass_transition_from_atmosphere(config, state, action, rng,
                                           periapsis_after_maneuver, density, heat_rate)
end

const _SPACEAGORA_RL_BACKEND_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const _SPACEAGORA_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const _SPACEAGORA_MODULE = Ref{Any}(nothing)
const _SPACEAGORA_LOAD_LOCK = ReentrantLock()
const _SPACEAGORA_GRAMSUITE_MODULE = Ref{Any}(nothing)
const _SPACEAGORA_GRAMSUITE_LOAD_LOCK = ReentrantLock()
const _SPACEAGORA_MARSGRAM_MODEL_CACHE = Dict{Tuple{DateTime, Float64, Int}, Any}()
const _SPACEAGORA_MARSGRAM_MODEL_LOCK = ReentrantLock()
const _SPACEAGORA_MARSGRAM_MODEL_CACHE_LIMIT = 64
const _SPACEAGORA_PHYSICS_GRAM_SURROGATE_FALLBACK_INSTALLED = Ref(false)
const _SPACEAGORA_PHYSICS_GRAM_SURROGATE_FALLBACK_LOCK = ReentrantLock()
const _MARS_ROTATION_RATE_RAD_S = 2pi / 88775.244

Base.@kwdef struct SpaceAGORAPhysicsGRAMSurrogateBase
    planet_name::String = "mars"
    temperature_k::Float64 = 200.0
end

function _spaceagora_repo_root()
    return _SPACEAGORA_REPO_ROOT
end

function _load_spaceagora!()
    loaded = _SPACEAGORA_MODULE[]
    loaded !== nothing && return loaded

    return lock(_SPACEAGORA_LOAD_LOCK) do
        loaded = _SPACEAGORA_MODULE[]
        loaded !== nothing && return loaded

        repo_root = _spaceagora_repo_root()
        isfile(joinpath(repo_root, "Project.toml")) || throw(ArgumentError(
            "SpaceAGORA physics backend requires the sibling SpaceAGORA project at $repo_root"
        ))
        repo_root in LOAD_PATH || pushfirst!(LOAD_PATH, repo_root)
        pkgid = Base.PkgId(Base.UUID("afbfb69f-5c0b-4832-b760-43725dff8540"), "SpaceAGORA")
        try
            loaded = Base.require(pkgid)
        catch err
            throw(ErrorException(
                "SpaceAGORA physics backend could not load the SpaceAGORA package from $repo_root. " *
                "Run it from the SpaceAGORA.jl project, or instantiate that project first. " *
                "Original load failure: $(sprint(showerror, err))"
            ))
        end
        _SPACEAGORA_MODULE[] = loaded
        return loaded
    end
end

function _spaceagora_gramsuite_path()
    return joinpath(_spaceagora_repo_root(), "data", "GRAMSuite.jl")
end

function _spaceagora_marsgram_root()
    configured = strip(get(ENV, "SPACEAGORA_GRAM_ROOT", get(ENV, "GRAM_ROOT", "")))
    return isempty(configured) ?
           joinpath(_spaceagora_gramsuite_path(), "GRAM Suite 2.0") :
           normpath(expanduser(configured))
end

function _gramsuite_missing_dependency_error(err)::Bool
    msg = sprint(showerror, err)
    return occursin("is required but does not seem to be installed", msg) ||
           occursin("Run `Pkg.instantiate()` to install all recorded dependencies", msg)
end

function _instantiate_vendored_gramsuite!(vendored_gramsuite::String)
    Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    previous_project = something(Base.active_project(), "")
    try
        Pkg.activate(vendored_gramsuite; io=devnull)
        Pkg.instantiate(; io=devnull)
    finally
        if !isempty(previous_project)
            Pkg.activate(dirname(previous_project); io=devnull)
        end
    end
    return nothing
end

function _load_spaceagora_gramsuite!()
    loaded = _SPACEAGORA_GRAMSUITE_MODULE[]
    loaded !== nothing && return loaded

    return lock(_SPACEAGORA_GRAMSUITE_LOAD_LOCK) do
        loaded = _SPACEAGORA_GRAMSUITE_MODULE[]
        loaded !== nothing && return loaded

        vendored_gramsuite = _spaceagora_gramsuite_path()
        isdir(vendored_gramsuite) || throw(ArgumentError(
            "SpaceAGORA MarsGRAM backend requires vendored GRAMSuite.jl at $vendored_gramsuite"
        ))
        vendored_gramsuite in LOAD_PATH || pushfirst!(LOAD_PATH, vendored_gramsuite)

        pkgid = Base.PkgId(Base.UUID("b50455af-6a46-4eae-bf92-8039261dd674"), "GRAMSuite")
        try
            loaded = Base.require(pkgid)
        catch err
            if _gramsuite_missing_dependency_error(err)
                _instantiate_vendored_gramsuite!(vendored_gramsuite)
                loaded = Base.require(pkgid)
            else
                throw(ErrorException(
                    "SpaceAGORA MarsGRAM backend could not load GRAMSuite.jl. " *
                    "Original load failure: $(sprint(showerror, err))"
                ))
            end
        end
        _SPACEAGORA_GRAMSUITE_MODULE[] = loaded
        return loaded
    end
end

function _marsgram_initial_time(gramsuite, epoch::DateTime)
    initial_time_ctor = Base.invokelatest(getproperty, gramsuite, :InitialTime)
    return Base.invokelatest(
        initial_time_ctor;
        year = Int32(Dates.year(epoch)),
        month = Int16(Dates.month(epoch)),
        day = Int16(Dates.day(epoch)),
        hour = Int16(Dates.hour(epoch)),
        minute = Int16(Dates.minute(epoch)),
        second = Float32(Dates.second(epoch) + Dates.millisecond(epoch) / 1000),
    )
end

function _marsgram_live_model(config::AerobrakingScenarioConfig, seed::Integer)
    gramsuite = _load_spaceagora_gramsuite!()
    perturbation_scale = config.randomization_config.process_noise ?
                         config.randomization_config.marsgram_perturbation_scale : 0.0
    key = (config.nominal_epoch, Float64(perturbation_scale), Int(seed))
    lock(_SPACEAGORA_MARSGRAM_MODEL_LOCK) do
        if haskey(_SPACEAGORA_MARSGRAM_MODEL_CACHE, key)
            return _SPACEAGORA_MARSGRAM_MODEL_CACHE[key]
        end
        if length(_SPACEAGORA_MARSGRAM_MODEL_CACHE) >= _SPACEAGORA_MARSGRAM_MODEL_CACHE_LIMIT
            empty!(_SPACEAGORA_MARSGRAM_MODEL_CACHE)
        end
        gram_root = _spaceagora_marsgram_root()
        model_ctor = Base.invokelatest(getproperty, gramsuite, :GRAMAtmosphereModel)
        model = try
            Base.invokelatest(
                model_ctor;
                gram_root_directory = gram_root,
                gram_data_directory = gram_root,
                spice_directory = joinpath(gram_root, "SPICE"),
                planet_name = "mars",
                seed = Int(seed),
                initial_time = _marsgram_initial_time(gramsuite, config.nominal_epoch),
                gram_perturbation_scales = perturbation_scale > 0 ? perturbation_scale : nothing,
                mars_dust_storm = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            )
        catch err
            throw(ErrorException(
                "SpaceAGORA live MarsGRAM backend is not ready. " *
                "Expected a host-native GRAM shared library under $gram_root/Build/lib. " *
                "Run `julia --project=SpaceAGORA.jl SpaceAGORA.jl/scripts/ensure_gram_native.jl` " *
                "after installing build tools such as `make`. Original error: $(sprint(showerror, err))"
            ))
        end
        _SPACEAGORA_MARSGRAM_MODEL_CACHE[key] = model
        return model
    end
end

function _periapsis_lat_lon_rad(config::AerobrakingScenarioConfig, state::AerobrakingDecisionState)
    cosΩ, sinΩ = cos(state.raan_rad), sin(state.raan_rad)
    cosω, sinω = cos(state.argument_of_periapsis_rad), sin(state.argument_of_periapsis_rad)
    cosi, sini = cos(state.inclination_rad), sin(state.inclination_rad)
    x = cosΩ * cosω - sinΩ * sinω * cosi
    y = sinΩ * cosω + cosΩ * sinω * cosi
    z = sinω * sini
    lat = asin(clamp(z, -1.0, 1.0))
    epoch_elapsed_s = Dates.value(state.epoch - config.nominal_epoch) / 1000
    lon = mod(atan(y, x) - _MARS_ROTATION_RATE_RAD_S * epoch_elapsed_s, 2pi)
    return lat, lon
end

function _marsgram_density_state(config::AerobrakingScenarioConfig,
                                 state::AerobrakingDecisionState,
                                 altitude_m::Real,
                                 mode::Symbol;
                                 prediction::Bool=false)
    gramsuite = _load_spaceagora_gramsuite!()
    lat, lon = _periapsis_lat_lon_rad(config, state)
    elapsed_s = Dates.value(state.epoch - config.nominal_epoch) / 1000
    h = Float64(altitude_m)
    if mode == :spaceagora_marsgram
        seed = prediction ? state.marsgram_prediction_seed : state.marsgram_seed
        model = _marsgram_live_model(config, seed)
        density_state = Base.invokelatest(getproperty, gramsuite, :density_state)
        return Base.invokelatest(
            density_state,
            model,
            h,
            lat,
            lon,
            elapsed_s,
            true;
            vacuum_temperature = config.temperature_k,
        )
    end
    throw(ArgumentError("unsupported MarsGRAM backend mode: $mode"))
end

Base.@kwdef mutable struct SpaceAGORAPhysicsPassStats
    max_density_kg_m3::Float64 = 0.0
    max_heat_rate_w_cm2::Float64 = 0.0
    min_altitude_m::Float64 = Inf
    in_atmosphere::Bool = false
    entry_time_s::Float64 = NaN
    exit_time_s::Float64 = NaN
end

function _spaceagora_physics_initial_time(spaceagora, epoch::DateTime)
    initial_time_ctor = Base.invokelatest(getproperty, getproperty(spaceagora, :SimulationModel), :InitialTime)
    return Base.invokelatest(
        initial_time_ctor;
        year = Int32(Dates.year(epoch)),
        month = Int16(Dates.month(epoch)),
        day = Int16(Dates.day(epoch)),
        hour = Int16(Dates.hour(epoch)),
        minute = Int16(Dates.minute(epoch)),
        second = Float32(Dates.second(epoch) + Dates.millisecond(epoch) / 1000),
    )
end

function _ensure_spaceagora_physics_gram_surrogate_fallback!(spaceagora)
    _SPACEAGORA_PHYSICS_GRAM_SURROGATE_FALLBACK_INSTALLED[] && return nothing
    return lock(_SPACEAGORA_PHYSICS_GRAM_SURROGATE_FALLBACK_LOCK) do
        _SPACEAGORA_PHYSICS_GRAM_SURROGATE_FALLBACK_INSTALLED[] && return nothing
        em = getproperty(getproperty(spaceagora, :SimulationModel), :EnvironmentModels)
        fallback_type = SpaceAGORAPhysicsGRAMSurrogateBase
        Core.eval(em, quote
            function _gram_point_density(model::$fallback_type,
                                         h::Float64,
                                         lat::Float64,
                                         lon::Float64,
                                         el_time::Float64,
                                         wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
                return 0.0, model.temperature_k, SVector{3, Float64}(0.0, 0.0, 0.0)
            end
        end)
        _SPACEAGORA_PHYSICS_GRAM_SURROGATE_FALLBACK_INSTALLED[] = true
        return nothing
    end
end

function _spaceagora_physics_density_model(spaceagora,
                                           config::AerobrakingScenarioConfig,
                                           state::AerobrakingDecisionState,
                                           initial_time;
                                           prediction::Bool=false)
    gramsuite = _load_spaceagora_gramsuite!()
    _ensure_spaceagora_physics_gram_surrogate_fallback!(spaceagora)
    sm = getproperty(spaceagora, :SimulationModel)
    gram_root = _spaceagora_marsgram_root()
    surrogate_file = Base.invokelatest(
        Base.invokelatest(getproperty, gramsuite, :gram_default_surrogate_file),
        "mars";
        gram_root = gram_root,
    )
    isfile(surrogate_file) || throw(ArgumentError("GRAM surrogate payload not found: $surrogate_file"))
    model_ctor = Base.invokelatest(getproperty, sm, :GRAMAtmosphereModelSurrogate)
    base_model = SpaceAGORAPhysicsGRAMSurrogateBase(
        planet_name = "mars",
        temperature_k = config.temperature_k,
    )
    try
        return Base.invokelatest(
            model_ctor,
            base_model,
            surrogate_file,
            nothing,
        )
    catch err
        throw(ErrorException(
            "SpaceAGORA physics backend could not initialize the GRAM surrogate density model. " *
            "Expected a prebuilt Mars surrogate payload at $surrogate_file. " *
            "If the file is missing, ensure GRAMSuite LFS assets are present or rebuild the offline surrogate. " *
            "Original error: $(sprint(showerror, err))"
        ))
    end
end

function _ensure_spaceagora_scaled_aero_effector!(spaceagora)
    sm = getproperty(spaceagora, :SimulationModel)
    if !isdefined(sm, :RLScaledAerodynamicCoefficientfM)
        Core.eval(sm, quote
            struct RLScaledAerodynamicCoefficientfM <: AbstractForceTorqueModel
                model::AerodynamicCoefficientfM
                drag_scale::Float64
                lift_scale::Float64

                function RLScaledAerodynamicCoefficientfM(model::AerodynamicCoefficientfM,
                                                          drag_scale::Real,
                                                          lift_scale::Real)
                    drag = Float64(drag_scale)
                    lift = Float64(lift_scale)
                    drag > 0.0 || throw(ArgumentError("drag_scale must be > 0.0, got $drag"))
                    lift > 0.0 || throw(ArgumentError("lift_scale must be > 0.0, got $lift"))
                    return new(model, drag, lift)
                end
            end

            @inline EffectorSampling.environment_requirements(model::RLScaledAerodynamicCoefficientfM) =
                EffectorSampling.environment_requirements(model.model)
            @inline EffectorSampling.solver_partition(model::RLScaledAerodynamicCoefficientfM) =
                EffectorSampling.solver_partition(model.model)

            @inline function _rl_scaled_aero_cache_value(cache, idx::Int)
                if idx <= length(cache) && isassigned(cache, idx)
                    return cache[idx]
                end
                return SVector{3, Float64}(0.0, 0.0, 0.0)
            end

            @inline function _rl_scaled_aero_store_cache!(cache, idx::Int, value::SVector{3, Float64})
                idx <= length(cache) || resize!(cache, idx)
                cache[idx] = value
                return nothing
            end

            function DynamicEffectors.calcForceTorque(model::RLScaledAerodynamicCoefficientfM,
                                                      x::AbstractVector{Float64},
                                                      param::ODEParams,
                                                      i::Int64)
                force, torque = DynamicEffectors.calcForceTorque(model.model, x, param, i)
                drag = _rl_scaled_aero_cache_value(param.save_cache.drag_cache, i)
                lift = _rl_scaled_aero_cache_value(param.save_cache.lift_cache, i)
                cross_force = _rl_scaled_aero_cache_value(param.save_cache.cross_cache, i)
                if norm(drag) + norm(lift) + norm(cross_force) <= eps(Float64)
                    return model.drag_scale * force, model.drag_scale * torque
                end

                scaled_drag = model.drag_scale * drag
                scaled_lift = model.lift_scale * lift
                scaled_cross = model.drag_scale * cross_force
                _rl_scaled_aero_store_cache!(param.save_cache.drag_cache, i, scaled_drag)
                _rl_scaled_aero_store_cache!(param.save_cache.lift_cache, i, scaled_lift)
                _rl_scaled_aero_store_cache!(param.save_cache.cross_cache, i, scaled_cross)
                return scaled_drag + scaled_lift + scaled_cross, model.drag_scale * torque
            end
        end)
    end

    se = getproperty(spaceagora, :SimulationEngine)
    scaled_type = Base.invokelatest(getproperty, sm, :RLScaledAerodynamicCoefficientfM)
    Core.eval(se, :(@inline _dynamic_effector_threadsafe(::$scaled_type)::Bool = true))
    return nothing
end

function _spaceagora_scaled_aero_effector(spaceagora, state::AerobrakingDecisionState)
    sm = getproperty(spaceagora, :SimulationModel)
    _ensure_spaceagora_scaled_aero_effector!(spaceagora)
    aero_base = Base.invokelatest(Base.invokelatest(getproperty, sm, :AerodynamicCoefficientfM))
    scaled_ctor = Base.invokelatest(getproperty, sm, :RLScaledAerodynamicCoefficientfM)
    return Base.invokelatest(
        scaled_ctor,
        aero_base,
        state.drag_coefficient_scale,
        state.lift_coefficient_scale,
    )
end

function _spaceagora_physics_spacecraft(spaceagora, ic)
    make_spacecraft = Base.invokelatest(
        getproperty,
        Base.invokelatest(getproperty, spaceagora, :TelemetryVerification),
        :make_three_body_spacecraft,
    )
    return Base.invokelatest(
        make_spacecraft;
        bus_dims = (2.2, 2.6, 1.7),
        panel_dims = (0.01, 5.5 / 1.35, 2.6),
        bus_mass = 391.0,
        panel_mass_each = 10.0,
        panel_offset_y = 2.6 / 2.0 + 5.5 / 4.0,
        ic = ic,
        reflection_coefficient = 0.9,
        prop_mass = 50.0,
        id = 100,
    )
end

function _spaceagora_physics_dynamic_effectors(spaceagora, state::AerobrakingDecisionState,
                                               _planet, _spacecraft)
    sm = getproperty(spaceagora, :SimulationModel)
    gravity = Base.invokelatest(Base.invokelatest(getproperty, sm, :InverseSquaredJ2GravityModel))
    aero = _spaceagora_scaled_aero_effector(spaceagora, state)
    return (gravity, aero)
end

function _spaceagora_physics_simulation_configuration(config::AerobrakingScenarioConfig,
                                                      state::AerobrakingDecisionState,
                                                      action::AerobrakingAction;
                                                      prediction::Bool=false,
                                                      campaign_max_passes::Union{Nothing,Int}=nothing)
    spaceagora = _load_spaceagora!()
    sm = getproperty(spaceagora, :SimulationModel)
    gram_root = _spaceagora_marsgram_root()
    spice_path = joinpath(gram_root, "SPICE")
    initial_time = _spaceagora_physics_initial_time(spaceagora, state.epoch)
    planet = Base.invokelatest(Base.invokelatest(getproperty, sm, :Mars))
    density_model = _spaceagora_physics_density_model(
        spaceagora,
        config,
        state,
        initial_time;
        prediction=prediction,
    )

    periapsis_after_maneuver = periapsis_after_action_m(config, state, action)
    planet_radius = Float64(getproperty(planet, :Rp_e))
    rp_radius = planet_radius + periapsis_after_maneuver
    ra_radius = max(state.apoapsis_radius_m, rp_radius + 1.0)
    ic = Base.invokelatest(
        Base.invokelatest(getproperty, sm, :InitialCondition);
        ra = ra_radius,
        rp = rp_radius,
        i = rad2deg(state.inclination_rad),
        ω = rad2deg(state.argument_of_periapsis_rad),
        Ω = rad2deg(state.raan_rad),
        # Start just after apoapsis so the orbit-end callback records the next apoapsis.
        ν = 180.0 + 1.0e-3,
    )
    spacecraft = _spaceagora_physics_spacecraft(spaceagora, ic)
    dynamic_effectors = _spaceagora_physics_dynamic_effectors(spaceagora, state, planet, spacecraft)
    orbital_period = orbital_period_s(config, ra_radius, periapsis_after_maneuver)
    pass_mission_time = max(1.5 * orbital_period + 600.0, orbital_period + 1800.0)
    mission_time = if campaign_max_passes === nothing
        pass_mission_time
    else
        max(1, campaign_max_passes) * pass_mission_time
    end

    file_paths = Base.invokelatest(
        Base.invokelatest(getproperty, sm, :FilePaths);
        GRAM = gram_root,
        SPICE = spice_path,
        gravity_harmonics = joinpath(_spaceagora_repo_root(), "data", "Gravity_harmonics_data"),
    )
    simulation_settings = Base.invokelatest(
        Base.invokelatest(getproperty, sm, :SimulationSettings);
        results = false,
        verbose = false,
        results_directory = joinpath(tempdir(), "spaceagora_rl_physics"),
        generate_plots = false,
        normalize = false,
        save_csv = false,
    )
    mission_configuration = Base.invokelatest(
        Base.invokelatest(getproperty, sm, :MissionConfiguration);
        mission_type = campaign_max_passes === nothing ?
                       Base.invokelatest(getproperty, sm, :MissionOrbits) :
                       Base.invokelatest(getproperty, sm, :MissionTime),
        keplerian = false,
        number_of_orbits = 1,
        mission_time = mission_time,
        orientation_sim = false,
        num_steps_to_save = 1000,
        data_rate = 10.0,
    )
    environment_model = Base.invokelatest(
        Base.invokelatest(getproperty, sm, :EnvironmentModel);
        planet = planet,
        EI = config.termination_config.out_of_passage_periapsis_altitude_m / 1000,
        density_model = density_model,
        ephemerides_model = Base.invokelatest(Base.invokelatest(getproperty, sm, :SimpleEphemeridesModel)),
        topography = false,
        wind = true,
        thermal_model = Base.invokelatest(
            Base.invokelatest(getproperty, sm, :MaxwellianHeat);
            thermal_accomodation_factor = 1.0,
            planet = planet,
        ),
    )
    args = Base.invokelatest(
        Base.invokelatest(getproperty, sm, :SimulationConfiguration);
        file_paths = file_paths,
        simulation_settings = simulation_settings,
        mission_configuration = mission_configuration,
        environment_model = environment_model,
        dynamics_model = Base.invokelatest(
            Base.invokelatest(getproperty, sm, :DynamicsModel),
            [spacecraft],
            dynamic_effectors,
        ),
        guidance_model = Base.invokelatest(
            Base.invokelatest(getproperty, sm, :GuidanceModel);
            guidance_effectors = (),
            guidance_rates = Float64[],
        ),
        navigation_model = Base.invokelatest(
            Base.invokelatest(getproperty, sm, :NavigationModel);
            navigation_effectors = (),
            navigation_rates = Float64[],
        ),
        control_model = Base.invokelatest(
            Base.invokelatest(getproperty, sm, :ControlModel);
            control_effectors = (),
            control_rates = Float64[],
        ),
        initial_time = initial_time,
        integration_tolerances = Base.invokelatest(
            Base.invokelatest(getproperty, sm, :IntegrationTolerances);
            reltol_orbit = 1e-8,
            abstol_orbit = 1e-8,
            dt_max_orbit = 30.0,
            reltol_atmosphere = 1e-8,
            abstol_atmosphere = 1e-8,
            dt_max_atmosphere = 5.0,
        ),
        solver_config = Base.invokelatest(
            Base.invokelatest(getproperty, sm, :SolverConfig);
            solver_mode = :split_imex,
            split_imex_solver = :kencarp4,
        ),
    )
    return args, periapsis_after_maneuver
end

@inline function _spaceagora_solution_satellite_state(u)
    return hasproperty(u, :sc) ? u.sc[1] : u
end

function _record_spaceagora_physics_sample!(spaceagora,
                                            stats::SpaceAGORAPhysicsPassStats,
                                            u,
                                            p,
                                            t::Float64)
    engine = getproperty(spaceagora, :SimulationEngine)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    sc = _spaceagora_solution_satellite_state(u)
    planet_frame = Base.invokelatest(
        Base.invokelatest(getproperty, engine, :sample_planet_frame),
        sc,
        p,
        1,
        t,
    )
    altitude = Float64(getproperty(planet_frame, :alt_m))
    stats.min_altitude_m = min(stats.min_altitude_m, altitude)
    atmosphere_limit_m = Float64(getproperty(getproperty(p, :args), :environment_model).EI) * 1000
    altitude <= atmosphere_limit_m || return nothing

    if !stats.in_atmosphere
        stats.in_atmosphere = true
        stats.entry_time_s = t
    end
    stats.exit_time_s = t

    atmosphere = Base.invokelatest(
        Base.invokelatest(getproperty, engine, :sample_atmosphere),
        sc,
        p,
        1,
        t;
        write_buffers = false,
    )
    stats.max_density_kg_m3 = max(
        stats.max_density_kg_m3,
        max(0.0, Float64(getproperty(atmosphere, :rho_kg_m3))),
    )

    heat_rates = try
        getproperty(getproperty(p, :shared_buffers), :heat_rates)[1]
    catch
        Float64[]
    end
    if isempty(heat_rates) || maximum(heat_rates) <= 0.0
        heat_func = Base.invokelatest(getproperty, callbacks, :_compute_stage_heat_rates!)
        heat_rates = Base.invokelatest(
            heat_func,
            p,
            sc,
            1,
            t;
            use_buffered_density = false,
        )
    end
    for heat_rate in heat_rates
        if isfinite(heat_rate)
            stats.max_heat_rate_w_cm2 = max(stats.max_heat_rate_w_cm2, max(0.0, Float64(heat_rate)))
        end
    end
    return nothing
end

function _record_spaceagora_physics_sample!(spaceagora,
                                            stats::SpaceAGORAPhysicsPassStats,
                                            integrator)
    return _record_spaceagora_physics_sample!(
        spaceagora,
        stats,
        getproperty(integrator, :u),
        getproperty(integrator, :p),
        Float64(getproperty(integrator, :t)),
    )
end

function _collect_spaceagora_physics_solution_stats!(spaceagora,
                                                     stats::SpaceAGORAPhysicsPassStats,
                                                     sol)
    p = getproperty(getproperty(sol, :prob), :p)
    sample_count = max(length(sol.u), 1000)
    t_start = Float64(first(sol.t))
    t_stop = Float64(last(sol.t))
    for t in range(t_start, t_stop; length=sample_count)
        u = Base.invokelatest(sol, t)
        _record_spaceagora_physics_sample!(spaceagora, stats, u, p, Float64(t))
    end
    return stats
end

function _spaceagora_physics_stats_callback(spaceagora, stats::SpaceAGORAPhysicsPassStats)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    discrete_callback = Base.invokelatest(getproperty, callbacks, :DiscreteCallback)
    condition(u, t, integrator) = true
    affect!(integrator) = _record_spaceagora_physics_sample!(spaceagora, stats, integrator)
    initialize = (cb, u, t, integrator) -> affect!(integrator)
    return Base.invokelatest(discrete_callback, condition, affect!; initialize=initialize)
end

function _spaceagora_physics_next_state_from_u(spaceagora,
                                               config::AerobrakingScenarioConfig,
                                               state::AerobrakingDecisionState,
                                               action::AerobrakingAction,
                                               args,
                                               stats::SpaceAGORAPhysicsPassStats,
                                               final_u,
                                               elapsed_s::Real,
                                               periapsis_after_maneuver::Real)
    engine = getproperty(spaceagora, :SimulationEngine)
    planet = getproperty(getproperty(args, :environment_model), :planet)
    pos = Base.invokelatest(Base.invokelatest(getproperty, engine, :_state_position_ii), final_u, 1)
    vel = Base.invokelatest(Base.invokelatest(getproperty, engine, :_state_velocity_ii), final_u, 1)
    oe = getfield(
        Base.invokelatest(Base.invokelatest(getproperty, engine, :rvtoorbitalelement), pos, vel, planet),
        :data,
    )
    a = Float64(oe[1])
    e = Float64(oe[2])
    planet_radius = Float64(getproperty(planet, :Rp_e))
    periapsis_radius = (isfinite(a) && isfinite(e)) ? a * (1 - e) :
                       planet_radius + min(Float64(periapsis_after_maneuver), stats.min_altitude_m)
    apoapsis_radius = (isfinite(a) && isfinite(e)) ? a * (1 + e) :
                      max(state.apoapsis_radius_m, periapsis_radius + 1.0)
    next_apoapsis = max(apoapsis_radius, periapsis_radius + 1.0)
    next_periapsis_altitude = periapsis_radius - planet_radius
    drag_time = if isfinite(stats.entry_time_s) && isfinite(stats.exit_time_s)
        max(0.0, stats.exit_time_s - stats.entry_time_s)
    else
        drag_passage_duration_s(periapsis_after_maneuver)
    end

    total_delta_v = state.total_delta_v_mps + action.magnitude_mps
    maneuver_count = state.maneuver_count + (action.magnitude_mps > 1e-9 ? 1 : 0)
    elapsed = Float64(elapsed_s)
    next_epoch = state.epoch + Millisecond(round(Int, elapsed * 1000))
    return AerobrakingDecisionState(
        pass_index = state.pass_index + 1,
        apoapsis_radius_m = next_apoapsis,
        periapsis_altitude_m = next_periapsis_altitude,
        inclination_rad = isfinite(oe[3]) ? Float64(oe[3]) : state.inclination_rad,
        raan_rad = isfinite(oe[4]) ? Float64(oe[4]) : state.raan_rad,
        argument_of_periapsis_rad = isfinite(oe[5]) ? Float64(oe[5]) : state.argument_of_periapsis_rad,
        epoch = next_epoch,
        mission_elapsed_s = state.mission_elapsed_s + elapsed,
        total_delta_v_mps = total_delta_v,
        maneuver_count = maneuver_count,
        previous_density_kg_m3 = stats.max_density_kg_m3,
        previous_heat_rate_w_cm2 = stats.max_heat_rate_w_cm2,
        last_drag_passage_time_s = drag_time,
        drag_coefficient_scale = state.drag_coefficient_scale,
        lift_coefficient_scale = state.lift_coefficient_scale,
        marsgram_seed = state.marsgram_seed,
        marsgram_prediction_seed = state.marsgram_prediction_seed,
    )
end

function _spaceagora_physics_next_state(spaceagora,
                                        config::AerobrakingScenarioConfig,
                                        state::AerobrakingDecisionState,
                                        action::AerobrakingAction,
                                        args,
                                        stats::SpaceAGORAPhysicsPassStats,
                                        sol,
                                        periapsis_after_maneuver::Real)
    return _spaceagora_physics_next_state_from_u(
        spaceagora,
        config,
        state,
        action,
        args,
        stats,
        sol.u[end],
        Float64(sol.t[end]),
        periapsis_after_maneuver,
    )
end

function _spaceagora_physics_propagate(config::AerobrakingScenarioConfig,
                                       state::AerobrakingDecisionState,
                                       action::AerobrakingAction;
                                       prediction::Bool=false)
    spaceagora = _load_spaceagora!()
    args, periapsis_after_maneuver =
        _spaceagora_physics_simulation_configuration(config, state, action; prediction=prediction)
    stats = SpaceAGORAPhysicsPassStats()
    stats_callback = _spaceagora_physics_stats_callback(spaceagora, stats)
    run_simulation_fn = Base.invokelatest(getproperty, spaceagora, :run_simulation)
    sol = try
        Base.invokelatest(
            run_simulation_fn,
            args;
            return_solution = true,
            extra_callbacks = (stats_callback,),
        )
    catch err
        throw(ErrorException(
            "SpaceAGORA physics backend failed while propagating one aerobraking pass. " *
            "No simplified transition or surrogate fallback was used. " *
            "Original error: $(sprint(showerror, err))"
        ))
    end
    _collect_spaceagora_physics_solution_stats!(spaceagora, stats, sol)
    next_state = _spaceagora_physics_next_state(
        spaceagora,
        config,
        state,
        action,
        args,
        stats,
        sol,
        periapsis_after_maneuver,
    )
    return next_state, stats
end

function spaceagora_physics_pass(config::AerobrakingScenarioConfig,
                                 state::AerobrakingDecisionState,
                                 action::AerobrakingAction,
                                 rng::AbstractRNG)
    next_state, _ = _spaceagora_physics_propagate(config, state, action; prediction=false)
    return next_state
end

function spaceagora_physics_predicted_heat_rate(config::AerobrakingScenarioConfig,
                                                state::AerobrakingDecisionState,
                                                action::AerobrakingAction)
    _, stats = _spaceagora_physics_propagate(config, state, action; prediction=true)
    return stats.max_heat_rate_w_cm2
end

function density_heat_for_spaceagora_marsgram(config::AerobrakingScenarioConfig,
                                              state::AerobrakingDecisionState,
                                              action::AerobrakingAction,
                                              rng::AbstractRNG,
                                              mode::Symbol;
                                              apply_process_noise::Bool=true)
    periapsis_after_maneuver = periapsis_after_action_m(config, state, action)
    periapsis_velocity = periapsis_velocity_mps(config, state, periapsis_after_maneuver)
    rho, _, _ = _marsgram_density_state(config, state, periapsis_after_maneuver, mode)
    density = max(0.0, Float64(rho))
    heat_rate = paper_heat_rate_from_density_w_cm2(config, density, periapsis_velocity)
    if apply_process_noise && mode != :spaceagora_marsgram
        density, heat_rate = apply_density_process_noise(config, density, heat_rate, rng)
    end
    return periapsis_after_maneuver, density, heat_rate
end

function deterministic_predicted_heat_rate(config::AerobrakingScenarioConfig,
                                           state::AerobrakingDecisionState,
                                           action::AerobrakingAction)
    periapsis_after_maneuver = periapsis_after_action_m(config, state, action)
    periapsis_velocity = periapsis_velocity_mps(config, state, periapsis_after_maneuver)
    if config.backend_mode == :spaceagora_physics || config.backend_mode == :spaceagora_full_physics
        return spaceagora_physics_predicted_heat_rate(config, state, action)
    end
    if config.backend_mode == :spaceagora_marsgram
        rho, _, _ = _marsgram_density_state(config, state, periapsis_after_maneuver,
                                            config.backend_mode; prediction=true)
        return paper_heat_rate_from_density_w_cm2(config, max(0.0, Float64(rho)), periapsis_velocity)
    end
    return paper_heat_rate_w_cm2(config, periapsis_after_maneuver, periapsis_velocity)
end

function spaceagora_marsgram_pass(config::AerobrakingScenarioConfig, state::AerobrakingDecisionState,
                                  action::AerobrakingAction, rng::AbstractRNG, mode::Symbol)
    periapsis_after_maneuver, density, heat_rate =
        density_heat_for_spaceagora_marsgram(config, state, action, rng, mode)
    return pass_transition_from_atmosphere(config, state, action, rng,
                                           periapsis_after_maneuver, density, heat_rate)
end

function step_scenario(config::AerobrakingScenarioConfig, state::AerobrakingDecisionState,
                       action::AerobrakingAction, rng::AbstractRNG)
    sim_config = build_simulation_configuration(config, state, action)
    adapter = SpaceAGORACoreAdapter(config.backend_mode)
    next_state = propagate_pass(adapter, config, state, action, rng)
    obs = observe_state(config, next_state)
    flags = classify_termination(obs, config; training=config.training, pass_count=next_state.pass_index)
    reward = paper_reward(obs, config, action, flags, config.reward_config)
    normalized = normalize_observation(obs, config.normalization_bounds)
    metrics = pass_metrics_from_state(next_state)
    return AerobrakingStepResult(next_state, action, obs, normalized, reward, flags, metrics, sim_config)
end

step_scenario(config::AerobrakingScenarioConfig, state::AerobrakingDecisionState,
              action_index::Integer, rng::AbstractRNG) =
    step_scenario(config, state, action_from_index(action_index), rng)

step_scenario(backend::SpaceAGORAAerobrakingBackend, state::AerobrakingDecisionState,
              action, rng::AbstractRNG) =
    step_scenario(backend.config, state, action isa AerobrakingAction ? action : action_from_index(action), rng)
