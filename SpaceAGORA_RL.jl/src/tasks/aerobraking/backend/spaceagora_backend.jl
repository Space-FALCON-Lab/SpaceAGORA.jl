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
                                    randomization_config::Union{Nothing,AerobrakingRandomizationConfig}=nothing)
    constants = mars_odyssey_phase_constants(phase)
    term = termination_config === nothing ? TerminationConfig(max_passes=max_passes) : termination_config
    rand_cfg = randomization_config === nothing ? AerobrakingRandomizationConfig(nominal=nominal) : randomization_config
    return AerobrakingScenarioConfig(
        constants.phase,
        constants.final_apoapsis_radius_m,
        constants.r_norm_m,
        constants.initial_apoapsis_radius_m,
        constants.nominal_periapsis_altitude_m,
        deg2rad(constants.nominal_inclination_deg),
        deg2rad(constants.nominal_raan_deg),
        deg2rad(constants.nominal_argument_of_periapsis_deg),
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
        hp = config.nominal_periapsis_altitude_m + uniform_jitter(rng, 1.5e3)
        inc = deg2rad(rand_cfg.nonnominal_inclination_low_deg +
                      rand(rng) * (rand_cfg.nonnominal_inclination_high_deg -
                                   rand_cfg.nonnominal_inclination_low_deg))
        raan = deg2rad(rand_cfg.nonnominal_angle_low_deg +
                       rand(rng) * (rand_cfg.nonnominal_angle_high_deg -
                                    rand_cfg.nonnominal_angle_low_deg))
        omega = config.nominal_argument_of_periapsis_rad + deg2rad(uniform_jitter(rng, 15.0))
        epoch = config.nominal_epoch + Day(rand(rng, 0:27)) + Hour(rand(rng, 0:23))
    end
    return AerobrakingDecisionState(
        apoapsis_radius_m = ra,
        periapsis_altitude_m = hp,
        inclination_rad = inc,
        raan_rad = raan,
        argument_of_periapsis_rad = omega,
        epoch = epoch,
    )
end

reset_scenario(backend::SpaceAGORAAerobrakingBackend, rng::AbstractRNG) =
    reset_scenario(backend.config, rng)

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

function paper_surrogate_pass(config::AerobrakingScenarioConfig, state::AerobrakingDecisionState,
                              action::AerobrakingAction, rng::AbstractRNG)
    periapsis_after_maneuver = apply_apoapsis_maneuver(config, state, action)
    periapsis_after_maneuver = clamp(periapsis_after_maneuver, 50e3, 180e3)

    rp = config.mars_radius_m + periapsis_after_maneuver
    a_before_drag = (state.apoapsis_radius_m + rp) / 2
    periapsis_velocity = sqrt(config.mu_m3_s2 * (2 / rp - 1 / a_before_drag))
    density = exponential_mars_density(config, periapsis_after_maneuver)
    heat_rate = paper_heat_rate_w_cm2(config, periapsis_after_maneuver, periapsis_velocity)
    drag_time = drag_passage_duration_s(periapsis_after_maneuver)

    heat_scale = clamp((max(heat_rate, 1e-6) / config.heat_nominal_w_cm2)^1.15, 0.05, 4.0)
    apoapsis_drop = config.base_apoapsis_decay_m + config.heat_decay_gain_m * heat_scale
    if config.randomization_config.process_noise
        apoapsis_drop *= 1 + config.randomization_config.process_noise_scale * randn(rng)
    end
    next_apoapsis = state.apoapsis_radius_m - max(0.0, apoapsis_drop)

    period = orbital_period_s(config, max(next_apoapsis, config.mars_radius_m + periapsis_after_maneuver + 1), periapsis_after_maneuver)
    elapsed_s = max(period, drag_time)
    raan, omega = j2_angle_update(config, state, max(next_apoapsis, rp + 1), periapsis_after_maneuver, elapsed_s)
    next_epoch = state.epoch + Millisecond(round(Int, elapsed_s * 1000))

    total_delta_v = state.total_delta_v_mps + action.magnitude_mps
    maneuver_count = state.maneuver_count + (action.magnitude_mps > 1e-9 ? 1 : 0)

    next_state = AerobrakingDecisionState(
        pass_index = state.pass_index + 1,
        apoapsis_radius_m = next_apoapsis,
        periapsis_altitude_m = periapsis_after_maneuver,
        inclination_rad = state.inclination_rad,
        raan_rad = raan,
        argument_of_periapsis_rad = omega,
        epoch = next_epoch,
        mission_elapsed_s = state.mission_elapsed_s + elapsed_s,
        total_delta_v_mps = total_delta_v,
        maneuver_count = maneuver_count,
        previous_density_kg_m3 = density,
        previous_heat_rate_w_cm2 = heat_rate,
        last_drag_passage_time_s = drag_time,
    )
    return next_state
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
