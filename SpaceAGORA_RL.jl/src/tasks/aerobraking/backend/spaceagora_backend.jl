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

function _base_mars_model()
    # Keep config parsing lightweight; live SpaceAGORA is loaded only when a
    # physics simulation is actually constructed.
    return (
        Rp_e=3.396190e6,
        Rp_p=3.3762e6,
        Rp_m=3.3895e6,
        μ=0.4282837285418775e5 * 1e9,
        J2=1.96045e-3,
        g_ref=3.72076,
        ρ_ref=8.7489231e-7,
        h_ref=90.0e3,
        H=6.308278108e3,
        R=188.92,
        T_ref=150.0,
        γ=1.33,
    )
end

function mars_odyssey_phase_constants(phase::AbstractString)
    key = lowercase(String(phase))
    mars_radius = _base_mars_model().Rp_e
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
    spaceagora_atmosphere_model::Symbol
    spaceagora_tabulated_flight_file::String
    spaceagora_tabulated_flight_sigma::Float64
    spaceagora_gravity_harmonics_degree::Int
    spaceagora_gravity_harmonics_order::Int
    spaceagora_gravity_harmonics_file::String
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
                                    spaceagora_atmosphere_model::Symbol=:gram,
                                    spaceagora_tabulated_flight_file::AbstractString="data/telemetry/Odyssey/odyssey_accelerometer_density.feather",
                                    spaceagora_tabulated_flight_sigma::Real=0.0,
                                    spaceagora_gravity_harmonics_degree::Int=50,
                                    spaceagora_gravity_harmonics_order::Int=50,
                                    spaceagora_gravity_harmonics_file::AbstractString="data/Gravity_harmonics_data/Mars50c.csv",
                                    reward_config::RewardConfig=RewardConfig(),
                                    termination_config::Union{Nothing,TerminationConfig}=nothing,
                                    randomization_config::Union{Nothing,AerobrakingRandomizationConfig}=nothing)
    constants = mars_odyssey_phase_constants(phase)
    mars = _base_mars_model()
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
        mars.J2,
        mars.μ,
        mars.ρ_ref,
        mars.h_ref,
        mars.H,
        mars.R,
        mars.γ,
        mars.T_ref,
        mars.g_ref,
        20e3,
        130e3,
        92e3,
        3450.0,
        0.15,
        backend_mode,
        spaceagora_atmosphere_model,
        String(spaceagora_tabulated_flight_file),
        Float64(spaceagora_tabulated_flight_sigma),
        spaceagora_gravity_harmonics_degree,
        spaceagora_gravity_harmonics_order,
        String(spaceagora_gravity_harmonics_file),
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
    true_anomaly_rad::Float64 = pi
    epoch::DateTime
    mission_elapsed_s::Float64 = 0.0
    total_delta_v_mps::Float64 = 0.0
    maneuver_count::Int = 0
    previous_density_kg_m3::Float64 = 0.0
    previous_heat_rate_w_cm2::Float64 = 0.0
    last_drag_passage_time_s::Float64 = 400.0
    gram_seed::Int = 1001
    aerodynamic_cd_scale::Float64 = 1.0
    aerodynamic_cl_scale::Float64 = 1.0
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

function _randomized_nonnominal_epoch(config::AerobrakingScenarioConfig,
                                      rand_cfg::AerobrakingRandomizationConfig,
                                      rng::AbstractRNG)
    if rand_cfg.initial_date_start !== nothing && rand_cfg.initial_date_days > 0
        date = rand_cfg.initial_date_start::Date
        epoch = DateTime(date + Day(rand(rng, 0:(rand_cfg.initial_date_days - 1))))
    else
        epoch = config.nominal_epoch + Day(rand(rng, 0:27))
    end
    if rand_cfg.randomize_initial_time_of_day
        epoch += Hour(rand(rng, 0:23)) + Minute(rand(rng, 0:59)) + Second(rand(rng, 0:59))
    end
    return epoch
end

function reset_scenario(config::AerobrakingScenarioConfig, rng::AbstractRNG)
    rand_cfg = config.randomization_config
    aero_cd_scale, aero_cl_scale = if rand_cfg.aerodynamic_coefficient_dispersion
        1.0 + uniform_jitter(rng, randomized_cd_span(rand_cfg)),
        1.0 + uniform_jitter(rng, randomized_cl_span(rand_cfg))
    else
        1.0, 1.0
    end
    gram_seed = rand_cfg.marsgram_seed_base + rand(rng, 0:1_000_000)
    true_anomaly = pi + deg2rad(uniform_jitter(rng, rand_cfg.initial_true_anomaly_jitter_deg))
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
        epoch = _randomized_nonnominal_epoch(config, rand_cfg, rng)
    end
    return AerobrakingDecisionState(
        apoapsis_radius_m = ra,
        periapsis_altitude_m = hp,
        inclination_rad = inc,
        raan_rad = raan,
        argument_of_periapsis_rad = omega,
        true_anomaly_rad = true_anomaly,
        epoch = epoch,
        gram_seed = gram_seed,
        aerodynamic_cd_scale = aero_cd_scale,
        aerodynamic_cl_scale = aero_cl_scale,
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
        gram_seed = state.gram_seed,
        aerodynamic_cd_scale = state.aerodynamic_cd_scale,
        aerodynamic_cl_scale = state.aerodynamic_cl_scale,
    )
    return next_state
end

Base.@kwdef mutable struct SpaceAGORAPhysicsPassStats
    max_density_kg_m3::Float64 = 0.0
    max_heat_rate_w_cm2::Float64 = 0.0
    drag_entry_time_s::Float64 = NaN
    drag_exit_time_s::Float64 = NaN
    last_sample_time_s::Float64 = NaN
    in_drag_passage::Bool = false
    min_altitude_m::Float64 = Inf
end

function _spaceagora_solution_satellite_state(u)
    if hasproperty(u, :sc)
        return u.sc[1]
    end
    return u
end

function _state_heat_rate_w_cm2(spaceagora, integrator, sat_idx::Int)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    try
        rates = Base.invokelatest(
            getproperty(callbacks, :_compute_stage_heat_rates!),
            integrator.p,
            getproperty(integrator.u, :sc)[sat_idx],
            sat_idx,
            Float64(integrator.t);
            use_buffered_density=false,
        )
        isempty(rates) && return 0.0
        return maximum(Float64.(rates)) / 1.0e4
    catch
        try
            heat_rates = integrator.p.shared_buffers.heat_rates
            sat_idx <= length(heat_rates) || return 0.0
            rates = heat_rates[sat_idx]
            isempty(rates) && return 0.0
            return maximum(Float64.(rates)) / 1.0e4
        catch
            return 0.0
        end
    end
end

function _record_spaceagora_physics_sample!(spaceagora, stats::SpaceAGORAPhysicsPassStats, integrator)
    engine = getproperty(spaceagora, :SimulationEngine)
    pos = Base.invokelatest(Base.invokelatest(getproperty, engine, :_state_position_ii), integrator.u, 1)
    planet = integrator.p.args.environment_model.planet
    altitude = norm(pos) - planet.Rp_e
    t = Float64(integrator.t)
    stats.min_altitude_m = min(stats.min_altitude_m, altitude)
    stats.last_sample_time_s = t

    rho = try
        densities = integrator.p.shared_buffers.densities
        isempty(densities) ? 0.0 : Float64(densities[1])
    catch
        0.0
    end
    heat_rate = _state_heat_rate_w_cm2(spaceagora, integrator, 1)
    stats.max_density_kg_m3 = max(stats.max_density_kg_m3, rho)
    stats.max_heat_rate_w_cm2 = max(stats.max_heat_rate_w_cm2, heat_rate)

    entry_altitude = integrator.p.args.environment_model.EI * 1e3
    if altitude <= entry_altitude
        if !stats.in_drag_passage
            stats.drag_entry_time_s = t
            stats.in_drag_passage = true
        end
        stats.drag_exit_time_s = t
    elseif stats.in_drag_passage
        stats.drag_exit_time_s = t
        stats.in_drag_passage = false
    end
    return nothing
end

function _spaceagora_orbital_elements(spaceagora, pos, vel, planet)
    control_hooks = getproperty(getproperty(spaceagora, :SimulationModel), :ControlHooks)
    return Base.invokelatest(getproperty(control_hooks, :rvtoorbitalelement), pos, vel, planet)
end

function _spaceagora_orbital_element(oe, idx::Int)
    if hasfield(typeof(oe), :data)
        data = getfield(oe, :data)
        1 <= idx <= length(data) || throw(BoundsError(oe, idx))
        return data[idx]
    end
    return oe[idx]
end

function _spaceagora_physics_next_state_from_u(spaceagora,
                                               config::AerobrakingScenarioConfig,
                                               state::AerobrakingDecisionState,
                                               action::AerobrakingAction,
                                               args,
                                               stats::SpaceAGORAPhysicsPassStats,
                                               final_u,
                                               elapsed_pass_s::Real,
                                               periapsis_after_maneuver::Real)
    engine = getproperty(spaceagora, :SimulationEngine)
    planet = args.environment_model.planet
    pos = Base.invokelatest(Base.invokelatest(getproperty, engine, :_state_position_ii), final_u, 1)
    vel = Base.invokelatest(Base.invokelatest(getproperty, engine, :_state_velocity_ii), final_u, 1)
    oe = _spaceagora_orbital_elements(spaceagora, pos, vel, planet)
    a = Float64(_spaceagora_orbital_element(oe, 1))
    e = Float64(_spaceagora_orbital_element(oe, 2))
    apoapsis_radius = a * (1.0 + e)
    periapsis_altitude = a * (1.0 - e) - planet.Rp_e
    drag_time = if isfinite(stats.drag_entry_time_s) && isfinite(stats.drag_exit_time_s)
        max(0.0, stats.drag_exit_time_s - stats.drag_entry_time_s)
    else
        drag_passage_duration_s(periapsis_after_maneuver)
    end
    elapsed_s = Float64(elapsed_pass_s)
    next_epoch = state.epoch + Millisecond(round(Int, elapsed_s * 1000))
    return AerobrakingDecisionState(
        pass_index = state.pass_index + 1,
        apoapsis_radius_m = apoapsis_radius,
        periapsis_altitude_m = periapsis_altitude,
        inclination_rad = Float64(_spaceagora_orbital_element(oe, 3)),
        raan_rad = Float64(_spaceagora_orbital_element(oe, 5)),
        argument_of_periapsis_rad = Float64(_spaceagora_orbital_element(oe, 4)),
        epoch = next_epoch,
        mission_elapsed_s = state.mission_elapsed_s + elapsed_s,
        total_delta_v_mps = state.total_delta_v_mps + action.magnitude_mps,
        maneuver_count = state.maneuver_count + (action.magnitude_mps > 1e-9 ? 1 : 0),
        previous_density_kg_m3 = stats.max_density_kg_m3,
        previous_heat_rate_w_cm2 = stats.max_heat_rate_w_cm2,
        last_drag_passage_time_s = drag_time,
        gram_seed = state.gram_seed,
        aerodynamic_cd_scale = state.aerodynamic_cd_scale,
        aerodynamic_cl_scale = state.aerodynamic_cl_scale,
    )
end

function _spaceagora_physics_propagate_single_pass(config::AerobrakingScenarioConfig,
                                                   state::AerobrakingDecisionState,
                                                   action::AerobrakingAction)
    spaceagora = _load_spaceagora!(; load_gramsuite=_spaceagora_live_needs_gramsuite(config))
    args, _ = _spaceagora_physics_simulation_configuration(
        config,
        state,
        action;
        prediction=false,
        campaign_max_passes=1,
    )
    stats = SpaceAGORAPhysicsPassStats()
    stats_callback = _spaceagora_single_pass_stats_callback(spaceagora, stats)
    run_simulation_fn = Base.invokelatest(getproperty, spaceagora, :run_simulation)
    sol = Base.invokelatest(
        run_simulation_fn,
        args;
        return_solution=true,
        extra_callbacks=(stats_callback,),
    )
    isempty(sol.u) && throw(ErrorException("SpaceAGORA physics propagation returned an empty solution."))
    return _spaceagora_physics_next_state_from_u(
        spaceagora,
        config,
        state,
        action,
        args,
        stats,
        sol.u[end],
        Float64(sol.t[end]),
        periapsis_after_action_m(config, state, action),
    )
end

function _spaceagora_single_pass_stats_callback(spaceagora, stats::SpaceAGORAPhysicsPassStats)
    callbacks = getproperty(getproperty(spaceagora, :SimulationModel), :SimulationCallbacks)
    discrete_callback = Base.invokelatest(getproperty, callbacks, :DiscreteCallback)
    condition(u, t, integrator) = true
    affect!(integrator) = _record_spaceagora_physics_sample!(spaceagora, stats, integrator)
    initialize = (cb, u, t, integrator) -> affect!(integrator)
    return Base.invokelatest(discrete_callback, condition, affect!; initialize=initialize)
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
