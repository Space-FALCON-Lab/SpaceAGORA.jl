using Dates
using Random
using Test
using SpaceAGORA_RL

mutable struct FakeEphemerisGrowthState
    enabled::Bool
    next_growth_t_s::Float64
end

@testset "actions" begin
    @test action_count() == 13
    @test zero_action_index() == 7
    @test action_from_index(1).lowers_periapsis
    @test action_from_index(13).raises_periapsis
    @test action_from_index(zero_action_index()).delta_v_mps == 0
    @test action_from_index(1).phi_deg == 0
    @test action_from_index(13).phi_deg == 180
    @test action_from_delta_v(-0.123).lowers_periapsis
    @test action_from_delta_v(0.123).raises_periapsis
    @test action_from_delta_v(0.123).index == nearest_action_index(0.123)
    continuous = continuous_action_from_delta_v(-2.375)
    @test continuous.delta_v_mps == -2.375
    @test continuous.magnitude_mps == 2.375
    @test continuous.lowers_periapsis
    @test continuous.index == nearest_action_index(-2.375)
end

@testset "normalization" begin
    bounds = paper_normalization_bounds("Main")
    low_obs = PaperObservation(
        400,
        bounds.lower[2],
        85e3,
        deg2rad(0.0),
        deg2rad(90.0),
        deg2rad(79.0),
        730486,
        0,
        0,
    )
    high_obs = PaperObservation(
        1200,
        bounds.upper[2],
        145e3,
        deg2rad(360.0),
        deg2rad(180.0),
        deg2rad(100.0),
        731216,
        1.5e-7,
        0.5,
    )
    @test normalize_observation(low_obs, bounds) ≈ zeros(Float32, 9)
    @test normalize_observation(high_obs, bounds) ≈ ones(Float32, 9)
end

@testset "pr_drl randomized initial geometry" begin
    raw = load_config(default_config_path())
    config = resolve_config(raw; source_path=default_config_path())
    @test config.training.algorithm == :pr_drl
    @test !config.scenario.randomization_config.nominal
    @test config.scenario.randomization_config.process_noise
    @test config.scenario.termination_config.out_of_passage_periapsis_altitude_m == 135e3
    @test config.scenario.randomization_config.periapsis_jitter_m == 2500.0
    @test config.scenario.randomization_config.initial_date_start == Date(2001, 12, 1)
    @test config.scenario.randomization_config.initial_date_days == 31
    @test config.scenario.randomization_config.initial_true_anomaly_jitter_deg == 0.025
    rng = MersenneTwister(4)
    state = reset_scenario(config.scenario, rng)
    @test 88.6 <= rad2deg(state.inclination_rad) <= 98.6
    @test 60.0 <= rad2deg(state.argument_of_periapsis_rad) <= 90.0
    @test 110.0 <= rad2deg(state.raan_rad) <= 120.0
    @test Date(2001, 12, 1) <= Date(state.epoch) <= Date(2001, 12, 31)
    @test abs(rad2deg(state.true_anomaly_rad) - 180.0) <= 0.025
    @test abs(state.periapsis_altitude_m - config.scenario.nominal_periapsis_altitude_m) <=
          config.scenario.randomization_config.periapsis_jitter_m
    @test 0.9 <= state.aerodynamic_cd_scale <= 1.1
    @test 0.9 <= state.aerodynamic_cl_scale <= 1.1
    @test config.scenario.randomization_config.marsgram_seed_base <= state.gram_seed <=
          config.scenario.randomization_config.marsgram_seed_base + 1_000_000

    flight_config = paper_odyssey_flight_evaluation_config(process_noise_scale=0.0)
    flight_state = reset_scenario(flight_config, MersenneTwister(4))
    @test rad2deg(flight_state.inclination_rad) ≈ 93.6
    @test rad2deg(flight_state.raan_rad) ≈ 115.0
    @test rad2deg(flight_state.argument_of_periapsis_rad) ≈ 89.0
end

@testset "successful case repetition scheduling" begin
    @test SpaceAGORA_RL._next_successful_case_repeat(17, 0, false, 9) === nothing
    @test SpaceAGORA_RL._next_successful_case_repeat(17, 0, true, 0) === nothing
    @test SpaceAGORA_RL._next_successful_case_repeat(17, 0, true, 9) ==
          (seed=17, repeat_index=1)
    @test SpaceAGORA_RL._next_successful_case_repeat(17, 1, false, 9) ==
          (seed=17, repeat_index=2)
    @test SpaceAGORA_RL._next_successful_case_repeat(17, 8, true, 9) ==
          (seed=17, repeat_index=9)
    @test SpaceAGORA_RL._next_successful_case_repeat(17, 9, true, 9) === nothing

    active_config_path = joinpath(
        dirname(default_config_path()),
        "pr_drl_spaceagora_physics.toml",
    )
    active_config = resolve_config(active_config_path)
    @test active_config.ddqn.target_update == 10_000
    @test active_config.epsilon.stop == 0.01
    @test active_config.training.successful_case_repetitions == 9

    marsgram_config_path = joinpath(
        dirname(default_config_path()),
        "pr_drl_spaceagora_physics_marsgram.toml",
    )
    marsgram_config = resolve_config(marsgram_config_path)
    @test marsgram_config.ddqn.target_update == 10_000
    @test marsgram_config.epsilon.stop == 0.01
    @test marsgram_config.training.successful_case_repetitions == 9
end

@testset "paper AADS uses continuous seeded bisection maneuver" begin
    config = default_aerobraking_config(phase="Main", training=false)
    state = AerobrakingDecisionState(
        apoapsis_radius_m = config.initial_apoapsis_radius_m,
        periapsis_altitude_m = 120e3,
        inclination_rad = config.nominal_inclination_rad,
        raan_rad = config.nominal_raan_rad,
        argument_of_periapsis_rad = config.nominal_argument_of_periapsis_rad,
        epoch = config.nominal_epoch,
    )
    policy = AADSHeuristicPolicy()
    action = policy_action_index(policy, config, state, observe_state(config, state), MersenneTwister(3))
    @test action isa AerobrakingAction
    @test -6.0 <= action.delta_v_mps < -1.0
    @test !(action.delta_v_mps in PAPER_ACTIONS_MPS)
    @test SpaceAGORA_RL.predicted_heat_rate(config, state, action) ≈
          policy.target_heat_rate_w_cm2 atol=0.01
    action_index, resolved_action = SpaceAGORA_RL._resolve_policy_action(action)
    @test action_index == action.index
    @test resolved_action.delta_v_mps == action.delta_v_mps

    evaluation_config = paper_odyssey_flight_evaluation_config(max_passes=1)
    evaluation = evaluate_policy(
        policy,
        evaluation_config;
        episodes=1,
        seed=4,
        paper_protocol=false,
        protected_initialization=ProtectedInitializationConfig(enabled=false),
    )
    executed_action = only(only(evaluation.summaries).action_trace)
    @test abs(executed_action) > 1.0
    @test !(executed_action in PAPER_ACTIONS_MPS)
end

@testset "thermal status and reward" begin
    config = default_aerobraking_config(phase="Main", training=true)
    far_error = 500e3
    @test thermal_status(0.01, far_error, config.reward_config) == thermal_low
    @test thermal_status(0.15, far_error, config.reward_config) == thermal_nominal
    @test thermal_status(0.27, far_error, config.reward_config) == thermal_high
    @test thermal_status(0.35, far_error, config.reward_config) == thermal_medium
    @test thermal_status(0.50, far_error, config.reward_config) == thermal_hard

    action = action_from_index(zero_action_index())
    nominal_obs = PaperObservation(800, config.final_apoapsis_radius_m + 500e3, 92e3, 1, 2, 1.5, 730700, 1e-8, 0.15)
    flags = TerminationFlags(false, false, false, false, false, false, false)
    @test paper_reward(nominal_obs, config, action, flags, config.reward_config) > -0.2

    low_obs = PaperObservation(800, config.final_apoapsis_radius_m + 500e3, 92e3, 1, 2, 1.5, 730700, 1e-8, 0.01)
    low_flags = TerminationFlags(false, false, false, false, true, true, false)
    @test paper_reward(low_obs, config, action, low_flags, config.reward_config) == -4

    medium_obs = PaperObservation(800, config.final_apoapsis_radius_m + 500e3, 92e3, 1, 2, 1.5, 730700, 1e-8, 0.35)
    medium_flags = TerminationFlags(false, false, false, false, true, true, false)
    @test paper_reward(medium_obs, config, action, medium_flags, config.reward_config) == -5

    hard_obs = PaperObservation(800, config.final_apoapsis_radius_m + 500e3, 92e3, 1, 2, 1.5, 730700, 1e-8, 0.50)
    hard_flags = TerminationFlags(false, false, false, false, true, true, false)
    @test paper_reward(hard_obs, config, action, hard_flags, config.reward_config) == -6

    success_obs = PaperObservation(800, config.final_apoapsis_radius_m + 2e3, 92e3, 1, 2, 1.5, 730700, 1e-8, 0.15)
    success_flags = TerminationFlags(true, false, false, false, false, true, false)
    @test paper_reward(success_obs, config, action, success_flags, config.reward_config) > 9
    @test config.reward_config.low_heat_action_cap_mps == 1.0

    rounded_success_obs = PaperObservation(
        800,
        config.final_apoapsis_radius_m + 2.6e3,
        92e3,
        1,
        2,
        1.5,
        730700,
        1e-8,
        0.50,
    )
    rounded_success_flags = TerminationFlags(true, false, false, false, true, true, false)
    @test paper_reward(
        rounded_success_obs,
        config,
        action,
        rounded_success_flags,
        config.reward_config,
    ) ≈ 2.8

    undershoot_obs = PaperObservation(800, config.final_apoapsis_radius_m - 20e3, 92e3, 1, 2, 1.5, 730700, 1e-8, 0.15)
    undershoot_flags = TerminationFlags(false, true, false, false, false, true, false)
    @test paper_reward(undershoot_obs, config, action, undershoot_flags, config.reward_config) == -4
end

@testset "paper metric delta-v accounting" begin
    config = default_aerobraking_config(phase="Main", training=false, max_passes=5)
    rng = MersenneTwister(9)
    state = reset_scenario(config, rng)
    summary = empty_episode_summary()
    for action in (13, 7)
        result = step_scenario(config, state, action, rng)
        summary = SpaceAGORA_RL.update_episode_summary(summary, result)
        state = result.state
    end
    final_summary = SpaceAGORA_RL.finalize_episode_summary(summary, config)
    @test final_summary.abm_delta_v_mps ≈ 1.0
    @test final_summary.total_delta_v_mps ≈ final_summary.abm_delta_v_mps +
                                             final_summary.apoapsis_correction_delta_v_mps +
                                             final_summary.periapsis_raise_delta_v_mps
    @test final_summary.total_mission_delta_v_mps == final_summary.total_delta_v_mps
    @test final_summary.total_delta_v_mps >= final_summary.abm_delta_v_mps
end

@testset "SpaceAGORA heat rates retain W/cm2 units" begin
    @test SpaceAGORA_RL._max_heat_rate_w_cm2([0.04, 0.18, 0.11]) == 0.18
    @test SpaceAGORA_RL._max_heat_rate_w_cm2([NaN, -1.0, 0.25]) == 0.25
end

@testset "passage minimum altitude is distinct from osculating periapsis" begin
    config = default_aerobraking_config(training=false)
    state = AerobrakingDecisionState(
        apoapsis_radius_m=config.initial_apoapsis_radius_m,
        periapsis_altitude_m=84.9e3,
        inclination_rad=config.nominal_inclination_rad,
        raan_rad=config.nominal_raan_rad,
        argument_of_periapsis_rad=config.nominal_argument_of_periapsis_rad,
        epoch=config.nominal_epoch,
        previous_pass_minimum_altitude_m=91.8e3,
    )

    observation = observe_state(config, state)
    flags = classify_termination(observation, config; training=false)
    metrics = SpaceAGORA_RL.pass_metrics_from_state(state)

    @test observation.periapsis_altitude_m == 91.8e3
    @test !flags.impact
    @test metrics.periapsis_altitude_m == 91.8e3
    @test state.periapsis_altitude_m == 84.9e3
end

@testset "SpaceAGORA propagated orbital-element ordering" begin
    semimajor_axis_m = 4.0e6
    eccentricity = 0.1
    inclination_rad = 0.3
    raan_rad = 1.2
    argument_of_periapsis_rad = 2.4
    true_anomaly_rad = 0.5
    orbital_elements = [
        semimajor_axis_m,
        eccentricity,
        inclination_rad,
        raan_rad,
        argument_of_periapsis_rad,
        true_anomaly_rad,
    ]
    fake_spaceagora = (
        SimulationEngine = (
            _state_position_ii = (state, spacecraft_index) -> :position,
            _state_velocity_ii = (state, spacecraft_index) -> :velocity,
        ),
        SimulationModel = (
            ControlHooks = (
                rvtoorbitalelement = (position, velocity, planet) -> orbital_elements,
            ),
        ),
    )
    planet = (Rp_e = 3.4e6,)
    args = (environment_model = (planet = planet,),)
    config = default_aerobraking_config()
    state = reset_scenario(config, MersenneTwister(12))
    stats = SpaceAGORA_RL.SpaceAGORAPhysicsPassStats(
        max_density_kg_m3 = 1.0e-8,
        max_heat_rate_w_cm2 = 0.1,
        min_altitude_m = 91.8e3,
    )

    next_state = SpaceAGORA_RL._spaceagora_physics_next_state_from_u(
        fake_spaceagora,
        config,
        state,
        action_from_index(zero_action_index()),
        args,
        stats,
        nothing,
        100.0,
        state.periapsis_altitude_m,
    )

    @test next_state.inclination_rad == inclination_rad
    @test next_state.raan_rad == raan_rad
    @test next_state.argument_of_periapsis_rad == argument_of_periapsis_rad
    @test next_state.periapsis_altitude_m == 200e3
    @test next_state.previous_pass_minimum_altitude_m == 91.8e3
    @test observe_state(config, next_state).periapsis_altitude_m == 91.8e3
end

@testset "continuous physics campaign evaluation policy selection" begin
    config = default_aerobraking_config(
        backend_mode=:spaceagora_physics,
        spaceagora_atmosphere_model=:marsgram,
        spaceagora_gram_once_per_step=true,
        training=false,
        max_passes=3,
    )
    state = reset_scenario(config, MersenneTwister(12))
    rollout = SpaceAGORA_RL.SpaceAGORAPhysicsCampaignRollout(
        config = config,
        schedule = EpsilonSchedule(),
        ddqn_config = DDQNConfig(),
        policy_snapshot = nothing,
        evaluation_policy = NoManeuverPolicy(),
        rng = MersenneTwister(12),
        episode_index = 1,
        worker_id = 2,
        seed = 12,
        max_passes_per_campaign = 3,
        global_step_start = 0,
        train = false,
        state = state,
        norm_obs = normalize_observation(
            observe_state(config, state),
            config.normalization_bounds,
        ),
        action_index = 1,
        action = action_from_index(1),
        summary = empty_episode_summary(),
        transitions = Transition[],
    )

    selected = SpaceAGORA_RL._spaceagora_physics_campaign_select_action!(rollout)
    @test selected.index == zero_action_index()
    @test rollout.action_index == zero_action_index()
    continuous = continuous_action_from_delta_v(2.375)
    SpaceAGORA_RL._spaceagora_physics_campaign_set_action!(rollout, continuous)
    @test rollout.action.delta_v_mps == 2.375
    @test rollout.action_index == nearest_action_index(2.375)
    @test SpaceAGORA_RL._spaceagora_physics_campaign_decode_action_command(3) ==
          (3, false, 0)
    @test SpaceAGORA_RL._spaceagora_physics_campaign_decode_action_command(
        (action_index=4, protected=true, policy_version=7),
    ) == (4, true, 7)
    @test SpaceAGORA_RL._spaceagora_physics_gram_once_per_step(config)

    # A campaign initialized just before apoapsis must not count that first
    # high-altitude root as a completed atmospheric passage.
    @test !rollout.periapsis_seen_since_pass_start
    @test isnothing(
        SpaceAGORA_RL._spaceagora_physics_campaign_record_apoapsis!(
            nothing,
            rollout,
            (t=2.0,),
            1,
        ),
    )
    @test isempty(rollout.transitions)
    SpaceAGORA_RL._spaceagora_physics_campaign_mark_periapsis!(rollout, 2)
    @test !rollout.periapsis_seen_since_pass_start
    SpaceAGORA_RL._spaceagora_physics_campaign_mark_periapsis!(rollout, 1)
    @test rollout.periapsis_seen_since_pass_start

    surrogate_config = default_aerobraking_config(
        backend_mode=:spaceagora_physics,
        spaceagora_atmosphere_model=:marsgram_surrogate,
        spaceagora_gram_once_per_step=false,
    )
    @test !SpaceAGORA_RL._spaceagora_physics_gram_once_per_step(surrogate_config)
end

@testset "protected initialization accounting and exclusion" begin
    config = default_aerobraking_config(phase="Main", training=false, max_passes=3)
    rng = MersenneTwister(27)
    state = reset_scenario(config, rng)
    summary = empty_episode_summary()
    initial = SpaceAGORA_RL.run_protected_initializer(
        config,
        state,
        rng,
        summary;
        settings=ProtectedInitializationConfig(corridor_maneuver=false),
    )

    @test length(initial.results) == 1
    @test initial.summary.pass_count == 1
    @test initial.summary.protected_passes == 1
    @test initial.summary.episode_reward == 0.0
    @test initial.summary.thermal_violations == 0
    @test initial.summary.protected_trace == [true]
    @test isnan(only(initial.summary.reward_trace))
    @test length(initial.summary.heat_rate_trace) == 1

    training_config = default_aerobraking_config(
        phase="Main",
        training=true,
        max_passes=3,
    )
    tolerant_config = paper_evaluation_scenario(
        training_config;
        terminal_on_thermal_violation=false,
    )
    evaluation = evaluate_policy(
        NoManeuverPolicy(),
        tolerant_config;
        episodes=1,
        seed=31,
        paper_protocol=false,
        protected_initialization=ProtectedInitializationConfig(
            corridor_maneuver=false,
        ),
    )
    evaluated = only(evaluation.summaries)
    @test evaluated.pass_count == 3
    @test evaluated.protected_passes == 1
    @test length(evaluation.transitions) == 2
    @test count(row -> row.protected, evaluation.pass_rows) == 1
    @test isnan(first(evaluation.pass_rows).reward)
end

@testset "paper evaluation protocol defaults" begin
    @test PAPER_IID_EVALUATION_EPISODES == 40
    @test PAPER_GENERALIZATION_EVALUATION_EPISODES == 100
    @test PAPER_EVALUATION_MODES == ("conservative", "tolerant")
    config = default_aerobraking_config(training=true, max_passes=7)
    evaluation_config = paper_evaluation_scenario(config)
    @test !evaluation_config.training
    @test evaluation_config.termination_config.max_passes == 7
    mode_scenarios = paper_evaluation_mode_scenarios(evaluation_config)
    @test mode_scenarios["conservative"].termination_config.terminal_on_thermal_violation
    @test !mode_scenarios["tolerant"].termination_config.terminal_on_thermal_violation
    thermal_observation = PaperObservation(
        800,
        evaluation_config.final_apoapsis_radius_m + 500e3,
        92e3,
        1,
        2,
        1.5,
        730700,
        1e-8,
        0.50,
    )
    conservative_flags = classify_termination(
        thermal_observation,
        mode_scenarios["conservative"];
        training=false,
    )
    tolerant_flags = classify_termination(
        thermal_observation,
        mode_scenarios["tolerant"];
        training=false,
    )
    @test conservative_flags.thermal_violation
    @test conservative_flags.terminated
    @test tolerant_flags.thermal_violation
    @test !tolerant_flags.terminated
    suites = generalization_suite_configs(evaluation_config)
    @test all(!scenario.training for scenario in values(suites))
    @test suites["nominal"].randomization_config.nominal
    @test !suites["iid_randomized"].randomization_config.nominal

    spaceagora_training = default_aerobraking_config(
        backend_mode=:spaceagora_physics,
        spaceagora_atmosphere_model=:marsgram,
        spaceagora_gram_wind_mode=:perturbed,
        training=true,
        randomization_config=AerobrakingRandomizationConfig(
            nominal=false,
            marsgram_perturbation_scale=1.25,
            process_noise=true,
            process_noise_scale=0.4,
            aerodynamic_coefficient_dispersion=true,
        ),
    )
    generalization = generalization_evaluation_suite(spaceagora_training)
    @test first.(generalization) == [
        GENERALIZATION_EVALUATION_REFERENCE_CASE,
        collect(GENERALIZATION_EVALUATION_CASES)...,
    ]
    cases = Dict(generalization)
    @test !cases["iid_reference"].training
    @test all(!case.termination_config.terminal_on_thermal_violation for case in values(cases))
    @test !cases["iid_reference"].randomization_config.nominal
    @test cases["iid_reference"].randomization_config.process_noise
    @test all(case.backend_mode == :spaceagora_physics for case in values(cases))
    @test all(case.spaceagora_gram_wind_mode == :perturbed for case in values(cases))
    @test cases["nominal"].randomization_config.nominal
    @test cases["nominal"].randomization_config.apoapsis_jitter_m == 0.0
    @test !cases["nominal"].randomization_config.process_noise
    @test cases["exponential_density"].spaceagora_atmosphere_model == :exponential
    @test cases["aggressive_atmosphere"].randomization_config.marsgram_perturbation_scale == 2.5
    @test cases["aggressive_atmosphere"].spaceagora_mars_mgcm_dust_levels ==
          (0.3, 0.0, 0.0)
    @test cases["aggressive_atmosphere"].spaceagora_mars_dust_storm[3] == 3.0
    @test cases["short_campaign"].phase == "Walkout"
    @test cases["long_campaign"].phase == "Campaign"
    @test rad2deg(cases["long_campaign"].nominal_argument_of_periapsis_rad) ≈ 109.0
    accurate = cases["high_accuracy_spaceagora"].spaceagora_integration_config
    nominal_integration = cases["nominal"].spaceagora_integration_config
    @test accurate.solver_mode == nominal_integration.solver_mode
    @test accurate.split_imex_solver == nominal_integration.split_imex_solver
    @test accurate.reltol_orbit < nominal_integration.reltol_orbit
    @test accurate.dt_max_atmosphere_s < nominal_integration.dt_max_atmosphere_s
    @test_throws ArgumentError generalization_evaluation_suite(evaluation_config)
    @test_throws ArgumentError default_aerobraking_config(
        spaceagora_integration_config=SpaceAGORAIntegrationConfig(dt_max_orbit_s=0.0),
    )

    paper_training = resolve_config(
        joinpath(dirname(default_config_path()), "pr_drl_paper_surrogate.toml"),
    )
    @test paper_training.training.protected_first_pass
    @test paper_training.training.protected_initial_corridor_maneuver
    @test paper_training.scenario.termination_config.max_passes == 1000
end

@testset "epsilon schedule floors" begin
    schedule = EpsilonSchedule(
        start=1.0,
        stop=0.0,
        decay_steps=5_000,
        decay_start_step=10_000,
    )
    @test epsilon_value(schedule, 10_000) == 1.0
    @test epsilon_value(schedule, 12_500) == 0.5
    @test epsilon_value(schedule, 15_000) == 0.0
    @test epsilon_value(schedule, 1_100_000) == 0.0

    config_directory = dirname(default_config_path())
    active_physics = resolve_config(
        joinpath(config_directory, "pr_drl_spaceagora_physics.toml"),
    )
    @test active_physics.epsilon.stop == 0.01
    @test epsilon_value(active_physics.epsilon, 15_000) ≈ 0.01
    @test epsilon_value(active_physics.epsilon, 1_100_000) ≈ 0.01

    marsgram_physics = resolve_config(
        joinpath(config_directory, "pr_drl_spaceagora_physics_marsgram.toml"),
    )
    @test marsgram_physics.epsilon.stop == 0.01
    @test epsilon_value(marsgram_physics.epsilon, 15_000) ≈ 0.01
    @test epsilon_value(marsgram_physics.epsilon, 1_100_000) ≈ 0.01

    marsgram = resolve_config(
        joinpath(config_directory, "pr_drl_spaceagora_marsgram.toml"),
    )
    @test marsgram.epsilon.stop == 0.0
    @test marsgram.training.validate_checkpoints
    @test marsgram.training.validation_episodes == PAPER_IID_EVALUATION_EPISODES
    @test marsgram.training.validation_checkpoint_stride == 5
end

@testset "validation thermal attribution and checkpoint selection" begin
    config = default_aerobraking_config(training=false)
    target = config.final_apoapsis_radius_m
    summary = EpisodeSummary(
        success=false,
        target_error_m=500e3,
        thermal_violations=3,
        heat_rate_trace=[0.01, 0.04, 0.27, 0.35, 0.50],
        apoapsis_trace_m=fill(target + 500e3, 5),
        protected_trace=[true, false, false, false, false],
    )
    counts = thermal_violation_breakdown(summary, config)
    @test counts == (low=1, high=1, medium=1, hard=1, total=4)
    @test terminal_thermal_violation_type(summary, config) == thermal_hard
    episode = episode_thermal_violation_metrics(summary, config)
    @test episode.soft_thermal_violations == 1
    @test episode.terminal_thermal_violation_type == "hard"

    aggregate = aggregate_thermal_violation_metrics([summary], config)
    @test aggregate.mean_low_thermal_violations == 1.0
    @test aggregate.mean_soft_thermal_violations == 1.0
    @test aggregate.terminal_hard_violation_rate == 1.0
    @test aggregate.thermal_terminal_failure_rate == 1.0

    records = [
        (
            checkpoint="checkpoint_5000.jls",
            mode="conservative",
            success_rate=0.6,
            thermal_terminal_failure_rate=0.3,
            mean_target_error_km=100.0,
            mean_thermal_violations=0.3,
        ),
        (
            checkpoint="checkpoint_10000.jls",
            mode="conservative",
            success_rate=0.6,
            thermal_terminal_failure_rate=0.2,
            mean_target_error_km=120.0,
            mean_thermal_violations=0.2,
        ),
        (
            checkpoint="checkpoint_10000.jls",
            mode="tolerant",
            success_rate=0.8,
            thermal_terminal_failure_rate=0.0,
            mean_target_error_km=20.0,
            mean_thermal_violations=0.5,
        ),
    ]
    best = select_best_validation_checkpoint(records)
    @test best.checkpoint == "checkpoint_10000.jls"
    @test best.mode == "conservative"

    checkpoint_metrics = [
        (
            episode_reward=1.0,
            target_error_km=5.0,
            thermal_violations=0,
            total_delta_v_mps=20.0,
            maneuver_count=4,
            pass_count=100,
            mission_duration_days=10.0,
            success=true,
            impact=false,
            out_of_drag_passage=false,
        ),
        (
            episode_reward=-1.0,
            target_error_km=-20.0,
            thermal_violations=2,
            total_delta_v_mps=30.0,
            maneuver_count=6,
            pass_count=150,
            mission_duration_days=15.0,
            success=false,
            impact=false,
            out_of_drag_passage=false,
        ),
    ]
    checkpoint_result = (
        metrics=checkpoint_metrics,
        aggregate=(
            episodes=2,
            success_rate=0.5,
            mean_target_error_km=12.5,
            mean_thermal_violations=1.0,
            thermal_terminal_failure_rate=0.0,
        ),
    )
    checkpoint_record = SpaceAGORA_RL._checkpoint_validation_record(
        "checkpoint_5000.jls",
        "tolerant",
        checkpoint_result,
        Dict(
            :global_step => 5_000,
            :mean_training_loss => 2.0,
            :training_loss_sum => 20.0,
            :training_loss_count => 10,
        ),
    )
    @test checkpoint_record.global_step == 5_000
    @test checkpoint_record.success_percent == 50.0
    @test checkpoint_record.surpassed_target_percent == 50.0
    @test checkpoint_record.mean_target_error_km == 12.5
    @test checkpoint_record.mean_maneuver_count == 5.0
    @test checkpoint_record.training_loss_count == 10
end

@testset "paper finite 4 N burn duration" begin
    @test SpaceAGORA_RL.PAPER_ABM_THRUST_N == 4.0
    @test SpaceAGORA_RL.paper_finite_burn_duration_s(461.0, 1.0) ≈ 115.25
    @test SpaceAGORA_RL.paper_finite_burn_duration_s(461.0, -0.5) ≈ 57.625
    @test SpaceAGORA_RL.paper_finite_burn_duration_s(461.0, 0.0) == 0.0
end

@testset "frozen checkpoint discovery" begin
    mktempdir() do directory
        for name in ("checkpoint_10000.jls", "checkpoint_5000.jls", "checkpoint_final.jls",
                     "manifest.toml")
            touch(joinpath(directory, name))
        end
        @test basename.(frozen_checkpoint_paths(directory)) == [
            "checkpoint_5000.jls",
            "checkpoint_10000.jls",
            "checkpoint_final.jls",
        ]
    end
end

@testset "aads prediction uses a distinct campaign seed" begin
    config = paper_pr_drl_evaluation_config(process_noise_scale=0.4)
    rng = MersenneTwister(12)
    state = reset_scenario(config, rng)
    prediction_state = SpaceAGORA_RL._aads_prediction_state(state)
    @test prediction_state.gram_seed != state.gram_seed
    @test prediction_state.gram_seed == SpaceAGORA_RL._aads_prediction_state(state).gram_seed
    @test prediction_state.apoapsis_radius_m == state.apoapsis_radius_m
    predicted = SpaceAGORA_RL.predicted_heat_rate(config, state, zero_action_index())
    result = step_scenario(config, state, zero_action_index(), rng)
    deterministic_config = paper_pr_drl_evaluation_config(process_noise_scale=0.0)
    deterministic_result = step_scenario(deterministic_config, state, zero_action_index(), MersenneTwister(12))
    @test predicted ≈ result.metrics.max_heat_rate_w_cm2
    @test !isapprox(result.metrics.apoapsis_radius_m,
                    deterministic_result.metrics.apoapsis_radius_m;
                    rtol=1e-10, atol=1e-6)
end

@testset "spaceagora marsgram backend selection" begin
    geometry = SpaceAGORA_RL._paper_odyssey_spacecraft_geometry()
    panel_area_each_m2 = geometry.panel_dims[2] * geometry.panel_dims[3]
    bus_area_m2 = geometry.bus_dims[1] * geometry.bus_dims[3]
    @test geometry.bus_dims == (2.2, 2.6, 1.7)
    @test geometry.panel_dims == (0.01, 1.9, 1.91)
    @test geometry.panel_offset_y_m == 2.25
    @test 2.0 * panel_area_each_m2 ≈ 7.258
    @test bus_area_m2 + 2.0 * panel_area_each_m2 ≈ 10.998

    config = paper_pr_drl_marsgram_evaluation_config(process_noise_scale=0.0)
    @test config.backend_mode == :spaceagora_marsgram
    @test SpaceAGORA_RL._spaceagora_gram_wind_enabled(config)
    physics_config = paper_pr_drl_physics_evaluation_config(process_noise_scale=0.0)
    @test physics_config.backend_mode == :spaceagora_physics
    live = resolve_config(joinpath(dirname(default_config_path()), "pr_drl_spaceagora_marsgram.toml"))
    @test live.scenario.randomization_config.marsgram_perturbation_scale == 1.0
    @test live.scenario.randomization_config.process_noise_scale == 1.0
    physics = resolve_config(joinpath(dirname(default_config_path()), "pr_drl_spaceagora_physics.toml"))
    @test physics.scenario.backend_mode == :spaceagora_physics
    @test physics.scenario.spaceagora_atmosphere_model == :tabulated_flight
    @test physics.scenario.spaceagora_tabulated_flight_file ==
          "data/telemetry/Odyssey/odyssey_accelerometer_density.feather"
    @test physics.scenario.spaceagora_tabulated_flight_sigma == 0.0
    @test physics.scenario.spaceagora_gravity_harmonics_degree == 20
    @test physics.scenario.spaceagora_gravity_harmonics_order == 20
    @test physics.scenario.spaceagora_gravity_harmonics_file ==
          "data/Gravity_harmonics_data/Mars50c.csv"
    @test physics.training.protected_first_pass
    @test physics.training.protected_initial_corridor_maneuver
    @test physics.training.protected_first_pass_suppress_thermal_terminal
    @test physics.training.validate_checkpoints
    @test physics.training.validation_checkpoint_stride == 5
    @test physics.training.worker_backend == :threads
    @test !physics.scenario.spaceagora_gram_once_per_step
    @test !SpaceAGORA_RL._spaceagora_gram_wind_enabled(physics.scenario)
    native_gram = resolve_config(
        joinpath(dirname(default_config_path()), "pr_drl_spaceagora_physics_marsgram.toml"),
    )
    @test native_gram.scenario.spaceagora_atmosphere_model == :marsgram
    @test native_gram.scenario.spaceagora_gram_wind_mode == :perturbed
    @test native_gram.scenario.spaceagora_gram_once_per_step
    @test native_gram.training.worker_backend == :processes
    @test native_gram.training.validate_checkpoints
    @test native_gram.training.validation_checkpoint_stride == 5
    @test native_gram.scenario.randomization_config.marsgram_perturbation_scale == 1.0
    @test SpaceAGORA_RL._spaceagora_gram_wind_enabled(native_gram.scenario)
    gram_evaluation = paper_evaluation_scenario(native_gram.scenario)
    @test gram_evaluation.spaceagora_gram_wind_mode == :perturbed
    @test SpaceAGORA_RL._spaceagora_gram_wind_enabled(gram_evaluation)
    zero_wind = resolve_config(
        joinpath(dirname(default_config_path()), "pr_drl_spaceagora_physics_marsgram.toml");
        gram_wind_mode=:zero,
    )
    @test zero_wind.scenario.spaceagora_gram_wind_mode == :zero
    @test !SpaceAGORA_RL._spaceagora_gram_wind_enabled(zero_wind.scenario)
    @test SpaceAGORA_RL._spaceagora_live_needs_gramsuite(zero_wind.scenario)
    nominal_wind = resolve_config(
        joinpath(dirname(default_config_path()), "pr_drl_spaceagora_physics_marsgram.toml");
        gram_wind_mode="nominal",
    )
    @test nominal_wind.scenario.spaceagora_gram_wind_mode == :nominal
    @test SpaceAGORA_RL._spaceagora_gram_wind_enabled(nominal_wind.scenario)
    @test_throws ArgumentError resolve_config(
        joinpath(dirname(default_config_path()), "pr_drl_spaceagora_physics_marsgram.toml");
        gram_wind_mode="variable",
    )
    invalid = deepcopy(load_config(default_config_path()))
    invalid["training"]["worker_backend"] = "tasks"
    @test_throws ArgumentError resolve_config(invalid)
    invalid = deepcopy(load_config(default_config_path()))
    invalid["spaceagora_physics"] = Dict(
        "atmosphere_model" => "tabulated_flight",
        "gram_once_per_step" => true,
    )
    @test_throws ArgumentError resolve_config(invalid)
    @test_throws ArgumentError SpaceAGORA_RL.propagate_pass(
        SpaceAGORA_RL.SpaceAGORACoreAdapter(:unsupported_backend),
        config,
        reset_scenario(config, MersenneTwister(21)),
        action_from_index(zero_action_index()),
        MersenneTwister(22),
    )
end

@testset "spaceagora physics solver maxiters budget" begin
    @test SpaceAGORA_RL._spaceagora_physics_solver_maxiters(1) == 5_000_000
    @test SpaceAGORA_RL._spaceagora_physics_solver_maxiters(80) == 20_000_000
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "1234",
            "SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS" => nothing) do
        @test SpaceAGORA_RL._spaceagora_physics_solver_maxiters(80) == 1234
    end
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "1234",
            "SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS" => "5678") do
        @test SpaceAGORA_RL._spaceagora_physics_solver_maxiters(80) == 5678
    end
    withenv("SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS" => "0") do
        @test_throws ArgumentError SpaceAGORA_RL._spaceagora_physics_solver_maxiters(80)
    end
end

@testset "RL shared ephemeris coverage and injection" begin
    physics = resolve_config(
        joinpath(dirname(default_config_path()), "pr_drl_spaceagora_physics.toml"),
    )
    coverage = SpaceAGORA_RL._spaceagora_rl_ephemeris_coverage(
        physics.scenario,
        physics.training.max_passes_per_campaign,
    )
    @test coverage.start == DateTime(2001, 12, 1)
    @test coverage.latest_start == DateTime(2002, 1, 1)
    @test coverage.pass_cap == 1000
    @test coverage.dt_s == 30.0
    @test coverage.sample_count == max(2, ceil(Int, coverage.total_span_s / coverage.dt_s) + 1)
    withenv("SPACEAGORA_RL_SHARED_EPHEMERIS_MAX_SAMPLES" => nothing) do
        @test SpaceAGORA_RL._spaceagora_rl_shared_ephemeris_max_samples() == 1_000_000
        @test coverage.sample_count <=
              SpaceAGORA_RL._spaceagora_rl_shared_ephemeris_max_samples()
    end
    withenv("SPACEAGORA_RL_SHARED_EPHEMERIS_MAX_SAMPLES" => "500000") do
        @test SpaceAGORA_RL._spaceagora_rl_shared_ephemeris_max_samples() == 500_000
    end
    @test coverage.total_span_s ==
          Dates.value(coverage.latest_start - coverage.start) / 1.0e3 +
          coverage.mission_duration_s

    nominal = default_aerobraking_config(phase="Main", nominal=true, max_passes=5)
    nominal_coverage = SpaceAGORA_RL._spaceagora_rl_ephemeris_coverage(nominal, 20)
    @test nominal_coverage.start == nominal.nominal_epoch
    @test nominal_coverage.latest_start == nominal.nominal_epoch
    @test nominal_coverage.pass_cap == 5
    @test SpaceAGORA_RL._spaceagora_rl_shared_ephemeris_cache(
        nominal,
        5;
        build_if_missing=false,
    ) === nothing

    cache = SpaceAGORA_RL.SpaceAGORARLSharedEphemerisCache(
        :nbody,
        :srp,
        :frame,
        100.0,
        200.0,
        30.0,
        4,
    )
    growth = FakeEphemerisGrowthState(true, 150.0)
    shared_buffers = (
        et_start=Ref(120.0),
        nbody_ephemeris_cache=Ref{Any}(nothing),
        srp_sun_ephemeris_cache=Ref{Any}(nothing),
        planet_frame_ephemeris_cache=Ref{Any}(nothing),
        ephemeris_cache_growth_state=growth,
    )
    params = (
        shared_buffers=shared_buffers,
        args=(mission_configuration=(mission_time=50.0,),),
    )
    SpaceAGORA_RL._install_spaceagora_rl_shared_ephemeris_cache!(params, cache)
    @test shared_buffers.nbody_ephemeris_cache[] == :nbody
    @test shared_buffers.srp_sun_ephemeris_cache[] == :srp
    @test shared_buffers.planet_frame_ephemeris_cache[] == :frame
    @test !growth.enabled
    @test growth.next_growth_t_s == Inf

    out_of_range = (
        shared_buffers=merge(shared_buffers, (et_start=Ref(175.0),)),
        args=(mission_configuration=(mission_time=50.0,),),
    )
    @test_throws ArgumentError SpaceAGORA_RL._install_spaceagora_rl_shared_ephemeris_cache!(
        out_of_range,
        cache,
    )
end

@testset "spaceagora physics outer parallel routing" begin
    withenv("SPACEAGORA_RL_SHARED_EPHEMERIS_CACHE" => "1",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
            "SPACEAGORA_INNER_THREAD_BUDGET" => nothing,
            "SPACEAGORA_NBODY_EPHEMERIS_CACHE" => nothing,
            "SPACEAGORA_SRP_EPHEMERIS_CACHE" => nothing,
            "SPACEAGORA_PLANET_FRAME_CACHE" => nothing) do
        observed = SpaceAGORA_RL._with_spaceagora_physics_outer_parallelism(2) do
            (
                outer = ENV["SPACEAGORA_OUTER_PARALLEL_ACTIVE"],
                inner = ENV["SPACEAGORA_INNER_THREAD_BUDGET"],
                nbody_cache = ENV["SPACEAGORA_NBODY_EPHEMERIS_CACHE"],
                srp_cache = ENV["SPACEAGORA_SRP_EPHEMERIS_CACHE"],
                planet_frame_cache = ENV["SPACEAGORA_PLANET_FRAME_CACHE"],
            )
        end
        @test observed.outer == "1"
        @test parse(Int, observed.inner) == max(1, Threads.nthreads() ÷ 2)
        @test observed.nbody_cache == "0"
        @test observed.srp_cache == "0"
        @test observed.planet_frame_cache == "0"
        @test !haskey(ENV, "SPACEAGORA_OUTER_PARALLEL_ACTIVE")
        @test !haskey(ENV, "SPACEAGORA_INNER_THREAD_BUDGET")
        @test !haskey(ENV, "SPACEAGORA_NBODY_EPHEMERIS_CACHE")
        @test !haskey(ENV, "SPACEAGORA_SRP_EPHEMERIS_CACHE")
        @test !haskey(ENV, "SPACEAGORA_PLANET_FRAME_CACHE")
    end

    withenv("SPACEAGORA_RL_SHARED_EPHEMERIS_CACHE" => "1",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
            "SPACEAGORA_INNER_THREAD_BUDGET" => "3") do
        observed = SpaceAGORA_RL._with_spaceagora_physics_outer_parallelism(2) do
            (
                outer = ENV["SPACEAGORA_OUTER_PARALLEL_ACTIVE"],
                inner = ENV["SPACEAGORA_INNER_THREAD_BUDGET"],
            )
        end
        @test observed == (outer = "1", inner = "3")
        @test ENV["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "0"
        @test ENV["SPACEAGORA_INNER_THREAD_BUDGET"] == "3"
    end

    withenv("SPACEAGORA_RL_SHARED_EPHEMERIS_CACHE" => "0",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        observed = SpaceAGORA_RL._with_spaceagora_physics_outer_parallelism(1) do
            haskey(ENV, "SPACEAGORA_OUTER_PARALLEL_ACTIVE")
        end
        @test !observed
    end

    gram_config = default_aerobraking_config(
        backend_mode=:spaceagora_physics,
        spaceagora_atmosphere_model=:marsgram,
        spaceagora_gram_once_per_step=true,
    )
    withenv("SPACEAGORA_GRAM_ONCE_PER_STEP" => nothing,
            "SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
        observed = SpaceAGORA_RL._with_spaceagora_physics_outer_parallelism(
            2,
            gram_config,
            :processes,
        ) do
            (
                once_per_step=ENV["SPACEAGORA_GRAM_ONCE_PER_STEP"],
                inner=ENV["SPACEAGORA_INNER_THREAD_BUDGET"],
            )
        end
        @test observed == (once_per_step="1", inner="1")
        @test !haskey(ENV, "SPACEAGORA_GRAM_ONCE_PER_STEP")
    end
end

@testset "replay buffer wraparound" begin
    buffer = ReplayBuffer(9, 3)
    for i in 1:5
        obs = fill(Float32(i), 9)
        push!(buffer, Transition(obs, 1, Float32(i), obs .+ 1, false, false, i))
    end
    @test length(buffer) == 3
    batch = sample_batch(buffer, 2, MersenneTwister(1))
    @test size(batch.observations) == (9, 2)
    @test length(batch.actions) == 2
end

@testset "ddqn target" begin
    online_q = Float32[1 4; 3 2]
    target_q = Float32[10 40; 30 20]
    rewards = Float32[1, 2]
    terminated = Bool[false, true]
    truncated = Bool[false, false]
    targets = compute_ddqn_targets(online_q, target_q, rewards, terminated, truncated, 0.5)
    @test targets ≈ Float32[16, 2]
end

@testset "ddqn learner updates on replay batch" begin
    rng = MersenneTwister(5)
    config = DDQNConfig(hidden_dim=8, obs_dim=3, action_dim=4, batch_size=2,
                        replay_size=8, train_start=0, train_frequency=1)
    learner = DDQNLearner(rng, config)
    push!(learner.replay, Transition(Float32[0.1, 0.2, 0.3], 1, 0.5f0,
                                     Float32[0.2, 0.3, 0.4], false, false, 1))
    push!(learner.replay, Transition(Float32[0.3, 0.2, 0.1], 3, -0.2f0,
                                     Float32[0.4, 0.3, 0.2], true, false, 2))
    before = copy(learner.online.W1)
    loss = SpaceAGORA_RL.train_step!(learner, rng)
    @test isfinite(loss)
    @test learner.train_steps == 1
    @test SpaceAGORA_RL.mean_training_loss(learner) == loss
    @test learner.online.W1 != before
end

@testset "a2c discounted returns respect campaign boundaries" begin
    rewards = Float32[1 1 1; 2 2 2]
    done = Bool[false true false; false false false]
    valid = trues(2, 3)
    bootstraps = Float32[10, 5]
    returns = compute_discounted_returns(rewards, done, valid, bootstraps, 0.5)
    @test returns[1, :] ≈ Float32[1.5, 1.0, 6.0]
    @test returns[2, :] ≈ Float32[4.125, 4.25, 4.5]
end

@testset "a2c bootstraps truncations without crossing episode resets" begin
    rewards = reshape(Float32[1, 2, 3], 1, :)
    episode_end = reshape(Bool[false, true, false], 1, :)
    terminated = reshape(Bool[false, false, false], 1, :)
    valid = trues(1, 3)
    next_values = reshape(Float32[10, 20, 30], 1, :)
    returns = compute_discounted_returns(
        rewards,
        episode_end,
        terminated,
        valid,
        next_values,
        0.5,
    )
    @test returns[1, :] ≈ Float32[7, 12, 18]

    terminated[1, 2] = true
    terminal_returns = compute_discounted_returns(
        rewards,
        episode_end,
        terminated,
        valid,
        next_values,
        0.5,
    )
    @test terminal_returns[1, :] ≈ Float32[2, 2, 18]
    @test_throws ArgumentError compute_discounted_returns(
        rewards,
        falses(1, 3),
        terminated,
        valid,
        next_values,
        0.5,
    )
end

@testset "a2c learner updates on rollout batch" begin
    rng = MersenneTwister(7)
    config = A2CConfig(hidden_dim=8, obs_dim=3, action_dim=4, segment_length=2)
    learner = A2CLearner(rng, config)
    batch = A2CRolloutBatch(
        Float32[0.1 0.2 0.3; 0.0 0.2 0.4; 1.0 0.5 0.0],
        [1, 3, 4],
        Float32[0.5, -0.2, 0.1],
    )
    before = copy(learner.actor.W1)
    loss = SpaceAGORA_RL.train_step!(learner, batch)
    @test isfinite(loss)
    @test learner.train_steps == 1
    @test SpaceAGORA_RL.mean_training_loss(learner) == loss
    @test learner.actor.W1 != before
    @test learner.policy_version == 1
    @test isfinite(learner.last_actor_gradient_norm)
    @test isfinite(learner.last_critic_gradient_norm)
    @test_throws ArgumentError SpaceAGORA_RL.train_step!(learner, batch)
end

@testset "A2C config matches PR-DRL live physics setup" begin
    config_directory = dirname(default_config_path())
    a2c_path = joinpath(config_directory, "a2c_spaceagora_physics_marsgram.toml")
    pr_drl_path = joinpath(config_directory, "pr_drl_spaceagora_physics_marsgram.toml")
    a2c_config = resolve_config(a2c_path)
    @test a2c_config.training.algorithm == :a2c
    @test a2c_config.training.worker_backend == :processes
    @test a2c_config.a2c.normalize_advantages
    @test a2c_config.actor_critic_action.mode == :continuous
    @test a2c_config.scenario.backend_mode == :spaceagora_physics
    @test a2c_config.scenario.spaceagora_atmosphere_model == :marsgram
    @test a2c_config.scenario.spaceagora_gram_wind_mode == :perturbed
    @test a2c_config.scenario.spaceagora_gram_once_per_step

    a2c_raw = load_config(a2c_path)
    pr_drl_raw = load_config(pr_drl_path)
    for section in ("scenario", "reward", "termination", "randomization", "spaceagora_physics")
        @test a2c_raw[section] == pr_drl_raw[section]
    end
    for key in (
        "seed",
        "device",
        "global_steps",
        "episodes",
        "max_passes_per_campaign",
        "n_workers",
        "worker_backend",
        "protected_first_pass",
        "protected_initial_corridor_maneuver",
        "protected_first_pass_suppress_thermal_terminal",
        "protected_corridor_low_w_cm2",
        "protected_corridor_high_w_cm2",
        "successful_case_repetitions",
    )
        @test a2c_raw["training"][key] == pr_drl_raw["training"][key]
    end

    pr_drl_raw = load_config(default_config_path())
    pr_drl_raw["a2c"]["segment_length"] = 0
    @test resolve_config(pr_drl_raw).training.algorithm == :pr_drl

    a2c_raw = deepcopy(pr_drl_raw)
    a2c_raw["training"]["algorithm"] = "a2c"
    @test_throws ArgumentError resolve_config(a2c_raw)
end

@testset "parallel a2c uses exact on-policy transition batches" begin
    mktempdir() do output_dir
        raw = load_config(default_config_path())
        raw["scenario"]["backend_mode"] = "paper_surrogate"
        raw["training"]["algorithm"] = "a2c"
        raw["training"]["global_steps"] = 3
        raw["training"]["episodes"] = 20
        raw["training"]["n_workers"] = 2
        raw["training"]["checkpoint_frequency"] = 0
        raw["training"]["progress_frequency"] = 0
        raw["training"]["output_dir"] = output_dir
        raw["training"]["protected_first_pass"] = true
        raw["a2c"]["hidden_dim"] = 8
        raw["a2c"]["segment_length"] = 2
        raw["a2c"]["normalize_advantages"] = true
        config = resolve_config(raw)
        session = build_training_session(config; run_id="a2c-parallel-test")
        result = train_parallel!(session)
        @test result.global_step == 3
        @test length(result.transitions) == 3
        @test session.learner.train_steps == session.learner.policy_version
        @test session.learner.train_steps > 0
        @test all(summary -> summary.protected_passes >= 1, result.summaries)
        payload = load_checkpoint(joinpath(result.output_dir, "checkpoint_final.jls"))
        @test payload[:algorithm] == :a2c
        @test payload[:policy_version] == session.learner.policy_version
        @test load_trained_a2c_policy(joinpath(result.output_dir, "checkpoint_final.jls")) isa
              GreedyA2CPolicy
    end
end

@testset "run manifest title distinguishes config algorithm" begin
    raw = load_config(default_config_path())
    pr_drl_config = resolve_config(deepcopy(raw); source_path="configs/aerobraking/paper_replication.toml")
    ddqn_raw = deepcopy(raw)
    ddqn_raw["training"]["algorithm"] = "ddqn"
    ddqn_config = resolve_config(ddqn_raw; source_path="configs/aerobraking/paper_replication.toml")
    a2c_raw = deepcopy(raw)
    a2c_raw["training"]["algorithm"] = "a2c"
    a2c_config = resolve_config(a2c_raw; source_path="configs/aerobraking/paper_replication.toml")

    pr_drl_manifest = RunManifest(pr_drl_config)
    ddqn_manifest = RunManifest(ddqn_config)
    a2c_manifest = RunManifest(a2c_config)

    @test SpaceAGORA_RL.run_title(pr_drl_config) == "paper_replication-pr_drl"
    @test SpaceAGORA_RL.run_title(ddqn_config) == "paper_replication-ddqn"
    @test SpaceAGORA_RL.run_title(a2c_config) == "paper_replication-a2c"
    @test endswith(pr_drl_manifest.run_id, "_paper_replication-pr_drl")
    @test endswith(ddqn_manifest.run_id, "_paper_replication-ddqn")
    @test endswith(a2c_manifest.run_id, "_paper_replication-a2c")
    @test SpaceAGORA_RL.manifest_dict(pr_drl_manifest)["algorithm"] == "pr_drl"
    @test SpaceAGORA_RL.manifest_dict(
        pr_drl_manifest;
        gram_wind_mode=:zero,
    )["gram_wind_mode"] == "zero"
    @test SpaceAGORA_RL.manifest_dict(a2c_manifest)["title"] == "paper_replication-a2c"
    @test RunManifest(a2c_config; run_id="custom-run").run_id == "custom-run"
end

@testset "backend determinism" begin
    config = default_aerobraking_config(phase="Main", training=false, max_passes=10)
    actions = [7, 5, 9, 7]

    function rollout(seed)
        rng = MersenneTwister(seed)
        state = reset_scenario(config, rng)
        metrics = AerobrakingPassMetrics[]
        for action in actions
            result = step_scenario(config, state, action, rng)
            push!(metrics, result.metrics)
            state = result.state
        end
        return metrics
    end

    a = rollout(11)
    b = rollout(11)
    @test getfield.(a, :apoapsis_radius_m) == getfield.(b, :apoapsis_radius_m)
    @test getfield.(a, :max_heat_rate_w_cm2) == getfield.(b, :max_heat_rate_w_cm2)
    @test getfield.(a, :periapsis_altitude_m) == getfield.(b, :periapsis_altitude_m)
end

@testset "training history stays disk backed" begin
    mktempdir() do directory
        summaries = SpaceAGORA_RL.DiskBackedHistory(
            EpisodeSummary,
            joinpath(directory, "summaries.jls"),
        )
        first_summary = EpisodeSummary(episode_index=1, episode_reward=2.5)
        second_summary = EpisodeSummary(episode_index=2, episode_reward=-1.0)
        push!(summaries, first_summary)
        push!(summaries, second_summary)

        @test length(summaries) == 2
        @test summaries[1].episode_reward == 2.5
        @test getfield.(summaries[1:2], :episode_index) == [1, 2]
        @test getfield.(collect(summaries), :episode_reward) == [2.5, -1.0]

        metrics = SpaceAGORA_RL.MappedHistory(summaries, summary -> summary.episode_reward)
        @test collect(metrics) == [2.5, -1.0]
        SpaceAGORA_RL.close_history!(summaries)
        @test summaries[2].episode_index == 2
    end
end

@testset "recent training stats use paper outcome metrics" begin
    summaries = [
        EpisodeSummary(
            episode_index=1,
            episode_reward=-2.0,
            thermal_violations=4,
            pass_count=10,
            target_error_m=-50e3,
        ),
        EpisodeSummary(
            episode_index=2,
            episode_reward=6.0,
            success=false, # Intentionally inconsistent: the final distance drives this metric.
            thermal_violations=0,
            pass_count=8,
            target_error_m=5e3,
        ),
        EpisodeSummary(
            episode_index=3,
            episode_reward=2.0,
            success=true,
            thermal_violations=3,
            pass_count=6,
            target_error_m=15e3,
        ),
    ]

    stats = SpaceAGORA_RL._recent_training_stats(summaries, 2)
    @test stats.mean_reward == 4.0
    @test stats.reached_goal_percent == 100.0
    @test stats.mean_thermal_violations == 1.5
    @test stats.mean_passes_to_end == 7.0
    @test stats.mean_end_distance_km == 5.0
    @test SpaceAGORA_RL._recent_training_stats(
        summaries,
        2;
        target_tolerance_m=4e3,
    ).reached_goal_percent == 0.0
    thermal_only_stats = SpaceAGORA_RL._recent_training_stats(summaries, 1)
    @test isnan(thermal_only_stats.reached_goal_percent)
    @test isnan(thermal_only_stats.mean_end_distance_km)

    empty_stats = SpaceAGORA_RL._recent_training_stats(EpisodeSummary[], 100)
    @test isnan(empty_stats.reached_goal_percent)
    @test isnan(empty_stats.mean_thermal_violations)
    @test isnan(empty_stats.mean_passes_to_end)
    @test isnan(empty_stats.mean_end_distance_km)
end

@testset "episode traces update in place" begin
    config = default_aerobraking_config(phase="Main", training=false, max_passes=2)
    rng = MersenneTwister(123)
    state = reset_scenario(config, rng)
    result = step_scenario(config, state, zero_action_index(), rng)
    summary = empty_episode_summary()
    heat_trace = summary.heat_rate_trace

    updated = SpaceAGORA_RL.update_episode_summary(summary, result)
    @test updated === summary
    @test updated.heat_rate_trace === heat_trace
    @test updated.heat_rate_trace == [result.metrics.max_heat_rate_w_cm2]
    @test updated.protected_trace == [false]
end

@testset "float32 network forward reuses its input" begin
    rng = MersenneTwister(8)
    network = init_q_network(rng; input_dim=3, hidden_dim=4, output_dim=2)
    observations = rand(rng, Float32, 3, 5)
    cache = SpaceAGORA_RL.forward_cache(network, observations)
    @test cache[1] === observations
    @test predict_q(network, observations) == cache[end]

    default_network = init_q_network(MersenneTwister(19); input_dim=3,
                                     hidden_dim=4, output_dim=2)
    explicit_network = init_q_network(MersenneTwister(19); input_dim=3,
                                      hidden_dim=4, output_dim=2,
                                      output_gain=sqrt(2))
    @test default_network.W3 == explicit_network.W3
end

include("a3c_tests.jl")
include("continuous_actor_critic_tests.jl")
include("td3_tests.jl")
include("evaluate_rl_run_checkpoint_selection.jl")
include("rpo_hypr_rl_tests.jl")
