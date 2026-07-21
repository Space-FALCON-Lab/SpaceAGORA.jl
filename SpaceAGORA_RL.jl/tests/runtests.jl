using Random
using Test
using SpaceAGORA_RL

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
    @test action_from_delta_v(0.123).index == 0
end

@testset "normalization" begin
    bounds = paper_normalization_bounds("Main")
    low_obs = PaperObservation(400, 3.3962e6, 85e3, 0, 1.57, 1.39, 730486, 0, 0)
    high_obs = PaperObservation(1200, 10100e3, 145e3, 2pi, 3.14, 1.75, 731216, 1.5e-7, 0.5)
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
    rng = MersenneTwister(4)
    state = reset_scenario(config.scenario, rng)
    @test 88.6 <= rad2deg(state.inclination_rad) <= 98.6
    @test 60.0 <= rad2deg(state.argument_of_periapsis_rad) <= 90.0
    @test 110.0 <= rad2deg(state.raan_rad) <= 120.0
    @test abs(state.periapsis_altitude_m - config.scenario.nominal_periapsis_altitude_m) <=
          config.scenario.randomization_config.periapsis_jitter_m
    @test 0.9 <= state.drag_coefficient_scale <= 1.1
    @test 0.9 <= state.lift_coefficient_scale <= 1.1
    @test state.marsgram_seed != state.marsgram_prediction_seed

    flight_config = paper_odyssey_flight_evaluation_config(process_noise_scale=0.0)
    flight_state = reset_scenario(flight_config, MersenneTwister(4))
    @test rad2deg(flight_state.inclination_rad) ≈ 93.6
    @test rad2deg(flight_state.argument_of_periapsis_rad) ≈ 115.0
    @test rad2deg(flight_state.raan_rad) ≈ 89.0
end

@testset "paper AADS uses bisection maneuver" begin
    config = default_aerobraking_config(phase="Main", training=false)
    state = AerobrakingDecisionState(
        apoapsis_radius_m = config.initial_apoapsis_radius_m,
        periapsis_altitude_m = 120e3,
        inclination_rad = config.nominal_inclination_rad,
        raan_rad = config.nominal_raan_rad,
        argument_of_periapsis_rad = config.nominal_argument_of_periapsis_rad,
        epoch = config.nominal_epoch,
    )
    action = policy_action_index(AADSHeuristicPolicy(), config, state, observe_state(config, state),
                                 MersenneTwister(3))
    @test action isa AerobrakingAction
    @test minimum(PAPER_ACTIONS_MPS) <= action.delta_v_mps <= maximum(PAPER_ACTIONS_MPS)
    @test action.index == 0
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

@testset "aads prediction is nominal under randomized actual pass" begin
    config = paper_pr_drl_evaluation_config(process_noise_scale=0.4)
    rng = MersenneTwister(12)
    state = reset_scenario(config, rng)
    predicted = SpaceAGORA_RL.predicted_heat_rate(config, state, zero_action_index())
    result = step_scenario(config, state, zero_action_index(), rng)
    @test !isapprox(predicted, result.metrics.max_heat_rate_w_cm2; rtol=1e-10, atol=1e-12)
end

@testset "spaceagora marsgram backend selection" begin
    config = paper_pr_drl_marsgram_evaluation_config(process_noise_scale=0.0)
    @test config.backend_mode == :spaceagora_marsgram
    physics_config = paper_pr_drl_physics_evaluation_config(process_noise_scale=0.0)
    @test physics_config.backend_mode == :spaceagora_physics
    live = resolve_config(joinpath(dirname(default_config_path()), "pr_drl_spaceagora_marsgram.toml"))
    @test live.scenario.randomization_config.marsgram_perturbation_scale == 1.0
    @test live.scenario.randomization_config.process_noise_scale == 1.0
    physics = resolve_config(joinpath(dirname(default_config_path()), "pr_drl_spaceagora_physics.toml"))
    @test physics.scenario.backend_mode == :spaceagora_physics
    @test_throws ArgumentError SpaceAGORA_RL.propagate_pass(
        SpaceAGORA_RL.SpaceAGORACoreAdapter(:unsupported_backend),
        config,
        reset_scenario(config, MersenneTwister(21)),
        action_from_index(zero_action_index()),
        MersenneTwister(22),
    )
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
    @test learner.actor.W1 != before
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
