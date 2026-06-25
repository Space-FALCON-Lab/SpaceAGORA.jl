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
end

@testset "normalization" begin
    bounds = paper_normalization_bounds("Main")
    low_obs = PaperObservation(400, 3.3962e6, 85e3, 0, 1.57, 1.39, 730486, 0, 0)
    high_obs = PaperObservation(1200, 10100e3, 145e3, pi, 3.14, 1.75, 731216, 1.5e-7, 0.5)
    @test normalize_observation(low_obs, bounds) ≈ zeros(Float32, 9)
    @test normalize_observation(high_obs, bounds) ≈ ones(Float32, 9)
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
