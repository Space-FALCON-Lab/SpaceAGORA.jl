@testset "TD3 configuration and MarsGRAM training case" begin
    config_directory = dirname(default_config_path())
    td3_path = joinpath(config_directory, "td3_spaceagora_physics_marsgram.toml")
    a2c_path = joinpath(config_directory, "a2c_spaceagora_physics_marsgram.toml")
    config = resolve_config(td3_path)
    @test config.training.algorithm == :td3
    @test config.training.worker_backend == :processes
    @test config.scenario.backend_mode == :spaceagora_physics
    @test config.scenario.spaceagora_atmosphere_model == :marsgram
    @test config.scenario.spaceagora_gram_wind_mode == :perturbed
    @test config.td3.action_dim == 1
    @test config.td3.policy_delay == 2
    @test config.td3.tau == 0.005

    td3_raw = load_config(td3_path)
    a2c_raw = load_config(a2c_path)
    for section in ("scenario", "reward", "termination", "randomization", "spaceagora_physics")
        @test td3_raw[section] == a2c_raw[section]
    end
    manifest = SpaceAGORA_RL.manifest_dict(RunManifest(config))
    @test manifest["algorithm"] == "td3"
    @test manifest["td3"]["policy_delay"] == 2
    @test manifest["continuous_action"]["low_mps"] == -1.0
    @test manifest["continuous_action"]["high_mps"] == 1.0

    invalid = deepcopy(td3_raw)
    invalid["td3"]["tau"] = 0.0
    @test_throws ArgumentError resolve_config(invalid)
    invalid = deepcopy(td3_raw)
    invalid["td3"]["replay_size"] = 10
    @test_throws ArgumentError resolve_config(invalid)
end

@testset "continuous replay preserves exact actions" begin
    buffer = ContinuousReplayBuffer(3, 1, 2)
    for (index, action) in enumerate((0.237f0, -0.413f0, 0.619f0))
        observation = fill(Float32(index), 3)
        push!(buffer, ContinuousTransition(
            observation,
            action,
            Float32(index),
            observation .+ 1f0,
            false,
            false,
            index,
        ))
    end
    @test length(buffer) == 2
    @test sort(vec(buffer.actions[:, 1:2])) == sort(Float32[-0.413, 0.619])
    batch = sample_batch(buffer, 2, MersenneTwister(1))
    @test size(batch.actions) == (1, 2)
    @test all(action -> action in Float32[-0.413, 0.619], vec(batch.actions))

    discrete = Transition(Float32[0, 0, 0], nearest_action_index(0.237), 1f0,
                          Float32[1, 1, 1], false, false, 1)
    exact = SpaceAGORA_RL.continuous_transition(
        discrete,
        continuous_action_from_delta_v(0.237),
    )
    @test only(exact.action) == 0.237f0
    @test only(exact.action) != Float32(PAPER_ACTIONS_MPS[discrete.action_index])
end

@testset "network input gradients and Polyak updates" begin
    rng = MersenneTwister(211)
    network = init_q_network(rng; input_dim=3, hidden_dim=5, output_dim=1)
    observations = rand(rng, Float32, 3, 2)
    output_delta = ones(Float32, 1, 2)
    _, input_delta = SpaceAGORA_RL.network_gradients_and_input_delta(
        network,
        observations,
        output_delta,
    )
    original = observations[2, 1]
    step = 1f-3
    observations[2, 1] = original + step
    plus = sum(predict_q(network, observations))
    observations[2, 1] = original - step
    minus = sum(predict_q(network, observations))
    observations[2, 1] = original
    finite_difference = (plus - minus) / (2 * step)
    @test isapprox(input_delta[2, 1], finite_difference; atol=2e-3, rtol=2e-2)

    source = copy(network)
    target = copy(network)
    fill!(source.W1, 2f0)
    fill!(target.W1, 0f0)
    SpaceAGORA_RL.polyak_update!(target, source, 0.25)
    @test all(==(0.5f0), target.W1)
    @test_throws ArgumentError SpaceAGORA_RL.polyak_update!(target, source, 1.1)
end

@testset "TD3 clipped twin targets and delayed actor update" begin
    config = TD3Config(
        hidden_dim=8,
        obs_dim=3,
        action_dim=1,
        batch_size=2,
        replay_size=8,
        train_start=0,
        random_steps=0,
        target_policy_noise=0.2,
        target_noise_clip=0.1,
        discount=0.5,
        policy_delay=2,
    )
    learner = TD3Learner(MersenneTwister(223), config)
    for network in (learner.target_actor, learner.target_critic1, learner.target_critic2)
        for field in (:W1, :b1, :W2, :b2, :W3, :b3)
            fill!(getfield(network, field), 0f0)
        end
    end
    learner.target_critic1.b3[1] = 4f0
    learner.target_critic2.b3[1] = 2f0
    batch = ContinuousReplayBatch(
        Float32[0.1 0.2; 0.2 0.3; 0.3 0.4],
        reshape(Float32[0.237, -0.413], 1, :),
        Float32[1, 3],
        Float32[0.2 0.3; 0.3 0.4; 0.4 0.5],
        Bool[false, true],
        Bool[false, false],
        [1, 2],
    )
    targets = td3_targets(
        learner,
        batch,
        MersenneTwister(1);
        target_noise=fill(1f0, 1, 2),
    )
    @test targets == Float32[2, 3]

    for index in 1:2
        push!(learner.replay, ContinuousTransition(
            batch.observations[:, index],
            batch.actions[1, index],
            batch.rewards[index],
            batch.next_observations[:, index],
            batch.terminated[index],
            batch.truncated[index],
            index,
        ))
    end
    actor_before = copy(learner.actor.W1)
    critic_before = copy(learner.critic1.W1)
    SpaceAGORA_RL.train_step!(learner, MersenneTwister(2))
    @test learner.critic1.W1 != critic_before
    @test learner.actor.W1 == actor_before
    @test learner.actor_updates == 0
    SpaceAGORA_RL.train_step!(learner, MersenneTwister(3))
    @test learner.actor.W1 != actor_before
    @test learner.actor_updates == 1
    @test learner.train_steps == 2
    @test isfinite(learner.last_actor_loss)
    @test isfinite(learner.last_critic1_loss)
    @test isfinite(learner.last_critic2_loss)
end

@testset "TD3 bounded policy, checkpoint, and surrogate training smoke" begin
    learner = TD3Learner(
        MersenneTwister(227),
        TD3Config(hidden_dim=8, obs_dim=3, action_dim=1, random_steps=0),
    )
    observation = Float32[0.1, 0.2, 0.3]
    deterministic_a = select_action(learner, observation; rng=MersenneTwister(1), test=true)
    deterministic_b = select_action(learner, observation; rng=MersenneTwister(2), test=true)
    @test deterministic_a.delta_v_mps == deterministic_b.delta_v_mps
    @test -1.0 <= deterministic_a.delta_v_mps <= 1.0

    mktempdir() do output_dir
        checkpoint_path = joinpath(output_dir, "td3.jls")
        save_checkpoint(checkpoint_path, learner)
        payload = load_checkpoint(checkpoint_path)
        @test payload[:algorithm] == :td3
        @test payload[:action_bounds_mps] == (-1.0, 1.0)
        @test load_trained_td3_policy(checkpoint_path) isa GreedyTD3Policy
    end

    mktempdir() do output_dir
        raw = load_config(default_config_path())
        raw["scenario"]["backend_mode"] = "paper_surrogate"
        raw["training"]["algorithm"] = "td3"
        raw["training"]["device"] = "cpu"
        raw["training"]["global_steps"] = 4
        raw["training"]["episodes"] = 20
        raw["training"]["n_workers"] = 1
        raw["training"]["checkpoint_frequency"] = 0
        raw["training"]["progress_frequency"] = 0
        raw["training"]["protected_first_pass"] = false
        raw["training"]["output_dir"] = output_dir
        raw["td3"] = Dict{String,Any}(
            "hidden_dim" => 8,
            "batch_size" => 2,
            "replay_size" => 32,
            "train_start" => 0,
            "random_steps" => 0,
            "policy_delay" => 2,
            "exploration_noise" => 0.1,
        )
        session = build_training_session(
            resolve_config(raw);
            run_id="td3-surrogate-test",
        )
        result = train_parallel!(session)
        @test result.global_step == 4
        @test length(result.transitions) == 4
        @test all(transition -> transition isa ContinuousTransition, result.transitions)
        @test session.learner.train_steps == 3
        @test session.learner.actor_updates == 1
        checkpoint_path = joinpath(result.output_dir, "checkpoint_final.jls")
        @test load_trained_td3_policy(checkpoint_path) isa GreedyTD3Policy
        evaluated = evaluate_policy(
            load_trained_td3_policy(checkpoint_path),
            session.config.scenario;
            episodes=1,
            paper_protocol=false,
        )
        @test length(evaluated.summaries) == 1
    end
end
