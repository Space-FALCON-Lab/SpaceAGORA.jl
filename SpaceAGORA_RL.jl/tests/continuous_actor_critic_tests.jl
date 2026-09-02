@testset "continuous actor-critic configuration" begin
    config_directory = dirname(default_config_path())
    a2c_config = resolve_config(joinpath(
        config_directory,
        "a2c_spaceagora_physics_marsgram.toml",
    ))
    a3c_config = resolve_config(joinpath(
        config_directory,
        "a3c_spaceagora_physics_marsgram.toml",
    ))
    @test a2c_config.actor_critic_action.mode == :continuous
    @test a3c_config.actor_critic_action.mode == :continuous
    a2c_manifest = SpaceAGORA_RL.manifest_dict(RunManifest(a2c_config))
    a3c_manifest = SpaceAGORA_RL.manifest_dict(RunManifest(a3c_config))
    @test a2c_manifest["a2c"]["action_space"] == "continuous"
    @test a3c_manifest["a3c"]["action_space"] == "continuous"
    @test a2c_manifest["continuous_action"]["low_mps"] == -1.0
    @test a3c_manifest["continuous_action"]["high_mps"] == 1.0
    @test CONTINUOUS_ACTION_LOW_MPS == minimum(PAPER_ACTIONS_MPS) == -1.0
    @test CONTINUOUS_ACTION_HIGH_MPS == maximum(PAPER_ACTIONS_MPS) == 1.0

    default_raw = load_config(default_config_path())
    default_raw["training"]["algorithm"] = "a2c"
    @test resolve_config(default_raw).actor_critic_action.mode == :discrete

    invalid = deepcopy(default_raw)
    invalid["a2c"]["action_space"] = "hybrid"
    @test_throws ArgumentError resolve_config(invalid)
    invalid = deepcopy(default_raw)
    invalid["a2c"]["action_space"] = "continuous"
    invalid["a2c"]["initial_log_std"] = 2.0
    invalid["a2c"]["log_std_max"] = 1.0
    @test_throws ArgumentError resolve_config(invalid)
end

@testset "bounded Gaussian policy and output gradients" begin
    action_config = ActorCriticActionConfig(mode=:continuous, initial_log_std=-1.0)
    learner = A2CLearner(
        MersenneTwister(101),
        A2CConfig(hidden_dim=8, obs_dim=3, action_dim=4, normalize_advantages=false);
        action_config=action_config,
    )
    observation = Float32[0.1, -0.2, 0.3]
    sampled = [
        select_action(learner, observation; rng=MersenneTwister(seed))
        for seed in 1:32
    ]
    @test all(action -> action isa AerobrakingAction, sampled)
    @test all(action -> CONTINUOUS_ACTION_LOW_MPS <= action.delta_v_mps <=
                        CONTINUOUS_ACTION_HIGH_MPS, sampled)
    @test any(action -> !(action.delta_v_mps in PAPER_ACTIONS_MPS), sampled)
    deterministic_a = select_action(learner, observation; rng=MersenneTwister(1), test=true)
    deterministic_b = select_action(learner, observation; rng=MersenneTwister(2), test=true)
    @test deterministic_a.delta_v_mps == deterministic_b.delta_v_mps

    observations = Float32[0.1 0.2; -0.2 0.1; 0.3 0.4]
    actions_mps = Float32[-0.25, 0.37]
    advantages = Float32[0.7, -0.4]
    entropy_coef = 0.02
    policy_loss, entropy, output_delta =
        SpaceAGORA_RL.continuous_policy_loss_and_output_delta(
            learner.actor,
            action_config,
            observations,
            actions_mps,
            advantages,
            entropy_coef,
        )
    gradients = SpaceAGORA_RL.network_gradients_from_output_delta(
        learner.actor,
        observations,
        output_delta,
    )
    objective() = begin
        loss, ent, _ = SpaceAGORA_RL.continuous_policy_loss_and_output_delta(
            learner.actor,
            action_config,
            observations,
            actions_mps,
            advantages,
            entropy_coef,
        )
        loss - entropy_coef * ent
    end
    for index in 1:2
        original = learner.actor.b3[index]
        step = 1f-3
        learner.actor.b3[index] = original + step
        plus = objective()
        learner.actor.b3[index] = original - step
        minus = objective()
        learner.actor.b3[index] = original
        finite_difference = (plus - minus) / (2 * step)
        @test isapprox(gradients.b3[index], finite_difference; atol=2e-3, rtol=2e-2)
    end
    @test isfinite(policy_loss)
    @test isfinite(entropy)
end

@testset "continuous A2C and A3C updates" begin
    action_config = ActorCriticActionConfig(mode=:continuous)
    a2c = A2CLearner(
        MersenneTwister(103),
        A2CConfig(hidden_dim=8, obs_dim=3, action_dim=4, normalize_advantages=false);
        action_config=action_config,
    )
    a2c_batch = ContinuousA2CRolloutBatch(
        Float32[0.1 0.2; 0.2 0.3; 0.3 0.4],
        Float32[-0.25, 0.37],
        Float32[0.5, -0.2],
    )
    @test isfinite(SpaceAGORA_RL.train_step!(a2c, a2c_batch))
    @test a2c.train_steps == 1
    @test size(a2c.actor.W3, 1) == 2

    a3c = A3CLearner(
        MersenneTwister(107),
        A3CConfig(hidden_dim=8, obs_dim=3, action_dim=4, t_max=1);
        action_config=action_config,
    )
    local_model = A3CLocalModel(a3c)
    rollout = A3CRollout()
    transition = Transition(
        Float32[0.1, 0.2, 0.3],
        nearest_action_index(0.23),
        0.5f0,
        Float32[0.2, 0.3, 0.4],
        false,
        true,
        1,
    )
    SpaceAGORA_RL.push_a3c_transition!(
        rollout,
        transition,
        continuous_action_from_delta_v(0.23),
    )
    batch = SpaceAGORA_RL.a3c_rollout_batch(rollout, local_model, a3c.config)
    @test batch isa ContinuousA3CRolloutBatch
    @test batch.actions_mps == Float32[0.23]
    @test isfinite(SpaceAGORA_RL.a3c_update!(a3c, local_model, batch))
    @test a3c.train_steps == 1
    @test size(a3c.actor.W3, 1) == 2
end

@testset "continuous action worker command" begin
    action = continuous_action_from_delta_v(0.237)
    selected, protected, policy_version =
        SpaceAGORA_RL._spaceagora_physics_campaign_decode_action_command(
            (action=action, protected=false, policy_version=7),
        )
    @test selected == action
    @test !protected
    @test policy_version == 7
    @test selected.delta_v_mps == 0.237
end

@testset "continuous A2C and A3C parallel smokes" begin
    for algorithm in ("a2c", "a3c")
        mktempdir() do output_dir
            raw = load_config(default_config_path())
            raw["scenario"]["backend_mode"] = "paper_surrogate"
            raw["training"]["algorithm"] = algorithm
            raw["training"]["device"] = "cpu"
            raw["training"]["global_steps"] = 4
            raw["training"]["episodes"] = 20
            raw["training"]["n_workers"] = algorithm == "a2c" ? 1 : 2
            raw["training"]["checkpoint_frequency"] = 0
            raw["training"]["progress_frequency"] = 0
            raw["training"]["output_dir"] = output_dir
            algorithm_table = get!(raw, algorithm, Dict{String,Any}())
            algorithm_table["action_space"] = "continuous"
            algorithm_table["hidden_dim"] = 8
            algorithm_table[algorithm == "a2c" ? "segment_length" : "t_max"] = 2
            session = build_training_session(
                resolve_config(raw);
                run_id="continuous-$algorithm-test",
            )
            result = train_parallel!(session)
            @test result.global_step == 4
            @test session.learner.train_steps > 0
            @test size(session.learner.actor.W3, 1) == 2
            checkpoint_path = joinpath(result.output_dir, "checkpoint_final.jls")
            payload = load_checkpoint(checkpoint_path)
            @test payload[:action_space] == :continuous
            @test payload[:action_bounds_mps] == (-1.0, 1.0)
            policy = algorithm == "a2c" ?
                     load_trained_a2c_policy(checkpoint_path) :
                     load_trained_a3c_policy(checkpoint_path)
            @test policy.action_config.mode == :continuous
        end
    end
end
