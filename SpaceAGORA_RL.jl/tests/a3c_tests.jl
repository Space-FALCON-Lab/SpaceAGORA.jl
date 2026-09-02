@testset "A3C configuration and live-physics parity" begin
    config_directory = dirname(default_config_path())
    a3c_path = joinpath(config_directory, "a3c_spaceagora_physics_marsgram.toml")
    a2c_path = joinpath(config_directory, "a2c_spaceagora_physics_marsgram.toml")
    config = resolve_config(a3c_path)
    @test config.training.algorithm == :a3c
    @test config.training.worker_backend == :processes
    @test config.a3c.t_max == 10
    @test config.a3c.max_policy_lag == -1
    @test config.a3c.adam_beta1 == 0.9
    @test config.a3c.adam_beta2 == 0.99
    @test config.actor_critic_action.mode == :continuous
    @test config.scenario.backend_mode == :spaceagora_physics
    @test config.scenario.spaceagora_atmosphere_model == :marsgram

    a3c_raw = load_config(a3c_path)
    a2c_raw = load_config(a2c_path)
    @test !haskey(SpaceAGORA_RL.manifest_dict(RunManifest(resolve_config(a2c_path))), "a3c")
    for section in ("scenario", "reward", "termination", "randomization", "spaceagora_physics")
        @test a3c_raw[section] == a2c_raw[section]
    end

    invalid = deepcopy(a3c_raw)
    invalid["a3c"]["t_max"] = 0
    @test_throws ArgumentError resolve_config(invalid)
    invalid = deepcopy(a3c_raw)
    invalid["a3c"]["max_policy_lag"] = -2
    @test_throws ArgumentError resolve_config(invalid)
end

@testset "A3C terminal returns and local update" begin
    rng = MersenneTwister(31)
    config = A3CConfig(
        hidden_dim=8,
        obs_dim=3,
        action_dim=4,
        discount=0.5,
        t_max=2,
    )
    learner = A3CLearner(rng, config)
    local_model = A3CLocalModel(learner)
    rollout = A3CRollout([
        Transition(Float32[0.1, 0.2, 0.3], 1, 1f0,
                   Float32[0.2, 0.3, 0.4], false, false, 1),
        Transition(Float32[0.2, 0.3, 0.4], 2, 2f0,
                   Float32[0.3, 0.4, 0.5], true, false, 2),
    ])
    batch = SpaceAGORA_RL.a3c_rollout_batch(rollout, local_model, config)
    @test batch.returns ≈ Float32[2, 2]
    before = copy(learner.actor.W1)
    loss = SpaceAGORA_RL.a3c_update!(learner, local_model, batch)
    @test isfinite(loss)
    @test learner.actor.W1 != before
    @test learner.train_steps == 1
    @test learner.policy_version == 1
    @test local_model.policy_version == learner.policy_version
    @test local_model.actor.W1 == learner.actor.W1
end

@testset "A3C accepts and records asynchronous policy lag" begin
    rng = MersenneTwister(37)
    config = A3CConfig(hidden_dim=8, obs_dim=3, action_dim=4, t_max=1)
    learner = A3CLearner(rng, config)
    first_local = A3CLocalModel(learner)
    stale_local = A3CLocalModel(learner)
    batch = A3CRolloutBatch(
        reshape(Float32[0.1, 0.2, 0.3], 3, 1),
        [1],
        Float32[0.5],
        0,
    )
    SpaceAGORA_RL.a3c_update!(learner, first_local, batch)
    SpaceAGORA_RL.a3c_update!(learner, stale_local, batch)
    @test learner.policy_version == 2
    @test learner.train_steps == 2
    @test learner.max_observed_policy_lag == 1
    @test SpaceAGORA_RL.mean_policy_lag(learner) == 0.5
    @test stale_local.policy_version == 2

    guarded = A3CLearner(
        MersenneTwister(37),
        A3CConfig(hidden_dim=8, obs_dim=3, action_dim=4, t_max=1, max_policy_lag=0),
    )
    current = A3CLocalModel(guarded)
    stale = A3CLocalModel(guarded)
    SpaceAGORA_RL.a3c_update!(guarded, current, batch)
    @test SpaceAGORA_RL.a3c_update!(guarded, stale, batch) === nothing
    @test guarded.policy_version == 1
    @test guarded.dropped_stale_updates == 1
    @test stale.policy_version == 1
end

@testset "parallel A3C trains without the A2C rollout barrier" begin
    mktempdir() do output_dir
        raw = load_config(default_config_path())
        raw["scenario"]["backend_mode"] = "paper_surrogate"
        raw["training"]["algorithm"] = "a3c"
        raw["training"]["device"] = "cpu"
        raw["training"]["global_steps"] = 7
        raw["training"]["episodes"] = 20
        raw["training"]["n_workers"] = 2
        raw["training"]["checkpoint_frequency"] = 0
        raw["training"]["progress_frequency"] = 0
        raw["training"]["output_dir"] = output_dir
        raw["training"]["protected_first_pass"] = true
        raw["a3c"] = Dict{String,Any}(
            "hidden_dim" => 8,
            "t_max" => 2,
            "normalize_advantages" => false,
        )
        config = resolve_config(raw)
        session = build_training_session(config; run_id="a3c-parallel-test")
        result = train_parallel!(session)
        @test result.global_step == 7
        @test length(result.transitions) == 7
        @test session.learner.train_steps >= 4
        @test session.learner.train_steps == session.learner.policy_version
        @test session.learner.gradient_lag_count == session.learner.train_steps
        payload = load_checkpoint(joinpath(result.output_dir, "checkpoint_final.jls"))
        @test payload[:algorithm] == :a3c
        @test payload[:global_step] == 7
        manifest = SpaceAGORA_RL.manifest_dict(session.manifest)
        @test haskey(manifest, "a3c")
        @test manifest["algorithm"] == "a3c"
        @test load_trained_a3c_policy(joinpath(result.output_dir, "checkpoint_final.jls")) isa
              GreedyA3CPolicy
    end
end
