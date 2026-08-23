using LinearAlgebra

function _fake_rpo_evaluator(scenario, config, points, progress, quaternions)
    path_length = sum(
        norm(points[:, index + 1] - points[:, index])
        for index in 1:(size(points, 2) - 1)
    )
    waypoint_penalty = 0.01 * (size(points, 2) + size(quaternions, 2))
    attitude_penalty = 0.01 * sum(abs2, quaternions[1:3, :])
    objective = path_length + waypoint_penalty + attitude_penalty
    t_ref = [0.0, 1.0]
    r_ref = hcat(points[:, 1], points[:, end])
    v_ref = zeros(3, 2)
    q_ref = SpaceAGORA_RL._attitude_reference(progress, quaternions, [0.0, 1.0])
    return RPOHyPRRLEvaluation(
        true, objective, path_length, 1.0e-4 * objective, 1.0,
        0.0, 0.0, 0.0, 0.0, 10.0, 0.0,
        t_ref, r_ref, v_ref, q_ref, zeros(6), (fake=true,),
    )
end

function _fake_rpo_mdp(; max_edits=6)
    config = RPOHyPRRLConfig(
        max_translation_waypoints=2,
        max_attitude_waypoints=2,
        max_edits=max_edits,
        safe_distance_m=0.0,
    )
    path = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    scenario = RPOHyPRRLScenario(
        start_rtn=[0.0, 0.0, 0.0],
        goal_rtn=[1.0, 0.0, 0.0],
        geometry=nothing,
        seed_path_rtn=path,
    )
    return RPOHyPRRLMDP(config, scenario; evaluator=_fake_rpo_evaluator)
end

@testset "masked DDQN targets and replay" begin
    online = reshape(Float32[10.0, 2.0, 8.0], 3, 1)
    target = reshape(Float32[1.0, 5.0, 20.0], 3, 1)
    mask = reshape(Bool[false, true, false], 3, 1)
    targets = compute_ddqn_targets(
        online, target, Float32[1.0], Bool[false], Bool[false], 0.5;
        next_action_masks=mask,
    )
    @test targets == Float32[3.5]

    buffer = MaskedReplayBuffer(2, 3, 2)
    transition = MaskedTransition(
        Float32[0.0, 1.0], BitVector([true, false, true]), 3, 1.0f0,
        Float32[1.0, 0.0], BitVector([false, true, false]),
        false, false, 7,
    )
    push!(buffer, transition)
    batch = sample_batch(buffer, 1, MersenneTwister(2))
    @test batch.action_masks[:, 1] == transition.action_mask
    @test batch.next_action_masks[:, 1] == transition.next_action_mask

    ddqn = DDQNConfig(obs_dim=2, action_dim=3, hidden_dim=4,
                      replay_size=10, batch_size=1, train_start=0)
    learner = DDQNLearner(MersenneTwister(3), ddqn)
    @test select_action(
        learner, Float32[0.0, 0.0]; test=true,
        action_mask=BitVector([false, true, false]),
    ) == 2
end

@testset "HyPR-RL sequential editor and dynamic waypoints" begin
    mdp = _fake_rpo_mdp()
    state = reset_scenario(mdp, MersenneTwister(4))
    @test size(state.control_points_rtn, 2) == 2
    @test length(observe_state(mdp, state)) == rpo_hypr_rl_observation_dim(mdp.config)
    actions = rpo_editor_actions(mdp.config)

    insert_translation = findfirst(
        action -> action.kind == :insert_translation && action.slot == 1, actions,
    )
    @test valid_action_mask(mdp, state)[insert_translation]
    inserted = step_scenario(mdp, state, insert_translation, MersenneTwister(5))
    @test inserted.accepted
    @test size(inserted.state.control_points_rtn, 2) == 3

    insert_attitude = findfirst(
        action -> action.kind == :insert_attitude && action.slot == 1, actions,
    )
    attitude_inserted = step_scenario(
        mdp, inserted.state, insert_attitude, MersenneTwister(6),
    )
    @test attitude_inserted.accepted
    @test length(attitude_inserted.state.attitude_progress) == 3

    rotate = findfirst(
        action -> action.kind == :rotate_attitude && action.slot == 1, actions,
    )
    @test valid_action_mask(mdp, attitude_inserted.state)[rotate]
    rotated = step_scenario(mdp, attitude_inserted.state, rotate, MersenneTwister(7))
    @test rotated.accepted
    @test rotated.state.attitude_quaternions[:, 2] !=
          attitude_inserted.state.attitude_quaternions[:, 2]

    stop = findfirst(action -> action.kind == :stop, actions)
    plan = hypr_rl_plan(mdp, (observation, mask, planner_state) -> stop;
                        rng=MersenneTwister(8))
    @test plan.valid
    @test plan.stopped
    @test plan.cost <= plan.seed_cost + 1.0e-12

    convenience_plan = rpo_hypr_rl_plan_path(
        mdp.scenario.start_rtn, mdp.scenario.goal_rtn, nothing,
        (observation, mask, planner_state) -> stop;
        config=mdp.config,
        seed_path_rtn=mdp.scenario.seed_path_rtn,
        evaluator=_fake_rpo_evaluator,
        rng=MersenneTwister(10),
    )
    @test convenience_plan.valid
    @test convenience_plan.stopped
end

@testset "attitude changes six-thruster propellant accounting" begin
    config = RPOHyPRRLConfig(
        max_translation_waypoints=0,
        max_attitude_waypoints=0,
        thruster_min_firing_time_s=0.0,
        reaction_wheel_kp=0.0,
        reaction_wheel_kd=0.0,
    )
    command = reshape([0.005, 0.0, 0.0], 3, 1)
    identity = [0.0, 0.0, 0.0, 1.0]
    tilted = [0.0, 0.0, sin(pi / 8), cos(pi / 8)]
    identity_scenario = RPOHyPRRLScenario(
        start_rtn=zeros(3), goal_rtn=ones(3), geometry=nothing,
        seed_path_rtn=hcat(zeros(3), ones(3)),
        initial_attitude_rtn_to_body=identity,
    )
    tilted_scenario = RPOHyPRRLScenario(
        start_rtn=zeros(3), goal_rtn=ones(3), geometry=nothing,
        seed_path_rtn=hcat(zeros(3), ones(3)),
        initial_attitude_rtn_to_body=tilted,
    )
    identity_result = SpaceAGORA_RL._account_thruster_and_attitude_history(
        command, reshape(identity, 4, 1), 1.0, identity_scenario, config,
    )
    tilted_result = SpaceAGORA_RL._account_thruster_and_attitude_history(
        command, reshape(tilted, 4, 1), 1.0, tilted_scenario, config,
    )
    @test tilted_result.propellant_used_kg > 1.3 * identity_result.propellant_used_kg
end

@testset "HyPR-RL masked training smoke" begin
    mdp = _fake_rpo_mdp(max_edits=2)
    ddqn = rpo_hypr_rl_ddqn_config(
        mdp.config; hidden_dim=8, batch_size=2, train_start=1,
        replay_size=32, target_update=2,
    )
    training = RPOHyPRRLTrainingConfig(
        episodes=6, seed=9, checkpoint_every_episodes=0,
    )
    result = train_hypr_rl!(mdp; training=training, ddqn_config=ddqn)
    @test length(result.episode_returns) == 6
    @test length(result.replay) > 0
    @test result.learner.train_steps > 0
end
