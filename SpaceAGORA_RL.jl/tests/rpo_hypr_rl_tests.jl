using LinearAlgebra
using TOML

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

@testset "HyPR-RL episode epsilon schedule and ETA" begin
    training = RPOHyPRRLTrainingConfig(
        episodes=100_000,
        epsilon_start=1.0,
        epsilon_stop=0.01,
        epsilon_decay_start_episode=10_000,
        epsilon_decay_end_episode=20_000,
    )
    schedule = rpo_hypr_rl_epsilon_schedule(training)
    @test epsilon_value(schedule, 10_000) == 1.0
    @test epsilon_value(schedule, 15_000) ≈ 0.505
    @test epsilon_value(schedule, 20_000) ≈ 0.01
    @test epsilon_value(schedule, 100_000) ≈ 0.01
    @test SpaceAGORA_RL._hypr_rl_eta_seconds(100, 500, 20.0) ≈ 80.0
    @test RPOHyPRRLConfig().max_translation_waypoints == 20
    @test RPOHyPRRLConfig().max_edits == 64
    @test RPOHyPRRLTrainingConfig().episodes == 100_000
    @test RPOHyPRRLTrainingConfig().epsilon_decay_end_episode == 20_000
end

@testset "HyPR-RL bounded fuel and wheel reward" begin
    config = RPOHyPRRLConfig()
    components = (fuel_ref=0.002, total=1.0e6, J_len=1.0e6)
    wheel_reference = sum(
        config.reaction_wheel_max_momentum_nms^2 / (2.0 * inertia)
        for inertia in config.reaction_wheel_inertia_kgm2
    )
    fuel_only = SpaceAGORA_RL._rpo_fuel_wheel_objective(
        components, 0.001, 0.0, config,
    )
    with_wheel = SpaceAGORA_RL._rpo_fuel_wheel_objective(
        components, 0.001, wheel_reference, config,
    )
    @test fuel_only.total ≈ 0.25
    @test with_wheel.total - fuel_only.total ≈ config.wheel_weight

    evaluation = objective -> RPOHyPRRLEvaluation(
        true, objective, objective, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, Float64[], zeros(3, 0), zeros(3, 0), zeros(4, 0),
        zeros(6), (fake=true,),
    )
    infeasible = RPOHyPRRLEvaluation(
        false, Inf, Inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        -Inf, Inf, Float64[], zeros(3, 0), zeros(3, 0), zeros(4, 0),
        zeros(6), (fake=true,),
    )
    @test SpaceAGORA_RL._rpo_terminal_reward(
        evaluation(1.0e9), evaluation(0.0), config,
    ) == config.completion_bonus
    @test SpaceAGORA_RL._rpo_terminal_reward(
        evaluation(1.0), evaluation(2.0), config,
    ) == 0.0
    @test SpaceAGORA_RL._rpo_terminal_reward(
        infeasible, evaluation(1.0), config,
    ) == config.infeasible_edit_penalty
    @test SpaceAGORA_RL._rpo_terminal_reward(
        evaluation(1.0), infeasible, config,
    ) == -config.infeasible_edit_penalty
end

@testset "HyPR-RL canonical endpoint distribution" begin
    raw = TOML.parsefile(joinpath(
        @__DIR__, "..", "configs", "rpo", "hypr_rl.toml",
    ))
    @test raw["training"]["episodes"] == 100_000
    @test raw["training"]["epsilon_decay_start_episode"] == 10_000
    @test raw["training"]["epsilon_decay_end_episode"] == 20_000
    @test raw["task"]["max_edits"] == 64
    @test raw["task"]["fuel_weight"] == 1.0
    @test raw["task"]["wheel_weight"] == 0.01
    @test raw["task"]["duration_weight"] == 0.0
    @test raw["task"]["allocation_error_weight"] == 0.0
    @test !haskey(raw["scenario"], "position_randomization_m")
    @test !haskey(raw["scenario"], "start_rtn")
    @test !haskey(raw["scenario"], "goal_rtn")

    base = build_rpo_hypr_rl_scenario(station_points=800, station_seed=740)
    sampler = build_rpo_hypr_rl_endpoint_sampler(
        base;
        safe_distance_m=raw["task"]["safe_distance_m"],
        endpoint_clearance_margin_m=
            raw["scenario"]["endpoint_clearance_margin_m"],
        endpoint_max_clearance_m=
            raw["scenario"]["endpoint_max_clearance_m"],
        min_separation_m=raw["scenario"]["min_separation_m"],
        surrounded_max_distance_m=
            raw["scenario"]["surrounded_max_distance_m"],
        max_sampling_tries=raw["scenario"]["max_sampling_tries"],
    )
    first = sample_rpo_hypr_rl_scenario(sampler, MersenneTwister(741))
    repeated = sample_rpo_hypr_rl_scenario(sampler, MersenneTwister(741))
    @test first.start_rtn == repeated.start_rtn
    @test first.goal_rtn == repeated.goal_rtn
    @test norm(first.goal_rtn - first.start_rtn) >= 1.5
    @test SpaceAGORA_RL._rpo_endpoint_clearance_m(
        first.start_rtn, first.geometry,
    ) >= 0.55
    @test SpaceAGORA_RL._rpo_endpoint_clearance_m(
        first.goal_rtn, first.geometry,
    ) >= 0.55
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
    @test plan.diagnostics.final_position_error_m == 0.0

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

@testset "HyPR-RL planner terminal evaluator" begin
    terminal_calls = Ref(0)
    terminal_evaluator = function (scenario, config, points, progress, quaternions)
        terminal_calls[] += 1
        evaluation = _fake_rpo_evaluator(
            scenario, config, points, progress, quaternions,
        )
        return RPOHyPRRLEvaluation(
            evaluation.feasible, evaluation.objective + 2.0,
            evaluation.path_cost, evaluation.propellant_used_kg,
            evaluation.duration_s, evaluation.allocation_error_impulse_mps,
            evaluation.thruster_saturation_fraction, evaluation.wheel_energy_j,
            evaluation.wheel_peak_momentum_nms, evaluation.min_clearance_m,
            evaluation.final_position_error_m, evaluation.t_ref_s,
            evaluation.r_ref_rtn, evaluation.v_ref_rtn,
            evaluation.q_ref_rtn_to_body, evaluation.thruster_impulse_ns,
            (fake=true, evaluator_mode=:terminal),
        )
    end
    mdp = _fake_rpo_mdp()
    stop = findfirst(
        action -> action.kind == :stop, rpo_editor_actions(mdp.config),
    )
    plan = hypr_rl_plan(
        mdp, (observation, mask, planner_state) -> stop;
        rng=MersenneTwister(81), terminal_evaluator=terminal_evaluator,
    )
    @test terminal_calls[] == 1
    @test plan.cost ≈ plan.seed_cost + 2.0
    @test plan.diagnostics.edit_cost ≈ plan.seed_cost
    @test plan.diagnostics.evaluator.evaluator_mode == :terminal
end

@testset "HyPR-RL infeasible stop reward remains finite" begin
    failed_evaluator = function (scenario, config, points, progress, quaternions)
        return SpaceAGORA_RL._failed_rpo_evaluation(reason=:test_infeasible)
    end
    base_mdp = _fake_rpo_mdp(max_edits=2)
    mdp = RPOHyPRRLMDP(
        base_mdp.config, base_mdp.scenario; evaluator=failed_evaluator,
    )
    state = reset_scenario(mdp, MersenneTwister(31))
    stop = findfirst(
        action -> action.kind == :stop, rpo_editor_actions(mdp.config),
    )
    result = step_scenario(mdp, state, stop, MersenneTwister(32))
    @test isfinite(result.reward)
    @test result.reward == -mdp.config.infeasible_edit_penalty
end

@testset "HyPR-RL incrementally repairs an infeasible seed" begin
    config = RPOHyPRRLConfig(
        max_translation_waypoints=1,
        max_attitude_waypoints=0,
        max_edits=4,
        safe_distance_m=0.8,
    )
    seed_path = [0.0 0.5 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    scenario = RPOHyPRRLScenario(
        start_rtn=[0.0, 0.0, 0.0],
        goal_rtn=[1.0, 0.0, 0.0],
        geometry=nothing,
        seed_path_rtn=seed_path,
    )
    evaluator = function (scenario, config, points, progress, quaternions)
        clearance = points[2, 2]
        deficit = max(0.0, config.safe_distance_m - clearance)
        feasible = deficit <= 1.0e-12
        components = (
            total=1.0 + deficit,
            J_obs=deficit,
            violation_count=feasible ? 0 : 1,
        )
        return RPOHyPRRLEvaluation(
            feasible, components.total, components.total,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, clearance, 0.0,
            Float64[], zeros(3, 0), zeros(3, 0), zeros(4, 0), zeros(6),
            (path_components=components,),
        )
    end
    mdp = RPOHyPRRLMDP(config, scenario; evaluator=evaluator)
    state = reset_scenario(mdp, MersenneTwister(33))
    move_outward = findfirst(
        action -> action.kind == :move_translation && action.slot == 1 &&
                  action.axis == 2 && action.direction == 1 && action.scale == 2,
        rpo_editor_actions(config),
    )

    first_edit = step_scenario(mdp, state, move_outward, MersenneTwister(34))
    @test first_edit.accepted
    @test !first_edit.state.evaluation.feasible
    @test first_edit.info.reason == :accepted_infeasible
    @test first_edit.reward > 0.0

    second_edit = step_scenario(
        mdp, first_edit.state, move_outward, MersenneTwister(35),
    )
    @test second_edit.accepted
    @test second_edit.state.evaluation.feasible
    @test second_edit.state.best_evaluation.feasible
    @test second_edit.reward > 0.0
end

@testset "HyPR-RL uses the original HyPR clearance penalty" begin
    config = RPOHyPRRLConfig(safe_distance_m=0.5)
    components = (
        total=42.0,
        J_obs=1.0e-6,
        violation_count=7,
        min_clearance=0.1,
    )
    evaluation = RPOHyPRRLEvaluation(
        false, 42.0, 42.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.1, 0.0, Float64[], zeros(3, 0), zeros(3, 0), zeros(4, 0),
        zeros(6), (path_components=components, clearance_penalty=42.0),
    )
    @test SpaceAGORA_RL._rpo_infeasibility_score(evaluation, config) == 42.0

    pso_components = SpaceAGORA_RL._rpo_clearance_penalty_components(
        components, 7,
    )
    @test pso_components.total == 42.0
    @test pso_components.clearance_penalty == 42.0
    @test pso_components.violation_count == 7
    @test !pso_components.retimed_feasible
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

@testset "HyPR-RL two-level terminal evaluation" begin
    commands = SpaceAGORA_RL._rpo_feedforward_acceleration_history(
        [1.0 1.0; 0.0 0.0; 0.0 0.0], zeros(3, 2), 0.1, 1.0,
    )
    @test commands[:, 1] ≈ [-0.03, 0.0, 0.0]

    edit_calls = Ref(0)
    terminal_calls = Ref(0)
    edit_evaluator = function (scenario, config, points, progress, quaternions)
        edit_calls[] += 1
        return _fake_rpo_evaluator(scenario, config, points, progress, quaternions)
    end
    terminal_evaluator = function (scenario, config, points, progress, quaternions)
        terminal_calls[] += 1
        evaluation = _fake_rpo_evaluator(
            scenario, config, points, progress, quaternions,
        )
        return RPOHyPRRLEvaluation(
            evaluation.feasible, evaluation.objective + 2.0,
            evaluation.path_cost, evaluation.propellant_used_kg,
            evaluation.duration_s, evaluation.allocation_error_impulse_mps,
            evaluation.thruster_saturation_fraction, evaluation.wheel_energy_j,
            evaluation.wheel_peak_momentum_nms, evaluation.min_clearance_m,
            evaluation.final_position_error_m, evaluation.t_ref_s,
            evaluation.r_ref_rtn, evaluation.v_ref_rtn,
            evaluation.q_ref_rtn_to_body, evaluation.thruster_impulse_ns,
            (fake=true, evaluator_mode=:terminal),
        )
    end
    base_mdp = _fake_rpo_mdp(max_edits=2)
    mdp = RPOHyPRRLMDP(
        base_mdp.config, base_mdp.scenario; evaluator=edit_evaluator,
    )
    ddqn = rpo_hypr_rl_ddqn_config(
        mdp.config; hidden_dim=8, batch_size=1, train_start=0,
        replay_size=8,
    )
    learner = DDQNLearner(
        MersenneTwister(12), ddqn;
        schedule=EpsilonSchedule(start=0.0, stop=0.0),
    )
    snapshot = SpaceAGORA_RL.cpu_network(learner.online)
    fill!(snapshot.W1, 0.0f0)
    fill!(snapshot.b1, 0.0f0)
    fill!(snapshot.W2, 0.0f0)
    fill!(snapshot.b2, 0.0f0)
    fill!(snapshot.W3, 0.0f0)
    fill!(snapshot.b3, 0.0f0)
    stop = findfirst(action -> action.kind == :stop,
                     rpo_editor_actions(mdp.config))
    snapshot.b3[stop] = 1.0f0

    rollout = SpaceAGORA_RL._run_hypr_rl_worker_episode(
        mdp, learner.schedule, learner.config, snapshot,
        1, 1, 101, 202, 0, terminal_evaluator,
    )
    @test edit_calls[] == 1
    @test terminal_calls[] == 1
    @test rollout.terminal_feasible
    @test rollout.terminal_reward_correction ≈ 0.0
    @test rollout.episode_return ≈ 0.0
    @test length(rollout.transitions) == 1
    @test rollout.transitions[end].reward ≈ 0.0f0
    @test rollout.transitions[end].terminated
end

@testset "HyPR-RL PR-DRL parallel training smoke" begin
    mdp = _fake_rpo_mdp(max_edits=2)
    ddqn = rpo_hypr_rl_ddqn_config(
        mdp.config; hidden_dim=8, batch_size=2, train_start=1,
        replay_size=32, target_update=2,
    )
    training = RPOHyPRRLTrainingConfig(
        episodes=6, seed=9, n_workers=2, worker_backend=:threads,
        successful_case_repetitions=0,
        progress_every_episodes=0, checkpoint_every_episodes=0,
    )
    result = train_hypr_rl!(mdp; training=training, ddqn_config=ddqn)
    @test length(result.episode_returns) == 6
    @test length(result.replay) > 0
    @test result.learner.train_steps > 0
    @test result.worker_backend == :threads
    @test result.active_workers == min(2, Threads.nthreads())
    @test all(result.terminal_feasible)
    @test all(iszero, result.terminal_reward_corrections)

    snapshot = SpaceAGORA_RL.cpu_network(result.learner.online)
    first_rollout = SpaceAGORA_RL._run_hypr_rl_worker_episode(
        mdp, result.learner.schedule, result.learner.config, snapshot,
        1, 1, 101, 202, 0,
    )
    repeated_rollout = SpaceAGORA_RL._run_hypr_rl_worker_episode(
        mdp, result.learner.schedule, result.learner.config, snapshot,
        1, 1, 101, 202, 0,
    )
    @test getfield.(first_rollout.transitions, :action_index) ==
          getfield.(repeated_rollout.transitions, :action_index)
    @test getfield.(first_rollout.transitions, :reward) ==
          getfield.(repeated_rollout.transitions, :reward)
    @test all(transition -> transition.action_mask[transition.action_index],
              first_rollout.transitions)

    mktempdir() do directory
        checkpoint = save_hypr_rl_checkpoint(
            joinpath(directory, "hypr_rl.jls"), result.learner, mdp.config,
        )
        payload = load_checkpoint(checkpoint)
        @test payload[:algorithm] == :pr_drl
        @test payload[:task] == :rpo_hypr_rl
    end
end


@testset "HyPR-RL successful scenario replay" begin
    mdp = _fake_rpo_mdp(max_edits=2)
    ddqn = rpo_hypr_rl_ddqn_config(
        mdp.config; hidden_dim=8, batch_size=2, train_start=100,
        replay_size=64, target_update=10,
    )
    training = RPOHyPRRLTrainingConfig(
        episodes=5, seed=19, n_workers=1, worker_backend=:threads,
        successful_case_repetitions=3,
        progress_every_episodes=0, checkpoint_every_episodes=0,
    )
    result = train_hypr_rl!(mdp; training=training, ddqn_config=ddqn)
    @test result.successful_case_repeat_indices == [0, 1, 2, 3, 0]
    @test result.scenario_seeds[1:4] == fill(20, 4)
    @test result.scenario_seeds[5] == 21
    @test all(result.terminal_feasible)
end
