#!/usr/bin/env julia

using Random
using TOML
using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ?
        joinpath(@__DIR__, "..", "..", "configs", "rpo", "hypr_rl.toml") : args[1]
    raw = TOML.parsefile(config_path)
    task = raw["task"]
    scenario_config = raw["scenario"]
    training_config = raw["training"]
    Symbol(training_config["algorithm"]) == :pr_drl ||
        throw(ArgumentError("HyPR-RL training currently requires algorithm = \"pr_drl\""))
    Symbol(get(task, "curve_type", "bezier")) == :bezier ||
        throw(ArgumentError("HyPR-RL translation waypoints must use curve_type = \"bezier\""))
    config = RPOHyPRRLConfig(
        safe_distance_m=task["safe_distance_m"],
        max_translation_waypoints=task["max_translation_waypoints"],
        max_attitude_waypoints=task["max_attitude_waypoints"],
        max_edits=task["max_edits"],
        coordinate_scale_m=task["coordinate_scale_m"],
        fuel_weight=task["fuel_weight"],
        duration_weight=task["duration_weight"],
        allocation_error_weight=task["allocation_error_weight"],
        wheel_weight=task["wheel_weight"],
    )
    base_scenario = build_rpo_hypr_rl_scenario(
        station_asset=Symbol(scenario_config["station_asset"]),
        station_points=scenario_config["station_points"],
        station_seed=scenario_config["station_seed"],
        station_keepout_radius_m=scenario_config["station_keepout_radius_m"],
    )
    scenario_sampler = build_rpo_hypr_rl_endpoint_sampler(
        base_scenario;
        station_asset=Symbol(scenario_config["station_asset"]),
        safe_distance_m=config.safe_distance_m,
        endpoint_clearance_margin_m=
            scenario_config["endpoint_clearance_margin_m"],
        endpoint_max_clearance_m=scenario_config["endpoint_max_clearance_m"],
        min_separation_m=scenario_config["min_separation_m"],
        surrounded_max_distance_m=scenario_config["surrounded_max_distance_m"],
        max_sampling_tries=scenario_config["max_sampling_tries"],
    )
    scenario = sample_rpo_hypr_rl_scenario(
        scenario_sampler, MersenneTwister(training_config["seed"] + 1),
    )
    mdp = RPOHyPRRLMDP(
        config, scenario; evaluator=evaluate_rpo_training_candidate,
    )
    ddqn = rpo_hypr_rl_ddqn_config(
        config;
        hidden_dim=training_config["hidden_dim"],
        learning_rate=training_config["learning_rate"],
        discount=training_config["discount"],
        batch_size=training_config["batch_size"],
        train_start=training_config["train_start"],
        replay_size=training_config["replay_size"],
        target_update=training_config["target_update"],
    )
    training = RPOHyPRRLTrainingConfig(
        episodes=training_config["episodes"],
        seed=training_config["seed"],
        n_workers=training_config["n_workers"],
        worker_backend=Symbol(training_config["worker_backend"]),
        epsilon_start=training_config["epsilon_start"],
        epsilon_stop=training_config["epsilon_stop"],
        epsilon_decay_start_episode=
            training_config["epsilon_decay_start_episode"],
        epsilon_decay_end_episode=training_config["epsilon_decay_end_episode"],
        successful_case_repetitions=
            training_config["successful_case_repetitions"],
        progress_every_episodes=training_config["progress_every_episodes"],
        checkpoint_every_episodes=training_config["checkpoint_every_episodes"],
        checkpoint_directory=training_config["checkpoint_directory"],
    )
    println("HyPR-RL path representation=bezier_control_waypoints " *
            "evaluators edit=retimed_feedforward terminal=full_lqmpc")
    println(
        "HyPR-RL endpoints=canonical_hypr_surface_shell " *
        "clearance_m=$(scenario_sampler.endpoint_min_clearance_m)-" *
        "$(scenario_sampler.endpoint_max_clearance_m) " *
        "min_separation_m=$(scenario_sampler.min_separation_m) " *
        "surrounded_max_distance_m=$(scenario_sampler.surrounded_max_distance_m)",
    )
    result = train_hypr_rl!(
        mdp;
        training=training,
        ddqn_config=ddqn,
        scenario_sampler=scenario_sampler,
        terminal_evaluator=evaluate_rpo_candidate,
    )
    final_path = joinpath(training.checkpoint_directory, "hypr_rl_final.jls")
    save_hypr_rl_checkpoint(
        final_path, result.learner, config;
        training_metadata=(
            episodes=training.episodes,
            curve_type=:bezier,
            max_bezier_waypoints=config.max_translation_waypoints,
            epsilon_start=training.epsilon_start,
            epsilon_stop=training.epsilon_stop,
            epsilon_decay_start_episode=training.epsilon_decay_start_episode,
            epsilon_decay_end_episode=training.epsilon_decay_end_episode,
            successful_case_repetitions=training.successful_case_repetitions,
            endpoint_distribution=:canonical_hypr_surface_shell,
            endpoint_min_clearance_m=scenario_sampler.endpoint_min_clearance_m,
            endpoint_max_clearance_m=scenario_sampler.endpoint_max_clearance_m,
            min_separation_m=scenario_sampler.min_separation_m,
            surrounded_max_distance_m=
                scenario_sampler.surrounded_max_distance_m,
            edit_evaluator=:retimed_feedforward,
            terminal_evaluator=:full_lqmpc,
        ),
    )
    println("wrote HyPR-RL checkpoint to ", final_path)
    println("final 100-episode mean return: ",
            sum(result.episode_returns[max(1, end - 99):end]) /
            min(100, length(result.episode_returns)))
    return result
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
