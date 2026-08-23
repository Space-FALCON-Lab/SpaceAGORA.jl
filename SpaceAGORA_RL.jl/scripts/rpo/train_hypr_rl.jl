#!/usr/bin/env julia

using Random
using TOML
using SpaceAGORA_RL

function _rpo_training_scenario(base::RPOHyPRRLScenario, sigma_m::Real,
                                rng::AbstractRNG)
    start = base.start_rtn .+ Float64(sigma_m) .* randn(rng, 3)
    goal = base.goal_rtn .+ Float64(sigma_m) .* randn(rng, 3)
    return RPOHyPRRLScenario(
        start_rtn=start,
        goal_rtn=goal,
        geometry=base.geometry,
        pso_config=base.pso_config,
        tracking_settings=base.tracking_settings,
        rrt_settings=base.rrt_settings,
        initial_attitude_rtn_to_body=base.initial_attitude_rtn_to_body,
        final_attitude_rtn_to_body=base.final_attitude_rtn_to_body,
    )
end

function main(args=ARGS)
    config_path = isempty(args) ?
        joinpath(@__DIR__, "..", "..", "configs", "rpo", "hypr_rl.toml") : args[1]
    raw = TOML.parsefile(config_path)
    task = raw["task"]
    scenario_config = raw["scenario"]
    training_config = raw["training"]
    config = RPOHyPRRLConfig(
        safe_distance_m=task["safe_distance_m"],
        max_translation_waypoints=task["max_translation_waypoints"],
        max_attitude_waypoints=task["max_attitude_waypoints"],
        max_edits=task["max_edits"],
        coordinate_scale_m=task["coordinate_scale_m"],
        fuel_weight=task["fuel_weight"],
        duration_weight=task["duration_weight"],
        allocation_error_weight=task["allocation_error_weight"],
    )
    scenario = build_rpo_hypr_rl_scenario(
        start_rtn=scenario_config["start_rtn"],
        goal_rtn=scenario_config["goal_rtn"],
        station_asset=Symbol(scenario_config["station_asset"]),
        station_points=scenario_config["station_points"],
        station_seed=scenario_config["station_seed"],
        station_keepout_radius_m=scenario_config["station_keepout_radius_m"],
    )
    mdp = RPOHyPRRLMDP(config, scenario)
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
        checkpoint_every_episodes=training_config["checkpoint_every_episodes"],
        checkpoint_directory=training_config["checkpoint_directory"],
    )
    sigma_m = scenario_config["position_randomization_m"]
    result = train_hypr_rl!(
        mdp;
        training=training,
        ddqn_config=ddqn,
        scenario_sampler=(rng, _) -> _rpo_training_scenario(scenario, sigma_m, rng),
    )
    final_path = joinpath(training.checkpoint_directory, "hypr_rl_final.jls")
    save_hypr_rl_checkpoint(
        final_path, result.learner, config;
        training_metadata=(episodes=training.episodes,),
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
