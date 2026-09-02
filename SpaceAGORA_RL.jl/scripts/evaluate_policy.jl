#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    checkpoint_path = length(args) >= 2 ? args[2] : nothing
    episodes = length(args) >= 3 ?
               parse(Int, args[3]) :
               PAPER_IID_EVALUATION_EPISODES
    algorithm = checkpoint_path === nothing ? nothing :
                Symbol(get(load_checkpoint(checkpoint_path), :algorithm, :pr_drl))
    policy = if checkpoint_path === nothing
        AADSHeuristicPolicy()
    elseif algorithm == :td3
        load_trained_td3_policy(checkpoint_path)
    elseif algorithm == :a2c
        load_trained_a2c_policy(checkpoint_path)
    elseif algorithm == :a3c
        load_trained_a3c_policy(checkpoint_path)
    else
        load_trained_pr_drl_policy(checkpoint_path)
    end
    policy_name = checkpoint_path === nothing ? "aads_heuristic" : string(algorithm)
    scenario = paper_evaluation_scenario(
        config.scenario;
        max_passes=max(1000, config.scenario.termination_config.max_passes),
    )
    result = evaluate_policy(policy, scenario;
                             episodes=episodes,
                             seed=config.training.seed,
                             policy_name=policy_name)
    println(result.aggregate)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
