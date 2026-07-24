#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    checkpoint_path = length(args) >= 2 ? args[2] : nothing
    episodes = length(args) >= 3 ?
               parse(Int, args[3]) :
               PAPER_IID_EVALUATION_EPISODES
    policy = checkpoint_path === nothing ?
             AADSHeuristicPolicy() :
             load_trained_pr_drl_policy(checkpoint_path)
    policy_name = checkpoint_path === nothing ? "aads_heuristic" : "pr_drl"
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
