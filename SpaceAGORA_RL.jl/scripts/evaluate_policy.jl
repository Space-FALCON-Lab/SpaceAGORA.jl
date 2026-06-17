#!/usr/bin/env julia

using SpaceAGORA_RL

function main(args=ARGS)
    config_path = isempty(args) ? default_config_path() : args[1]
    config = resolve_config(config_path)
    result = evaluate_policy(AADSHeuristicPolicy(), config.scenario;
                             episodes=config.training.episodes,
                             seed=config.training.seed,
                             policy_name="aads_heuristic")
    println(result.aggregate)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
